"""Surface extraction and finishing with MeshLib."""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass

import meshlib.mrmeshnumpy as mrmeshnumpy
import meshlib.mrmeshpy as mrmeshpy
import numpy as np
import SimpleITK as sitk

from .defaults import SUPPORTED_MESH_EXTENSIONS
from .outputs import extension


@dataclass(frozen=True)
class SmoothingSafeguard:
    requested_iterations: int
    accepted: bool
    initial_self_intersecting_faces: int
    remaining_self_intersecting_faces: int
    initial_disoriented_faces: int
    remaining_disoriented_faces: int
    protected_vertices: int
    protected_vertex_fraction: float
    expanded_rings: int
    rms_displacement_mm: float
    max_displacement_mm: float
    rejection_reason: str | None


@dataclass(frozen=True)
class DecimationSafeguard:
    simplify_error_mm: float
    input_faces: int
    actual_faces: int
    attempted: bool
    accepted: bool
    error_introduced_mm: float
    initial_self_intersecting_faces: int
    remaining_self_intersecting_faces: int
    initial_disoriented_faces: int
    remaining_disoriented_faces: int
    repair_attempts: int
    protected_input_faces: int
    protected_input_face_fraction: float
    rejection_reason: str | None


_DECIMATION_PROTECTION_RINGS = 4
_MAX_DECIMATION_REPAIR_ATTEMPTS = 8


def index_to_physical(image: sitk.Image) -> np.ndarray:
    """Affine mapping voxel indices to the image's physical millimetres."""
    direction = np.asarray(image.GetDirection(), dtype=float).reshape(3, 3)
    spacing = np.asarray(image.GetSpacing(), dtype=float)
    origin = np.asarray(image.GetOrigin(), dtype=float)
    affine = np.eye(4, dtype=float)
    affine[:3, :3] = direction @ np.diag(spacing)
    affine[:3, 3] = origin
    return affine


def marching_cubes(image: sitk.Image, isovalue: float = 0.5) -> mrmeshpy.Mesh:
    """Extract an isosurface in voxel-index coordinates with MeshLib."""
    values_zyx = sitk.GetArrayViewFromImage(image)
    values_xyz = np.ascontiguousarray(values_zyx.transpose(2, 1, 0), dtype=np.float32)
    volume = mrmeshnumpy.simpleVolumeFrom3Darray(values_xyz)
    volume.voxelSize = mrmeshpy.Vector3f(1.0, 1.0, 1.0)
    params = mrmeshpy.MarchingCubesParams()
    params.iso = float(isovalue)
    params.lessInside = False
    mesh = mrmeshpy.marchingCubes(volume, params)

    # MeshLib places samples at voxel-cell centers [i, i+1], while SimpleITK
    # treats array values as samples at integer indices. The half-voxel shift
    # preserves the established convention before the physical-space affine.
    vertices, faces = to_arrays(mesh)
    vertices -= 0.5
    return from_arrays(vertices, faces)


def transform(mesh: mrmeshpy.Mesh, matrix: np.ndarray) -> mrmeshpy.Mesh:
    vertices, faces = to_arrays(mesh)
    linear = np.asarray(matrix[:3, :3], dtype=float)
    transformed = vertices @ linear.T
    transformed += np.asarray(matrix[:3, 3], dtype=float)
    if np.linalg.det(linear) < 0:
        faces = faces[:, ::-1]
    return from_arrays(transformed, faces)


def smooth(mesh: mrmeshpy.Mesh, iterations: int, force: float) -> mrmeshpy.Mesh:
    """Apply MeshLib relaxation while approximately preserving volume."""
    if iterations <= 0:
        return mesh
    result = mrmeshpy.Mesh(mesh)
    params = mrmeshpy.MeshRelaxParams()
    params.iterations = int(iterations)
    params.force = float(force)
    if not mrmeshpy.relaxKeepVolume(result, params):
        raise RuntimeError("MeshLib smoothing was interrupted")
    return result


def selected_self_intersecting_faces(mesh: mrmeshpy.Mesh) -> np.ndarray:
    """Mark both faces in each MeshLib self-collision pair."""
    faces = to_arrays(mesh)[1]
    pairs = mrmeshpy.findSelfCollidingTriangles(mrmeshpy.MeshPart(mesh))
    selected = np.zeros(len(faces), dtype=bool)
    for pair in pairs:
        selected[int(pair.aFace)] = True
        selected[int(pair.bFace)] = True
    return selected


def selected_disoriented_faces(mesh: mrmeshpy.Mesh) -> np.ndarray:
    """Mark faces whose local orientation is inconsistent."""
    selected = mrmeshpy.findDisorientedFaces(mesh)
    return np.asarray(mrmeshnumpy.getNumpyBitSet(selected), dtype=bool)


def _smoothing_displacement(
    original_vertices: np.ndarray,
    result_vertices: np.ndarray,
) -> tuple[float, float]:
    distances = np.linalg.norm(result_vertices - original_vertices, axis=1)
    if not len(distances):
        return 0.0, 0.0
    return float(np.sqrt(np.mean(distances**2))), float(np.max(distances))


def _smoothing_stats(
    *,
    requested_iterations: int,
    accepted: bool,
    initial_self_intersections: int,
    remaining_self_intersections: int,
    initial_disoriented: int,
    remaining_disoriented: int,
    protected: np.ndarray,
    expanded_rings: int,
    original_vertices: np.ndarray,
    result_vertices: np.ndarray,
    rejection_reason: str | None = None,
) -> SmoothingSafeguard:
    protected_count = int(np.count_nonzero(protected))
    rms, maximum = _smoothing_displacement(original_vertices, result_vertices)
    return SmoothingSafeguard(
        requested_iterations=int(requested_iterations),
        accepted=accepted,
        initial_self_intersecting_faces=initial_self_intersections,
        remaining_self_intersecting_faces=remaining_self_intersections,
        initial_disoriented_faces=initial_disoriented,
        remaining_disoriented_faces=remaining_disoriented,
        protected_vertices=protected_count,
        protected_vertex_fraction=(
            float(protected_count / len(original_vertices))
            if len(original_vertices)
            else 0.0
        ),
        expanded_rings=expanded_rings,
        rms_displacement_mm=rms,
        max_displacement_mm=maximum,
        rejection_reason=rejection_reason,
    )


def protect_smoothed_surface(
    original: mrmeshpy.Mesh,
    smoothed: mrmeshpy.Mesh,
    requested_iterations: int,
) -> tuple[mrmeshpy.Mesh, SmoothingSafeguard]:
    """Keep smoothing where it preserves intersections and face orientation."""
    original_verts, original_faces = to_arrays(original)
    smoothed_verts, smoothed_faces = to_arrays(smoothed)
    if not np.array_equal(original_faces, smoothed_faces):
        raise ValueError("smoothing changed mesh connectivity")

    colliding = selected_self_intersecting_faces(smoothed)
    disoriented = selected_disoriented_faces(smoothed)
    initial_intersections = int(np.count_nonzero(colliding))
    initial_disoriented = int(np.count_nonzero(disoriented))
    protected = np.zeros(len(original_verts), dtype=bool)
    if not initial_intersections and not initial_disoriented:
        return smoothed, _smoothing_stats(
            requested_iterations=requested_iterations,
            accepted=True,
            initial_self_intersections=0,
            remaining_self_intersections=0,
            initial_disoriented=0,
            remaining_disoriented=0,
            protected=protected,
            expanded_rings=0,
            original_vertices=original_verts,
            result_vertices=smoothed_verts,
        )

    unsafe = np.logical_or(colliding, disoriented)
    protected[np.unique(original_faces[unsafe])] = True
    rings = 0
    while True:
        candidate_verts = smoothed_verts.copy()
        candidate_verts[protected] = original_verts[protected]
        candidate = from_arrays(candidate_verts, original_faces)
        remaining_colliding = selected_self_intersecting_faces(candidate)
        remaining_disoriented_faces = selected_disoriented_faces(candidate)
        remaining_intersections = int(np.count_nonzero(remaining_colliding))
        remaining_disoriented = int(np.count_nonzero(remaining_disoriented_faces))
        if not remaining_intersections and not remaining_disoriented:
            return candidate, _smoothing_stats(
                requested_iterations=requested_iterations,
                accepted=True,
                initial_self_intersections=initial_intersections,
                remaining_self_intersections=0,
                initial_disoriented=initial_disoriented,
                remaining_disoriented=0,
                protected=protected,
                expanded_rings=rings,
                original_vertices=original_verts,
                result_vertices=candidate_verts,
            )

        incident_faces = protected[original_faces].any(axis=1)
        expanded = protected.copy()
        expanded[np.unique(original_faces[incident_faces])] = True
        if np.array_equal(expanded, protected):
            reasons = []
            if remaining_intersections:
                reasons.append(
                    "%d self-intersecting face(s) remained" % remaining_intersections
                )
            if remaining_disoriented:
                reasons.append(
                    "%d disoriented face(s) remained" % remaining_disoriented
                )
            return original, _smoothing_stats(
                requested_iterations=requested_iterations,
                accepted=False,
                initial_self_intersections=initial_intersections,
                remaining_self_intersections=remaining_intersections,
                initial_disoriented=initial_disoriented,
                remaining_disoriented=remaining_disoriented,
                protected=protected,
                expanded_rings=rings,
                original_vertices=original_verts,
                result_vertices=original_verts,
                rejection_reason=" and ".join(reasons),
            )
        protected = expanded
        rings += 1


def smooth_safely(
    mesh: mrmeshpy.Mesh,
    iterations: int,
    force: float,
) -> tuple[mrmeshpy.Mesh, SmoothingSafeguard]:
    if iterations <= 0:
        vertices = to_arrays(mesh)[0]
        return mesh, _smoothing_stats(
            requested_iterations=iterations,
            accepted=True,
            initial_self_intersections=0,
            remaining_self_intersections=0,
            initial_disoriented=0,
            remaining_disoriented=0,
            protected=np.zeros(len(vertices), dtype=bool),
            expanded_rings=0,
            original_vertices=vertices,
            result_vertices=vertices,
        )
    return protect_smoothed_surface(mesh, smooth(mesh, iterations, force), iterations)


def largest_component(mesh: mrmeshpy.Mesh) -> tuple[mrmeshpy.Mesh, int]:
    components = mrmeshpy.getAllComponents(mrmeshpy.MeshPart(mesh))
    count = len(components)
    if not count:
        return mesh, 0
    largest = max(components, key=lambda region: region.count())
    return mesh.cloneRegion(largest), count


def component_count(mesh: mrmeshpy.Mesh) -> int:
    """Count disconnected surface shells without changing the mesh."""
    return len(mrmeshpy.getAllComponents(mrmeshpy.MeshPart(mesh)))


def count_defects(mesh: mrmeshpy.Mesh) -> tuple[int, int]:
    """Return boundary-edge and hole counts from MeshLib topology."""
    topology = mesh.topology
    return int(topology.findLeftBdEdges().count()), int(topology.findNumHoles())


def to_arrays(mesh: mrmeshpy.Mesh) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.asarray(mrmeshnumpy.getNumpyVerts(mesh), dtype=np.float64).copy(),
        np.asarray(mrmeshnumpy.getNumpyFaces(mesh.topology), dtype=np.int32).copy(),
    )


def from_arrays(vertices: np.ndarray, faces: np.ndarray) -> mrmeshpy.Mesh:
    return mrmeshnumpy.meshFromFacesVerts(
        np.ascontiguousarray(faces, dtype=np.int32),
        np.ascontiguousarray(vertices, dtype=np.float64),
    )


def _meshlib_topology_signature(mesh: mrmeshpy.Mesh) -> tuple[int, int, int]:
    topology = mesh.topology
    euler = (
        topology.numValidVerts()
        - topology.computeNotLoneUndirectedEdges()
        + topology.numValidFaces()
    )
    components = len(mrmeshpy.getAllComponents(mrmeshpy.MeshPart(mesh)))
    return int(components), int(topology.findNumHoles()), int(euler)


def _decimate_candidate(
    mesh: mrmeshpy.Mesh,
    simplify_error_mm: float,
    protected_faces: mrmeshpy.FaceBitSet | None = None,
) -> tuple[mrmeshpy.Mesh, tuple[int, int, int], tuple[int, int, int], float]:
    candidate = mrmeshpy.Mesh(mesh)
    before = _meshlib_topology_signature(candidate)
    settings = mrmeshpy.DecimateSettings()
    settings.strategy = mrmeshpy.DecimateStrategy.MinimizeError
    settings.maxDeletedFaces = 2**31 - 1
    settings.maxDeletedVertices = 2**31 - 1
    settings.maxError = float(simplify_error_mm)
    settings.subdivideParts = 1
    settings.packMesh = True
    if protected_faces is not None and protected_faces.any():
        region = candidate.topology.getValidFaces()
        region.subtract(protected_faces, 0)
        settings.region = region
    result = mrmeshpy.decimateMesh(candidate, settings)
    after = _meshlib_topology_signature(candidate)
    return candidate, before, after, float(result.errorIntroduced)


def _meshlib_self_intersecting_face_ids(mesh: mrmeshpy.Mesh) -> np.ndarray:
    return np.flatnonzero(selected_self_intersecting_faces(mesh))


def _meshlib_disoriented_face_ids(mesh: mrmeshpy.Mesh) -> np.ndarray:
    return np.flatnonzero(selected_disoriented_faces(mesh))


def _protect_decimation_face_neighborhood(
    reference_mesh: mrmeshpy.Mesh,
    protected_faces: mrmeshpy.FaceBitSet,
    candidate: mrmeshpy.Mesh,
    unsafe_face_ids: np.ndarray,
) -> None:
    candidate_verts, candidate_faces = to_arrays(candidate)
    reference = mrmeshpy.MeshPart(reference_mesh)
    for face_id in unsafe_face_ids:
        triangle = candidate_verts[candidate_faces[int(face_id)]]
        for point in (triangle.mean(axis=0), *triangle):
            projection = mrmeshpy.findProjection(
                mrmeshpy.Vector3f(*map(float, point)),
                reference,
            )
            if projection.valid():
                protected_faces.set(projection.proj.face)
    mrmeshpy.expand(
        reference_mesh.topology, protected_faces, _DECIMATION_PROTECTION_RINGS
    )


def decimate_safely(
    mesh: mrmeshpy.Mesh,
    simplify_error_mm: float,
) -> tuple[mrmeshpy.Mesh, DecimationSafeguard]:
    current = int(mesh.topology.numValidFaces())
    if simplify_error_mm <= 0:
        return mesh, DecimationSafeguard(
            simplify_error_mm=float(simplify_error_mm),
            input_faces=current,
            actual_faces=current,
            attempted=False,
            accepted=True,
            error_introduced_mm=0.0,
            initial_self_intersecting_faces=0,
            remaining_self_intersecting_faces=0,
            initial_disoriented_faces=0,
            remaining_disoriented_faces=0,
            repair_attempts=0,
            protected_input_faces=0,
            protected_input_face_fraction=0.0,
            rejection_reason=None,
        )

    before_defects = count_defects(mesh)
    reference_mesh = mrmeshpy.Mesh(mesh)
    protected_faces = mrmeshpy.FaceBitSet(current)
    initial_intersections = 0
    remaining_intersections = 0
    initial_disoriented = 0
    remaining_disoriented = 0
    repair_attempts = 0
    error_introduced_mm = 0.0
    reasons: list[str] = []
    candidate = mesh

    for attempt in range(_MAX_DECIMATION_REPAIR_ATTEMPTS + 1):
        candidate, topology_before, topology_after, error_introduced_mm = (
            _decimate_candidate(mesh, simplify_error_mm, protected_faces)
        )
        after_defects = count_defects(candidate)
        reasons = []
        if topology_after != topology_before:
            reasons.append(
                "topology changed from components/holes/Euler %s to %s"
                % (topology_before, topology_after)
            )
        if any(after > before for before, after in zip(before_defects, after_defects)):
            reasons.append(
                "boundary/hole defects changed from %s to %s"
                % (before_defects, after_defects)
            )
        if reasons:
            break

        colliding_ids = _meshlib_self_intersecting_face_ids(candidate)
        disoriented_ids = _meshlib_disoriented_face_ids(candidate)
        remaining_intersections = int(len(colliding_ids))
        remaining_disoriented = int(len(disoriented_ids))
        if attempt == 0:
            initial_intersections = remaining_intersections
            initial_disoriented = remaining_disoriented
        if not remaining_intersections and not remaining_disoriented:
            break
        if attempt == _MAX_DECIMATION_REPAIR_ATTEMPTS:
            if remaining_intersections:
                reasons.append(
                    "%d self-intersecting face(s) remained after local protection"
                    % remaining_intersections
                )
            if remaining_disoriented:
                reasons.append(
                    "%d disoriented face(s) remained after local protection"
                    % remaining_disoriented
                )
            break

        unsafe_ids = np.union1d(colliding_ids, disoriented_ids)
        protected_before = protected_faces.count()
        _protect_decimation_face_neighborhood(
            reference_mesh, protected_faces, candidate, unsafe_ids
        )
        if protected_faces.count() == protected_before:
            if remaining_intersections:
                reasons.append(
                    "%d self-intersecting face(s) remained because their source "
                    "neighborhood could not be expanded" % remaining_intersections
                )
            if remaining_disoriented:
                reasons.append(
                    "%d disoriented face(s) remained because their source "
                    "neighborhood could not be expanded" % remaining_disoriented
                )
            break
        repair_attempts += 1

    accepted = (
        not reasons and remaining_intersections == 0 and remaining_disoriented == 0
    )
    result = candidate if accepted else mesh
    protected_count = int(protected_faces.count())
    return result, DecimationSafeguard(
        simplify_error_mm=float(simplify_error_mm),
        input_faces=current,
        actual_faces=int(result.topology.numValidFaces()),
        attempted=True,
        accepted=accepted,
        error_introduced_mm=error_introduced_mm if accepted else 0.0,
        initial_self_intersecting_faces=initial_intersections,
        remaining_self_intersecting_faces=remaining_intersections,
        initial_disoriented_faces=initial_disoriented,
        remaining_disoriented_faces=remaining_disoriented,
        repair_attempts=repair_attempts,
        protected_input_faces=protected_count,
        protected_input_face_fraction=float(protected_count / current),
        rejection_reason="; ".join(reasons) if reasons else None,
    )


def decimate(mesh: mrmeshpy.Mesh, simplify_error_mm: float) -> mrmeshpy.Mesh:
    return decimate_safely(mesh, simplify_error_mm)[0]


def vertex_normals(mesh: mrmeshpy.Mesh) -> tuple[np.ndarray, np.ndarray]:
    """Return packed vertices and area-weighted unit vertex normals."""
    vertices, faces = to_arrays(mesh)
    return vertices, vertex_normals_from_arrays(vertices, faces)


def vertex_normals_from_arrays(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Return area-weighted unit vertex normals for indexed triangles."""
    triangles = vertices[faces]
    face_normals = np.cross(
        triangles[:, 1] - triangles[:, 0],
        triangles[:, 2] - triangles[:, 0],
    )
    normals = np.zeros_like(vertices)
    for corner in range(3):
        np.add.at(normals, faces[:, corner], face_normals)
    lengths = np.linalg.norm(normals, axis=1)
    normals[lengths > 0] /= lengths[lengths > 0, None]
    return normals


def validate_output_path(path: str) -> None:
    """Reject unsupported mesh destinations before expensive processing."""
    ext = extension(path)
    if ext not in SUPPORTED_MESH_EXTENSIONS:
        raise ValueError(
            "unsupported output extension %r; supported: %s"
            % (ext, ", ".join(SUPPORTED_MESH_EXTENSIONS))
        )


def write(mesh: mrmeshpy.Mesh, path: str) -> None:
    validate_output_path(path)
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    mrmeshpy.saveMesh(mesh, path)


def write_validated(mesh: mrmeshpy.Mesh, path: str) -> dict:
    from . import validate

    validate_output_path(path)
    vertices, faces = to_arrays(mesh)
    in_memory = validate.validate_arrays(vertices, faces)
    if not in_memory["valid"]:
        raise ValueError(
            "in-memory output mesh is invalid: %s" % "; ".join(in_memory["problems"])
        )

    destination = os.path.abspath(path)
    parent = os.path.dirname(destination)
    os.makedirs(parent, exist_ok=True)
    ext = os.path.splitext(destination)[1]
    fd, temporary = tempfile.mkstemp(
        prefix=".%s." % os.path.basename(destination), suffix=ext, dir=parent
    )
    os.close(fd)
    try:
        write(mesh, temporary)
        report = validate.validate(temporary)
        if not report["valid"]:
            raise ValueError(
                "serialized output mesh is invalid: %s" % "; ".join(report["problems"])
            )
        os.replace(temporary, destination)
        report["file"] = os.path.basename(destination)
        return report
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def bounds_mm(mesh: mrmeshpy.Mesh) -> tuple[float, ...]:
    bounds = mesh.computeBoundingBox()
    return (
        float(bounds.min.x),
        float(bounds.max.x),
        float(bounds.min.y),
        float(bounds.max.y),
        float(bounds.min.z),
        float(bounds.max.z),
    )
