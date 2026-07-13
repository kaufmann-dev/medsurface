"""Surface extraction and finishing with MeshLib."""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass

import meshlib.mrmeshnumpy as mrmeshnumpy
import meshlib.mrmeshpy as mrmeshpy
import numpy as np
import SimpleITK as sitk

_SUPPORTED_EXTENSIONS = (".obj", ".ply", ".stl")


@dataclass(frozen=True)
class SmoothingSafeguard:
    requested_iterations: int
    initial_self_intersecting_faces: int
    protected_vertices: int
    protected_vertex_fraction: float
    expanded_rings: int


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
    transformed = vertices @ np.asarray(matrix[:3, :3], dtype=float).T
    transformed += np.asarray(matrix[:3, 3], dtype=float)
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


def protect_smoothed_surface(
    original: mrmeshpy.Mesh,
    smoothed: mrmeshpy.Mesh,
    requested_iterations: int,
) -> tuple[mrmeshpy.Mesh, SmoothingSafeguard]:
    """Keep full smoothing except where it makes non-adjacent faces collide."""
    original_verts, original_faces = to_arrays(original)
    smoothed_verts, smoothed_faces = to_arrays(smoothed)
    if not np.array_equal(original_faces, smoothed_faces):
        raise ValueError("smoothing changed mesh connectivity")

    selected = selected_self_intersecting_faces(smoothed)
    initial = int(np.count_nonzero(selected))
    if initial == 0:
        return smoothed, SmoothingSafeguard(
            requested_iterations=int(requested_iterations),
            initial_self_intersecting_faces=0,
            protected_vertices=0,
            protected_vertex_fraction=0.0,
            expanded_rings=0,
        )

    protected = np.zeros(len(original_verts), dtype=bool)
    protected[np.unique(original_faces[selected])] = True
    rings = 0
    while True:
        candidate_verts = smoothed_verts.copy()
        candidate_verts[protected] = original_verts[protected]
        candidate = from_arrays(candidate_verts, original_faces)
        if not selected_self_intersecting_faces(candidate).any():
            count = int(np.count_nonzero(protected))
            return candidate, SmoothingSafeguard(
                requested_iterations=int(requested_iterations),
                initial_self_intersecting_faces=initial,
                protected_vertices=count,
                protected_vertex_fraction=float(count / len(original_verts)),
                expanded_rings=rings,
            )

        incident_faces = protected[original_faces].any(axis=1)
        expanded = protected.copy()
        expanded[np.unique(original_faces[incident_faces])] = True
        if np.array_equal(expanded, protected):
            raise ValueError("could not produce an intersection-free smoothed surface")
        protected = expanded
        rings += 1


def smooth_safely(
    mesh: mrmeshpy.Mesh,
    iterations: int,
    force: float,
) -> tuple[mrmeshpy.Mesh, SmoothingSafeguard]:
    if iterations <= 0:
        return mesh, SmoothingSafeguard(
            requested_iterations=int(iterations),
            initial_self_intersecting_faces=0,
            protected_vertices=0,
            protected_vertex_fraction=0.0,
            expanded_rings=0,
        )
    return protect_smoothed_surface(mesh, smooth(mesh, iterations, force), iterations)


def largest_component(mesh: mrmeshpy.Mesh) -> tuple[mrmeshpy.Mesh, int]:
    components = mrmeshpy.getAllComponents(mrmeshpy.MeshPart(mesh))
    count = len(components)
    if not count:
        return mesh, 0
    largest = max(components, key=lambda region: region.count())
    return mesh.cloneRegion(largest), count


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


def _protect_decimation_collision_neighborhood(
    reference_mesh: mrmeshpy.Mesh,
    protected_faces: mrmeshpy.FaceBitSet,
    candidate: mrmeshpy.Mesh,
    colliding_face_ids: np.ndarray,
) -> None:
    candidate_verts, candidate_faces = to_arrays(candidate)
    reference = mrmeshpy.MeshPart(reference_mesh)
    for face_id in colliding_face_ids:
        triangle = candidate_verts[candidate_faces[int(face_id)]]
        for point in (triangle.mean(axis=0), *triangle):
            projection = mrmeshpy.findProjection(
                mrmeshpy.Vector3f(*map(float, point)),
                reference,
            )
            if projection.valid():
                protected_faces.set(projection.proj.face)
    mrmeshpy.expand(reference_mesh.topology, protected_faces, _DECIMATION_PROTECTION_RINGS)


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
    repair_attempts = 0
    error_introduced_mm = 0.0
    reasons: list[str] = []
    candidate = mesh

    for attempt in range(_MAX_DECIMATION_REPAIR_ATTEMPTS + 1):
        candidate, topology_before, topology_after, error_introduced_mm = _decimate_candidate(
            mesh, simplify_error_mm, protected_faces
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
        remaining_intersections = int(len(colliding_ids))
        if attempt == 0:
            initial_intersections = remaining_intersections
        if not remaining_intersections:
            break
        if attempt == _MAX_DECIMATION_REPAIR_ATTEMPTS:
            reasons.append(
                "%d self-intersecting face(s) remained after local protection"
                % remaining_intersections
            )
            break

        protected_before = protected_faces.count()
        _protect_decimation_collision_neighborhood(
            reference_mesh, protected_faces, candidate, colliding_ids
        )
        if protected_faces.count() == protected_before:
            reasons.append(
                "%d self-intersecting face(s) remained because their source "
                "neighborhood could not be expanded" % remaining_intersections
            )
            break
        repair_attempts += 1

    accepted = not reasons and remaining_intersections == 0
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
    return vertices, normals


def write(mesh: mrmeshpy.Mesh, path: str) -> None:
    ext = os.path.splitext(path)[1].lower()
    if ext not in _SUPPORTED_EXTENSIONS:
        raise ValueError(
            "unsupported output extension %r; supported: %s"
            % (ext, ", ".join(_SUPPORTED_EXTENSIONS))
        )
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    mrmeshpy.saveMesh(mesh, path)


def write_validated(mesh: mrmeshpy.Mesh, path: str) -> dict:
    from . import validate

    vertices, faces = to_arrays(mesh)
    in_memory = validate.validate_arrays(vertices, faces)
    if not in_memory["valid"]:
        raise ValueError("in-memory output mesh is invalid: %s" % "; ".join(in_memory["problems"]))

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
        float(bounds.min.x), float(bounds.max.x),
        float(bounds.min.y), float(bounds.max.y),
        float(bounds.min.z), float(bounds.max.z),
    )
