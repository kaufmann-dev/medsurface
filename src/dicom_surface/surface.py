"""Surface extraction: label volume -> triangle mesh, in DICOM patient space."""

from __future__ import annotations

import os
from dataclasses import dataclass

import meshlib.mrmeshnumpy as mrmeshnumpy
import meshlib.mrmeshpy as mrmeshpy
import numpy as np
import pymeshlab
import SimpleITK as sitk
import vtk
from vtk.util import numpy_support  # noqa: N813

_EXT_WRITERS = {
    ".stl": vtk.vtkSTLWriter,
    ".ply": vtk.vtkPLYWriter,
    ".obj": vtk.vtkOBJWriter,
    ".vtp": vtk.vtkXMLPolyDataWriter,
}


@dataclass(frozen=True)
class SmoothingSafeguard:
    requested_iterations: int
    initial_self_intersecting_faces: int
    protected_vertices: int
    protected_vertex_fraction: float
    expanded_rings: int


@dataclass(frozen=True)
class DecimationSafeguard:
    requested_faces: int
    input_faces: int
    actual_faces: int
    attempted: bool
    accepted: bool
    initial_self_intersecting_faces: int
    remaining_self_intersecting_faces: int
    repair_attempts: int
    protected_input_faces: int
    protected_input_face_fraction: float
    rejection_reason: str | None


_DECIMATION_PROTECTION_RINGS = 4
_MAX_DECIMATION_REPAIR_ATTEMPTS = 8


def to_vtk_image(image: sitk.Image) -> vtk.vtkImageData:
    """Copy a SimpleITK image into vtkImageData in *index* space.

    Spacing and origin are deliberately left at identity: the image's true
    physical placement (including any oblique orientation) is applied later as a
    single affine, which handles direction cosines that vtkImageData cannot
    represent.

    Integer volumes are carried as uint8 (label volumes), floating-point volumes
    as float32 (distance fields).
    """
    arr = sitk.GetArrayFromImage(image)  # (z, y, x)
    if np.issubdtype(arr.dtype, np.floating):
        arr = np.ascontiguousarray(arr, dtype=np.float32)
        vtk_type = vtk.VTK_FLOAT
    else:
        arr = np.ascontiguousarray(arr, dtype=np.uint8)
        vtk_type = vtk.VTK_UNSIGNED_CHAR

    vtk_arr = numpy_support.numpy_to_vtk(arr.ravel(order="C"), deep=True, array_type=vtk_type)
    img = vtk.vtkImageData()
    sx, sy, sz = image.GetSize()
    img.SetDimensions(sx, sy, sz)
    img.SetSpacing(1.0, 1.0, 1.0)
    img.SetOrigin(0.0, 0.0, 0.0)
    img.GetPointData().SetScalars(vtk_arr)
    return img


def index_to_physical(image: sitk.Image) -> vtk.vtkMatrix4x4:
    """Affine mapping voxel index -> physical (LPS) millimetres.

    ``p = origin + Direction @ diag(spacing) @ index``
    """
    d = np.asarray(image.GetDirection(), dtype=float).reshape(3, 3)
    s = np.asarray(image.GetSpacing(), dtype=float)
    o = np.asarray(image.GetOrigin(), dtype=float)

    m = vtk.vtkMatrix4x4()
    m.Identity()
    linear = d @ np.diag(s)
    for i in range(3):
        for j in range(3):
            m.SetElement(i, j, float(linear[i, j]))
        m.SetElement(i, 3, float(o[i]))
    return m


def marching_cubes(img: vtk.vtkImageData, isovalue: float = 0.5) -> vtk.vtkPolyData:
    fe = vtk.vtkFlyingEdges3D()
    fe.SetInputData(img)
    fe.SetValue(0, isovalue)
    fe.ComputeNormalsOff()
    fe.ComputeGradientsOff()
    fe.Update()
    return fe.GetOutput()


def transform(poly: vtk.vtkPolyData, matrix: vtk.vtkMatrix4x4) -> vtk.vtkPolyData:
    t = vtk.vtkTransform()
    t.SetMatrix(matrix)
    f = vtk.vtkTransformPolyDataFilter()
    f.SetInputData(poly)
    f.SetTransform(t)
    f.Update()
    return f.GetOutput()


def smooth(poly: vtk.vtkPolyData, iterations: int, passband: float) -> vtk.vtkPolyData:
    """Windowed-sinc (Taubin-family) smoothing: no volumetric shrinkage.

    A plain Laplacian smoother would contract a closed surface toward its
    centroid, quietly shrinking the anatomy with every iteration.
    """
    if iterations <= 0:
        return poly
    f = vtk.vtkWindowedSincPolyDataFilter()
    f.SetInputData(poly)
    f.SetNumberOfIterations(int(iterations))
    f.SetPassBand(float(passband))
    f.BoundarySmoothingOn()
    f.NonManifoldSmoothingOn()
    f.NormalizeCoordinatesOn()
    f.FeatureEdgeSmoothingOff()
    f.Update()
    return f.GetOutput()


def selected_self_intersecting_faces(poly: vtk.vtkPolyData) -> np.ndarray:
    """Boolean face selection using the same MeshLab predicate as validation."""
    verts, faces = to_arrays(poly)
    ms = pymeshlab.MeshSet()
    ms.add_mesh(pymeshlab.Mesh(vertex_matrix=verts, face_matrix=faces))
    ms.apply_filter("compute_selection_by_self_intersections_per_face")
    return np.asarray(ms.current_mesh().face_selection_array(), dtype=bool)


def protect_smoothed_surface(
    original: vtk.vtkPolyData,
    smoothed: vtk.vtkPolyData,
    requested_iterations: int,
) -> tuple[vtk.vtkPolyData, SmoothingSafeguard]:
    """Keep full smoothing except where it makes non-adjacent faces collide.

    Intersecting vertices return to their known-valid pre-smooth positions. If
    that is not enough, the protected set grows one topological ring at a time.
    The process is deterministic and must terminate at the original geometry.
    """
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
    poly: vtk.vtkPolyData,
    iterations: int,
    passband: float,
) -> tuple[vtk.vtkPolyData, SmoothingSafeguard]:
    """Run every requested smoothing iteration, then protect collision patches."""
    if iterations <= 0:
        return poly, SmoothingSafeguard(
            requested_iterations=int(iterations),
            initial_self_intersecting_faces=0,
            protected_vertices=0,
            protected_vertex_fraction=0.0,
            expanded_rings=0,
        )
    return protect_smoothed_surface(poly, smooth(poly, iterations, passband), iterations)


def largest_component(poly: vtk.vtkPolyData) -> tuple[vtk.vtkPolyData, int]:
    """Keep only the biggest connected surface.

    This is what removes enclosed internal cavities (sinuses, trabecular air
    cells, marrow space). Each cavity is a closed shell disconnected from the
    outer surface, so it survives every labelmap-level cleanup and only
    disappears here.
    """
    conn = vtk.vtkPolyDataConnectivityFilter()
    conn.SetInputData(poly)
    conn.SetExtractionModeToAllRegions()
    conn.Update()
    n = conn.GetNumberOfExtractedRegions()

    conn.SetExtractionModeToLargestRegion()
    conn.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputConnection(conn.GetOutputPort())
    clean.Update()
    return clean.GetOutput(), n


def count_defects(poly: vtk.vtkPolyData) -> tuple[int, int]:
    """(boundary edges, non-manifold edges). Cheap: no file IO, no welding."""
    counts = []
    for boundary, nonmanifold in ((True, False), (False, True)):
        fe = vtk.vtkFeatureEdges()
        fe.SetInputData(poly)
        fe.SetBoundaryEdges(boundary)
        fe.SetNonManifoldEdges(nonmanifold)
        fe.FeatureEdgesOff()
        fe.ManifoldEdgesOff()
        fe.Update()
        counts.append(int(fe.GetOutput().GetNumberOfCells()))
    return counts[0], counts[1]


def to_arrays(poly: vtk.vtkPolyData) -> tuple[np.ndarray, np.ndarray]:
    """(vertices, triangles) as numpy arrays. Input must be triangulated."""
    verts = numpy_support.vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float64)
    conn = numpy_support.vtk_to_numpy(poly.GetPolys().GetConnectivityArray())
    faces = conn.reshape(-1, 3).astype(np.int32)
    return verts, faces


def from_arrays(verts: np.ndarray, faces: np.ndarray) -> vtk.vtkPolyData:
    points = vtk.vtkPoints()
    points.SetData(numpy_support.numpy_to_vtk(np.ascontiguousarray(verts, dtype=np.float64),
                                              deep=True))
    n = len(faces)
    offsets = np.arange(0, 3 * (n + 1), 3, dtype=np.int64)
    connectivity = np.ascontiguousarray(faces, dtype=np.int64).ravel()

    cells = vtk.vtkCellArray()
    cells.SetData(
        numpy_support.numpy_to_vtkIdTypeArray(offsets, deep=True),
        numpy_support.numpy_to_vtkIdTypeArray(connectivity, deep=True),
    )

    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    poly.SetPolys(cells)
    return poly


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
    poly: vtk.vtkPolyData,
    target_faces: int,
    protected_faces: mrmeshpy.FaceBitSet | None = None,
) -> tuple[vtk.vtkPolyData, tuple[int, int, int], tuple[int, int, int]]:
    tri = vtk.vtkTriangleFilter()
    tri.SetInputData(poly)
    tri.Update()

    verts, faces = to_arrays(tri.GetOutput())
    mesh = mrmeshnumpy.meshFromFacesVerts(faces, verts)
    before = _meshlib_topology_signature(mesh)

    settings = mrmeshpy.DecimateSettings()
    settings.maxDeletedFaces = int(len(faces) - target_faces)
    settings.maxError = float("inf")
    # Partitioned decimation trades quality for speed and can make the result
    # depend on CPU count. One part is already much faster than the old backend.
    settings.subdivideParts = 1
    settings.packMesh = True
    region = None
    if protected_faces is not None and protected_faces.any():
        region = mesh.topology.getValidFaces()
        region.subtract(protected_faces, 0)
        settings.region = region
    mrmeshpy.decimateMesh(mesh, settings)

    after = _meshlib_topology_signature(mesh)
    out_verts = mrmeshnumpy.getNumpyVerts(mesh)
    out_faces = mrmeshnumpy.getNumpyFaces(mesh.topology)
    return from_arrays(out_verts, out_faces), before, after


def _meshlib_self_intersecting_face_ids(poly: vtk.vtkPolyData) -> np.ndarray:
    """Return colliding face ids quickly for iterative decimation repair."""
    verts, faces = to_arrays(poly)
    mesh = mrmeshnumpy.meshFromFacesVerts(faces, verts)
    selected = mrmeshpy.localFindSelfIntersections(mesh)
    return np.fromiter(
        (int(selected.nthSetBit(i)) for i in range(selected.count())),
        dtype=np.int64,
        count=selected.count(),
    )


def _protect_decimation_collision_neighborhood(
    reference_mesh: mrmeshpy.Mesh,
    protected_faces: mrmeshpy.FaceBitSet,
    candidate: vtk.vtkPolyData,
    colliding_face_ids: np.ndarray,
) -> None:
    """Map collision patches back to source faces and exclude their collapse."""
    candidate_verts, candidate_faces = to_arrays(candidate)
    reference = mrmeshpy.MeshPart(reference_mesh)

    for face_id in colliding_face_ids:
        triangle = candidate_verts[candidate_faces[int(face_id)]]
        sample_points = (triangle.mean(axis=0), *triangle)
        for point in sample_points:
            projection = mrmeshpy.findProjection(
                mrmeshpy.Vector3f(*map(float, point)),
                reference,
            )
            if projection.valid():
                protected_faces.set(projection.proj.face)

    mrmeshpy.expand(
        reference_mesh.topology,
        protected_faces,
        _DECIMATION_PROTECTION_RINGS,
    )


def decimate_safely(
    poly: vtk.vtkPolyData,
    target_faces: int,
) -> tuple[vtk.vtkPolyData, DecimationSafeguard]:
    """Simplify without accepting changed topology or intersecting faces."""
    current = poly.GetNumberOfPolys()
    if target_faces <= 0 or current <= target_faces:
        return poly, DecimationSafeguard(
            requested_faces=int(target_faces),
            input_faces=int(current),
            actual_faces=int(current),
            attempted=False,
            accepted=True,
            initial_self_intersecting_faces=0,
            remaining_self_intersecting_faces=0,
            repair_attempts=0,
            protected_input_faces=0,
            protected_input_face_fraction=0.0,
            rejection_reason=None,
        )

    before_defects = count_defects(poly)
    source_verts, source_faces = to_arrays(poly)
    reference_mesh = mrmeshnumpy.meshFromFacesVerts(source_faces, source_verts)
    protected_faces = mrmeshpy.FaceBitSet(len(source_faces))
    initial_intersections = 0
    remaining_intersections = 0
    repair_attempts = 0
    reasons: list[str] = []
    candidate = poly

    for attempt in range(_MAX_DECIMATION_REPAIR_ATTEMPTS + 1):
        candidate, topology_before, topology_after = _decimate_candidate(
            poly,
            target_faces,
            protected_faces,
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
                "boundary/non-manifold defects changed from %s to %s"
                % (before_defects, after_defects)
            )
        if candidate.GetNumberOfPolys() > target_faces + 2:
            reasons.append(
                "could only reduce to %d face(s)" % candidate.GetNumberOfPolys()
            )
        if reasons:
            break

        colliding_ids = _meshlib_self_intersecting_face_ids(candidate)
        # MeshLab is the independent final authority. Usually the fast MeshLib
        # detector finds the same faces; only pay for both once MeshLib is clean.
        if not len(colliding_ids):
            mesh_lab_selection = selected_self_intersecting_faces(candidate)
            colliding_ids = np.flatnonzero(mesh_lab_selection)
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
            reference_mesh,
            protected_faces,
            candidate,
            colliding_ids,
        )
        if protected_faces.count() == protected_before:
            reasons.append(
                "%d self-intersecting face(s) remained because their source "
                "neighborhood could not be expanded" % remaining_intersections
            )
            break
        repair_attempts += 1

    accepted = not reasons and remaining_intersections == 0
    result = candidate if accepted else poly
    protected_count = int(protected_faces.count())
    return result, DecimationSafeguard(
        requested_faces=int(target_faces),
        input_faces=int(current),
        actual_faces=int(result.GetNumberOfPolys()),
        attempted=True,
        accepted=accepted,
        initial_self_intersecting_faces=initial_intersections,
        remaining_self_intersecting_faces=remaining_intersections,
        repair_attempts=repair_attempts,
        protected_input_faces=protected_count,
        protected_input_face_fraction=float(protected_count / current),
        rejection_reason="; ".join(reasons) if reasons else None,
    )


def decimate(poly: vtk.vtkPolyData, target_faces: int) -> vtk.vtkPolyData:
    """Return the safely decimated mesh, or the valid input if safeguards fail."""
    return decimate_safely(poly, target_faces)[0]


def compute_normals(poly: vtk.vtkPolyData) -> vtk.vtkPolyData:
    n = vtk.vtkPolyDataNormals()
    n.SetInputData(poly)
    n.ConsistencyOn()
    n.SplittingOff()
    n.AutoOrientNormalsOn()
    n.Update()
    return n.GetOutput()


def write(poly: vtk.vtkPolyData, path: str) -> None:
    ext = os.path.splitext(path)[1].lower()
    try:
        writer_cls = _EXT_WRITERS[ext]
    except KeyError:
        raise ValueError(
            "unsupported output extension %r; supported: %s"
            % (ext, ", ".join(sorted(_EXT_WRITERS)))
        ) from None

    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)

    w = writer_cls()
    w.SetFileName(path)
    w.SetInputData(poly)
    if hasattr(w, "SetFileTypeToBinary"):
        w.SetFileTypeToBinary()
    if not w.Write():
        raise IOError("failed to write %s" % path)


def bounds_mm(poly: vtk.vtkPolyData) -> tuple[float, ...]:
    return tuple(poly.GetBounds())
