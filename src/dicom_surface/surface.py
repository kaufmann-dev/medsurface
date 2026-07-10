"""Surface extraction: label volume -> triangle mesh, in DICOM patient space."""

from __future__ import annotations

import os

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


def decimate(poly: vtk.vtkPolyData, target_faces: int) -> vtk.vtkPolyData:
    """Topology-constrained quadric edge-collapse toward a face target.

    PyMeshLab exposes the boundary, normal, and topology constraints used by the
    pipeline. Those constraints can prevent a mesh from reaching the requested
    face count.
    """
    current = poly.GetNumberOfPolys()
    if target_faces <= 0 or current <= target_faces:
        return poly

    tri = vtk.vtkTriangleFilter()
    tri.SetInputData(poly)
    tri.Update()

    verts, faces = to_arrays(tri.GetOutput())
    ms = pymeshlab.MeshSet()
    ms.add_mesh(pymeshlab.Mesh(vertex_matrix=verts, face_matrix=faces))
    ms.apply_filter(
        "meshing_decimation_quadric_edge_collapse",
        targetfacenum=int(target_faces),
        preserveboundary=True,
        preservenormal=True,
        preservetopology=True,
        planarquadric=True,
        qualitythr=0.3,
    )
    m = ms.current_mesh()
    return from_arrays(m.vertex_matrix(), m.face_matrix())


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
