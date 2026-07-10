"""Mesh quality metrics.

STL stores triangle facets without an explicit shared-vertex graph. "Watertight",
"manifold", and "hole" therefore depend on welding coincident positions first.
Trimesh processing performs an approximate, tolerance-derived weld before the
metrics below are calculated.
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import trimesh
import vtk
import meshlib.mrmeshnumpy as mrmeshnumpy
import meshlib.mrmeshpy as mrmeshpy
from vtk.util import numpy_support


def _load_mesh(path: str) -> trimesh.Trimesh:
    """Load a triangle mesh without limiting validation to Trimesh formats."""
    if os.path.splitext(path)[1].lower() != ".vtp":
        return trimesh.load_mesh(path, process=True)

    reader = vtk.vtkXMLPolyDataReader()
    if not reader.CanReadFile(path):
        raise ValueError("invalid VTP file: %s" % path)
    reader.SetFileName(path)

    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(reader.GetOutputPort())
    triangles.PassLinesOff()
    triangles.PassVertsOff()
    triangles.Update()
    poly = triangles.GetOutput()
    if poly.GetPoints() is None or poly.GetNumberOfPolys() == 0:
        raise ValueError("VTP file contains no surface triangles: %s" % path)

    vertices = numpy_support.vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float64)
    connectivity = numpy_support.vtk_to_numpy(poly.GetPolys().GetConnectivityArray())
    faces = connectivity.reshape(-1, 3).astype(np.int64)
    return trimesh.Trimesh(vertices=vertices, faces=faces, process=True)


def _self_intersections(_path: str | None, mesh: trimesh.Trimesh) -> int | str:
    """Count unique faces in MeshLib self-collision pairs."""
    try:
        vertices = np.ascontiguousarray(mesh.vertices, dtype=np.float64)
        faces = np.ascontiguousarray(mesh.faces, dtype=np.int32)
        mr_mesh = mrmeshnumpy.meshFromFacesVerts(faces, vertices)
        pairs = mrmeshpy.findSelfCollidingTriangles(mrmeshpy.MeshPart(mr_mesh))
        selected = set()
        for pair in pairs:
            selected.add(int(pair.aFace))
            selected.add(int(pair.bFace))
        return len(selected)
    except Exception as exc:  # noqa: BLE001
        return "error: %s" % type(exc).__name__


def _component_count(mesh: trimesh.Trimesh) -> int:
    """Number of connected surface shells.

    Counted straight from face adjacency. ``Trimesh.split()`` would be the obvious
    call, but it builds a submesh per component with ``repair=True``, i.e. it
    *fills holes while measuring them*. A validator must not mutate what it
    reports on.

    Faces are adjacent only across edges shared by exactly two faces, so a
    non-manifold edge severs adjacency and inflates this count. Treat it as an
    upper bound on a non-manifold mesh.
    """
    n_faces = len(mesh.faces)
    if n_faces == 0:
        return 0
    components = trimesh.graph.connected_components(
        mesh.face_adjacency, nodes=np.arange(n_faces)
    )
    return int(len(components))


def _problems(report: dict[str, Any]) -> list[str]:
    """Reasons a measured mesh does not satisfy the CLI validity contract."""
    problems = []
    if not report["watertight"]:
        problems.append("mesh is not watertight")
    if not report["winding_consistent"]:
        problems.append("face winding is inconsistent")
    if not report["is_volume"]:
        problems.append("mesh is not a valid enclosed volume")
    if report["boundary_edges"]:
        problems.append("%d boundary edge(s)" % report["boundary_edges"])
    if report["nonmanifold_edge_uses"]:
        problems.append("%d non-manifold edge use(s)" % report["nonmanifold_edge_uses"])
    if report["degenerate_faces"]:
        problems.append("%d degenerate face(s)" % report["degenerate_faces"])

    intersections = report["self_intersecting_faces"]
    if isinstance(intersections, int):
        if intersections:
            problems.append("%d self-intersecting face(s)" % intersections)
    else:
        problems.append("self-intersection check did not complete (%s)" % intersections)
    return problems


def _validate_mesh(mesh: trimesh.Trimesh, *, file: str, bytes: int) -> dict[str, Any]:
    """Full quality report and a single validity result for one triangle mesh."""

    edges = np.sort(mesh.edges_sorted, axis=1)
    _uniq, inverse, counts = np.unique(edges, axis=0, return_inverse=True, return_counts=True)
    per_edge = counts[inverse]

    areas = mesh.area_faces
    watertight = bool(mesh.is_watertight)

    report: dict[str, Any] = {
        "file": file,
        "bytes": bytes,
        "triangles": int(len(mesh.faces)),
        "vertices": int(len(mesh.vertices)),
        "components": _component_count(mesh),
        "watertight": watertight,
        "winding_consistent": bool(mesh.is_winding_consistent),
        "is_volume": bool(mesh.is_volume),
        "euler_number": int(mesh.euler_number),
        "boundary_edges": int(np.count_nonzero(per_edge == 1)),
        "nonmanifold_edge_uses": int(np.count_nonzero(per_edge > 2)),
        "degenerate_faces": int(np.count_nonzero(areas <= 1e-12)),
        "area_mm2": float(areas.sum()),
        "volume_mm3": float(mesh.volume) if watertight else None,
        "bbox_min": [float(v) for v in mesh.bounds[0]],
        "bbox_max": [float(v) for v in mesh.bounds[1]],
        "bbox_extents_mm": [float(v) for v in mesh.extents],
    }

    if watertight and report["components"] == 1:
        # genus = (2 - euler) / 2 for a closed orientable surface
        report["genus"] = int((2 - report["euler_number"]) // 2)

    report["self_intersecting_faces"] = _self_intersections(None, mesh)
    report["problems"] = _problems(report)
    report["valid"] = not report["problems"]

    return report


def validate_arrays(vertices: np.ndarray, faces: np.ndarray) -> dict[str, Any]:
    """Apply the complete validation contract before serialization."""
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
    return _validate_mesh(mesh, file="<in-memory>", bytes=0)


def validate(path: str) -> dict[str, Any]:
    """Full quality report and a single validity result for a mesh file."""
    return _validate_mesh(
        _load_mesh(path),
        file=os.path.basename(path),
        bytes=os.path.getsize(path),
    )


def summarise(report: dict[str, Any]) -> str:
    ok = "yes" if report["watertight"] else "NO"
    lines = [
        "  valid               %s" % ("yes" if report["valid"] else "NO"),
        "  triangles           %s" % f"{report['triangles']:,}",
        "  vertices            %s" % f"{report['vertices']:,}",
        "  components          %s" % f"{report['components']:,}",
        "  watertight          %s" % ok,
        "  winding consistent  %s" % ("yes" if report["winding_consistent"] else "NO"),
        "  boundary edges      %s" % f"{report['boundary_edges']:,}",
        "  non-manifold edges  %s" % f"{report['nonmanifold_edge_uses']:,}",
        "  degenerate faces    %s" % f"{report['degenerate_faces']:,}",
    ]
    if "genus" in report:
        lines.append("  genus               %d" % report["genus"])
    if report.get("volume_mm3") is not None:
        lines.append("  volume              %.0f mm3" % report["volume_mm3"])
    else:
        lines.append("  volume              undefined (mesh is not closed)")
    si = report.get("self_intersecting_faces")
    if si is not None:
        lines.append("  self-intersections  %s" % (f"{si:,}" if isinstance(si, int) else si))
    e = report["bbox_extents_mm"]
    lines.append("  bounding box        %.1f x %.1f x %.1f mm" % (e[0], e[1], e[2]))
    lines.append("  size                %.1f MB" % (report["bytes"] / 1048576.0))
    if report["problems"]:
        lines.append("  problems")
        lines.extend("    - %s" % problem for problem in report["problems"])
    return "\n".join(lines)
