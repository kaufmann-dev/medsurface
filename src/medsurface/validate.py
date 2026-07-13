"""Mesh quality metrics measured with MeshLib."""

from __future__ import annotations

import os
from typing import Any

import meshlib.mrmeshnumpy as mrmeshnumpy
import meshlib.mrmeshpy as mrmeshpy
import numpy as np

from .defaults import SUPPORTED_MESH_EXTENSIONS


def _load_mesh(path: str) -> mrmeshpy.Mesh:
    """Load a supported triangle mesh with MeshLib."""
    ext = os.path.splitext(path)[1].lower()
    if ext not in SUPPORTED_MESH_EXTENSIONS:
        raise ValueError(
            "unsupported mesh extension %r; supported: %s"
            % (ext, ", ".join(SUPPORTED_MESH_EXTENSIONS))
        )
    return mrmeshpy.loadMesh(path)


def _self_intersections(mesh: mrmeshpy.Mesh) -> int | str:
    """Count unique faces in MeshLib self-collision pairs."""
    try:
        pairs = mrmeshpy.findSelfCollidingTriangles(mrmeshpy.MeshPart(mesh))
        selected = set()
        for pair in pairs:
            selected.add(int(pair.aFace))
            selected.add(int(pair.bFace))
        return len(selected)
    except Exception as exc:  # noqa: BLE001
        return "error: %s" % type(exc).__name__


def _problems(report: dict[str, Any]) -> list[str]:
    """Reasons a MeshLib-loaded mesh does not satisfy the validity contract."""
    problems = []
    if not report["watertight"]:
        problems.append("mesh is not watertight")
    if report["holes"]:
        problems.append("%d hole(s)" % report["holes"])
    if report["boundary_edges"]:
        problems.append("%d boundary edge(s)" % report["boundary_edges"])
    if not report["winding_consistent"]:
        problems.append("%d disoriented face(s)" % report["disoriented_faces"])
    if not report["is_volume"]:
        problems.append("mesh is not a valid enclosed volume")

    intersections = report["self_intersecting_faces"]
    if isinstance(intersections, int):
        if intersections:
            problems.append("%d self-intersecting face(s)" % intersections)
    else:
        problems.append("self-intersection check did not complete (%s)" % intersections)
    return problems


def _vector(vector: mrmeshpy.Vector3f) -> list[float]:
    return [float(vector.x), float(vector.y), float(vector.z)]


def _validate_mesh(mesh: mrmeshpy.Mesh, *, file: str, bytes: int) -> dict[str, Any]:
    """Return MeshLib-native metrics for one imported triangle mesh."""
    topology = mesh.topology
    triangles = int(topology.numValidFaces())
    vertices = int(topology.numValidVerts())
    components = int(len(mrmeshpy.getAllComponents(mrmeshpy.MeshPart(mesh))))
    holes = int(topology.findNumHoles())
    boundary_edges = int(topology.findLeftBdEdges().count())
    disoriented_faces = int(mrmeshpy.findDisorientedFaces(mesh).count())
    watertight = bool(topology.isClosed())
    winding_consistent = disoriented_faces == 0
    volume = float(mesh.volume())
    is_volume = bool(watertight and winding_consistent and abs(volume) > 1e-12)
    euler = int(
        vertices - topology.computeNotLoneUndirectedEdges() + triangles
    )
    bounds = mesh.computeBoundingBox()
    bbox_min = _vector(bounds.min)
    bbox_max = _vector(bounds.max)

    report: dict[str, Any] = {
        "file": file,
        "bytes": bytes,
        "triangles": triangles,
        "vertices": vertices,
        "components": components,
        "watertight": watertight,
        "winding_consistent": winding_consistent,
        "disoriented_faces": disoriented_faces,
        "is_volume": is_volume,
        "euler_number": euler,
        "holes": holes,
        "boundary_edges": boundary_edges,
        "area_mm2": float(mesh.area()),
        "volume_mm3": abs(volume) if watertight else None,
        "bbox_min": bbox_min,
        "bbox_max": bbox_max,
        "bbox_extents_mm": [hi - lo for lo, hi in zip(bbox_min, bbox_max)],
    }

    if watertight and components == 1:
        report["genus"] = int((2 - euler) // 2)

    report["self_intersecting_faces"] = _self_intersections(mesh)
    report["problems"] = _problems(report)
    report["valid"] = not report["problems"]
    return report


def validate_arrays(vertices: np.ndarray, faces: np.ndarray) -> dict[str, Any]:
    """Apply the complete MeshLib validation contract before serialization."""
    mesh = mrmeshnumpy.meshFromFacesVerts(
        np.ascontiguousarray(faces, dtype=np.int32),
        np.ascontiguousarray(vertices, dtype=np.float64),
    )
    return _validate_mesh(mesh, file="<in-memory>", bytes=0)


def validate(path: str) -> dict[str, Any]:
    """Full quality report and a single validity result for a mesh file."""
    return _validate_mesh(
        _load_mesh(path),
        file=os.path.basename(path),
        bytes=os.path.getsize(path),
    )


def summarise(report: dict[str, Any]) -> str:
    lines = [
        "  valid               %s" % ("yes" if report["valid"] else "NO"),
        "  triangles           %s" % f"{report['triangles']:,}",
        "  vertices            %s" % f"{report['vertices']:,}",
        "  components          %s" % f"{report['components']:,}",
        "  watertight          %s" % ("yes" if report["watertight"] else "NO"),
        "  winding consistent  %s" % ("yes" if report["winding_consistent"] else "NO"),
        "  disoriented faces   %s" % f"{report['disoriented_faces']:,}",
        "  holes               %s" % f"{report['holes']:,}",
        "  boundary edges      %s" % f"{report['boundary_edges']:,}",
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
