"""Topological repair for meshes that are not watertight.

The extraction pipeline normally produces a closed, manifold surface. This
module repairs meshes from elsewhere or an output that still contains a defect.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import meshlib.mrmeshpy as mm

from . import surface
from .paths import same_file
from .stages import Logger, ProgressSink, stage_runner


@dataclass
class RepairResult:
    stats: dict[str, int]
    quality: dict[str, Any]


def repair(
    in_path: str,
    out_path: str,
    log: Logger | None = None,
    progress: ProgressSink | None = None,
) -> RepairResult:
    """Weld, fix multiple edges, collapse degeneracies, fill every hole."""
    if same_file(in_path, out_path):
        raise ValueError("input and output must be different files; repair is not in-place")
    surface.validate_output_path(out_path)

    def say(msg: str) -> None:
        if log:
            log(msg)

    step = stage_runner(say, progress, 34)

    mesh = step("load mesh", lambda: mm.loadMesh(in_path))
    stats = {
        "faces_in": int(mesh.topology.numValidFaces()),
        "holes_in": int(mesh.topology.findNumHoles()),
    }
    say("loaded %d faces, %d holes" % (stats["faces_in"], stats["holes_in"]))

    united = step(
        "unite close vertices", lambda: mm.uniteCloseVertices(mesh, 1e-6, False)
    )
    say("united %d close vertices" % united)

    step("fix multiple edges", lambda: mm.fixMultipleEdges(mesh))

    params = mm.FixMeshDegeneraciesParams()
    # Avoid the default remeshing mode, which can subdivide the entire mesh.
    params.mode = mm.FixMeshDegeneraciesParams.Mode.Decimate
    params.maxDeviation = 1e-5
    params.tinyEdgeLength = 1e-4
    step("fix mesh degeneracies", lambda: mm.fixMeshDegeneracies(mesh, params))
    say("fixed degeneracies -> %d faces" % mesh.topology.numValidFaces())

    def fill_holes() -> tuple[int, int]:
        holes = mesh.topology.findHoleRepresentiveEdges()
        fill = mm.FillHoleParams()
        fill.metric = mm.getUniversalMetric(mesh)
        filled = 0
        for i in range(holes.size()):
            try:
                mm.fillHole(mesh, holes[i], fill)
                filled += 1
            except Exception:  # noqa: BLE001,PERF203
                pass
        return filled, int(holes.size())

    filled, found = step("fill boundary holes", fill_holes)
    say("filled %d/%d holes" % (filled, found))

    stats.update(
        faces_out=int(mesh.topology.numValidFaces()),
        holes_out=int(mesh.topology.findNumHoles()),
        holes_filled=filled,
        vertices_united=int(united),
    )
    quality = step(
        "validate and publish repaired mesh",
        lambda: surface.write_validated(mesh, out_path),
    )
    return RepairResult(stats=stats, quality=quality)
