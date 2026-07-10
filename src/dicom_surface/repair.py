"""Topological repair for meshes that are not watertight.

The extraction pipeline normally produces a closed, manifold surface. This
module repairs meshes from elsewhere or an output that still contains a defect.
"""

from __future__ import annotations

import os

import meshlib.mrmeshpy as mm


def _same_file(left: str, right: str) -> bool:
    """Compare paths safely even through symlinks or hard links."""
    if os.path.realpath(os.path.abspath(left)) == os.path.realpath(os.path.abspath(right)):
        return True
    try:
        return os.path.samefile(left, right)
    except OSError:
        return False


def repair(in_path: str, out_path: str, log=None) -> dict:
    """Weld, fix multiple edges, collapse degeneracies, fill every hole."""
    if _same_file(in_path, out_path):
        raise ValueError("input and output must be different files; repair is not in-place")

    def say(msg):
        if log:
            log(msg)

    mesh = mm.loadMesh(in_path)
    stats = {
        "faces_in": int(mesh.topology.numValidFaces()),
        "holes_in": int(mesh.topology.findNumHoles()),
    }
    say("loaded %d faces, %d holes" % (stats["faces_in"], stats["holes_in"]))

    united = mm.uniteCloseVertices(mesh, 1e-6, False)
    say("united %d close vertices" % united)

    mm.fixMultipleEdges(mesh)

    params = mm.FixMeshDegeneraciesParams()
    # Avoid the default remeshing mode, which can subdivide the entire mesh.
    params.mode = mm.FixMeshDegeneraciesParams.Mode.Decimate
    params.maxDeviation = 1e-5
    params.tinyEdgeLength = 1e-4
    mm.fixMeshDegeneracies(mesh, params)
    say("fixed degeneracies -> %d faces" % mesh.topology.numValidFaces())

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
    say("filled %d/%d holes" % (filled, holes.size()))

    mm.saveMesh(mesh, out_path)

    stats.update(
        faces_out=int(mesh.topology.numValidFaces()),
        holes_out=int(mesh.topology.findNumHoles()),
        holes_filled=filled,
        vertices_united=int(united),
    )
    return stats
