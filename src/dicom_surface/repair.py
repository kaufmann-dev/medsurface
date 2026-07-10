"""Topological repair for meshes that are not watertight.

The extraction pipeline in this package normally produces a closed, manifold
surface directly. This module exists for meshes from elsewhere, or for the rare
case where aggressive morphology leaves a defect.

MeshLib was chosen over the alternatives after a head-to-head comparison on a
7.5M-triangle CT bone surface: it was the only library that returned a watertight,
2-manifold, hole-free result, and it did so in ~9 seconds.
"""

from __future__ import annotations

import meshlib.mrmeshpy as mm


def repair(in_path: str, out_path: str, log=None) -> dict:
    """Weld, fix multiple edges, collapse degeneracies, fill every hole."""
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
    # The default is Mode.Remesh, which SUBDIVIDES: on a 7.5M-face mesh it
    # produced 30M faces and a 1.5 GB file while reporting success.
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
