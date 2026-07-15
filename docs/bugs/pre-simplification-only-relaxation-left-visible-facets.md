# Pre-simplification-only relaxation left visible facets

Fixed: 2026-07-15

Baseline commit: `954db0af47a6029a4ce9506e4ead03f0bb3cd970`

## Symptom

Default intensity conversion again preserved the skull anatomy after mask
smoothing was disabled, but the result looked visibly rougher than the earlier
valid surface. The difference was most apparent as small polygonal facets over
the calvarium and broad facial-bone surfaces.

## Confirmed root cause

The physical-coordinate finishing change had removed the final mesh-relaxation
stage and moved its iterations before simplification. Bone and auto therefore
ran 60 relaxation iterations before quadric reduction instead of 20 before and
40 after. Skin similarly changed from 25 plus 10 to 35 before. These totals are
not visually equivalent: simplification changes the triangulation and can
create larger planar facets after the only smoothing pass has finished.

## Fix

The shared finishing pipeline now has explicit guarded relaxation on both sides
of simplification. Bone and auto default to 20 pre plus 40 post iterations; skin
uses 25 plus 10; teeth stays at 10 plus 0. Labelmaps retain 20 plus 0 because
their 0.8 mm physical mask smoothing remains the appropriate default for voxel
terracing. All four commands expose the final pass as
`--post-mesh-smooth-iters`.

Both relaxation stages now reject self-intersections and inconsistent winding.
Unsafe vertices are restored from the stage input and protection expands by
topological rings; if no clean result is possible, the valid input surface is
retained. JSON provenance records pre-smoothing, simplification, and
post-smoothing separately, including RMS and maximum post-pass displacement.
The CLI reports that displacement because MeshLib's simplification estimate
does not include later smoothing.

Regressions cover preset and labelmap defaults, stage order, all four command
paths, local winding protection, full-pass fallback, provenance, and the
separate displacement measurement.
