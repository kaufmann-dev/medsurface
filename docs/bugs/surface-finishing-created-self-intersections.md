# Surface finishing created self-intersections

Fixed: 2026-07-10 17:45:53 CEST (+0200)

Baseline commit: `e63a97caad415f3eeb06d11f4aa3e5f65a48887b`

## Symptom

Converting the 2024 skull CT with the default bone preset wrote a watertight,
consistently wound 600,000-triangle STL that failed validation with 79
self-intersecting faces. Reducing global smoothing avoided some collisions but
made anatomical surface quality depend on relaxing the preset. Merge used a
separate finishing sequence and had the same class of risk.

Stage-by-stage measurement found no intersections after marching cubes, two
after initial smoothing, four after the old decimator, and 79 after final
smoothing. The individual stages could therefore turn a valid input surface
into an invalid output even though boundary and non-manifold edge checks passed.

## Cause

Windowed-sinc smoothing and aggressive quadric decimation preserve connectivity
but do not guarantee an embedding without triangle intersections. The pipeline
applied both operations without checking their in-memory result. Its decimation
guard initially rejected an invalid candidate wholesale, which guaranteed a
valid file but retained millions of source triangles instead of meeting the
requested face budget. Conversion and merge also duplicated their finishing
logic, so a safeguard added to one path could drift from the other.

## Fix

Conversion and merge now share one intersection-safe finishing function.
Smoothing still runs every requested iteration. If it introduces collisions,
only vertices in the affected patches return to their pre-smooth positions; the
protected set expands by topological rings until the result is collision-free.

Decimation now uses MeshLib and validates the component/hole/Euler signature,
face budget, boundary edges, non-manifold edges, and self-intersections. When a
candidate intersects itself, its collision triangles are projected onto the
valid source mesh. Four-ring source neighborhoods are excluded from collapse
and decimation restarts, allowing simplification elsewhere to retain the exact
target. MeshLib performs fast per-pass detection and PyMeshLab independently
checks a candidate before acceptance. The valid source is retained only if the
topology checks fail or eight local-protection passes cannot produce a clean
candidate. Safeguard statistics are included in JSON provenance and human
warnings.

Regression tests cover unchanged clean smoothing, local smoothing protection,
unsafe decimation rejection, and successful source-patch protection. The two
supplied real CT conversions and their merge all completed with exit status 0,
600,000 triangles, zero self-intersections, zero boundary edges, and zero
non-manifold edges. The formerly failing 29-face standalone candidate protected
2,355 of 5,459,072 source faces; the merged candidate protected 3,967 of
6,006,220 source faces.
