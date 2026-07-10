# Surface finishing created self-intersections

Fixed: 2026-07-11 01:15:44 CEST (+0200)

Baseline commit: `248e35f084d5ee5e96a08f1f5ccf3f6a081de79e`

## Symptom

The default bone conversion of the 2024 skull CT could be watertight and
consistently wound yet contain self-intersecting faces after surface finishing.
The former finishing path also required PyMeshLab as an independent final
authority, despite its selection disagreeing with MeshLib on real candidates.

## Confirmed cause

Windowed-sinc smoothing and quadric simplification can preserve connectivity
while causing non-adjacent triangles to collide. The protection architecture in
`248e35f0` correctly repaired these local patches, but final acceptance used
two different predicates. PyMeshLab's selection was a false-positive source for
the MeshLib-validated candidates, and it added an OpenGL runtime dependency.

## Fix

The restored `248e35f0` smoothing and source-neighborhood protection flow is
retained. MeshLib `findSelfCollidingTriangles` is now the sole collision
predicate in smoothing, simplification, file validation, and VTP's in-memory
path; both faces of every colliding pair are counted.

MeshLib simplification now uses `MinimizeError` with an estimated QEM
surface-deviation limit in millimetres rather than an exact face budget. It
preserves topology and defect checks, retries collision protection up to eight
times, and keeps the valid pre-simplification mesh if no candidate succeeds.
The configured limit, introduced-error estimate, resulting face count, retry
count, and protected source faces are recorded in provenance and warnings.

Conversion and merge validate the in-memory mesh, serialize to a temporary file
in the destination directory, validate the file, then atomically replace the
requested destination. An invalid result cannot overwrite an existing output.

Regression coverage verifies collision-pair counting, clean and VTP validation,
error-limited simplification, zero-error no-op behavior, protection retries,
and atomic failure handling.
