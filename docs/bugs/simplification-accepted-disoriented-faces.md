# Simplification accepted disoriented faces

Fixed: 2026-07-15 02:22:17 CEST (+0200)

Baseline commit: `0ec26c52e03682b70da182fb90e293224be4a678`

## Symptom

A default bone conversion simplified a valid 3.2-million-triangle surface to
about 318,000 triangles, then failed during final publication because the
simplified mesh contained one disoriented face. Disabling simplification
produced a valid surface, confirming that segmentation, physical smoothing,
marching cubes, and surface relaxation were not the source of the defect.

## Confirmed cause

Safe simplification accepted a candidate after checking its component, hole,
Euler, boundary, and self-intersection properties. It did not run the same
MeshLib disoriented-face check used by final validation. An edge collapse could
therefore introduce inconsistent winding without changing any property in the
simplification acceptance gate, and the defect was discovered only after all
geometry processing had finished.

## Fix

Simplification now finds both self-intersecting and disoriented candidate faces
before accepting a result. Their candidate triangles are projected onto the
pre-simplification mesh, four-ring source neighborhoods are protected from
collapse, and simplification retries through the existing bounded local repair
loop. If no clean candidate is possible, the valid pre-simplification surface
is retained instead of publishing an invalid mesh.

Safeguard provenance records initial and remaining counts for both defect
types, and warnings identify which defects caused source-face protection.
Regression coverage verifies the MeshLib disoriented-face detector, local
protection retry, and safe rejection when the defect cannot be removed. The
reported conversion now succeeds after one retry with zero disoriented or
self-intersecting faces.
