# Print profiles terraced surfaces without guaranteeing feature size

Fixed: 2026-07-11 10:22:00 CEST (+0200)

Baseline commit: `da194407f018d6f50fc1d484fed3a1c8298b5e38`

## Symptom

Applying the resin print profile to the 2024 skull CT roughly doubled the STL
triangle count and made the surface visibly more terraced than the standard
conversion. The profile's stated minimum-feature target was only a mask-space
heuristic and did not guarantee the corresponding final-mesh wall thickness.

## Confirmed cause

The 0.6 mm resin target resampled the complete 0.315 x 0.315 x 0.8 mm source
mask onto a 0.3 mm isotropic grid before surface extraction. Raw triangles grew
from about 4.0 million to 8.0 million, while the same topology-based smoothing
had less physical reach on the denser mesh and retained the source slice
terracing. The thin-region selector also had a known core-reach blind spot.

Mesh-space alternatives did not meet the combined correctness contract. A
local MeshLib offset produced 7.3 million faces and 35 components before its
repair path encountered nested self-intersections. A uniform signed-distance
offset changed skull volume by about 19%; normal simplification lost the
containment property, while strict simplification produced a 92.5 MB mesh that
failed winding validation. Printer minimums also depend on material, process,
orientation, support, and slicer behavior rather than one model-independent
number.

## Fix

Print profiles and their minimum-feature option were removed instead of
retaining a transformation that could neither preserve surface quality nor
deliver its advertised property. The associated profile registry, CLI options,
mask resampling, selective thickening, warnings, provenance, merge branch,
tests, and documentation were deleted. Tissue presets and explicit median,
opening, closing, island, smoothing, and simplification controls remain.
