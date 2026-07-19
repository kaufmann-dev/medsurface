# NIfTI affine rounding failed strict verification

Fixed: 2026-07-19 12:44:45 CEST (+0200)

Baseline commit: `e713ccecdde5686052188f47580a2107276a5a16`

## Symptom

Strict conversion of a real CT volume to `.nii.gz` wrote and read back identical
voxels, type, dimensions, spacing, and direction, but rejected the temporary
file because its origin changed during serialization. The same volume converted
successfully to NRRD.

## Confirmed root cause

The verifier applied one fixed `1e-5` absolute tolerance to every geometry
field. NIfTI-1 stores its affine values as float32. The source origin coordinate
`708.1 mm` consequently read back as `708.0999755859375 mm`, a difference of
`0.000024414 mm` (about 24 nanometres) and less than one float32 representable
step at that magnitude. This was harmless format precision, not a physical
geometry change.

## Fix

NRRD and MetaImage retain the existing `1e-5` absolute geometry tolerance.
NIfTI verification uses the greater of `1e-5` and one float32 ULP for each
spacing, origin, and direction component. Values outside float32 range receive
no relaxed tolerance. Pixel digest, type, components, dimensions, and
compression remain exact requirements.

Regression coverage verifies real NIfTI writer rounding at a 708.1 mm origin
and rejects an origin shift larger than one float32 ULP. The formerly rejected
real CT now publishes `.nii.gz` successfully, while the existing sheared-affine
test remains rejected.
