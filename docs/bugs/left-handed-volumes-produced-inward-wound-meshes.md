# Left-handed volumes produced inward-wound meshes

Fixed: 2026-07-13 20:50:55 CEST (+0200)

Baseline commit: `32d426835336dba6ad790b5272a54b6f93008cea`

## Symptom

Converting a NIfTI volume with a negative-determinant direction matrix produced
a closed mesh whose faces were all disoriented. Final validation rejected the
output even though the same voxels converted successfully after reorientation
to a right-handed coordinate system.

## Confirmed cause

Marching cubes created outward-wound faces in voxel-index space. The subsequent
index-to-physical affine reflected the vertices for a left-handed image but kept
the original face order. A reflection reverses orientation, so the complete
surface became inward-wound. A synthetic reflected sphere reproduced the
failure with every face disoriented; reversing each face restored validity.

## Fix

The shared mesh transform now checks the determinant of its linear component
and reverses triangle order for reflections. Positive-determinant transforms are
unchanged. Unit regressions cover both determinant signs, and a real
SimpleITK-written left-handed NIfTI converts to a valid, consistently wound STL.
