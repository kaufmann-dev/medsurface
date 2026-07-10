# Default CLI validation skipped self-intersections

Fixed: 2026-07-10 08:23:30 UTC (+0000)

Baseline commit: `c37b53f7469b6e205e9fb1fa7d70e0b8e4df9f3c`

## Symptom

The CLI called its mesh-quality check after conversion, merge, repair, and from
`validate`, but self-intersections were omitted unless a separate option was
given. A watertight mesh with intersecting closed shells therefore returned a
successful validation result.

## Cause

The Python validation function defaulted to measuring self-intersections, while
each CLI command explicitly overrode that default from an opt-in flag. Command
success then depended only on watertightness, even though winding, volume,
degenerate geometry, and intersections were already measured or available.

## Fix

All CLI validation uses one complete contract. It always measures
self-intersections and reports `valid` plus concrete `problems`. Every validating
command uses that result for its exit status, while `convert --no-validate` and
`merge --no-validate` remain the explicit way to skip the whole check. Multiple
independently closed shells remain valid.

In the development environment, complete validation of a 589,824-face
procedural torus took approximately 9.3 seconds. The cost depends on mesh size
and geometry.
