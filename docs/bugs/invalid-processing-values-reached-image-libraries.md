# Invalid processing values reached image and mesh libraries

Fixed: 2026-07-13 12:50:14 UTC (+0000)

Baseline commit: `260e7198a198f73001eb4695a548f0a1b2d2b729`

## Symptom

Values such as `--grid-mm 0`, negative morphology extents, non-finite spacing,
or excessive smoothing strength were accepted by the CLI. Some silently
disabled work; others produced division warnings, allocation failures, toolkit
exceptions, or Python tracebacks. Extremely small positive spacings could plan
impractically large arrays.

## Confirmed cause

Typed CLI conversion checked syntax but not semantic ranges or finiteness.
Programmatic conversion and merge entry points had no equivalent validation,
resampling allocated before enforcing a size limit, and fused-grid voxel counts
used a fixed-width NumPy product that could overflow.

## Fix

The CLI rejects invalid processing numbers before discovery with concise usage
errors. Presets and programmatic entry points enforce the same constraints.
Convert resampling and merge grid planning use non-overflowing counts, reject
non-finite or non-positive spacing, and refuse grids above 800 million voxels
before allocation. Regressions cover zero, negative, NaN, infinity, excessive
force, huge resampling plans, and tiny merge spacing.

