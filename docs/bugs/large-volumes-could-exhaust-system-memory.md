# Large volumes could exhaust system memory

Fixed: 2026-07-13 20:50:55 CEST (+0200)

Baseline commit: `32d426835336dba6ad790b5272a54b6f93008cea`

## Symptom

Loading and converting a native 2.7-billion-voxel volume exhausted system
memory instead of failing before pixel allocation. Conversion succeeded only
after the volume was resampled externally to a coarser spacing. Existing guards
covered newly planned resampling and merge grids, not native input loading.

## Confirmed cause

Discovery already retained header dimensions, but the loader did not inspect
their product before asking SimpleITK to read pixels. Separate resampled and
fused-grid guards used an undocumented 800-million-voxel ceiling and had no
shared expert override.

## Fix

One 500-million-voxel default ceiling now covers source volumes, conversion
resampling grids, and merge fused grids using non-overflowing integer products.
Source checks run before pixel reads and grid checks run before allocation.
`--allow-large-volume` on `convert` and `merge` bypasses every voxel-count guard,
is propagated through the programmatic APIs, and is recorded in provenance. It
does not promise sufficient memory or bypass other safety checks. Regressions
cover the boundary, pre-read refusal, allocation refusal, both CLI paths, and
the override.
