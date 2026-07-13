# Merge omitted shared mask and component controls

Fixed: 2026-07-13 12:50:14 UTC (+0000)

Baseline commit: `260e7198a198f73001eb4695a548f0a1b2d2b729`

## Symptom

Conversion exposed opening and explicit keep-all island/component controls, but
merge did not. A user could not reproduce the same preset override while fusing
volumes even though both workflows call the same mask and surface functions.

## Confirmed cause

The merge command's option list and preset override omitted `--opening-mm`,
`--all-islands`, and `--all-components`. The processing layer already supported
all three values.

## Fix

Merge now exposes the three options and applies them through the same preset
override used by conversion. CLI coverage passes every shared override and
asserts that the resolved preset reaches the merge engine unchanged.

