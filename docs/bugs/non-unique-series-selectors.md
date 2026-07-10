# DICOM series numbers were ambiguous selectors

Fixed: 2026-07-10 09:41:27 UTC (+0000)

Baseline commit: `c4d6d212a5992bf6419a6362519ca85f709585bd`

## Symptom

The `#` column and `--series` selector treated DICOM `SeriesNumber` as a unique
identifier. Recursive mixed-study trees can contain several unrelated
SeriesInstanceUIDs with the same number. Orientation splits added dotted values
such as `1021.1`, but a bare number or UID could still resolve implicitly to one
usable stack and obscure the collision.

## Cause

The display identifier was derived from optional, non-unique DICOM metadata.
Selection then tried that derived value before the complete UID and description,
including a special dotted syntax that was not an identifier in the source
data.

## Fix

Discovery now assigns every sorted row a unique 1-based integer ID. Selection
accepts that displayed ID, a complete SeriesInstanceUID, or a description
substring. DICOM `SeriesNumber` remains visible as `DICOM #` metadata but is no
longer a selector, and every orientation split receives its own row ID. UID and
description collisions fail with the matching row IDs instead of choosing one
implicitly. JSON and provenance keep UID, SeriesNumber, and orientation-part
metadata separate; row IDs remain local to a fresh discovery result.
