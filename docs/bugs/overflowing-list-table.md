# Series list overflowed narrow terminals

Fixed: 2026-07-10 09:41:27 UTC (+0000)

Baseline commit: `c4d6d212a5992bf6419a6362519ca85f709585bd`

## Symptom

`dicom-surface list` printed a manually padded, 114-character row for every
series. Narrow terminals wrapped each physical line without preserving column
boundaries, so long descriptions and status notes appeared under unrelated
headings. Modalities longer than the four-character `MOD` field, including
`RTDOSE` and `RTPLAN`, shifted every later value.

## Cause

The CLI formatted rows with fixed `%` widths and truncated only descriptions.
It did not measure terminal width, wrap flexible fields, or reserve the actual
width needed by modality values.

## Fix

The list uses a responsive Rich table with separate ID, DICOM number, modality,
description, slice count, voxel spacing, plane, and status columns. Numeric and
geometry fields remain unwrapped, while descriptions and statuses wrap within
their columns. The recommended row includes explicit `default` text as well as
styling. Rendering tests cover wide and narrow terminals, long modalities and
notes, unusable rows, literal markup-like text, and redirected ANSI-free output.
