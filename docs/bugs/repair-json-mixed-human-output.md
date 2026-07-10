# Repair JSON output contained human-readable text

Fixed: 2026-07-10 08:42:32 UTC (+0000)

Baseline commit: `e926407dfe22af1449524f41966b530eb8298613`

## Symptom

`dicom-surface repair broken.stl -o fixed.stl --json` exited successfully but
printed progress, a human-readable quality summary, and then JSON. The complete
stdout stream could not be parsed as JSON.

## Cause

The repair handler printed its normal output unconditionally and treated JSON
as an additional output block instead of a mutually exclusive format.

## Fix

JSON mode suppresses progress and human output and writes one JSON object to
stdout. The object includes the output path, repair statistics, and complete
quality report.
