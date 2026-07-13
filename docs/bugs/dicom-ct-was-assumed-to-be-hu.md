# DICOM CT modality alone was treated as verified HU

Fixed: 2026-07-13 12:50:14 UTC (+0000)

Baseline commit: `260e7198a198f73001eb4695a548f0a1b2d2b729`

## Symptom

Every DICOM candidate with `Modality=CT` suppressed HU calibration warnings and
recorded verified calibration. Derived, multienergy, inconsistent, or malformed
CT series could therefore receive preset HU thresholds without evidence that
their rescaled output units were Hounsfield units.

## Confirmed cause

HU classification checked only the candidate's modality. Discovery did not
retain `ImageType`, `RescaleType`, the rescale transform, or multienergy state,
and it did not compare calibration metadata across slices.

## Fix

Discovery now reads and compares calibration evidence for every instance. HU is
verified only with a consistent finite rescale slope/intercept plus either an
explicit `RescaleType=HU` or original, non-localizer, non-multienergy CT image
metadata. Derived and multienergy CT require the explicit HU unit. List JSON and
provenance expose the evidence and verification result, while ambiguous series
receive the existing calibration warning.

