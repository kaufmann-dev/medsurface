# Multi-valued DICOM kernels rendered inconsistently

Fixed: 2026-07-13 18:21:02 CEST (+0200)

Baseline commit: `ea2012ea75522ae81ae6f8d38cf0a918f1744b4c`

## Symptom

For a DICOM `ConvolutionKernel` containing the two values `Hr68f` and `1`, the
human `list` status displayed `Hr68f\1` while convert and merge warnings
displayed the Python-escaped form `Hr68f\\1`. The output also made those kernel
values easy to confuse with the separate `SeriesDescription` text `Hr68 A1`.

## Confirmed cause

Discovery flattened DICOM's multi-valued `ConvolutionKernel` into one
backslash-delimited string. The list status interpolated that string normally,
while the shared convert and merge warning rendered it with `repr`, escaping the
already synthetic separator. List JSON and conversion provenance also exposed
the flattened representation instead of the original value structure.

## Fix

Series metadata now preserves convolution kernels as a tuple of DICOM values.
Sharp-kernel classification checks each value independently, and one shared
formatter renders human output as `Hr68f, 1`. The list table labels that value
as `Kernel` under Metadata, while Status contains only the derived `sharp
kernel` classification. List JSON and conversion provenance serialize the
values as `["Hr68f", "1"]`. Regressions cover DICOM discovery, the list status,
list JSON, the shared convert/merge warning, and provenance; the real 231-slice
skull series verifies the same paths end to end.
