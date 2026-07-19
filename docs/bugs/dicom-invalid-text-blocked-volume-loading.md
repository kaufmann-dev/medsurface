# Invalid DICOM text blocked volume loading

Fixed: 2026-07-19 12:32:03 CEST (+0200)

Baseline commit: `e713ccecdde5686052188f47580a2107276a5a16`

## Symptom

Strict conversion of two real CT studies failed after pixel loading with a
Python traceback from `Image.SetMetaData`. Extraction and fusion used the same
loader, so the malformed metadata also prevented those pixel workflows from
starting. Failed conversions did not publish partial volume or JSON outputs.

## Confirmed root cause

One DICOM text value contained a byte that was invalid under the declared
character encoding. SimpleITK exposed it as a Python string containing the lone
surrogate `U+DCFC`. The series loader unconditionally copied every first-slice
metadata value onto the 3-D image, but SimpleITK's `SetMetaData` binding cannot
convert a lone surrogate to its required C++ string and raised `TypeError`.

A synthetic DICOM stack containing the same invalid byte reproduced the exact
surrogate and exception without using private medical data. Pixel values and
physical geometry were valid and unrelated to the failure.

## Fix

Pixel workflows no longer load or promote source metadata. Strict conversion
requests metadata only when `--strip-metadata` is absent. Its DICOM promotion
path validates each value as UTF-8, copies valid values exactly, and omits only
malformed values with a count-only warning; it never guesses replacement text.
Unexpected failures while setting valid metadata still propagate.

File inputs follow the same operation boundary: strict conversion preserves
their source metadata by default, while extraction and fusion erase it from the
working image. Regression coverage uses a synthetic invalid DICOM tag to verify
conversion, stripping, warning, metadata isolation, exact image preservation,
and a complete DICOM-to-mesh workflow.
