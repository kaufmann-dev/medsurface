# NumPy 2.5 emitted warnings from image-toolkit adapters

Fixed: 2026-07-10 07:56:43 UTC (+0000)

Baseline commit: `c2c7810394e0c8710853c346e285d77b2eaf4ca8`

## Symptom

The test suite passed but emitted hundreds of repeated `DeprecationWarning`
messages after resolving NumPy 2.5.

## Cause

NumPy 2.5 deprecated assigning directly to `ndarray.shape`. SimpleITK 2.5.5
does this in `GetArrayViewFromImage`, and VTK 9.6.2 does it in
`vtk_to_numpy`. Those were the latest available releases, so upgrading could
not resolve the compatibility gap.

## Fix

The project constrains NumPy to `>=1.24,<2.5` until compatible SimpleITK and VTK
releases are available. NumPy 2.4.6 was verified with both affected adapter
calls while treating deprecation warnings as errors, followed by the complete
project test suite.
