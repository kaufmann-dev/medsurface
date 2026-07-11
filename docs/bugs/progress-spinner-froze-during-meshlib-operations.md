# Progress spinner froze during MeshLib operations

Fixed: 2026-07-11 09:38:56 CEST (+0200)

Baseline commit: `f935a9211267a567e74a545009ba3d254007a34f`

## Symptom

Interactive conversion progress sometimes appeared frozen during long stages such
as component analysis, simplification, and validation. Both the spinner and its
elapsed-time display stopped updating until the stage completed, even though mesh
processing continued successfully.

## Confirmed cause

Rich refreshed the live progress display from a Python thread in the main CLI
process. Several MeshLib Python bindings hold the Python global interpreter lock
for their complete native operation, preventing Rich's refresh thread from
running. Measurement with a one-million-face mesh showed `decimateMesh` blocking
all Python thread progress for its complete 1.75-second call.

## Fix

Interactive progress rendering now runs in a small helper process with its own
Python interpreter. The processing process sends stage updates and log messages
over a pipe, while the helper independently refreshes the spinner and elapsed
time. Redirected and test consoles retain the existing plain-text or in-process
fallback. A pseudo-terminal regression test verifies that the independent
renderer advances elapsed time and receives stage updates.
