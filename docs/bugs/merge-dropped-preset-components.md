# Merge discarded components that a preset promised to preserve

Fixed: 2026-07-13 12:50:14 UTC (+0000)

Baseline commit: `260e7198a198f73001eb4695a548f0a1b2d2b729`

## Symptom

`merge --preset teeth` could discard every disconnected surface shell except
the largest even though the teeth preset promises to retain separate teeth.
When component selection was disabled, conversion and merge also reported one
surface component without measuring the actual count.

## Confirmed cause

The merge finishing call hardcoded `keep_largest_component=True` instead of
using the selected preset. The shared finishing function initialized its
component count to one and only measured components while selecting the largest.

## Fix

Merge passes the preset's component policy unchanged. Shared finishing now
counts shells before optional selection, reports whether it kept one or all, and
returns the measured count. A real MHA-to-NRRD merge regression verifies that
the teeth preset publishes and reports two valid disconnected components.

