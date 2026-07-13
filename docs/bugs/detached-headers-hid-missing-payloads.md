# Detached image headers hid missing payloads during listing

Fixed: 2026-07-13 12:50:14 UTC (+0000)

Baseline commit: `260e7198a198f73001eb4695a548f0a1b2d2b729`

## Symptom

An `.mhd` or `.nhdr` header could appear usable and become the automatic default
after its referenced payload had been removed. Conversion then failed only when
SimpleITK attempted to load pixels.

## Confirmed cause

Discovery asked SimpleITK to read image information, which can parse these
headers without opening the detached payload. The catalog did not resolve or
check the reference itself.

## Fix

Discovery resolves single-file, listed, and integer-sequence payload references
relative to the header. Missing or unreadable payloads make the candidate
unusable immediately, and payload paths remain attached to the source for load
identity and output-overwrite protection. Regressions cover both `.mhd` and
`.nhdr`.

