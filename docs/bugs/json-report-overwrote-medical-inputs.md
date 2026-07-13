# JSON report paths could overwrite medical image data

Fixed: 2026-07-13 12:50:14 UTC (+0000)

Baseline commit: `260e7198a198f73001eb4695a548f0a1b2d2b729`

## Symptom

`convert --json` and `merge --json` accepted a source image, DICOM instance,
detached payload, or mesh output as the report destination. The command wrote
the mesh first and then replaced the colliding file with JSON while returning a
successful result. Symlink and hard-link aliases behaved the same way.

## Confirmed cause

The CLI opened the report destination directly after processing. It had no
filesystem-identity check against discovered input files or the mesh output,
and report publication was not atomic.

## Fix

Catalog file sources now retain every detached payload path, and all candidates
expose the complete set of files required to load them. Convert and merge reject
mesh/report aliases against each other and every discovered medical input before
loading the processing engine. The shared comparison resolves lexical aliases,
symlinks, and hard links. JSON is written, flushed, and atomically replaced from
a temporary sibling; a failed write leaves an existing report unchanged.

