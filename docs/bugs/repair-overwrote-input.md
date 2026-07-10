# Repair allowed its input to be overwritten

Fixed: 2026-07-10 08:42:32 UTC (+0000)

Baseline commit: `e926407dfe22af1449524f41966b530eb8298613`

## Symptom

`dicom-surface repair model.stl -o model.stl` completed successfully and
replaced the original mesh, despite repair being documented as producing a new
file.

## Cause

The repair library loaded the input and then saved directly to the output path
without checking whether both paths identified the same file.

## Fix

Repair now rejects identical paths, resolved path aliases, symlinks, and hard
links before loading or modifying the mesh. The guard is in the repair library,
so direct Python callers receive the same protection as the CLI.
