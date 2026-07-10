# Root command loaded VTK before showing help

Fixed: 2026-07-10 09:12:51 UTC (+0000)

Baseline commit: `aadebd0eac26f1608c6b2d3b4dfb116244747b25`

## Symptom

Running `dicom-surface` with no arguments appeared to hang during startup and
could show a traceback from VTK if interrupted. After startup completed, the
command reported a missing-subcommand error instead of showing the available
commands.

## Cause

The CLI imported merge, conversion, and validation modules before parsing its
arguments. Those modules transitively imported VTK even when the requested
operation only needed command help. The root parser also required a subcommand,
so an empty invocation could not fall back to its help text.

## Fix

Command handlers now import their processing dependencies only when invoked. A
lightweight defaults module keeps parser construction independent of the merge
implementation. Running the root command with no arguments prints the same
command overview as `--help` and exits successfully, without importing VTK or
the processing pipeline. Invalid command names, missing arguments, unknown
flags, and invalid choices are likewise rejected by the parser before heavy
processing imports and produce concise usage errors without tracebacks.
