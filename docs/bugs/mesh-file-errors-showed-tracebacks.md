# Mesh file errors exposed Python tracebacks

Fixed: 2026-07-10 08:42:32 UTC (+0000)

Baseline commit: `e926407dfe22af1449524f41966b530eb8298613`

## Symptom

Validating or repairing a missing, unreadable, or malformed mesh printed an
internal Python traceback. Other commands already presented expected user input
errors as concise CLI messages.

## Cause

The `validate` and `repair` handlers called their mesh libraries without an
error boundary, so expected file, format, and runtime loader exceptions escaped
through the CLI entry point.

## Fix

Both handlers catch expected `OSError`, `ValueError`, and `RuntimeError`
failures, print one actionable error to stderr, and return exit status 1.
Unexpected exception types still propagate so programming defects are not
silently hidden.
