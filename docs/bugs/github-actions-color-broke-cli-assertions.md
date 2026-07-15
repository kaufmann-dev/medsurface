# GitHub Actions color broke CLI assertions

Fixed: 2026-07-15 04:15:52 CEST (+0200)

Baseline commit: `e12b6365989ffd04168ffaa62177c798655891c6`

## Symptom

The complete test suite passed locally, but every GitHub Actions Python matrix
job failed CLI tests that checked rejected option names such as
`--post-smooth-iters` and `--components`. The first affected revision had four
failures; later CLI coverage increased the count to 24.

## Confirmed root cause

GitHub Actions sets `CI=true`, which makes Typer force terminal styling even
when its test runner captures stderr. ANSI escape sequences were inserted
inside styled option names, so the rendered diagnostics looked correct but the
captured strings no longer contained each option as one plain substring. The
new assertions exposed this pre-existing environment difference; no production
CLI behavior or dependency had changed.

## Fix

The shared `CliRunner` now sets Typer's test-only
`_TYPER_FORCE_DISABLE_TERMINAL` environment variable. Captured command output is
therefore deterministic plain text on developer machines and GitHub Actions.
Explicit Rich console and pseudo-terminal tests continue to exercise colored,
interactive CLI rendering independently. The formerly failing 24-test subset
passes with both `CI=true` and `GITHUB_ACTIONS=true`.
