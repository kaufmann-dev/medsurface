# Progress renderer printed a traceback on Ctrl+C

Fixed: 2026-07-15 04:05:19 CEST (+0200)

Baseline commit: `92fc2f8b39180da4e560318679a395dc78d79ba9`

## Symptom

Pressing `Ctrl+C` during an interactive command printed a Python traceback from
`medsurface.progress_renderer`, ending in `KeyboardInterrupt`. The main command
also received the interrupt, so the terminal could show repeated interrupt
characters while the two processes shut down independently.

## Confirmed cause

The main CLI and its progress-renderer helper belonged to the terminal's same
foreground process group. A terminal `SIGINT` therefore reached both processes.
The main CLI converted its `KeyboardInterrupt` into the expected exit status,
but the renderer had no signal policy and exposed its unhandled exception.

## Fix

The renderer now ignores `SIGINT`, leaving the main CLI as the single owner of
terminal cancellation. The progress context stops and reaps the renderer, resets
its active state, and prints one `Cancelled.` diagnostic before allowing the
interrupt to produce exit status 130. Additional interrupts are ignored during
that cleanup window rather than bypassing it. Quiet commands use the same
cancellation message without starting a renderer. Subprocess regressions send
`SIGINT` to a real pseudo-terminal process group and to a quiet CLI process,
verifying clean status 130 exits, no traceback, no orphaned renderer, and no
published output.
