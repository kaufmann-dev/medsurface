# List command hints broke paths with spaces

Fixed: 2026-07-13 18:04:34 CEST (+0200)

Baseline commit: `22ad71de771337dc660bdfc6d5a063871b32cbcd`

## Symptom

After listing a directory whose path contained spaces, the suggested
`medsurface convert` command interpolated the path without shell quoting. On a
narrow terminal, Rich also inserted hard line breaks into long paths. Copying
the displayed hint therefore produced multiple arguments or multiple commands
instead of one usable conversion command.

## Confirmed cause

The list handler constructed the hint with percent formatting against the raw
`Path` and printed it with Rich's normal width-constrained wrapping. It did not
serialize an argument vector for the current platform or distinguish visual
terminal wrapping from newline insertion.

## Fix

Suggested commands are now formatted from an argument vector with POSIX shell
quoting on POSIX and command-line quoting on Windows. Rich soft wrapping leaves
the hint as one logical line while still allowing the terminal to display it at
its available width. Regressions cover both automatic-default and explicit-ID
hints, including spaces and shell metacharacters in an input path.
