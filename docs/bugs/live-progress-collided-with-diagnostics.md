# Live progress collided with warnings and errors

Fixed: 2026-07-13 17:50:33 CEST (+0200)

Baseline commit: `b80012689a55938e84257268e33f1d0890246621`

## Symptom

During interactive conversion and merging, a warning could be appended directly
to the spinner's elapsed-time field, producing output such as
`0:00:04Warning: ...`. Merge also displayed fixed-volume geometry warnings only
after the moving volume had loaded. Successful validation ended with a
terminal-wide `Problems` panel whose only content was `None`.

## Confirmed cause

The GIL-independent helper process rendered progress on stdout while the parent
process rendered warnings and errors independently on stderr. Both Rich consoles
wrote to the same terminal without coordinating the live line. Merge loaded both
volumes before asking either one for post-load geometry warnings. The quality
renderer unconditionally expanded a problems panel even for an empty problem
list.

## Fix

An active progress display now routes warnings and errors to the helper process.
The helper clears the live spinner, writes the diagnostic to stderr, and redraws
the current stage, preserving both clean terminal output and stream semantics.
The in-process fallback applies the same sequence. Merge emits each volume's
post-load warnings immediately after that volume loads. Valid quality reports
omit the empty problems panel; invalid reports retain a compact `Mesh problems`
panel.
