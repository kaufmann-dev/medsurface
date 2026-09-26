"""Timed processing stages shared by every workflow."""

from __future__ import annotations

import time
from typing import Any, Callable

Logger = Callable[[str], None]
#: Receives machine-readable progress events such as
#: ``{"event": "stage_start", "stage": "marching cubes"}``.
ProgressSink = Callable[[dict[str, Any]], None]


def stage_runner(
    say: Logger,
    progress: ProgressSink | None,
    width: int,
    **context: Any,
):
    """Return ``step(message, function)`` that logs and reports one timed stage."""

    def step(message: str, function):
        say("%s ..." % message)
        if progress is not None:
            progress({"event": "stage_start", "stage": message, **context})
        before = time.time()
        value = function()
        seconds = time.time() - before
        say("  %-*s %6.1fs" % (width, message, seconds))
        if progress is not None:
            progress(
                {"event": "stage_end", "stage": message, "seconds": seconds, **context}
            )
        return value

    return step
