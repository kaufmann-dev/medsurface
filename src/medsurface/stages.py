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


class Reporter:
    """Optional log and warning callbacks plus timed stages for one workflow.

    Every warning is streamed to ``warn`` as it happens and collected in
    :attr:`warnings` for the workflow's result.
    """

    def __init__(
        self,
        log: Logger | None,
        warn: Logger | None,
        progress: ProgressSink | None,
        width: int,
    ) -> None:
        self._log = log
        self._warn = warn
        self.progress = progress
        self.width = width
        self.warnings: list[str] = []
        self.step = stage_runner(self.say, progress, width)

    def say(self, message: str) -> None:
        if self._log is not None:
            self._log(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)
        if self._warn is not None:
            self._warn(message)

    def warn_all(self, messages: list[str]) -> None:
        for message in messages:
            self.warn(message)
