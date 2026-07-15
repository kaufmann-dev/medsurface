"""Render interactive CLI progress in a GIL-independent helper process."""

from __future__ import annotations

import json
import os
import queue
import signal
import sys
import threading
import time

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.text import Text


def main() -> None:
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    initial = sys.argv[1] if len(sys.argv) > 1 else "Working ..."
    if sys.stdout.isatty() and os.environ.get("TERM", "").casefold() in {"dumb", "unknown"}:
        os.environ["TERM"] = "xterm-256color"
    # Rich treats TERM=dumb as non-interactive even when stdout is a real PTY.
    # The parent only starts this helper for a terminal, so the descriptor is
    # the authoritative signal here.
    console = Console(
        highlight=False,
        markup=False,
        force_terminal=sys.stdout.isatty(),
        force_interactive=sys.stdout.isatty(),
    )
    diagnostic_console = Console(
        stderr=True,
        highlight=False,
        markup=False,
        force_terminal=sys.stderr.isatty(),
    )
    messages: queue.SimpleQueue[dict[str, str]] = queue.SimpleQueue()
    input_closed = threading.Event()

    def read_messages() -> None:
        try:
            for line in sys.stdin:
                messages.put(json.loads(line))
        finally:
            input_closed.set()

    threading.Thread(target=read_messages, daemon=True).start()
    progress = Progress(
        SpinnerColumn(style="cyan"),
        TextColumn("{task.description}", markup=False),
        TimeElapsedColumn(),
        console=console,
        transient=True,
        auto_refresh=False,
    )
    progress.start()
    task_id = progress.add_task(initial, total=None)
    stopping = False
    try:
        while not stopping:
            while True:
                try:
                    message = messages.get_nowait()
                except queue.Empty:
                    break
                kind = message.get("kind")
                text = message.get("message", "")
                if kind == "update":
                    progress.reset(task_id, description=text, total=None)
                elif kind == "log":
                    progress.console.print(Text(text, style="cyan"))
                elif kind in {"warning", "error"}:
                    prefix, style = {
                        "warning": ("Warning: ", "bold yellow"),
                        "error": ("Error: ", "bold red"),
                    }[kind]
                    diagnostic = Text()
                    diagnostic.append(prefix, style=style)
                    diagnostic.append(text)
                    progress.stop()
                    diagnostic_console.print(diagnostic)
                    progress.start()
                elif kind == "stop":
                    stopping = True
            progress.refresh()
            if input_closed.is_set() and messages.empty():
                stopping = True
            if not stopping:
                time.sleep(0.1)
    finally:
        progress.stop()


if __name__ == "__main__":
    main()
