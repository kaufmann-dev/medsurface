"""Render interactive CLI progress in a GIL-independent helper process."""

from __future__ import annotations

import json
import queue
import sys
import threading
import time

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.text import Text


def main() -> None:
    initial = sys.argv[1] if len(sys.argv) > 1 else "Working ..."
    console = Console(highlight=False, markup=False)
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
