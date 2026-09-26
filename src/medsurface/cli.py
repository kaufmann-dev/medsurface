"""Typer command line interface with Rich human-readable output."""

from __future__ import annotations

import json
import math
import os
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import warnings
from collections import Counter
from contextlib import contextmanager, nullcontext
from contextvars import ContextVar, Token
from enum import Enum
from pathlib import Path
from typing import Any

import typer
from click.core import ParameterSource
from rich import box
from rich.console import Console, Group
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TaskID, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich.text import Text

from . import __version__, defaults
from . import presets as presets_mod
from .outputs import VolumeOutput, extension, volume_output
from .presets import PRESETS


class PresetChoice(str, Enum):
    """Tissue presets accepted by ``--preset``."""

    AUTO = "auto"
    BONE = "bone"
    SKIN = "skin"
    TEETH = "teeth"


class ComponentChoice(str, Enum):
    """Surface-component retention accepted by ``--components``."""

    ALL = "all"
    LARGEST = "largest"


class DestepChoice(str, Enum):
    """Stair-step fairing regions accepted by ``--destep``."""

    AUTO = "auto"
    ALL = "all"
    BAND = "band"


class AxisChoice(str, Enum):
    """Model axes accepted by ``--destep-axis``."""

    X = "x"
    Y = "y"
    Z = "z"


def _keep_largest_component(
    components: ComponentChoice | None,
) -> bool | None:
    if components is None:
        return None
    return components is ComponentChoice.LARGEST


_DESTEP_OPTION = typer.Option(
    None,
    "--destep",
    help="Final masked fairing that removes broad stair-step ripples: auto-detected "
    "smooth regions, all vertices, or an axis band.",
    show_default="off",
)
_DESTEP_AXIS_OPTION = typer.Option(
    AxisChoice.Z, "--destep-axis", help="Axis of the --destep band mask."
)
_DESTEP_FULL_OPTION = typer.Option(
    None,
    "--destep-full-mm",
    help="Model coordinate in mm at and beyond which the band is fully faired.",
)
_DESTEP_FROZEN_OPTION = typer.Option(
    None,
    "--destep-frozen-mm",
    help="Model coordinate in mm at and beyond which the band never moves.",
)
_DESTEP_ITERS_OPTION = typer.Option(
    defaults.DEFAULT_DESTEP_ITERS,
    "--destep-iters",
    help="Taubin fairing iterations for --destep.",
)
_DESTEP_MAX_OPTION = typer.Option(
    defaults.DEFAULT_DESTEP_MAX_MM,
    "--destep-max-mm",
    help="Largest per-vertex displacement in mm that --destep may introduce.",
)


def _given(ctx: typer.Context, parameter: str) -> bool:
    """Whether the user passed an option rather than relying on its default."""
    return ctx.get_parameter_source(parameter) is not ParameterSource.DEFAULT


def _destep_settings(
    ctx: typer.Context,
    region: DestepChoice | None,
    axis: AxisChoice,
    full_mm: float | None,
    frozen_mm: float | None,
    iterations: int,
    max_mm: float,
) -> presets_mod.DestepSettings | None:
    """Translate stair-step options into settings or exit with a usage error."""
    band_options = [
        name
        for name, given in (
            ("--destep-axis", _given(ctx, "destep_axis")),
            ("--destep-full-mm", full_mm is not None),
            ("--destep-frozen-mm", frozen_mm is not None),
        )
        if given
    ]
    tuning_options = [
        name
        for name, parameter in (
            ("--destep-iters", "destep_iters"),
            ("--destep-max-mm", "destep_max_mm"),
        )
        if _given(ctx, parameter)
    ]

    def require(options: list[str], requirement: str) -> None:
        verb = "requires" if len(options) == 1 else "require"
        _error("%s %s %s" % (", ".join(options), verb, requirement))
        raise typer.Exit(2)

    if region is None:
        if band_options or tuning_options:
            require(band_options + tuning_options, "--destep")
        return None
    if region is not DestepChoice.BAND and band_options:
        require(band_options, "--destep band")
    if region is DestepChoice.BAND:
        if full_mm is None or frozen_mm is None:
            _error("--destep band requires --destep-full-mm and --destep-frozen-mm")
            raise typer.Exit(2)
        if not (math.isfinite(full_mm) and math.isfinite(frozen_mm)):
            _error("--destep-full-mm and --destep-frozen-mm must be finite")
            raise typer.Exit(2)
        if full_mm == frozen_mm:
            _error("--destep-full-mm and --destep-frozen-mm must differ")
            raise typer.Exit(2)
    if iterations < 1:
        _error("--destep-iters must be at least 1")
        raise typer.Exit(2)
    _validate_processing_numbers(nonnegative=[], positive=(("--destep-max-mm", max_mm),))
    return presets_mod.DestepSettings(
        region=region.value,
        iterations=iterations,
        max_displacement_mm=float(max_mm),
        axis=axis.value,
        full_mm=full_mm,
        frozen_mm=frozen_mm,
    )


stdout_console = Console(highlight=False, markup=False)
stderr_console = Console(stderr=True, highlight=False, markup=False)

app = typer.Typer(
    add_completion=False,
    help="Convert volumes, fuse scans or labelmaps, and extract surface meshes.",
    invoke_without_command=True,
    no_args_is_help=False,
    pretty_exceptions_show_locals=False,
    rich_markup_mode="rich",
)
labelmap_app = typer.Typer(
    add_completion=False,
    help="Fuse labelmaps or extract their surface meshes.",
    invoke_without_command=True,
    no_args_is_help=False,
    pretty_exceptions_show_locals=False,
    rich_markup_mode="rich",
)
app.add_typer(labelmap_app, name="labelmap")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo("medsurface %s" % __version__)
        raise typer.Exit()


@app.callback()
def root(
    ctx: typer.Context,
    version: bool = typer.Option(
        False,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    """Convert volumes, fuse scans or labelmaps, and extract surface meshes."""
    if ctx.invoked_subcommand is None:
        stdout_console.print(ctx.get_help())


@labelmap_app.callback()
def labelmap_root(ctx: typer.Context) -> None:
    """Fuse labelmaps or extract their surface meshes."""
    if ctx.invoked_subcommand is None:
        stdout_console.print(ctx.get_help())


def _styled_message(prefix: str, prefix_style: str, message: object) -> Text:
    text = Text()
    text.append(prefix, style=prefix_style)
    text.append(str(message))
    return text


def _log(message: str) -> None:
    stdout_console.print(Text(str(message), style="cyan"))


_active_progress: ContextVar[Any] = ContextVar("active_progress", default=None)


class _ProgressDisplay:
    """Indeterminate stage progress with a plain-text redirected fallback."""

    def __init__(
        self,
        enabled: bool,
        initial: str,
        console: Console = stdout_console,
        diagnostic_console: Console = stderr_console,
    ) -> None:
        self.enabled = enabled
        self.initial = initial
        self.console = console
        self.diagnostic_console = diagnostic_console
        self.interactive = enabled and console.is_terminal and console.is_interactive
        self.progress: Progress | None = None
        self.task_id: TaskID | None = None
        self.renderer: subprocess.Popen[str] | None = None
        self._context_token: Token[Any] | None = None

    def _start_renderer(self) -> bool:
        """Start a renderer process when the console has a real output descriptor."""
        try:
            self.console.file.fileno()
            self.diagnostic_console.file.fileno()
        except (AttributeError, OSError, ValueError):
            return False
        try:
            self.renderer = subprocess.Popen(
                [sys.executable, "-m", "medsurface.progress_renderer", self.initial],
                stdin=subprocess.PIPE,
                stdout=self.console.file,
                stderr=self.diagnostic_console.file,
                text=True,
                bufsize=1,
            )
        except OSError:
            self.renderer = None
            return False
        return True

    def _send(self, kind: str, message: str = "") -> bool:
        if self.renderer is None or self.renderer.stdin is None:
            return False
        try:
            self.renderer.stdin.write(
                json.dumps({"kind": kind, "message": message}) + "\n"
            )
            self.renderer.stdin.flush()
        except (BrokenPipeError, OSError):
            return False
        return True

    def __enter__(self):
        self._context_token = _active_progress.set(self)
        if not self.enabled:
            return self
        if self.interactive:
            if self._start_renderer():
                return self
            self.progress = Progress(
                SpinnerColumn(style="cyan"),
                TextColumn("{task.description}", markup=False),
                TimeElapsedColumn(),
                console=self.console,
                transient=True,
            )
            self.progress.start()
            self.task_id = self.progress.add_task(self.initial, total=None)
        else:
            self._print(self.initial)
        return self

    def __exit__(self, _exc_type, _exc, _traceback) -> None:
        cancelled = _exc_type is KeyboardInterrupt
        previous_sigint_handler = None
        if cancelled:
            previous_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        try:
            if self.renderer is not None:
                self._send("stop")
                if self.renderer.stdin is not None:
                    self.renderer.stdin.close()
                try:
                    self.renderer.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    self.renderer.terminate()
                    self.renderer.wait(timeout=2)
            if self.progress is not None:
                self.progress.stop()
        finally:
            try:
                if self._context_token is not None:
                    _active_progress.reset(self._context_token)
                    self._context_token = None
                if cancelled:
                    self.diagnostic_console.print(Text("Cancelled.", style="yellow"))
            finally:
                if previous_sigint_handler is not None:
                    signal.signal(signal.SIGINT, previous_sigint_handler)

    def _print(self, message: str) -> None:
        self.console.print(Text(message, style="cyan"))

    def update(self, message: str) -> None:
        """Set the active stage and restart its elapsed-time counter."""
        if not self.enabled:
            return
        message = str(message).strip()
        if self.renderer is not None:
            self._send("update", message)
            return
        if self.progress is not None and self.task_id is not None:
            self.progress.reset(self.task_id, description=message, total=None)
            self.progress.refresh()
        else:
            self._print(message)

    def log(self, message: str) -> None:
        """Translate processing start messages into live stage updates."""
        if not self.enabled:
            return
        rendered = str(message)
        if rendered.strip().endswith("..."):
            self.update(rendered)
        elif self.renderer is not None:
            self._send("log", rendered)
        elif self.progress is not None:
            self.progress.console.print(Text(rendered, style="cyan"))
        else:
            self._print(rendered)

    def diagnostic(self, kind: str, message: object) -> None:
        """Print a warning or error without colliding with the live display."""
        rendered = str(message)
        if self.renderer is not None and self._send(kind, rendered):
            return
        if self.progress is not None:
            self.progress.stop()
            try:
                _print_diagnostic(self.diagnostic_console, kind, rendered)
            finally:
                self.progress.start()
            return
        _print_diagnostic(self.diagnostic_console, kind, rendered)


def _success(message: object) -> None:
    stdout_console.print(_styled_message("Success: ", "bold green", message))


def _print_diagnostic(console: Console, kind: str, message: object) -> None:
    prefix, style = {
        "warning": ("Warning: ", "bold yellow"),
        "error": ("Error: ", "bold red"),
    }[kind]
    console.print(_styled_message(prefix, style, message))


def _warn(message: object) -> None:
    progress = _active_progress.get()
    if progress is not None:
        progress.diagnostic("warning", message)
    else:
        _print_diagnostic(stderr_console, "warning", message)


def _error(message: object) -> None:
    progress = _active_progress.get()
    if progress is not None:
        progress.diagnostic("error", message)
    else:
        _print_diagnostic(stderr_console, "error", message)


class ProgressChoice(str, Enum):
    """Progress formats accepted by ``--progress``."""

    HUMAN = "human"
    JSON = "json"


_PROGRESS_OPTION = typer.Option(
    ProgressChoice.HUMAN,
    "--progress",
    help="Progress format: human-readable stages, or JSON lines on stderr for "
    "other programs (status, stage, label, warning, and error events).",
)


class _JsonProgress(_ProgressDisplay):
    """JSON-lines progress on stderr for programs driving the CLI."""

    def __init__(self, initial: str) -> None:
        super().__init__(False, initial)

    def __enter__(self):
        super().__enter__()
        self.update(self.initial)
        return self

    @staticmethod
    def emit(event: dict[str, Any]) -> None:
        sys.stderr.write(json.dumps(event, default=str) + "\n")
        sys.stderr.flush()

    def update(self, message: str) -> None:
        self.emit({"event": "status", "message": str(message).strip()})

    def log(self, message: str) -> None:
        return None

    def diagnostic(self, kind: str, message: object) -> None:
        self.emit({"event": kind, "message": str(message)})


def _progress_display(choice: ProgressChoice, quiet: bool, initial: str) -> _ProgressDisplay:
    if choice is ProgressChoice.JSON:
        return _JsonProgress(initial)
    return _ProgressDisplay(not quiet, initial)


def _progress_sink(progress: _ProgressDisplay):
    return progress.emit if isinstance(progress, _JsonProgress) else None


def _quality_status(report: dict[str, Any] | None) -> int:
    """Shared exit status for every command that validates a mesh."""
    return 0 if report is None or report["valid"] else 1


def _exit_for_quality(report: dict[str, Any] | None) -> None:
    status = _quality_status(report)
    if status:
        raise typer.Exit(status)


def _warn_if_invalid(report: dict[str, Any], subject: str) -> None:
    if not report["valid"]:
        _warn("%s failed validation: %s" % (subject, "; ".join(report["problems"])))


def _validate_mesh_output(path: Path) -> None:
    output_extension = extension(path)
    if output_extension not in defaults.SUPPORTED_MESH_EXTENSIONS:
        _error(
            "unsupported output extension %r; supported: %s"
            % (output_extension, ", ".join(defaults.SUPPORTED_MESH_EXTENSIONS))
        )
        raise typer.Exit(2)


def _validate_mesh_input(path: Path) -> None:
    input_extension = extension(path)
    if input_extension not in defaults.SUPPORTED_MESH_EXTENSIONS:
        _error(
            "unsupported mesh extension %r; supported: %s"
            % (input_extension, ", ".join(defaults.SUPPORTED_MESH_EXTENSIONS))
        )
        raise typer.Exit(2)


def _validate_volume_output(path: Path) -> VolumeOutput:
    try:
        output = volume_output(path)
        if path.exists() and path.is_dir():
            raise ValueError("volume output path is a directory: %s" % path)
        return output
    except ValueError as exc:
        _error(exc)
        raise typer.Exit(2) from None


def _fusion_result_payload(result: Any) -> dict[str, Any]:
    return {
        "output": result.output_path,
        "format": result.output_format,
        "compression": result.compression,
        "pixel_type": "uint8",
        "components": 1,
        "grid_mm": result.grid_mm,
        "grid_size": list(result.grid_size),
        "grid_origin_mm": list(result.grid_origin_mm),
        "grid_direction": list(result.grid_direction),
        "foreground_fixed_voxels": result.foreground_fixed_voxels,
        "foreground_moving_voxels": result.foreground_moving_voxels,
        "foreground_fused_voxels": result.foreground_fused_voxels,
        "volume_fixed_mm3": result.volume_fixed_mm3,
        "volume_moving_mm3": result.volume_moving_mm3,
        "volume_fused_mm3": result.volume_fused_mm3,
        "seconds": result.seconds,
        "warnings": result.warnings,
    }


def _print_fusion_result(result: Any) -> None:
    _success("wrote %s" % result.output_path)
    _log(
        "%s voxels at %.3f mm   foreground %.0f cm3   %.1fs"
        % (
            "x".join(str(value) for value in result.grid_size),
            result.grid_mm,
            result.volume_fused_mm3 / 1000.0,
            result.seconds,
        )
    )
    _print_command_hint(
        "Extract a surface with:  ",
        [
            "medsurface",
            "labelmap",
            "extract",
            result.output_path,
            "-o",
            "MODEL.stl",
        ],
    )


def _emit_remaining_warnings(messages: list[str], emitted: list[str]) -> None:
    """Render result warnings not already streamed by the processing layer."""
    streamed = Counter(emitted)
    for message in messages:
        if streamed[message]:
            streamed[message] -= 1
        else:
            _warn(message)


def _parse_threshold(
    value: str | None, option: str = "--threshold"
) -> float | str | None:
    """Return a numeric threshold, the ``auto`` sentinel, or ``None``."""
    if value is None:
        return None
    if value == "auto":
        return value
    try:
        parsed = float(value)
    except ValueError:
        _error("%s must be a number or 'auto'" % option)
        raise typer.Exit(2) from None
    if not math.isfinite(parsed):
        _error("%s must be finite or 'auto'" % option)
        raise typer.Exit(2)
    return parsed


def _validate_processing_numbers(
    *,
    nonnegative: list[tuple[str, float | int | None]],
    positive: tuple[tuple[str, float | int | None], ...] = (),
    unit_interval: tuple[tuple[str, float | int | None], ...] = (),
) -> None:
    for option, value in nonnegative:
        if value is not None and (not math.isfinite(value) or value < 0):
            _error("%s must be finite and non-negative" % option)
            raise typer.Exit(2)
    for option, value in positive:
        if value is not None and (not math.isfinite(value) or value <= 0):
            _error("%s must be finite and greater than zero" % option)
            raise typer.Exit(2)
    for option, value in unit_interval:
        if value is not None and (not math.isfinite(value) or not 0 < value <= 1):
            _error("%s must be finite, greater than zero, and at most one" % option)
            raise typer.Exit(2)


def _emit_discovery_warnings(captured: list[warnings.WarningMessage]) -> None:
    messages = [str(item.message) for item in captured]
    counts = Counter(messages)
    emitted: set[str] = set()
    for message in messages:
        if message in emitted:
            continue
        emitted.add(message)
        count = counts[message]
        suffix = " (repeated %d times)" % count if count > 1 else ""
        _warn(message + suffix)


def _discover(root: Path):
    """Discover volumes while presenting reader warnings concisely."""
    from . import catalog

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        found = catalog.discover(root)
    _emit_discovery_warnings(captured)
    return found


def _select_labelmap(path: Path, role: str | None = None):
    """Resolve one direct self-describing image file as a labelmap candidate."""
    found = _discover(path)
    prefix = "%s input: " % role if role else ""
    if not found:
        _error("%sno supported NIfTI, NRRD, or MetaImage volume found" % prefix)
        raise typer.Exit(1)
    if len(found) != 1 or found[0].dicom is not None:
        _error("%slabelmap must be one NIfTI, NRRD, or MetaImage file" % prefix)
        raise typer.Exit(2)

    from . import catalog

    try:
        return catalog.select(found, None), found
    except ValueError as exc:
        _error("%s%s" % (prefix, exc))
        raise typer.Exit(1) from None


def _plain(value: object, style: str | None = None) -> Text:
    return Text(str(value), style=style or "")


def _detail_block(details: list[tuple[str, object | None]]) -> str:
    """Render only human-readable details that carry an actual value."""
    lines: list[str] = []
    for label, value in details:
        if value is None:
            continue
        text = str(value).strip()
        if not text:
            continue
        lines.append("%s: %s" % (label, text))
    return "\n".join(lines)


def _format_shell_command(arguments: list[object]) -> str:
    rendered = [str(argument) for argument in arguments]
    if os.name == "nt":
        return subprocess.list2cmdline(rendered)
    return shlex.join(rendered)


def _print_command_hint(prefix: str, arguments: list[object]) -> None:
    command = Text(prefix)
    command.append(_format_shell_command(arguments), style="bold")
    stdout_console.print(command, soft_wrap=True)


def _volume_status(candidate: Any, recommended: Any | None) -> str:
    notes: list[str] = []
    if candidate is recommended:
        notes.append("default")
    if candidate.unusable_reason:
        notes.append(candidate.unusable_reason)
    else:
        notes.append("usable")
    series = candidate.dicom
    if series is not None and series.n_parts > 1:
        notes.append("orientation %d of %d in this UID" % (series.part, series.n_parts))
    if candidate.sharp_kernel:
        notes.append("sharp kernel")
    return "; ".join(notes)


def _volume_table(found: list[Any], recommended: Any | None) -> Table:
    """Build the responsive, format-neutral table used by ``list``."""
    table = Table(
        box=box.SIMPLE_HEAVY,
        expand=True,
        padding=(0, 1),
        row_styles=("", "dim"),
    )
    table.add_column("ID", justify="right", no_wrap=True)
    table.add_column("Input", ratio=2, overflow="fold")
    table.add_column("Metadata", ratio=3, overflow="fold")
    table.add_column("Geometry", ratio=2, overflow="fold")
    table.add_column("Status", ratio=2, overflow="fold", min_width=16)

    for candidate in found:
        voxel = None
        if candidate.spacing is not None and len(candidate.spacing) == 3:
            voxel = "%.3g × %.3g × %.3g" % candidate.spacing
        row_style = "bold cyan" if candidate is recommended else None
        series = candidate.dicom
        kernel = (
            series.kernel_display
            if series is not None and series.kernel_values
            else None
        )
        modality = None if candidate.modality == "?" else candidate.modality
        plane = candidate.plane
        if plane is not None and plane.casefold() == "unknown":
            plane = None
        input_details = _detail_block(
            [("Format", candidate.format), ("Source", candidate.source_name)]
        )
        metadata = _detail_block(
            [
                ("DICOM #", candidate.series_number),
                ("Modality", modality),
                ("Description", candidate.description),
                ("Kernel", kernel),
            ]
        )
        geometry = _detail_block(
            [
                ("Slices", candidate.slices),
                ("Voxel (mm)", voxel),
                ("Plane", plane),
            ]
        )
        table.add_row(
            _plain(candidate.id),
            _plain(input_details),
            _plain(metadata),
            _plain(geometry),
            _plain(_volume_status(candidate, recommended)),
            style=row_style,
        )
    return table


def _quality_table(report: dict[str, Any]) -> Table:
    """Build the shared two-column mesh-quality table."""
    table = Table(title="Mesh quality", box=box.SIMPLE_HEAVY, show_header=False)
    table.add_column("Metric", style="bold", no_wrap=True)
    table.add_column("Value")

    valid = bool(report["valid"])
    table.add_row(
        _plain("Status"),
        _plain("valid" if valid else "invalid", "green" if valid else "red"),
    )
    table.add_row(_plain("Triangles"), _plain(f"{report['triangles']:,}"))
    table.add_row(_plain("Vertices"), _plain(f"{report['vertices']:,}"))
    table.add_row(_plain("Components"), _plain(f"{report['components']:,}"))
    table.add_row(_plain("Watertight"), _plain("yes" if report["watertight"] else "no"))
    table.add_row(
        _plain("Winding consistent"),
        _plain("yes" if report["winding_consistent"] else "no"),
    )
    table.add_row(_plain("Boundary edges"), _plain(f"{report['boundary_edges']:,}"))
    table.add_row(
        _plain("Holes"),
        _plain(f"{report['holes']:,}"),
    )
    table.add_row(
        _plain("Disoriented faces"),
        _plain(f"{report['disoriented_faces']:,}"),
    )
    if "genus" in report:
        table.add_row(_plain("Genus"), _plain(report["genus"]))
    if report.get("volume_mm3") is None:
        table.add_row(_plain("Volume"), _plain("undefined (mesh is not closed)"))
    else:
        table.add_row(_plain("Volume"), _plain("%.0f mm³" % report["volume_mm3"]))
    intersections = report.get("self_intersecting_faces")
    if intersections is not None:
        value = (
            f"{intersections:,}"
            if isinstance(intersections, int)
            else str(intersections)
        )
        table.add_row(_plain("Self-intersections"), _plain(value))
    extents = report["bbox_extents_mm"]
    table.add_row(
        _plain("Bounding box"),
        _plain("%.1f × %.1f × %.1f mm" % (extents[0], extents[1], extents[2])),
    )
    table.add_row(_plain("Size"), _plain("%.1f MB" % (report["bytes"] / 1048576.0)))
    return table


def _print_quality(report: dict[str, Any], console: Console = stdout_console) -> None:
    console.print(_quality_table(report))
    problems = report["problems"]
    if not problems:
        return
    lines = [_plain("• " + problem) for problem in problems]
    console.print(
        Panel(
            Group(*lines),
            title=_plain("Mesh problems"),
            border_style="red",
            expand=False,
        )
    )


def _write_json_file(path: Path, payload: object) -> None:
    destination = path.absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=".%s." % destination.name,
        dir=destination.parent,
        text=True,
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


@contextmanager
def _restore_output_on_failure(path: Path):
    """Roll back a primary output if its optional report cannot be written."""
    destination = path.absolute()
    backup: str | None = None
    existed = os.path.lexists(destination)
    original_stat: os.stat_result | None = None
    if existed:
        if not destination.is_file() and not destination.is_symlink():
            raise ValueError("output path is not a file: %s" % destination)
        original_stat = os.lstat(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        descriptor, backup = tempfile.mkstemp(
            prefix=".%s.rollback." % destination.name,
            dir=destination.parent,
        )
        os.close(descriptor)
        os.unlink(backup)
        try:
            os.link(destination, backup, follow_symlinks=False)
        except OSError:
            shutil.copy2(destination, backup, follow_symlinks=False)
    try:
        yield
    except BaseException:
        unchanged = False
        if original_stat is not None and os.path.lexists(destination):
            unchanged = os.path.samestat(os.lstat(destination), original_stat)
        if unchanged and backup is not None:
            os.unlink(backup)
        elif backup is not None and os.path.lexists(backup):
            if os.path.lexists(destination):
                os.unlink(destination)
            os.replace(backup, destination)
        elif os.path.lexists(destination):
            os.unlink(destination)
        raise
    else:
        if backup is not None:
            if os.path.lexists(destination):
                os.unlink(backup)
            else:
                os.replace(backup, destination)
    finally:
        if backup is not None and os.path.lexists(backup):
            os.unlink(backup)


def _protect_output_paths(
    candidates: list[Any], output: Path, json_file: Path | None
) -> None:
    from . import catalog
    from .paths import protect_outputs

    inputs = [
        path for candidate in candidates for path in catalog.source_paths(candidate)
    ]
    protect_outputs(inputs, output, json_file)


@app.command("list")
def list_volumes(
    input_path: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        help="Volume file or directory tree containing supported volumes.",
    ),
    json_file: Path | None = typer.Option(
        None, "--json", help="Write the volume list to this JSON file."
    ),
    progress_format: ProgressChoice = _PROGRESS_OPTION,
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress the volume table and progress output."
    ),
) -> None:
    """Show every supported volume under an input path."""
    progress = _progress_display(progress_format, quiet, "Discovering volumes ...")
    with progress:
        found = _discover(input_path)
        if not found:
            _error("no supported volumes found under %s" % input_path)
            raise typer.Exit(1)
        if json_file is not None:
            try:
                _protect_output_paths(found, json_file, None)
            except ValueError as exc:
                _error(exc)
                raise typer.Exit(2) from None

    from . import catalog

    recommended = catalog.recommended(found)
    if json_file is not None:
        payload = [
            {
                "id": candidate.id,
                "default": candidate is recommended,
                "format": candidate.format,
                "source": candidate.source_name,
                "modality": candidate.modality,
                "description": candidate.description,
                "size": list(candidate.size) if candidate.size else None,
                "spacing": list(candidate.spacing) if candidate.spacing else None,
                "origin": list(candidate.origin) if candidate.origin else None,
                "direction": list(candidate.direction) if candidate.direction else None,
                "pixel_type": candidate.pixel_type,
                "components": candidate.components,
                "plane": candidate.plane,
                "usable": candidate.usable,
                "unusable_reason": candidate.unusable_reason,
                "dicom": (
                    {
                        "uid": candidate.dicom.uid,
                        "part": candidate.dicom.part,
                        "n_parts": candidate.dicom.n_parts,
                        "series_number": candidate.dicom.series_number,
                        "kernel": list(candidate.dicom.kernel_values),
                        "sharp_kernel": candidate.dicom.sharp_kernel,
                        "spacing_uniform": candidate.dicom.spacing_uniform,
                        "spacing_spread_mm": candidate.dicom.spacing_spread_mm,
                        "localizer": candidate.dicom.is_localizer,
                        "image_type": list(candidate.dicom.image_type),
                        "rescale_type": candidate.dicom.rescale_type,
                        "multi_energy_ct_acquisition": (
                            candidate.dicom.multi_energy_ct_acquisition
                        ),
                        "hu_calibration_verified": candidate.dicom.has_calibrated_hu,
                    }
                    if candidate.dicom is not None
                    else None
                ),
            }
            for candidate in found
        ]
        try:
            _write_json_file(json_file, payload)
        except OSError as exc:
            _error("cannot write JSON report %s: %s" % (json_file, exc))
            raise typer.Exit(1) from None
    if quiet:
        return

    stdout_console.print(_volume_table(found, recommended))
    stdout_console.print(
        Text(
            "Row IDs are local to this discovery result; run list again after directory contents change.",
            style="dim",
        )
    )
    if recommended is not None:
        _print_command_hint(
            "Extract the default with:  ",
            ["medsurface", "extract", input_path, "-o", "out.stl"],
        )
    elif any(candidate.usable for candidate in found):
        ambiguity = catalog.selection_ambiguity(found)
        if ambiguity is not None:
            stdout_console.print(
                Text("No automatic default: %s." % ambiguity, style="yellow")
            )
        _print_command_hint(
            "Choose a volume with:  ",
            [
                "medsurface",
                "extract",
                input_path,
                "--volume",
                "ID",
                "-o",
                "out.stl",
            ],
        )
    if json_file is not None:
        _success("wrote %s" % json_file)


@app.command()
def presets() -> None:
    """Describe the built-in tissue presets."""
    tissue = Table(title="Tissue presets (--preset)", box=box.SIMPLE_HEAVY, expand=True)
    tissue.add_column("Preset", no_wrap=True)
    tissue.add_column("Modality", no_wrap=True)
    tissue.add_column("Threshold", no_wrap=True)
    tissue.add_column("Processing", ratio=2, overflow="fold")
    tissue.add_column("Description", ratio=2, overflow="fold")
    for name in sorted(PRESETS):
        preset = PRESETS[name]
        threshold = (
            preset.threshold
            if isinstance(preset.threshold, str)
            else "%g" % preset.threshold
        )
        simplify = (
            "off"
            if preset.simplify_error_mm == 0
            else "%.2f mm" % preset.simplify_error_mm
        )
        mask_smooth = (
            "off" if preset.mask_smooth_mm == 0 else "%.1f mm" % preset.mask_smooth_mm
        )
        pre_surface_smooth = (
            "off"
            if preset.surface_smooth_iters == 0
            else "%d iter" % preset.surface_smooth_iters
        )
        post_surface_smooth = (
            "off"
            if preset.post_surface_smooth_iters == 0
            else "%d iter" % preset.post_surface_smooth_iters
        )
        processing = (
            "median %.1f mm; closing %.1f mm; mask smooth %s; "
            "pre-mesh smooth %s; simplify %s; post-mesh smooth %s"
        ) % (
            preset.median_mm,
            preset.closing_mm,
            mask_smooth,
            pre_surface_smooth,
            simplify,
            post_surface_smooth,
        )
        tissue.add_row(
            _plain(name),
            _plain("CT/HU" if preset.threshold_unit == "HU" else "any"),
            _plain(threshold),
            _plain(processing),
            _plain(preset.description),
        )
    stdout_console.print(tissue)


@app.command()
def convert(
    input_path: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        help="Volume file or directory tree containing supported volumes.",
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Atomic .nii/.nii.gz/.nrrd/.mha output file.",
    ),
    volume_id: int | None = typer.Option(
        None,
        "--volume",
        min=1,
        help="Integer volume ID displayed by 'medsurface list'.",
        show_default="automatic",
    ),
    strip_metadata: bool = typer.Option(
        False,
        "--strip-metadata",
        help="Remove source metadata while retaining required format headers.",
    ),
    allow_large_volume: bool = typer.Option(
        False,
        "--allow-large-volume",
        help="Bypass the 500-million-voxel source limit; may exhaust memory.",
    ),
    json_file: Path | None = typer.Option(
        None, "--json", help="Write results and provenance to this JSON file."
    ),
    progress_format: ProgressChoice = _PROGRESS_OPTION,
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Preserve one image while changing only its storage format."""
    _validate_volume_output(output)
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _progress_display(progress_format, quiet, "Discovering volumes ...")
    with progress:
        found = _discover(input_path)
        if not found:
            _error("no supported volumes found under %s" % input_path)
            raise typer.Exit(1)

        from . import catalog

        try:
            chosen = catalog.select(found, volume_id)
            _protect_output_paths(found, output, json_file)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        progress.log(
            "volume ID %d  %s  %s" % (chosen.id, chosen.format, chosen.source_name)
        )
        progress.update("Loading volume conversion engine ...")
        from . import volume as volume_mod

        transaction = (
            _restore_output_on_failure(output)
            if json_file is not None
            else nullcontext()
        )
        try:
            with transaction:
                result = volume_mod.convert(
                    candidate=chosen,
                    output_path=str(output),
                    strip_metadata=strip_metadata,
                    allow_large_volume=allow_large_volume,
                    log=progress.log,
                    warn=emit_warning,
                    progress=_progress_sink(progress),
                )
                if json_file is not None:
                    payload = {
                        "result": {
                            "output": result.output_path,
                            "format": result.format,
                            "compression": result.compression,
                            "dimensions": list(result.dimensions),
                            "pixel_type": result.pixel_type,
                            "components": result.components,
                            "spacing_mm": list(result.spacing),
                            "origin_mm": list(result.origin),
                            "direction": list(result.direction),
                            "seconds": result.seconds,
                            "warnings": result.warnings,
                            "metadata_policy": result.metadata_policy,
                        },
                        "provenance": result.provenance,
                    }
                    progress.update("Writing JSON report ...")
                    try:
                        _write_json_file(json_file, payload)
                    except OSError as exc:
                        _error("cannot write JSON report %s: %s" % (json_file, exc))
                        raise typer.Exit(1) from None
        except typer.Exit:
            raise
        except (OSError, RuntimeError, ValueError) as exc:
            _error(exc)
            raise typer.Exit(1) from None

    _emit_remaining_warnings(result.warnings, emitted_warnings)
    if not quiet:
        _success("wrote %s" % result.output_path)
        _log(
            "%s voxels   %s   %.1fs"
            % (
                "x".join(str(value) for value in result.dimensions),
                result.pixel_type,
                result.seconds,
            )
        )
        if json_file is not None:
            _success("wrote %s" % json_file)


@app.command()
def extract(
    ctx: typer.Context,
    input_path: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        help="Volume file or directory tree containing supported volumes.",
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output .stl/.ply/.obj file."
    ),
    volume_id: int | None = typer.Option(
        None,
        "--volume",
        min=1,
        help="Integer volume ID displayed by 'medsurface list'.",
        show_default="automatic",
    ),
    preset: PresetChoice = typer.Option(PresetChoice.BONE, "--preset"),
    threshold: str | None = typer.Option(
        None,
        "--threshold",
        help="Intensity (HU for CT) or 'auto' for Otsu.",
        show_default="preset",
    ),
    median_mm: float | None = typer.Option(
        None, "--median-mm", help="Despeckle kernel extent, mm.",
        show_default="preset",
    ),
    closing_mm: float | None = typer.Option(
        None, "--closing-mm", help="Pore-sealing kernel extent, mm.",
        show_default="preset",
    ),
    opening_mm: float | None = typer.Option(
        None, "--opening-mm", help="Bridge-breaking kernel extent, mm.",
        show_default="preset",
    ),
    min_island_mm3: float | None = typer.Option(
        None, "--min-island-mm3", help="Drop blobs smaller than this.",
        show_default="preset",
    ),
    all_islands: bool = typer.Option(
        False, "--all-islands", help="Keep every labelmap island."
    ),
    components: ComponentChoice | None = typer.Option(
        None,
        "--components",
        help="Surface components to keep.",
        show_default="preset",
    ),
    resample_mm: float | None = typer.Option(
        None,
        "--resample-mm",
        help="Isotropic surface-grid voxel size in mm (0 = native).",
        show_default="preset",
    ),
    mask_smooth_mm: float | None = typer.Option(
        None,
        "--mask-smooth-mm",
        help="Gaussian sigma in physical mm applied to the segmented mask before meshing (0 = off).",
        show_default="preset",
    ),
    surface_smooth_iters: int | None = typer.Option(
        None,
        "--mesh-smooth-iters",
        help="Topology-preserving surface relaxation iterations after meshing (0 = off).",
        show_default="preset",
    ),
    simplify_error_mm: float | None = typer.Option(
        None,
        "--simplify-error-mm",
        help="MeshLib estimated surface-deviation/QEM limit in model mm, not a certified Hausdorff bound (0 = off).",
        show_default="preset",
    ),
    post_surface_smooth_iters: int | None = typer.Option(
        None,
        "--post-mesh-smooth-iters",
        help="Final topology-preserving surface relaxation iterations after simplification (0 = off).",
        show_default="preset",
    ),
    destep: DestepChoice | None = _DESTEP_OPTION,
    destep_axis: AxisChoice = _DESTEP_AXIS_OPTION,
    destep_full_mm: float | None = _DESTEP_FULL_OPTION,
    destep_frozen_mm: float | None = _DESTEP_FROZEN_OPTION,
    destep_iters: int = _DESTEP_ITERS_OPTION,
    destep_max_mm: float = _DESTEP_MAX_OPTION,
    no_cap: bool = typer.Option(
        False, "--no-cap", help="Do not close anatomy at the field-of-view boundary."
    ),
    allow_large_volume: bool = typer.Option(
        False,
        "--allow-large-volume",
        help="Bypass the 500-million-voxel source and processing-grid limits; may exhaust memory.",
    ),
    json_file: Path | None = typer.Option(
        None, "--json", help="Write results and provenance to this JSON file."
    ),
    progress_format: ProgressChoice = _PROGRESS_OPTION,
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Segment an intensity volume and extract its surface mesh."""
    _validate_mesh_output(output)
    threshold_value = _parse_threshold(threshold)
    _validate_processing_numbers(
        nonnegative=[
            ("--median-mm", median_mm),
            ("--closing-mm", closing_mm),
            ("--opening-mm", opening_mm),
            ("--min-island-mm3", min_island_mm3),
            ("--resample-mm", resample_mm),
            ("--mask-smooth-mm", mask_smooth_mm),
            ("--mesh-smooth-iters", surface_smooth_iters),
            ("--simplify-error-mm", simplify_error_mm),
            ("--post-mesh-smooth-iters", post_surface_smooth_iters),
        ],
    )
    destep_settings = _destep_settings(
        ctx,
        destep,
        destep_axis,
        destep_full_mm,
        destep_frozen_mm,
        destep_iters,
        destep_max_mm,
    )
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _progress_display(progress_format, quiet, "Discovering volumes ...")
    with progress:
        found = _discover(input_path)
        if not found:
            _error("no supported volumes found under %s" % input_path)
            raise typer.Exit(1)

        from . import catalog

        try:
            chosen = catalog.select(found, volume_id)
            _protect_output_paths(found, output, json_file)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        progress.log(
            "volume ID %d  %s  %s" % (chosen.id, chosen.format, chosen.source_name)
        )
        resolved_preset = presets_mod.get(preset.value)
        resolved_preset = presets_mod.override(
            resolved_preset,
            median_mm=median_mm,
            closing_mm=closing_mm,
            opening_mm=opening_mm,
            min_island_mm3=min_island_mm3,
            resample_mm=resample_mm,
            mask_smooth_mm=mask_smooth_mm,
            surface_smooth_iters=surface_smooth_iters,
            simplify_error_mm=simplify_error_mm,
            post_surface_smooth_iters=post_surface_smooth_iters,
            keep_largest_island=False if all_islands else None,
            keep_largest_component=_keep_largest_component(components),
            destep=destep_settings,
        )
        progress.update("Loading surface extraction engine ...")
        from . import pipeline

        try:
            result = pipeline.extract(
                candidate=chosen,
                preset=resolved_preset,
                output_path=str(output),
                threshold=threshold_value,
                cap_field_of_view=not no_cap,
                allow_large_volume=allow_large_volume,
                log=progress.log,
                warn=emit_warning,
                progress=_progress_sink(progress),
            )
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(1) from None

        report = result.quality

        if json_file is not None:
            payload = {
                "result": {
                    "output": result.output_path,
                    "triangles": result.triangles,
                    "vertices": result.vertices,
                    "bounds_mm": list(result.bounds_mm),
                    "seconds": result.seconds,
                    "capped_field_of_view": result.capped_field_of_view,
                    "labelmap_components": result.labelmap_components,
                    "surface_components": result.surface_components,
                    "warnings": result.warnings,
                },
                "provenance": result.provenance,
                "quality": report,
            }
            progress.update("Writing JSON report ...")
            try:
                _write_json_file(json_file, payload)
            except OSError as exc:
                _error("cannot write JSON report %s: %s" % (json_file, exc))
                raise typer.Exit(1) from None

    _emit_remaining_warnings(result.warnings, emitted_warnings)

    if not quiet:
        _success("wrote %s" % result.output_path)
        _log(
            "triangles %s   vertices %s   %.1fs"
            % (f"{result.triangles:,}", f"{result.vertices:,}", result.seconds)
        )
        _print_quality(report)
        if json_file is not None:
            _success("wrote %s" % json_file)

    _warn_if_invalid(report, "output")

    _exit_for_quality(report)


def _parse_label_selection(value: str | None) -> list[int] | None:
    """Parse ``--labels`` such as ``3,7,12`` or ``10-14,20``."""
    if value is None:
        return None
    labels: set[int] = set()
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        try:
            if "-" in part:
                low_text, high_text = part.split("-", 1)
                low, high = int(low_text), int(high_text)
                if high < low:
                    raise ValueError
                labels.update(range(low, high + 1))
            else:
                labels.add(int(part))
        except ValueError:
            _error("--labels expects positive label IDs such as 3,7,12 or 10-14")
            raise typer.Exit(2) from None
    if not labels or min(labels) <= 0:
        _error("--labels expects positive label IDs such as 3,7,12 or 10-14")
        raise typer.Exit(2)
    return sorted(labels)


def _load_label_names(path: Path | None) -> dict[int, str] | None:
    if path is None:
        return None
    from .labelnames import load_names_json

    try:
        return load_names_json(path)
    except (OSError, ValueError) as exc:
        _error("cannot read --label-names %s: %s" % (path, exc))
        raise typer.Exit(2) from None


@labelmap_app.command("extract")
def extract_labelmap(
    ctx: typer.Context,
    input_path: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Direct NIfTI, NRRD, or MetaImage labelmap file.",
    ),
    output: Path | None = typer.Option(
        None,
        "-o",
        "--output",
        help="Output .stl/.ply/.obj file. With --split it is optional and receives "
        "the combined model of every selected label.",
    ),
    labels: str | None = typer.Option(
        None,
        "--labels",
        help="Only use these label IDs, e.g. 3,7,12 or 10-14.",
        show_default="every nonzero label",
    ),
    split: Path | None = typer.Option(
        None,
        "--split",
        file_okay=False,
        dir_okay=True,
        help="Write one mesh per label into this directory as NNN_name.stl.",
    ),
    label_names: Path | None = typer.Option(
        None,
        "--label-names",
        exists=True,
        dir_okay=False,
        readable=True,
        help='JSON object {"5": "liver"} naming --split files.',
        show_default="label table embedded in a NIfTI input",
    ),
    resample_mm: float = typer.Option(
        0.0,
        "--resample-mm",
        help="Isotropic surface-grid voxel size in mm (0 = native).",
    ),
    mask_smooth_mm: float = typer.Option(
        defaults.DEFAULT_LABELMAP_MASK_SMOOTH_MM,
        "--mask-smooth-mm",
        help="Gaussian sigma in physical mm applied to the labelmap before meshing (0 = off).",
    ),
    surface_smooth_iters: int = typer.Option(
        defaults.DEFAULT_LABELMAP_SURFACE_SMOOTH_ITERS,
        "--mesh-smooth-iters",
        help="Topology-preserving surface relaxation iterations after meshing (0 = off).",
    ),
    simplify_error_mm: float = typer.Option(
        defaults.DEFAULT_LABELMAP_SIMPLIFY_ERROR_MM,
        "--simplify-error-mm",
        help="MeshLib estimated surface-deviation/QEM limit in model mm (0 = off).",
    ),
    post_surface_smooth_iters: int = typer.Option(
        defaults.DEFAULT_LABELMAP_POST_SURFACE_SMOOTH_ITERS,
        "--post-mesh-smooth-iters",
        help="Final topology-preserving surface relaxation iterations after simplification (0 = off).",
    ),
    destep: DestepChoice | None = _DESTEP_OPTION,
    destep_axis: AxisChoice = _DESTEP_AXIS_OPTION,
    destep_full_mm: float | None = _DESTEP_FULL_OPTION,
    destep_frozen_mm: float | None = _DESTEP_FROZEN_OPTION,
    destep_iters: int = _DESTEP_ITERS_OPTION,
    destep_max_mm: float = _DESTEP_MAX_OPTION,
    components: ComponentChoice = typer.Option(
        ComponentChoice.ALL,
        "--components",
        help="Surface components to keep.",
    ),
    no_cap: bool = typer.Option(
        False, "--no-cap", help="Do not close foreground at the labelmap boundary."
    ),
    allow_large_volume: bool = typer.Option(
        False,
        "--allow-large-volume",
        help="Bypass the 500-million-voxel limits; may exhaust memory.",
    ),
    json_file: Path | None = typer.Option(
        None, "--json", help="Write results and provenance to this JSON file."
    ),
    progress_format: ProgressChoice = _PROGRESS_OPTION,
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Create one mesh from every nonzero (or selected) voxel, or one mesh per label."""
    if output is None and split is None:
        _error("pass -o/--output, --split, or both")
        raise typer.Exit(2)
    if output is not None:
        _validate_mesh_output(output)
    selection = _parse_label_selection(labels)
    names = _load_label_names(label_names)
    if names is not None and split is None:
        _error("--label-names requires --split")
        raise typer.Exit(2)
    _validate_processing_numbers(
        nonnegative=[
            ("--resample-mm", resample_mm),
            ("--mask-smooth-mm", mask_smooth_mm),
            ("--mesh-smooth-iters", surface_smooth_iters),
            ("--simplify-error-mm", simplify_error_mm),
            ("--post-mesh-smooth-iters", post_surface_smooth_iters),
        ],
    )
    destep_settings = _destep_settings(
        ctx,
        destep,
        destep_axis,
        destep_full_mm,
        destep_frozen_mm,
        destep_iters,
        destep_max_mm,
    )
    if split is not None:
        _extract_labelmap_split(
            input_path=input_path,
            split=split,
            output=output,
            selection=selection,
            names=names,
            resample_mm=resample_mm,
            mask_smooth_mm=mask_smooth_mm,
            surface_smooth_iters=surface_smooth_iters,
            simplify_error_mm=simplify_error_mm,
            post_surface_smooth_iters=post_surface_smooth_iters,
            keep_largest_component=_keep_largest_component(components),
            destep_settings=destep_settings,
            no_cap=no_cap,
            allow_large_volume=allow_large_volume,
            json_file=json_file,
            progress_format=progress_format,
            quiet=quiet,
        )
        return
    assert output is not None
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _progress_display(progress_format, quiet, "Discovering labelmap ...")
    with progress:
        chosen, found = _select_labelmap(input_path)
        try:
            _protect_output_paths(found, output, json_file)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        progress.update("Loading labelmap extraction engine ...")
        from . import labelmap as labelmap_mod

        try:
            result = labelmap_mod.extract(
                candidate=chosen,
                output_path=str(output),
                resample_mm=resample_mm,
                mask_smooth_mm=mask_smooth_mm,
                surface_smooth_iters=surface_smooth_iters,
                simplify_error_mm=simplify_error_mm,
                post_surface_smooth_iters=post_surface_smooth_iters,
                keep_largest_component=_keep_largest_component(components),
                destep=destep_settings,
                cap_field_of_view=not no_cap,
                allow_large_volume=allow_large_volume,
                log=progress.log,
                warn=emit_warning,
                labels=selection,
                progress=_progress_sink(progress),
            )
        except (OSError, RuntimeError, ValueError) as exc:
            _error(exc)
            raise typer.Exit(1) from None

        report = result.quality
        if json_file is not None:
            payload = {
                "result": {
                    "output": result.output_path,
                    "triangles": result.triangles,
                    "vertices": result.vertices,
                    "bounds_mm": list(result.bounds_mm),
                    "seconds": result.seconds,
                    "capped_field_of_view": result.capped_field_of_view,
                    "labelmap_components": result.labelmap_components,
                    "surface_components": result.surface_components,
                    "warnings": result.warnings,
                },
                "provenance": result.provenance,
                "quality": report,
            }
            progress.update("Writing JSON report ...")
            try:
                _write_json_file(json_file, payload)
            except OSError as exc:
                _error("cannot write JSON report %s: %s" % (json_file, exc))
                raise typer.Exit(1) from None

    _emit_remaining_warnings(result.warnings, emitted_warnings)
    if not quiet:
        _success("wrote %s" % result.output_path)
        _log(
            "triangles %s   vertices %s   %.1fs"
            % (f"{result.triangles:,}", f"{result.vertices:,}", result.seconds)
        )
        _print_quality(report)
        if json_file is not None:
            _success("wrote %s" % json_file)
    _warn_if_invalid(report, "output")
    _exit_for_quality(report)


def _extract_labelmap_split(
    *,
    input_path: Path,
    split: Path,
    output: Path | None,
    selection: list[int] | None,
    names: dict[int, str] | None,
    resample_mm: float | None,
    mask_smooth_mm: float,
    surface_smooth_iters: int,
    simplify_error_mm: float | None,
    post_surface_smooth_iters: int,
    keep_largest_component: bool | None,
    destep_settings: presets_mod.DestepSettings | None,
    no_cap: bool,
    allow_large_volume: bool,
    json_file: Path | None,
    progress_format: ProgressChoice,
    quiet: bool,
) -> None:
    if split.exists() and not split.is_dir():
        _error("--split must name a directory: %s" % split)
        raise typer.Exit(2)
    mesh_format = extension(output).lstrip(".") if output is not None else "stl"
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _progress_display(progress_format, quiet, "Discovering labelmap ...")
    with progress:
        chosen, found = _select_labelmap(input_path)
        try:
            if output is not None:
                _protect_output_paths(found, output, json_file)
            elif json_file is not None:
                _protect_output_paths(found, json_file, None)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        progress.update("Loading labelmap extraction engine ...")
        from . import labelmap as labelmap_mod

        try:
            settings = labelmap_mod.resolve_surface_settings(
                resample_mm=resample_mm,
                mask_smooth_mm=mask_smooth_mm,
                surface_smooth_iters=surface_smooth_iters,
                simplify_error_mm=simplify_error_mm,
                post_surface_smooth_iters=post_surface_smooth_iters,
                keep_largest_component=keep_largest_component,
                destep=destep_settings,
            )
            result = labelmap_mod.extract_labels(
                chosen,
                str(split),
                labels=selection,
                names=names,
                settings=settings,
                combined_path=None if output is None else str(output),
                mesh_format=mesh_format,
                cap_field_of_view=not no_cap,
                allow_large_volume=allow_large_volume,
                log=progress.log,
                warn=emit_warning,
                progress=_progress_sink(progress),
            )
        except (OSError, RuntimeError, ValueError) as exc:
            _error(exc)
            raise typer.Exit(1) from None

        if json_file is not None:
            payload = {
                "result": {
                    "output_dir": result.output_dir,
                    "meshes": [mesh.payload() for mesh in result.meshes],
                    "combined": (
                        None if result.combined is None else result.combined.payload()
                    ),
                    "skipped": result.skipped,
                    "failed": result.failed,
                    "seconds": result.seconds,
                    "warnings": result.warnings,
                },
                "provenance": result.provenance,
            }
            progress.update("Writing JSON report ...")
            try:
                _write_json_file(json_file, payload)
            except OSError as exc:
                _error("cannot write JSON report %s: %s" % (json_file, exc))
                raise typer.Exit(1) from None

    _emit_remaining_warnings(result.warnings, emitted_warnings)
    meshes = [*result.meshes, *([result.combined] if result.combined else [])]
    if not quiet:
        for mesh in meshes:
            label = "combined" if mesh.label is None else "label %d" % mesh.label
            if mesh.name:
                label += " (%s)" % mesh.name
            _success("wrote %s   %s" % (mesh.output_path, label))
        _log(
            "%d mesh(es)   triangles %s   %.1fs"
            % (
                len(meshes),
                f"{sum(mesh.triangles for mesh in meshes):,}",
                result.seconds,
            )
        )
        if json_file is not None:
            _success("wrote %s" % json_file)
    for mesh in meshes:
        _warn_if_invalid(mesh.quality, mesh.output_path)
    for failure in result.failed:
        _error("label %s failed: %s" % (failure["label"], failure["error"]))
    if result.failed or not all(mesh.quality["valid"] for mesh in meshes):
        raise typer.Exit(1)


def _labelmap_fusion_payload(result: Any) -> dict[str, Any]:
    return {
        "output": result.output_path,
        "format": result.output_format,
        "compression": result.compression,
        "pixel_type": result.pixel_type,
        "components": 1,
        "labels_preserved": result.preserve_labels,
        "grid_mm": result.grid_mm,
        "grid_size": list(result.grid_size),
        "grid_origin_mm": list(result.grid_origin_mm),
        "grid_direction": list(result.grid_direction),
        "labels": {str(label): item for label, item in result.labels.items()},
        "volume_fused_mm3": result.volume_fused_mm3,
        "seconds": result.seconds,
        "warnings": result.warnings,
    }


@labelmap_app.command("fuse")
def fuse_labelmaps(
    fixed_input: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Fixed labelmap; defines the output coordinate frame.",
    ),
    moving_inputs: list[Path] = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        metavar="MOVING...",
        help="One or more moving labelmaps to register to the fixed labelmap.",
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Atomic .nii/.nii.gz/.nrrd/.mha labelmap (binary unless --preserve-labels).",
    ),
    preserve_labels: bool = typer.Option(
        False,
        "--preserve-labels",
        help="Keep label IDs instead of writing a binary union; touching labels "
        "are split by their fused occupancy. NIfTI outputs keep label names.",
    ),
    labels: str | None = typer.Option(
        None,
        "--labels",
        help="Only fuse these label IDs, e.g. 3,7,12 or 10-14.",
        show_default="every nonzero label",
    ),
    label_names: Path | None = typer.Option(
        None,
        "--label-names",
        exists=True,
        dir_okay=False,
        readable=True,
        help='JSON object {"5": "liver"} of names embedded in a NIfTI output '
        "with --preserve-labels.",
        show_default="label tables embedded in the inputs",
    ),
    grid_mm: float = typer.Option(
        defaults.DEFAULT_FUSION_GRID_MM,
        "--grid-mm",
        help="Isotropic fused-grid voxel size in mm.",
    ),
    force: bool = typer.Option(
        False, "--force", help="Override registration-quality gates."
    ),
    allow_large_volume: bool = typer.Option(
        False,
        "--allow-large-volume",
        help="Bypass the 500-million-voxel limits; may exhaust memory.",
    ),
    json_file: Path | None = typer.Option(
        None, "--json", help="Write results and provenance to this JSON file."
    ),
    progress_format: ProgressChoice = _PROGRESS_OPTION,
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Register and fuse two or more labelmaps into one binary or labelled labelmap."""
    _validate_volume_output(output)
    _validate_processing_numbers(
        nonnegative=[],
        positive=(("--grid-mm", grid_mm),),
    )
    selection = _parse_label_selection(labels)
    names = _load_label_names(label_names)
    if names is not None and not preserve_labels:
        _error("--label-names requires --preserve-labels")
        raise typer.Exit(2)
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _progress_display(progress_format, quiet, "Discovering fixed labelmap ...")
    with progress:
        fixed, fixed_found = _select_labelmap(fixed_input, "fixed")
        movings = []
        found = list(fixed_found)
        for index, moving_input in enumerate(moving_inputs, start=1):
            role = "moving" if len(moving_inputs) == 1 else "moving %d" % index
            progress.update("Discovering %s labelmap ..." % role)
            moving, moving_found = _select_labelmap(moving_input, role)
            movings.append(moving)
            found.extend(moving_found)
        try:
            _protect_output_paths(found, output, json_file)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        progress.update("Loading labelmap fusion engine ...")
        from . import fusion as fusion_mod
        from . import labelmap as labelmap_mod

        transaction = (
            _restore_output_on_failure(output)
            if json_file is not None
            else nullcontext()
        )
        try:
            with transaction:
                result = labelmap_mod.fuse_labels(
                    fixed,
                    movings,
                    str(output),
                    grid_mm=grid_mm,
                    labels=selection,
                    preserve_labels=preserve_labels,
                    names=names,
                    force=force,
                    allow_large_volume=allow_large_volume,
                    log=progress.log,
                    warn=emit_warning,
                    progress=_progress_sink(progress),
                )
                if json_file is not None:
                    payload = {
                        "result": _labelmap_fusion_payload(result),
                        "provenance": result.provenance,
                    }
                    progress.update("Writing JSON report ...")
                    try:
                        _write_json_file(json_file, payload)
                    except OSError as exc:
                        _error("cannot write JSON report %s: %s" % (json_file, exc))
                        raise typer.Exit(1) from None
        except typer.Exit:
            raise
        except fusion_mod.FusionError as exc:
            _error(exc)
            raise typer.Exit(3) from None
        except (OSError, RuntimeError, ValueError) as exc:
            _error(exc)
            raise typer.Exit(1) from None

    _emit_remaining_warnings(result.warnings, emitted_warnings)
    if not quiet:
        _success("wrote %s" % result.output_path)
        _log(
            "%s voxels at %.3f mm   %d label(s)   %.0f cm3   %.1fs"
            % (
                "x".join(str(value) for value in result.grid_size),
                result.grid_mm,
                len(result.labels),
                result.volume_fused_mm3 / 1000.0,
                result.seconds,
            )
        )
        if preserve_labels:
            _print_command_hint(
                "Extract surfaces with:  ",
                [
                    "medsurface", "labelmap", "extract", result.output_path,
                    "--split", "MODELS", "-o", "COMBINED.stl",
                ],
            )
        else:
            _print_command_hint(
                "Extract a surface with:  ",
                ["medsurface", "labelmap", "extract", result.output_path, "-o", "MODEL.stl"],
            )
        if json_file is not None:
            _success("wrote %s" % json_file)


@app.command()
def fuse(
    fixed_input: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        help="Fixed volume file or directory; defines the output coordinate frame.",
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Atomic .nii/.nii.gz/.nrrd/.mha binary labelmap.",
    ),
    moving_input: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        help="Moving volume file or directory.",
    ),
    fixed_volume: int | None = typer.Option(
        None,
        "--fixed-volume",
        min=1,
        help="Integer fixed-volume ID displayed by 'medsurface list'.",
        show_default="automatic",
    ),
    moving_volume: int | None = typer.Option(
        None,
        "--moving-volume",
        min=1,
        help="Integer moving-volume ID displayed by 'medsurface list'.",
        show_default="automatic",
    ),
    preset: PresetChoice = typer.Option(PresetChoice.BONE, "--preset"),
    fixed_threshold: str | None = typer.Option(
        None,
        "--fixed-threshold",
        help="Fixed stored intensity or 'auto' for Otsu.",
        show_default="preset",
    ),
    moving_threshold: str | None = typer.Option(
        None,
        "--moving-threshold",
        help="Moving stored intensity or 'auto' for Otsu.",
        show_default="preset",
    ),
    median_mm: float | None = typer.Option(
        None, "--median-mm", help="Despeckle kernel extent, mm.",
        show_default="preset",
    ),
    closing_mm: float | None = typer.Option(
        None, "--closing-mm", help="Pore-sealing kernel extent, mm.",
        show_default="preset",
    ),
    opening_mm: float | None = typer.Option(
        None, "--opening-mm", help="Bridge-breaking kernel extent, mm.",
        show_default="preset",
    ),
    min_island_mm3: float | None = typer.Option(
        None, "--min-island-mm3", help="Drop blobs smaller than this.",
        show_default="preset",
    ),
    all_islands: bool = typer.Option(
        False, "--all-islands", help="Keep every labelmap island."
    ),
    grid_mm: float = typer.Option(
        defaults.DEFAULT_FUSION_GRID_MM,
        "--grid-mm",
        help="Isotropic fused-grid voxel size in mm.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Override registration-quality gates.",
    ),
    allow_large_volume: bool = typer.Option(
        False,
        "--allow-large-volume",
        help="Bypass the 500-million-voxel source and processing-grid limits; may exhaust memory.",
    ),
    json_file: Path | None = typer.Option(
        None, "--json", help="Write results and provenance to this JSON file."
    ),
    progress_format: ProgressChoice = _PROGRESS_OPTION,
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Segment, register, and union two volumes into a binary labelmap."""
    _validate_volume_output(output)
    fixed_threshold_value = _parse_threshold(fixed_threshold, "--fixed-threshold")
    moving_threshold_value = _parse_threshold(moving_threshold, "--moving-threshold")
    _validate_processing_numbers(
        nonnegative=[
            ("--median-mm", median_mm),
            ("--closing-mm", closing_mm),
            ("--opening-mm", opening_mm),
            ("--min-island-mm3", min_island_mm3),
        ],
        positive=(("--grid-mm", grid_mm),),
    )
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _progress_display(progress_format, quiet, "Discovering fixed volumes ...")
    with progress:
        fixed_found = _discover(fixed_input)
        if fixed_input.resolve() == moving_input.resolve():
            moving_found = fixed_found
        else:
            progress.update("Discovering moving volumes ...")
            moving_found = _discover(moving_input)
        if not fixed_found or not moving_found:
            _error("no supported volumes found for one or both inputs")
            raise typer.Exit(1)

        from . import catalog

        try:
            fixed_candidate = catalog.select(fixed_found, fixed_volume)
        except ValueError as exc:
            _error("fixed input: %s" % exc)
            raise typer.Exit(2) from None
        try:
            moving_candidate = catalog.select(moving_found, moving_volume)
        except ValueError as exc:
            _error("moving input: %s" % exc)
            raise typer.Exit(2) from None
        try:
            _protect_output_paths([*fixed_found, *moving_found], output, json_file)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        resolved_preset = presets_mod.get(preset.value)
        resolved_preset = presets_mod.override(
            resolved_preset,
            median_mm=median_mm,
            closing_mm=closing_mm,
            opening_mm=opening_mm,
            min_island_mm3=min_island_mm3,
            keep_largest_island=False if all_islands else None,
        )
        progress.update("Loading fusion engine ...")
        from . import fusion as fusion_mod

        transaction = (
            _restore_output_on_failure(output)
            if json_file is not None
            else nullcontext()
        )
        try:
            with transaction:
                result = fusion_mod.fuse(
                    fixed=fixed_candidate,
                    moving=moving_candidate,
                    preset=resolved_preset,
                    output_path=str(output),
                    fixed_threshold=fixed_threshold_value,
                    moving_threshold=moving_threshold_value,
                    grid_mm=grid_mm,
                    force=force,
                    allow_large_volume=allow_large_volume,
                    log=progress.log,
                    warn=emit_warning,
                    progress=_progress_sink(progress),
                )
                if json_file is not None:
                    payload = {
                        "result": _fusion_result_payload(result),
                        "provenance": result.provenance,
                    }
                    progress.update("Writing JSON report ...")
                    try:
                        _write_json_file(json_file, payload)
                    except OSError as exc:
                        _error("cannot write JSON report %s: %s" % (json_file, exc))
                        raise typer.Exit(1) from None
        except typer.Exit:
            raise
        except fusion_mod.FusionError as exc:
            _error(exc)
            raise typer.Exit(3) from None
        except (OSError, RuntimeError, ValueError) as exc:
            _error(exc)
            raise typer.Exit(2) from None

    _emit_remaining_warnings(result.warnings, emitted_warnings)

    if not quiet:
        _print_fusion_result(result)
        if json_file is not None:
            _success("wrote %s" % json_file)


@app.command()
def validate(
    mesh: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Mesh file to inspect.",
    ),
    json_file: Path | None = typer.Option(
        None, "--json", help="Write the quality report to this JSON file."
    ),
    progress_format: ProgressChoice = _PROGRESS_OPTION,
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress the quality report and progress output."
    ),
) -> None:
    """Report mesh quality without changing the file."""
    _validate_mesh_input(mesh)
    if json_file is not None:
        from .paths import same_file

        if same_file(json_file, mesh):
            _error("JSON report must not overwrite input file %s" % mesh)
            raise typer.Exit(2)
    progress = _progress_display(progress_format, quiet, "Loading validation engine ...")
    with progress:
        from . import validate as validate_mod

        progress.update("Validating mesh structure and self-intersections ...")
        try:
            report = validate_mod.validate(str(mesh))
        except (OSError, RuntimeError, ValueError) as exc:
            _error("cannot validate %s: %s" % (mesh, exc))
            raise typer.Exit(1) from None
        if json_file is not None:
            progress.update("Writing JSON report ...")
            try:
                _write_json_file(json_file, report)
            except OSError as exc:
                _error("cannot write JSON report %s: %s" % (json_file, exc))
                raise typer.Exit(1) from None
    if not quiet:
        stdout_console.print(_plain(mesh, "bold"))
        _print_quality(report)
        if json_file is not None:
            _success("wrote %s" % json_file)
    _exit_for_quality(report)


@app.command()
def repair(
    mesh: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Mesh file to repair.",
    ),
    output: Path = typer.Option(..., "-o", "--output", help="New repaired mesh file."),
    json_file: Path | None = typer.Option(
        None, "--json", help="Write repair statistics and quality to this JSON file."
    ),
    progress_format: ProgressChoice = _PROGRESS_OPTION,
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Make a non-watertight mesh watertight."""
    _validate_mesh_input(mesh)
    _validate_mesh_output(output)
    from .paths import protect_outputs

    try:
        protect_outputs([mesh], output, json_file)
    except ValueError as exc:
        _error(exc)
        raise typer.Exit(2) from None
    progress = _progress_display(progress_format, quiet, "Loading repair engine ...")
    with progress:
        from . import repair as repair_mod

        transaction = (
            _restore_output_on_failure(output)
            if json_file is not None
            else nullcontext()
        )
        try:
            with transaction:
                result = repair_mod.repair(
                    str(mesh),
                    str(output),
                    log=progress.log,
                    progress=_progress_sink(progress),
                )
                if json_file is not None:
                    progress.update("Writing JSON report ...")
                    try:
                        _write_json_file(
                            json_file,
                            {
                                "output": str(output),
                                "repair": result.stats,
                                "quality": result.quality,
                            },
                        )
                    except OSError as exc:
                        _error("cannot write JSON report %s: %s" % (json_file, exc))
                        raise typer.Exit(1) from None
        except typer.Exit:
            raise
        except (OSError, RuntimeError, ValueError) as exc:
            _error("cannot repair %s: %s" % (mesh, exc))
            raise typer.Exit(1) from None

    report = result.quality
    if not quiet:
        _success("wrote %s" % output)
        _print_quality(report)
        if json_file is not None:
            _success("wrote %s" % json_file)
