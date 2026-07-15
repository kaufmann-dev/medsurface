"""Typer command line interface with Rich human-readable output."""

from __future__ import annotations

import json
import math
import os
import shlex
import subprocess
import sys
import tempfile
import warnings
from collections import Counter
from contextvars import ContextVar, Token
from enum import Enum
from pathlib import Path
from typing import Any

import typer
from rich import box
from rich.console import Console, Group
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TaskID, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich.text import Text

from . import __version__, defaults
from . import presets as presets_mod
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


def _keep_largest_component(
    components: ComponentChoice | None,
) -> bool | None:
    if components is None:
        return None
    return components is ComponentChoice.LARGEST


stdout_console = Console(highlight=False, markup=False)
stderr_console = Console(stderr=True, highlight=False, markup=False)

app = typer.Typer(
    add_completion=False,
    help="Turn a medical image volume into a watertight 3D surface mesh.",
    invoke_without_command=True,
    no_args_is_help=False,
    pretty_exceptions_show_locals=False,
    rich_markup_mode="rich",
)
labelmap_app = typer.Typer(
    add_completion=False,
    help="Create surfaces from externally segmented labelmaps.",
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
    """Turn a medical image volume into a watertight 3D surface mesh."""
    if ctx.invoked_subcommand is None:
        stdout_console.print(ctx.get_help())


@labelmap_app.callback()
def labelmap_root(ctx: typer.Context) -> None:
    """Create surfaces from externally segmented labelmaps."""
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
            if self._context_token is not None:
                _active_progress.reset(self._context_token)
                self._context_token = None

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
    extension = path.suffix.casefold()
    if extension not in defaults.SUPPORTED_MESH_EXTENSIONS:
        _error(
            "unsupported output extension %r; supported: %s"
            % (extension, ", ".join(defaults.SUPPORTED_MESH_EXTENSIONS))
        )
        raise typer.Exit(2)


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
    json_output: bool = typer.Option(
        False,
        "--json",
        help="Write one plain JSON array to stdout.",
    ),
) -> None:
    """Show every supported volume under an input path."""
    with _ProgressDisplay(not json_output, "Discovering volumes ..."):
        found = _discover(input_path)
    if not found:
        _error("no supported volumes found under %s" % input_path)
        raise typer.Exit(1)

    from . import catalog

    recommended = catalog.recommended(found)
    if json_output:
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
        print(json.dumps(payload, indent=2))
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
            "Convert the default with:  ",
            ["medsurface", "convert", input_path, "-o", "out.stl"],
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
                "convert",
                input_path,
                "--volume",
                "ID",
                "-o",
                "out.stl",
            ],
        )


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
        ..., "-o", "--output", help="Output .stl/.ply/.obj file."
    ),
    volume_id: int | None = typer.Option(
        None,
        "--volume",
        min=1,
        help="Integer volume ID displayed by 'medsurface list'.",
    ),
    preset: PresetChoice = typer.Option(PresetChoice.BONE, "--preset"),
    threshold: str | None = typer.Option(
        None,
        "--threshold",
        help="Intensity (HU for CT) or 'auto' for Otsu.",
    ),
    median_mm: float | None = typer.Option(
        None, "--median-mm", help="Despeckle kernel extent, mm."
    ),
    closing_mm: float | None = typer.Option(
        None, "--closing-mm", help="Pore-sealing kernel extent, mm."
    ),
    opening_mm: float | None = typer.Option(
        None, "--opening-mm", help="Bridge-breaking kernel extent, mm."
    ),
    min_island_mm3: float | None = typer.Option(
        None, "--min-island-mm3", help="Drop blobs smaller than this."
    ),
    all_islands: bool = typer.Option(
        False, "--all-islands", help="Keep every labelmap island."
    ),
    components: ComponentChoice | None = typer.Option(
        None,
        "--components",
        help="Surface components to keep (default: preset).",
    ),
    resample_mm: float | None = typer.Option(
        None,
        "--resample-mm",
        help="Isotropic surface-grid voxel size in mm (0 = native).",
    ),
    mask_smooth_mm: float | None = typer.Option(
        None,
        "--mask-smooth-mm",
        help="Gaussian sigma in physical mm applied to the segmented mask before meshing (0 = off).",
    ),
    surface_smooth_iters: int | None = typer.Option(
        None,
        "--mesh-smooth-iters",
        help="Topology-preserving surface relaxation iterations after meshing (0 = off).",
    ),
    simplify_error_mm: float | None = typer.Option(
        None,
        "--simplify-error-mm",
        help="MeshLib estimated surface-deviation/QEM limit in model mm, not a certified Hausdorff bound (0 = off).",
    ),
    post_surface_smooth_iters: int | None = typer.Option(
        None,
        "--post-mesh-smooth-iters",
        help="Final topology-preserving surface relaxation iterations after simplification (0 = off).",
    ),
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
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Extract a surface mesh from a medical image volume."""
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
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _ProgressDisplay(not quiet, "Discovering volumes ...")
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
        )
        progress.update("Loading conversion engine ...")
        from . import pipeline

        try:
            result = pipeline.convert(
                candidate=chosen,
                preset=resolved_preset,
                output_path=str(output),
                threshold=threshold_value,
                cap_field_of_view=not no_cap,
                allow_large_volume=allow_large_volume,
                log=progress.log,
                warn=emit_warning,
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


@labelmap_app.command("convert")
def convert_labelmap(
    input_path: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Direct NIfTI, NRRD, or MetaImage labelmap file.",
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output .stl/.ply/.obj file."
    ),
    resample_mm: float | None = typer.Option(
        None,
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
    simplify_error_mm: float | None = typer.Option(
        None,
        "--simplify-error-mm",
        help="MeshLib estimated surface-deviation/QEM limit in model mm (0 = off).",
    ),
    post_surface_smooth_iters: int = typer.Option(
        defaults.DEFAULT_LABELMAP_POST_SURFACE_SMOOTH_ITERS,
        "--post-mesh-smooth-iters",
        help="Final topology-preserving surface relaxation iterations after simplification (0 = off).",
    ),
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
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Create one mesh from every nonzero voxel in a labelmap."""
    _validate_mesh_output(output)
    _validate_processing_numbers(
        nonnegative=[
            ("--resample-mm", resample_mm),
            ("--mask-smooth-mm", mask_smooth_mm),
            ("--mesh-smooth-iters", surface_smooth_iters),
            ("--simplify-error-mm", simplify_error_mm),
            ("--post-mesh-smooth-iters", post_surface_smooth_iters),
        ],
    )
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _ProgressDisplay(not quiet, "Discovering labelmap ...")
    with progress:
        chosen, found = _select_labelmap(input_path)
        try:
            _protect_output_paths(found, output, json_file)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        progress.update("Loading labelmap conversion engine ...")
        from . import labelmap as labelmap_mod

        try:
            result = labelmap_mod.convert(
                candidate=chosen,
                output_path=str(output),
                resample_mm=resample_mm,
                mask_smooth_mm=mask_smooth_mm,
                surface_smooth_iters=surface_smooth_iters,
                simplify_error_mm=simplify_error_mm,
                post_surface_smooth_iters=post_surface_smooth_iters,
                keep_largest_component=_keep_largest_component(components),
                cap_field_of_view=not no_cap,
                allow_large_volume=allow_large_volume,
                log=progress.log,
                warn=emit_warning,
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


@labelmap_app.command("merge")
def merge_labelmaps(
    fixed_input: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Fixed labelmap; defines the output coordinate frame.",
    ),
    moving_input: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Moving labelmap to register to the fixed labelmap.",
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output .stl/.ply/.obj file."
    ),
    grid_mm: float = typer.Option(
        defaults.DEFAULT_MERGE_GRID_MM,
        "--grid-mm",
        help="Isotropic fused-grid voxel size in mm.",
    ),
    mask_smooth_mm: float = typer.Option(
        defaults.DEFAULT_LABELMAP_MASK_SMOOTH_MM,
        "--mask-smooth-mm",
        help="Gaussian sigma in physical mm applied after labelmap fusion (0 = off).",
    ),
    surface_smooth_iters: int = typer.Option(
        defaults.DEFAULT_LABELMAP_SURFACE_SMOOTH_ITERS,
        "--mesh-smooth-iters",
        help="Topology-preserving surface relaxation iterations after meshing (0 = off).",
    ),
    simplify_error_mm: float | None = typer.Option(
        None,
        "--simplify-error-mm",
        help="MeshLib estimated surface-deviation/QEM limit in model mm (0 = off).",
    ),
    post_surface_smooth_iters: int = typer.Option(
        defaults.DEFAULT_LABELMAP_POST_SURFACE_SMOOTH_ITERS,
        "--post-mesh-smooth-iters",
        help="Final topology-preserving surface relaxation iterations after simplification (0 = off).",
    ),
    components: ComponentChoice = typer.Option(
        ComponentChoice.ALL,
        "--components",
        help="Surface components to keep.",
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
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Rigidly register two labelmaps and fuse their nonzero foreground."""
    _validate_mesh_output(output)
    _validate_processing_numbers(
        nonnegative=[
            ("--mask-smooth-mm", mask_smooth_mm),
            ("--mesh-smooth-iters", surface_smooth_iters),
            ("--simplify-error-mm", simplify_error_mm),
            ("--post-mesh-smooth-iters", post_surface_smooth_iters),
        ],
        positive=(("--grid-mm", grid_mm),),
    )
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _ProgressDisplay(not quiet, "Discovering fixed labelmap ...")
    with progress:
        fixed, fixed_found = _select_labelmap(fixed_input, "fixed")
        if fixed_input.resolve() == moving_input.resolve():
            moving, moving_found = fixed, fixed_found
        else:
            progress.update("Discovering moving labelmap ...")
            moving, moving_found = _select_labelmap(moving_input, "moving")
        try:
            _protect_output_paths([*fixed_found, *moving_found], output, json_file)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        progress.update("Loading labelmap merge engine ...")
        from . import labelmap as labelmap_mod
        from . import merge as merge_mod

        try:
            result = labelmap_mod.merge(
                fixed=fixed,
                moving=moving,
                output_path=str(output),
                grid_mm=grid_mm,
                mask_smooth_mm=mask_smooth_mm,
                surface_smooth_iters=surface_smooth_iters,
                simplify_error_mm=simplify_error_mm,
                post_surface_smooth_iters=post_surface_smooth_iters,
                keep_largest_component=_keep_largest_component(components),
                force=force,
                allow_large_volume=allow_large_volume,
                log=progress.log,
                warn=emit_warning,
            )
        except merge_mod.MergeError as exc:
            _error(exc)
            raise typer.Exit(3) from None
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
                    "grid_mm": result.grid_mm,
                    "grid_size": list(result.grid_size),
                    "volume_fixed_mm3": result.volume_fixed_mm3,
                    "volume_moving_mm3": result.volume_moving_mm3,
                    "volume_fused_mm3": result.volume_union_mm3,
                    "surface_components": result.surface_components,
                    "seconds": result.seconds,
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
    _warn_if_invalid(report, "fused output")
    _exit_for_quality(report)


@app.command()
def merge(
    fixed_input: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        help="Fixed volume file or directory; defines the output coordinate frame.",
    ),
    output: Path = typer.Option(
        ..., "-o", "--output", help="Output .stl/.ply/.obj file."
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
    ),
    moving_volume: int | None = typer.Option(
        None,
        "--moving-volume",
        min=1,
        help="Integer moving-volume ID displayed by 'medsurface list'.",
    ),
    preset: PresetChoice = typer.Option(PresetChoice.BONE, "--preset"),
    fixed_threshold: str | None = typer.Option(
        None,
        "--fixed-threshold",
        help="Fixed stored intensity or 'auto' for Otsu.",
    ),
    moving_threshold: str | None = typer.Option(
        None,
        "--moving-threshold",
        help="Moving stored intensity or 'auto' for Otsu.",
    ),
    median_mm: float | None = typer.Option(
        None, "--median-mm", help="Despeckle kernel extent, mm."
    ),
    closing_mm: float | None = typer.Option(
        None, "--closing-mm", help="Pore-sealing kernel extent, mm."
    ),
    opening_mm: float | None = typer.Option(
        None, "--opening-mm", help="Bridge-breaking kernel extent, mm."
    ),
    min_island_mm3: float | None = typer.Option(
        None, "--min-island-mm3", help="Drop blobs smaller than this."
    ),
    all_islands: bool = typer.Option(
        False, "--all-islands", help="Keep every labelmap island."
    ),
    components: ComponentChoice | None = typer.Option(
        None,
        "--components",
        help="Surface components to keep (default: preset).",
    ),
    grid_mm: float = typer.Option(
        defaults.DEFAULT_MERGE_GRID_MM,
        "--grid-mm",
        help="Isotropic fused-grid voxel size in mm.",
    ),
    mask_smooth_mm: float | None = typer.Option(
        None,
        "--mask-smooth-mm",
        help="Gaussian sigma in physical mm applied to the fused mask before meshing (0 = off).",
    ),
    surface_smooth_iters: int | None = typer.Option(
        None,
        "--mesh-smooth-iters",
        help="Topology-preserving surface relaxation iterations after meshing (0 = off).",
    ),
    simplify_error_mm: float | None = typer.Option(
        None,
        "--simplify-error-mm",
        help="MeshLib estimated surface-deviation/QEM limit in model mm, not a certified Hausdorff bound (0 = off).",
    ),
    post_surface_smooth_iters: int | None = typer.Option(
        None,
        "--post-mesh-smooth-iters",
        help="Final topology-preserving surface relaxation iterations after simplification (0 = off).",
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
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Register two scans of the same anatomy and fuse their surfaces."""
    _validate_mesh_output(output)
    fixed_threshold_value = _parse_threshold(fixed_threshold, "--fixed-threshold")
    moving_threshold_value = _parse_threshold(moving_threshold, "--moving-threshold")
    _validate_processing_numbers(
        nonnegative=[
            ("--median-mm", median_mm),
            ("--closing-mm", closing_mm),
            ("--opening-mm", opening_mm),
            ("--min-island-mm3", min_island_mm3),
            ("--mask-smooth-mm", mask_smooth_mm),
            ("--mesh-smooth-iters", surface_smooth_iters),
            ("--simplify-error-mm", simplify_error_mm),
            ("--post-mesh-smooth-iters", post_surface_smooth_iters),
        ],
        positive=(("--grid-mm", grid_mm),),
    )
    emitted_warnings: list[str] = []

    def emit_warning(message: str) -> None:
        emitted_warnings.append(message)
        _warn(message)

    progress = _ProgressDisplay(not quiet, "Discovering fixed volumes ...")
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
            keep_largest_component=_keep_largest_component(components),
        )
        progress.update("Loading merge engine ...")
        from . import merge as merge_mod

        try:
            result = merge_mod.merge(
                fixed=fixed_candidate,
                moving=moving_candidate,
                preset=resolved_preset,
                output_path=str(output),
                fixed_threshold=fixed_threshold_value,
                moving_threshold=moving_threshold_value,
                grid_mm=grid_mm,
                mask_smooth_mm=mask_smooth_mm,
                surface_smooth_iters=surface_smooth_iters,
                simplify_error_mm=simplify_error_mm,
                post_surface_smooth_iters=post_surface_smooth_iters,
                force=force,
                allow_large_volume=allow_large_volume,
                log=progress.log,
                warn=emit_warning,
            )
        except merge_mod.MergeError as exc:
            _error(exc)
            raise typer.Exit(3) from None
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        report = result.quality

        if json_file is not None:
            payload = {
                "result": {
                    "output": result.output_path,
                    "triangles": result.triangles,
                    "vertices": result.vertices,
                    "bounds_mm": list(result.bounds_mm),
                    "grid_mm": result.grid_mm,
                    "grid_size": list(result.grid_size),
                    "volume_fixed_mm3": result.volume_fixed_mm3,
                    "volume_moving_mm3": result.volume_moving_mm3,
                    "volume_fused_mm3": result.volume_union_mm3,
                    "surface_components": result.surface_components,
                    "seconds": result.seconds,
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

    _warn_if_invalid(report, "fused output")

    _exit_for_quality(report)


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
    json_output: bool = typer.Option(
        False, "--json", help="Write one plain JSON object to stdout."
    ),
) -> None:
    """Report mesh quality without changing the file."""
    progress = _ProgressDisplay(not json_output, "Loading validation engine ...")
    with progress:
        from . import validate as validate_mod

        progress.update("Validating mesh structure and self-intersections ...")
        try:
            report = validate_mod.validate(str(mesh))
        except (OSError, RuntimeError, ValueError) as exc:
            _error("cannot validate %s: %s" % (mesh, exc))
            raise typer.Exit(1) from None
    if json_output:
        print(json.dumps(report, indent=2))
    else:
        stdout_console.print(_plain(mesh, "bold"))
        _print_quality(report)
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
    json_output: bool = typer.Option(
        False, "--json", help="Write one plain JSON object to stdout."
    ),
    quiet: bool = typer.Option(
        False, "-q", "--quiet", help="Suppress normal progress output."
    ),
) -> None:
    """Make a non-watertight mesh watertight."""
    _validate_mesh_output(output)
    progress = _ProgressDisplay(
        not quiet and not json_output,
        "Loading repair engine ...",
    )
    with progress:
        from . import repair as repair_mod

        progress.update("Loading mesh for repair ...")
        try:
            result = repair_mod.repair(str(mesh), str(output), log=progress.log)
        except (OSError, RuntimeError, ValueError) as exc:
            _error("cannot repair %s: %s" % (mesh, exc))
            raise typer.Exit(1) from None

        stats = result.stats
        report = result.quality

    if json_output:
        print(
            json.dumps(
                {"output": str(output), "repair": stats, "quality": report}, indent=2
            )
        )
    elif not quiet:
        _success("wrote %s" % output)
        _print_quality(report)
