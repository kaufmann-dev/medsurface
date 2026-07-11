"""Typer command line interface with Rich human-readable output."""

from __future__ import annotations

import json
import subprocess
import sys
import warnings
from collections import Counter
from enum import Enum
from pathlib import Path
from typing import Any

import typer
from rich import box
from rich.console import Console, Group
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich.text import Text

from . import defaults, presets as presets_mod
from .presets import PRESETS


class PresetChoice(str, Enum):
    """Tissue presets accepted by ``--preset``."""

    AUTO = "auto"
    BONE = "bone"
    BONE_DETAIL = "bone-detail"
    SKIN = "skin"
    TEETH = "teeth"


stdout_console = Console(highlight=False, markup=False)
stderr_console = Console(stderr=True, highlight=False, markup=False)

app = typer.Typer(
    add_completion=False,
    help="Turn a DICOM series into a watertight 3D surface mesh.",
    invoke_without_command=True,
    no_args_is_help=False,
    pretty_exceptions_show_locals=False,
    rich_markup_mode="rich",
)


@app.callback()
def root(ctx: typer.Context) -> None:
    """Turn a DICOM series into a watertight 3D surface mesh."""
    if ctx.invoked_subcommand is None:
        stdout_console.print(ctx.get_help())


def _styled_message(prefix: str, prefix_style: str, message: object) -> Text:
    text = Text()
    text.append(prefix, style=prefix_style)
    text.append(str(message))
    return text


def _log(message: str) -> None:
    stdout_console.print(Text(str(message), style="cyan"))


class _ProgressDisplay:
    """Indeterminate stage progress with a plain-text redirected fallback."""

    def __init__(
        self,
        enabled: bool,
        initial: str,
        console: Console = stdout_console,
    ) -> None:
        self.enabled = enabled
        self.initial = initial
        self.console = console
        self.interactive = enabled and console.is_terminal and console.is_interactive
        self.progress: Progress | None = None
        self.task_id: int | None = None
        self.renderer: subprocess.Popen[str] | None = None

    def _start_renderer(self) -> bool:
        """Start a renderer process when the console has a real output descriptor."""
        try:
            self.console.file.fileno()
        except (AttributeError, OSError, ValueError):
            return False
        try:
            self.renderer = subprocess.Popen(
                [sys.executable, "-m", "dicom_surface.progress_renderer", self.initial],
                stdin=subprocess.PIPE,
                stdout=self.console.file,
                text=True,
                bufsize=1,
            )
        except OSError:
            self.renderer = None
            return False
        return True

    def _send(self, kind: str, message: str = "") -> None:
        if self.renderer is None or self.renderer.stdin is None:
            return
        try:
            self.renderer.stdin.write(json.dumps({"kind": kind, "message": message}) + "\n")
            self.renderer.stdin.flush()
        except (BrokenPipeError, OSError):
            pass

    def __enter__(self):
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


def _success(message: object) -> None:
    stdout_console.print(_styled_message("Success: ", "bold green", message))


def _warn(message: object) -> None:
    stderr_console.print(_styled_message("Warning: ", "bold yellow", message))


def _error(message: object) -> None:
    stderr_console.print(_styled_message("Error: ", "bold red", message))


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


def _parse_threshold(value: str | None) -> tuple[float | None, bool]:
    """Return an explicit numeric threshold and whether Otsu was requested."""
    if value is None:
        return None, False
    if value == "auto":
        return None, True
    try:
        return float(value), False
    except ValueError:
        _error("--threshold must be a number or 'auto'")
        raise typer.Exit(2) from None


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
    """Discover series while presenting pydicom warnings as concise CLI warnings."""
    from . import series as series_mod

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        found = series_mod.discover(str(root))
    _emit_discovery_warnings(captured)
    return found


def _plain(value: object, style: str | None = None) -> Text:
    return Text(str(value), style=style)


def _series_status(series: Any, recommended: Any | None) -> str:
    notes: list[str] = []
    if series is recommended:
        notes.append("default")
    if series.unusable_reason:
        notes.append(series.unusable_reason)
    else:
        notes.append("usable")
    if series.n_parts > 1:
        notes.append("orientation %d of %d in this UID" % (series.part, series.n_parts))
    if series.sharp_kernel:
        notes.append("sharp kernel %s" % series.kernel)
    return "; ".join(notes)


def _series_table(found: list[Any], recommended: Any | None) -> Table:
    """Build the responsive table used by ``list`` and rendering tests."""
    table = Table(
        box=box.SIMPLE_HEAVY,
        expand=True,
        padding=(0, 1),
        row_styles=("", "dim"),
    )
    table.add_column("ID", justify="right", no_wrap=True)
    table.add_column("DICOM #", justify="right", no_wrap=True)
    table.add_column("Modality", no_wrap=True)
    table.add_column("Description", ratio=2, overflow="fold")
    table.add_column("Slices", justify="right", no_wrap=True)
    table.add_column("Voxel (mm)", no_wrap=True)
    table.add_column("Plane", no_wrap=True)
    table.add_column("Status", ratio=3, overflow="fold")

    for series in found:
        voxel = "-"
        if series.pixel_spacing and series.slice_spacing:
            voxel = "%.3g × %.3g × %.3g" % (
                series.pixel_spacing[0],
                series.pixel_spacing[1],
                series.slice_spacing,
            )
        row_style = "bold cyan" if series is recommended else None
        table.add_row(
            _plain(series.id),
            _plain(series.series_number if series.series_number is not None else "-"),
            _plain(series.modality),
            _plain(series.description or "(none)"),
            _plain(series.n_slices),
            _plain(voxel),
            _plain(series.plane if series.usable else "-"),
            _plain(_series_status(series, recommended)),
            style=row_style,
        )
    return table


def _quality_table(report: dict[str, Any]) -> Table:
    """Build the shared two-column mesh-quality table."""
    table = Table(title="Mesh quality", box=box.SIMPLE_HEAVY, show_header=False)
    table.add_column("Metric", style="bold", no_wrap=True)
    table.add_column("Value")

    valid = bool(report["valid"])
    table.add_row(_plain("Status"), _plain("valid" if valid else "invalid", "green" if valid else "red"))
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
        value = f"{intersections:,}" if isinstance(intersections, int) else str(intersections)
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
    lines = [_plain("• " + problem) for problem in problems] if problems else [_plain("None", "green")]
    console.print(
        Panel(
            Group(*lines),
            title=_plain("Problems"),
            border_style="red" if problems else "green",
        )
    )


def _write_json_file(path: Path, payload: object) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


@app.command("list")
def list_series(
    dicom_dir: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        help="Directory tree containing DICOM instances.",
    ),
    json_output: bool = typer.Option(
        False,
        "--json",
        help="Write one plain JSON array to stdout.",
    ),
) -> None:
    """Show every series in a DICOM directory."""
    with _ProgressDisplay(not json_output, "Discovering DICOM series ..."):
        found = _discover(dicom_dir)
    if not found:
        _error("no DICOM instances found under %s" % dicom_dir)
        raise typer.Exit(1)

    if json_output:
        payload = [
            {
                "id": series.id,
                "uid": series.uid,
                "part": series.part,
                "n_parts": series.n_parts,
                "series_number": series.series_number,
                "modality": series.modality,
                "description": series.description,
                "slices": series.n_slices,
                "rows": series.rows,
                "columns": series.columns,
                "pixel_spacing": list(series.pixel_spacing) if series.pixel_spacing else None,
                "slice_spacing": series.slice_spacing,
                "plane": series.plane,
                "kernel": series.kernel,
                "sharp_kernel": series.sharp_kernel,
                "spacing_uniform": series.spacing_uniform,
                "spacing_spread_mm": series.spacing_spread_mm,
                "localizer": series.is_localizer,
                "usable": series.usable,
                "unusable_reason": series.unusable_reason,
            }
            for series in found
        ]
        print(json.dumps(payload, indent=2))
        return

    from . import series as series_mod

    ranked = series_mod.rank([series for series in found if series.usable])
    recommended = ranked[0] if ranked else None
    stdout_console.print(_series_table(found, recommended))
    stdout_console.print(
        Text(
            "Row IDs are local to this discovery result; run list again after directory contents change.",
            style="dim",
        )
    )
    command = Text("Convert the default with:  ")
    command.append("dicom-surface convert %s -o out.stl" % dicom_dir, style="bold")
    stdout_console.print(command)


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
        modality = ", ".join(preset.modalities) if preset.modalities else "any"
        threshold = preset.threshold if isinstance(preset.threshold, str) else "%g" % preset.threshold
        simplify = "off" if preset.simplify_error_mm == 0 else "%.2f mm" % preset.simplify_error_mm
        processing = "median %.1f mm; closing %.1f mm; smooth %d @ %.2f; simplify %s" % (
            preset.median_mm,
            preset.closing_mm,
            preset.smooth_iters,
            preset.smooth_force,
            simplify,
        )
        tissue.add_row(
            _plain(name),
            _plain(modality),
            _plain(threshold),
            _plain(processing),
            _plain(preset.description),
        )
    stdout_console.print(tissue)

@app.command()
def convert(
    dicom_dir: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        help="Directory tree containing the source DICOM series.",
    ),
    output: Path = typer.Option(..., "-o", "--output", help="Output .stl/.ply/.obj file."),
    series: str | None = typer.Option(
        None,
        "--series",
        help="Displayed row ID, complete SeriesInstanceUID, or description substring.",
    ),
    preset: PresetChoice = typer.Option(PresetChoice.BONE, "--preset"),
    threshold: str | None = typer.Option(
        None,
        "--threshold",
        help="Intensity (HU for CT) or 'auto' for Otsu.",
    ),
    median_mm: float | None = typer.Option(None, "--median-mm", help="Despeckle kernel extent, mm."),
    closing_mm: float | None = typer.Option(None, "--closing-mm", help="Pore-sealing kernel extent, mm."),
    opening_mm: float | None = typer.Option(None, "--opening-mm", help="Bridge-breaking kernel extent, mm."),
    min_island_mm3: float | None = typer.Option(None, "--min-island-mm3", help="Drop blobs smaller than this."),
    all_islands: bool = typer.Option(False, "--all-islands", help="Keep every labelmap island."),
    all_components: bool = typer.Option(False, "--all-components", help="Keep every surface shell."),
    resample_mm: float | None = typer.Option(
        None,
        "--resample-mm",
        help="Isotropic surface-grid voxel size in mm (0 = native).",
    ),
    smooth_iters: int | None = typer.Option(None, "--smooth-iters", help="MeshLib relaxation iterations."),
    smooth_force: float | None = typer.Option(
        None, "--smooth-force", help="MeshLib relaxation strength per iteration."
    ),
    simplify_error_mm: float | None = typer.Option(
        None,
        "--simplify-error-mm",
        help="MeshLib estimated surface-deviation/QEM limit in model mm, not a certified Hausdorff bound (0 = off).",
    ),
    post_smooth_iters: int | None = typer.Option(None, "--post-smooth-iters", help="Smoothing after simplification."),
    no_cap: bool = typer.Option(False, "--no-cap", help="Do not close anatomy at the field-of-view boundary."),
    json_file: Path | None = typer.Option(None, "--json", help="Write results and provenance to this JSON file."),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress normal progress output."),
) -> None:
    """Extract a surface mesh from a DICOM series."""
    threshold_value, auto_threshold = _parse_threshold(threshold)
    progress = _ProgressDisplay(not quiet, "Discovering DICOM series ...")
    with progress:
        found = _discover(dicom_dir)
        if not found:
            _error("no DICOM instances found under %s" % dicom_dir)
            raise typer.Exit(1)

        from . import series as series_mod

        try:
            chosen = series_mod.select(found, series)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        progress.log("series ID %d  %s  (%d slices)" % (chosen.id, chosen.label(), chosen.n_slices))
        resolved_preset = presets_mod.get(preset.value)
        resolved_preset = presets_mod.override(
            resolved_preset,
            median_mm=median_mm,
            closing_mm=closing_mm,
            opening_mm=opening_mm,
            min_island_mm3=min_island_mm3,
            resample_mm=resample_mm,
            smooth_iters=smooth_iters,
            smooth_force=smooth_force,
            simplify_error_mm=simplify_error_mm,
            post_smooth_iters=post_smooth_iters,
            keep_largest_island=False if all_islands else None,
            keep_largest_component=False if all_components else None,
        )
        if auto_threshold:
            resolved_preset = presets_mod.override(resolved_preset, threshold="auto")

        progress.update("Loading conversion engine ...")
        from . import pipeline

        try:
            result = pipeline.convert(
                series=chosen,
                preset=resolved_preset,
                output_path=str(output),
                threshold=threshold_value,
                cap_field_of_view=not no_cap,
                log=progress.log,
            )
        except pipeline.ModalityMismatch as exc:
            _error(exc)
            raise typer.Exit(2) from None
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
            _write_json_file(json_file, payload)

    for message in result.warnings:
        _warn(message)

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


@app.command()
def merge(
    dicom_dir_a: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        help="Fixed scan; defines the output coordinate frame.",
    ),
    output: Path = typer.Option(..., "-o", "--output", help="Output .stl/.ply/.obj file."),
    dicom_dir_b: Path | None = typer.Argument(
        None,
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        help="Moving scan; omit to choose two series from the first directory.",
    ),
    series_a: str | None = typer.Option(
        None,
        "--series-a",
        help="Fixed scan row ID, complete UID, or description substring.",
    ),
    series_b: str | None = typer.Option(
        None,
        "--series-b",
        help="Moving scan row ID, complete UID, or description substring.",
    ),
    preset: PresetChoice = typer.Option(PresetChoice.BONE, "--preset"),
    threshold: str | None = typer.Option(
        None,
        "--threshold",
        help="Intensity (HU for CT) or 'auto'; applies to both scans.",
    ),
    median_mm: float | None = typer.Option(None, "--median-mm"),
    closing_mm: float | None = typer.Option(None, "--closing-mm"),
    min_island_mm3: float | None = typer.Option(None, "--min-island-mm3"),
    grid_mm: float = typer.Option(
        defaults.DEFAULT_MERGE_GRID_MM,
        "--grid-mm",
        help="Isotropic fused-grid voxel size in mm.",
    ),
    smooth_iters: int | None = typer.Option(None, "--smooth-iters", help="Default: from the preset."),
    smooth_force: float | None = typer.Option(
        None, "--smooth-force", help="Default: from the preset."
    ),
    simplify_error_mm: float | None = typer.Option(
        None,
        "--simplify-error-mm",
        help="MeshLib estimated surface-deviation/QEM limit in model mm, not a certified Hausdorff bound (0 = off).",
    ),
    post_smooth_iters: int | None = typer.Option(None, "--post-smooth-iters", help="Default: from the preset."),
    force: bool = typer.Option(
        False,
        "--force",
        help="Override patient, modality, and registration gates.",
    ),
    json_file: Path | None = typer.Option(None, "--json", help="Write results and provenance to this JSON file."),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress normal progress output."),
) -> None:
    """Register two scans of the same anatomy and fuse their surfaces."""
    threshold_value, auto_threshold = _parse_threshold(threshold)
    directory_b = dicom_dir_b or dicom_dir_a
    progress = _ProgressDisplay(not quiet, "Discovering fixed DICOM series ...")
    with progress:
        found_a = _discover(dicom_dir_a)
        if directory_b == dicom_dir_a:
            found_b = found_a
        else:
            progress.update("Discovering moving DICOM series ...")
            found_b = _discover(directory_b)
        if not found_a or not found_b:
            _error("no DICOM instances found")
            raise typer.Exit(1)

        from . import series as series_mod

        try:
            chosen_a = series_mod.select(found_a, series_a)
            chosen_b = series_mod.select(found_b, series_b)
        except ValueError as exc:
            _error(exc)
            raise typer.Exit(2) from None

        resolved_preset = presets_mod.get(preset.value)
        resolved_preset = presets_mod.override(
            resolved_preset,
            median_mm=median_mm,
            closing_mm=closing_mm,
            min_island_mm3=min_island_mm3,
        )
        if auto_threshold:
            resolved_preset = presets_mod.override(resolved_preset, threshold="auto")

        progress.update("Loading merge engine ...")
        from . import merge as merge_mod, pipeline

        try:
            result = merge_mod.merge(
                series_a=chosen_a,
                series_b=chosen_b,
                preset=resolved_preset,
                output_path=str(output),
                threshold=threshold_value,
                grid_mm=grid_mm,
                smooth_iters=smooth_iters,
                smooth_force=smooth_force,
                simplify_error_mm=simplify_error_mm,
                post_smooth_iters=post_smooth_iters,
                force=force,
                log=progress.log,
            )
        except merge_mod.MergeError as exc:
            _error(exc)
            raise typer.Exit(3) from None
        except (pipeline.ModalityMismatch, ValueError) as exc:
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
                    "volume_fixed_mm3": result.volume_a_mm3,
                    "volume_moving_mm3": result.volume_b_mm3,
                    "volume_fused_mm3": result.volume_union_mm3,
                    "seconds": result.seconds,
                    "warnings": result.warnings,
                },
                "provenance": result.provenance,
                "quality": report,
            }
            progress.update("Writing JSON report ...")
            _write_json_file(json_file, payload)

    for message in result.warnings:
        _warn(message)

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
    json_output: bool = typer.Option(False, "--json", help="Write one plain JSON object to stdout."),
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
    json_output: bool = typer.Option(False, "--json", help="Write one plain JSON object to stdout."),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress normal progress output."),
) -> None:
    """Make a non-watertight mesh watertight."""
    progress = _ProgressDisplay(
        not quiet and not json_output,
        "Loading repair engine ...",
    )
    with progress:
        from . import repair as repair_mod

        progress.update("Loading mesh for repair ...")
        try:
            stats = repair_mod.repair(str(mesh), str(output), log=progress.log)
        except (OSError, RuntimeError, ValueError) as exc:
            _error("cannot repair %s: %s" % (mesh, exc))
            raise typer.Exit(1) from None

        progress.update("Loading validation engine ...")
        from . import validate as validate_mod

        progress.update("Validating repaired mesh ...")
        try:
            report = validate_mod.validate(str(output))
        except (OSError, RuntimeError, ValueError) as exc:
            _error("repaired mesh could not be validated: %s" % exc)
            raise typer.Exit(1) from None

    if json_output:
        print(json.dumps({"output": str(output), "repair": stats, "quality": report}, indent=2))
    elif not quiet:
        _success("wrote %s" % output)
        _print_quality(report)
    _warn_if_invalid(report, "repaired output")
    _exit_for_quality(report)
