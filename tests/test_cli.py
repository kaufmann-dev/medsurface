"""Cross-command Typer and Rich CLI contracts."""

from __future__ import annotations

import io
import json
import os
import pty
import re
import select
import shlex
import signal
import subprocess
import sys
import time
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest
from rich.console import Console
from typer.testing import CliRunner

from medsurface import __version__, cli
from medsurface.catalog import DicomSource, FileSource, VolumeCandidate
from medsurface.series import Series

runner = CliRunner(env={"_TYPER_FORCE_DISABLE_TERMINAL": "1"})


def _expected_shell_command(arguments) -> str:
    rendered = [str(argument) for argument in arguments]
    if os.name == "nt":
        return subprocess.list2cmdline(rendered)
    return shlex.join(rendered)


def _plain_subprocess_env() -> dict[str, str]:
    """Keep captured Typer output plain even when GitHub Actions forces color."""
    env = os.environ.copy()
    env["_TYPER_FORCE_DISABLE_TERMINAL"] = "1"
    return env


def _series(
    *,
    row_id: int = 1,
    uid: str = "1.2.3",
    number: int = 6,
    modality: str = "CT",
    description: str = "axial bone",
    kernel_values: tuple[str, ...] = (),
) -> Series:
    value = Series(
        uid=uid,
        modality=modality,
        description=description,
        series_number=number,
        id=row_id,
        kernel_values=kernel_values,
    )
    value.files = ["slice-%d" % index for index in range(8)]
    value.pixel_spacing = (0.5, 0.5)
    value.slice_spacing = 1.0
    value.normal = (0.0, 0.0, 1.0)
    return value


def _candidate(**kwargs) -> VolumeCandidate:
    series = _series(**kwargs)
    return VolumeCandidate(
        id=series.id,
        source=DicomSource(Path("scans"), series),
        format="DICOM",
        source_name="scans",
        modality=series.modality,
        description=series.description,
        size=(512, 512, series.n_slices),
        spacing=(0.5, 0.5, 1.0),
        direction=None,
        origin=None,
        pixel_type=None,
        components=1,
        plane=series.plane,
        unusable_reason=series.unusable_reason,
    )


def _file_candidate(path: Path, *, row_id: int = 1) -> VolumeCandidate:
    return VolumeCandidate(
        id=row_id,
        source=FileSource(path, "MetaImage"),
        format="MetaImage",
        source_name=path.name,
        modality=None,
        description=None,
        size=(8, 8, 8),
        spacing=(1.0, 1.0, 1.0),
        direction=(1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0),
        origin=(0.0, 0.0, 0.0),
        pixel_type="16-bit signed integer",
        components=1,
        plane="axial",
    )


def _quality(*, valid: bool = True) -> dict:
    problems = [] if valid else ["mesh is not watertight", "3 boundary edge(s)"]
    return {
        "file": "mesh.stl",
        "bytes": 2_097_152,
        "triangles": 1200,
        "vertices": 602,
        "components": 1,
        "watertight": valid,
        "winding_consistent": True,
        "is_volume": valid,
        "euler_number": 2,
        "boundary_edges": 0 if valid else 3,
        "holes": 0 if valid else 1,
        "disoriented_faces": 0,
        "area_mm2": 100.0,
        "volume_mm3": 90.0 if valid else None,
        "bbox_min": [0.0, 0.0, 0.0],
        "bbox_max": [1.0, 2.0, 3.0],
        "bbox_extents_mm": [1.0, 2.0, 3.0],
        "genus": 0,
        "self_intersecting_faces": 0,
        "problems": problems,
        "valid": valid,
    }


def _extract_result(
    output: str, *, warnings_: list[str] | None = None, quality: dict | None = None
):
    return SimpleNamespace(
        output_path=output,
        triangles=1200,
        vertices=602,
        bounds_mm=(0.0, 1.0, 0.0, 2.0, 0.0, 3.0),
        seconds=1.25,
        capped_field_of_view=False,
        labelmap_components=2,
        surface_components=1,
        warnings=warnings_ or [],
        provenance={"input": {"format": "DICOM"}},
        quality=_quality() if quality is None else quality,
    )


def _conversion_result(output: str, *, warnings_: list[str] | None = None):
    return SimpleNamespace(
        output_path=output,
        format="NIfTI",
        compression="gzip",
        dimensions=(8, 9, 10),
        pixel_type="16-bit signed integer",
        components=1,
        spacing=(0.7, 0.8, 0.9),
        origin=(10.0, 20.0, 30.0),
        direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        seconds=1.5,
        warnings=warnings_ or [],
        metadata_policy="preserved-best-effort",
        provenance={"input": {"format": "DICOM"}},
    )


def _fusion_result(output: str, *, warnings_: list[str] | None = None):
    return SimpleNamespace(
        output_path=output,
        output_format="NRRD",
        compression="gzip",
        grid_mm=0.8,
        grid_size=(10, 20, 30),
        grid_origin_mm=(1.0, 2.0, 3.0),
        grid_direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        foreground_fixed_voxels=195,
        foreground_moving_voxels=215,
        foreground_fused_voxels=293,
        volume_fixed_mm3=100.0,
        volume_moving_mm3=110.0,
        volume_fused_mm3=150.0,
        seconds=2.5,
        warnings=warnings_ or [],
        provenance={
            "output": {"kind": "binary labelmap"},
            "registration": {"surface_overlap": 0.95},
        },
    )


def _label_fusion_result(output: str, *, warnings_: list[str] | None = None):
    return SimpleNamespace(
        output_path=output,
        output_format="NRRD",
        compression="gzip",
        pixel_type="uint8",
        preserve_labels=False,
        grid_mm=0.8,
        grid_size=(10, 20, 30),
        grid_origin_mm=(1.0, 2.0, 3.0),
        grid_direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        labels={1: {"name": None, "voxels": 293, "volume_mm3": 150.0}},
        volume_fused_mm3=150.0,
        seconds=2.5,
        warnings=warnings_ or [],
        provenance={"output": {"kind": "binary labelmap"}},
    )


def _run_cli_in_clean_interpreter(argv):
    code = """
import json
import sys
from typer.testing import CliRunner
from medsurface import cli

heavy_modules = {
    "medsurface.labelmap",
    "medsurface.fusion",
    "medsurface.pipeline",
    "medsurface.registration",
    "medsurface.repair",
    "medsurface.surface",
    "medsurface.validate",
    "medsurface.volume",
    "meshlib",
}
assert heavy_modules.isdisjoint(sys.modules)
result = CliRunner().invoke(
    cli.app,
    json.loads(sys.argv[1]),
    prog_name="medsurface",
)
assert heavy_modules.isdisjoint(sys.modules)
sys.stdout.write(result.stdout)
sys.stderr.write(result.stderr)
raise SystemExit(result.exit_code)
"""
    return subprocess.run(
        [sys.executable, "-c", code, json.dumps(argv)],
        capture_output=True,
        text=True,
        check=False,
        env=_plain_subprocess_env(),
    )


def test_root_command_without_arguments_is_lightweight_help():
    result = _run_cli_in_clean_interpreter([])

    assert result.returncode == 0
    assert result.stderr == ""
    assert "\x1b" not in result.stdout
    assert "Usage: medsurface [OPTIONS] [COMMAND]" in result.stdout
    for command in (
        "list",
        "presets",
        "convert",
        "extract",
        "fuse",
        "validate",
        "repair",
        "labelmap",
    ):
        assert command in result.stdout


def test_version_is_lightweight_and_public():
    result = _run_cli_in_clean_interpreter(["--version"])

    assert result.returncode == 0
    assert result.stdout.strip() == "medsurface %s" % __version__
    assert result.stderr == ""


@pytest.mark.parametrize(
    "argv",
    [
        ["--help"],
        ["list", "--help"],
        ["presets", "--help"],
        ["convert", "--help"],
        ["extract", "--help"],
        ["fuse", "--help"],
        ["validate", "--help"],
        ["repair", "--help"],
        ["labelmap", "--help"],
        ["labelmap", "extract", "--help"],
        ["labelmap", "fuse", "--help"],
    ],
)
def test_every_help_surface_is_lightweight(argv):
    result = _run_cli_in_clean_interpreter(argv)

    assert result.returncode == 0
    assert "Usage: medsurface" in result.stdout
    assert result.stderr == ""
    assert "\x1b" not in result.stdout
    assert "Traceback" not in result.stdout


@pytest.mark.parametrize(
    "argv",
    [
        ["extract", "INPUT", "-o", "out.stl", "--post-smooth-iters", "1"],
        [
            "labelmap",
            "extract",
            "MASK",
            "-o",
            "out.stl",
            "--post-smooth-iters",
            "1",
        ],
    ],
)
def test_post_smoothing_option_is_removed(tmp_path, argv):
    mask = tmp_path / "mask.nii.gz"
    mask.write_bytes(b"mask")
    resolved = [
        str(tmp_path) if value == "INPUT" else str(mask) if value == "MASK" else value
        for value in argv
    ]

    result = runner.invoke(cli.app, resolved, prog_name="medsurface")

    assert result.exit_code == 2
    assert "--post-smooth-iters" in result.stderr
    assert "No such option" in result.stderr


@pytest.mark.parametrize(
    "option",
    ["--smooth-mm", "--smooth-iters", "--smooth-force"],
)
@pytest.mark.parametrize(
    "argv",
    [
        ["extract", "INPUT", "-o", "out.stl"],
        ["labelmap", "extract", "INPUT", "-o", "out.stl"],
    ],
)
def test_legacy_smoothing_options_are_removed(tmp_path, argv, option):
    mask = tmp_path / "mask.nii.gz"
    moving = tmp_path / "moving.nii.gz"
    mask.write_bytes(b"mask")
    moving.write_bytes(b"moving")
    resolved = [
        str(mask) if value == "INPUT" else str(moving) if value == "MOVING" else value
        for value in argv
    ]
    resolved.extend([option, "1"])

    result = runner.invoke(cli.app, resolved, prog_name="medsurface")

    assert result.exit_code == 2
    assert option in result.stderr
    assert "No such option" in result.stderr


def test_surface_help_describes_both_independent_smoothing_stages():
    for command in (
        ["extract", "--help"],
        ["labelmap", "extract", "--help"],
    ):
        result = _run_cli_in_clean_interpreter(command)
        normalized = " ".join(result.stdout.replace("│", " ").split())

        assert result.returncode == 0
        assert "--mask-smooth-mm" in normalized
        assert "Gaussian sigma" in normalized
        assert "--mesh-smooth-iters" in normalized
        assert "surface relaxation iterations" in normalized
        assert "--post-mesh-smooth-" in normalized
        assert "after simplification" in normalized
        assert "--smooth-mm" not in normalized
        assert "--smooth-iters" not in normalized
        assert "--smooth-force" not in normalized


def test_component_help_is_available_on_every_surface_command():
    for command in (
        ["extract", "--help"],
        ["labelmap", "extract", "--help"],
    ):
        result = _run_cli_in_clean_interpreter(command)
        normalized = " ".join(result.stdout.replace("│", " ").split())

        assert result.returncode == 0
        assert "--components" in normalized
        assert "[all|largest]" in normalized
        assert "--all-components" not in normalized


@pytest.mark.parametrize(
    "argv",
    [
        ["extract", "INPUT", "-o", "out.stl", "--all-components"],
        ["labelmap", "extract", "INPUT", "-o", "out.stl", "--all-components"],
    ],
)
def test_all_components_option_is_removed(tmp_path, argv):
    source = tmp_path / "input.nii.gz"
    moving = tmp_path / "moving.nii.gz"
    source.write_bytes(b"input")
    moving.write_bytes(b"moving")
    resolved = [
        str(source) if value == "INPUT" else str(moving) if value == "MOVING" else value
        for value in argv
    ]

    result = runner.invoke(cli.app, resolved, prog_name="medsurface")

    assert result.exit_code == 2
    assert "--all-components" in result.stderr
    assert "No such option" in result.stderr


@pytest.mark.parametrize(
    "argv",
    [
        ["extract", "INPUT", "-o", "out.stl", "--components", "first"],
        [
            "labelmap",
            "extract",
            "INPUT",
            "-o",
            "out.stl",
            "--components",
            "first",
        ],
    ],
)
def test_invalid_component_choice_fails_before_discovery(tmp_path, monkeypatch, argv):
    source = tmp_path / "input.nii.gz"
    moving = tmp_path / "moving.nii.gz"
    source.write_bytes(b"input")
    moving.write_bytes(b"moving")
    resolved = [
        str(source) if value == "INPUT" else str(moving) if value == "MOVING" else value
        for value in argv
    ]
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda *_args: pytest.fail("discovery must not start"),
    )
    monkeypatch.setattr(
        cli,
        "_select_labelmap",
        lambda *_args: pytest.fail("labelmap discovery must not start"),
    )

    result = runner.invoke(cli.app, resolved, prog_name="medsurface")

    assert result.exit_code == 2
    assert "--components" in result.stderr
    assert "all" in result.stderr
    assert "largest" in result.stderr


def test_fuse_help_describes_segmentation_options_without_surface_controls():
    result = _run_cli_in_clean_interpreter(["fuse", "--help"])

    assert result.returncode == 0
    normalized = " ".join(result.stdout.replace("│", " ").split())
    for option, description_start in (
        ("--median-mm", "Despeckle kernel"),
        ("--closing-mm", "Pore-sealing kernel"),
        ("--opening-mm", "Bridge-breaking kernel"),
        ("--min-island-mm3", "Drop blobs smaller"),
    ):
        assert option in normalized
        assert description_start in normalized
    for option in (
        "--components",
        "--mask-smooth-mm",
        "--mesh-smooth-iters",
        "--simplify-error-mm",
        "--post-mesh-smooth-iters",
        "--no-cap",
    ):
        assert option not in normalized


def test_strict_convert_help_exposes_only_storage_controls():
    result = _run_cli_in_clean_interpreter(["convert", "--help"])

    assert result.returncode == 0
    normalized = " ".join(result.stdout.replace("│", " ").split())
    for option in ("--volume", "--strip-metadata", "--allow-large-volume", "--json"):
        assert option in normalized
    for option in (
        "--preset",
        "--threshold",
        "--median-mm",
        "--components",
        "--mask-smooth-mm",
        "--mesh-smooth-iters",
        "--simplify-error-mm",
        "--no-cap",
    ):
        assert option not in normalized


def test_legacy_command_names_are_absent():
    root_help = _run_cli_in_clean_interpreter(["--help"])
    labelmap_help = _run_cli_in_clean_interpreter(["labelmap", "--help"])

    assert root_help.returncode == 0
    assert labelmap_help.returncode == 0
    assert " merge " not in " ".join(root_help.stdout.split())
    assert " convert " not in " ".join(labelmap_help.stdout.split())
    assert " merge " not in " ".join(labelmap_help.stdout.split())
    for argv in (["merge"], ["labelmap", "convert"], ["labelmap", "merge"]):
        result = _run_cli_in_clean_interpreter(argv)
        assert result.returncode == 2
        assert "No such command" in result.stderr


@pytest.mark.parametrize(
    "argv",
    [
        ["unknown"],
        ["list"],
        ["presets", "unexpected"],
        ["convert", "scan"],
        ["extract", "scan", "-o", "out.stl", "--preset", "unknown"],
        ["fuse"],
        ["validate"],
        ["repair", "broken.stl"],
        ["labelmap", "extract"],
        ["labelmap", "fuse"],
    ],
)
def test_malformed_invocations_fail_before_heavy_imports(argv):
    result = _run_cli_in_clean_interpreter(argv)

    assert result.returncode == 2
    assert result.stdout == ""
    assert result.stderr.startswith("Usage: medsurface")
    assert "\x1b" not in result.stderr
    assert "Error" in result.stderr
    assert "Traceback" not in result.stderr


def test_python_module_uses_public_program_name():
    result = subprocess.run(
        [sys.executable, "-m", "medsurface", "--help"],
        capture_output=True,
        text=True,
        check=False,
        env=_plain_subprocess_env(),
    )

    assert result.returncode == 0
    assert "Usage: medsurface" in result.stdout
    assert "\x1b" not in result.stdout


def test_registry_enums_match_presets_exactly():
    assert {choice.value for choice in cli.PresetChoice} == set(cli.PRESETS)


def test_removed_bone_detail_preset_is_rejected(tmp_path):
    result = runner.invoke(
        cli.app,
        ["extract", str(tmp_path), "-o", "out.stl", "--preset", "bone-detail"],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "bone-detail" in result.stderr


def test_every_validating_command_uses_the_same_quality_status():
    assert cli._quality_status(None) == 0
    assert cli._quality_status({"valid": True}) == 0
    assert cli._quality_status({"valid": False}) == 1


def test_progress_display_has_plain_redirected_fallback_and_live_elapsed_time():
    plain_stream = io.StringIO()
    plain_console = Console(
        file=plain_stream,
        width=100,
        color_system=None,
        force_terminal=False,
        highlight=False,
        markup=False,
    )
    with cli._ProgressDisplay(
        True, "Starting [literal] ...", plain_console
    ) as progress:
        progress.log("segment [literal] ...")
        progress.log("  segment  1.2s")

    plain = plain_stream.getvalue()
    assert "Starting [literal] ..." in plain
    assert "segment [literal] ..." in plain
    assert "segment  1.2s" in plain
    assert "\x1b" not in plain

    live_stream = io.StringIO()
    live_console = Console(
        file=live_stream,
        width=100,
        color_system="standard",
        force_terminal=True,
        force_interactive=True,
        highlight=False,
        markup=False,
    )
    with cli._ProgressDisplay(True, "Starting ...", live_console) as progress:
        progress.update("marching cubes ...")
        assert progress.interactive
        assert progress.progress is not None
        assert any(
            isinstance(column, cli.TimeElapsedColumn)
            for column in progress.progress.columns
        )
        assert progress.progress.tasks[0].description == "marching cubes ..."


def test_interactive_progress_renderer_runs_in_an_independent_process():
    master, slave = pty.openpty()
    output = os.fdopen(slave, "w", buffering=1, closefd=True)
    terminal_console = Console(
        file=output,
        width=100,
        color_system="standard",
        force_terminal=True,
        force_interactive=True,
        highlight=False,
        markup=False,
    )
    try:
        with cli._ProgressDisplay(True, "Starting ...", terminal_console) as progress:
            assert progress.renderer is not None
            progress.update("simplify ...")
            time.sleep(1.2)
        rendered = bytearray()
        while select.select([master], [], [], 0.1)[0]:
            rendered.extend(os.read(master, 65536))
    finally:
        output.close()
        os.close(master)

    text = rendered.decode(errors="replace")
    assert "simplify ..." in text
    assert "0:00:01" in text


def test_interactive_diagnostic_clears_live_progress_line():
    master, slave = pty.openpty()
    output = os.fdopen(slave, "w", buffering=1, closefd=True)
    terminal_console = Console(
        file=output,
        width=100,
        color_system="standard",
        force_terminal=True,
        force_interactive=True,
        highlight=False,
        markup=False,
    )
    try:
        with cli._ProgressDisplay(
            True,
            "load volume ...",
            terminal_console,
            terminal_console,
        ):
            time.sleep(0.2)
            cli._warn("voxels are strongly anisotropic")
            time.sleep(0.2)
        rendered = bytearray()
        while select.select([master], [], [], 0.1)[0]:
            rendered.extend(os.read(master, 65536))
    finally:
        output.close()
        os.close(master)

    text = rendered.decode(errors="replace")
    warning_at = text.index("Warning:")
    active_line = text[text.rfind("\r", 0, warning_at) + 1 : warning_at]
    assert "0:00:" not in active_line
    assert "voxels are strongly anisotropic" in text


def _interruptible_extract_script(input_path: Path, output_path: Path) -> str:
    return f"""
import time
from types import SimpleNamespace

from medsurface import catalog, cli, pipeline

chosen = SimpleNamespace(id=1, format="NIfTI", source_name="fixture.nii.gz")
cli._discover = lambda _path: [chosen]
catalog.select = lambda _found, _volume_id: chosen
cli._protect_output_paths = lambda *_paths: None

def fake_extract(**_kwargs):
    progress = cli._active_progress.get()
    renderer_pid = progress.renderer.pid if progress.renderer is not None else "none"
    print(f"READY renderer={{renderer_pid}}", flush=True)
    time.sleep(30)

pipeline.extract = fake_extract
cli.app(
    prog_name="medsurface",
    args=["extract", {str(input_path)!r}, "-o", {str(output_path)!r}],
)
"""


def _read_pty_until(master: int, needle: bytes, timeout: float) -> bytearray:
    rendered = bytearray()
    deadline = time.monotonic() + timeout
    while needle not in rendered and time.monotonic() < deadline:
        if not select.select([master], [], [], 0.1)[0]:
            continue
        try:
            rendered.extend(os.read(master, 65536))
        except OSError:
            break
    return rendered


def test_ctrl_c_is_owned_by_parent_without_renderer_traceback(tmp_path):
    output_path = tmp_path / "cancelled.stl"
    master, slave = pty.openpty()
    process = subprocess.Popen(
        [
            sys.executable,
            "-c",
            _interruptible_extract_script(tmp_path, output_path),
        ],
        stdin=slave,
        stdout=slave,
        stderr=slave,
        close_fds=True,
        start_new_session=True,
        env=_plain_subprocess_env(),
    )
    os.close(slave)
    rendered = bytearray()
    try:
        rendered.extend(_read_pty_until(master, b"READY renderer=", 5))
        assert b"READY renderer=" in rendered
        time.sleep(0.3)
        os.killpg(process.pid, signal.SIGINT)
        assert process.wait(timeout=5) == 130
        rendered.extend(_read_pty_until(master, b"Cancelled.", 2))
        rendered.extend(_read_pty_until(master, b"unused sentinel", 0.2))
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        os.close(master)

    text = rendered.decode(errors="replace")
    renderer_match = re.search(r"READY renderer=(\d+)", text)
    assert renderer_match is not None
    assert text.count("Cancelled.") == 1
    assert "Traceback" not in text
    assert "KeyboardInterrupt" not in text
    assert not output_path.exists()
    with pytest.raises(ProcessLookupError):
        os.kill(int(renderer_match.group(1)), 0)


def test_quiet_ctrl_c_exits_cleanly_without_renderer(tmp_path):
    output_path = tmp_path / "cancelled-quiet.stl"
    script = _interruptible_extract_script(tmp_path, output_path).replace(
        '"-o",', '"--quiet", "-o",'
    )
    process = subprocess.Popen(
        [sys.executable, "-c", script],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
        env=_plain_subprocess_env(),
    )
    try:
        assert process.stdout is not None
        assert select.select([process.stdout], [], [], 5)[0]
        assert process.stdout.readline().strip() == "READY renderer=none"
        process.send_signal(signal.SIGINT)
        stdout, stderr = process.communicate(timeout=5)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)

    assert process.returncode == 130
    assert stdout == ""
    assert stderr == "Cancelled.\n"
    assert not output_path.exists()


def test_repeated_ctrl_c_does_not_interrupt_cancellation_cleanup():
    class FakePipe:
        def write(self, _message):
            return None

        def flush(self):
            return None

        def close(self):
            return None

    class InterruptingRenderer:
        stdin = FakePipe()
        waited = False

        def wait(self, *, timeout):
            os.kill(os.getpid(), signal.SIGINT)
            self.waited = True

    diagnostic_stream = io.StringIO()
    diagnostic_console = Console(
        file=diagnostic_stream,
        color_system=None,
        force_terminal=False,
        highlight=False,
        markup=False,
    )
    display = cli._ProgressDisplay(
        False,
        "Starting ...",
        diagnostic_console=diagnostic_console,
    )
    renderer = InterruptingRenderer()

    with pytest.raises(KeyboardInterrupt):
        with display:
            display.renderer = renderer  # type: ignore[assignment]
            raise KeyboardInterrupt

    assert renderer.waited
    assert diagnostic_stream.getvalue() == "Cancelled.\n"


@pytest.mark.parametrize(
    "argv",
    [
        ["extract", ".", "-o", "out.stl", "--self-intersections"],
        ["fuse", ".", ".", "-o", "out.nrrd", "--self-intersections"],
        ["validate", __file__, "--self-intersections"],
        ["repair", __file__, "-o", "fixed.stl", "--self-intersections"],
        ["extract", ".", "-o", "out.stl", "--target-faces", "100"],
        ["fuse", ".", ".", "-o", "out.nrrd", "--target-faces", "100"],
        ["extract", ".", "-o", "out.stl", "--no-validate"],
        ["fuse", ".", ".", "-o", "out.nrrd", "--no-validate"],
        ["extract", ".", "-o", "out.stl", "--passband", "0.1"],
        ["fuse", ".", ".", "-o", "out.nrrd", "--passband", "0.1"],
        ["extract", ".", "-o", "out.stl", "--print-profile", "resin"],
        ["fuse", ".", ".", "-o", "out.nrrd", "--print-profile", "fdm"],
        ["extract", ".", "-o", "out.stl", "--min-feature-mm", "0.6"],
        ["fuse", ".", ".", "-o", "out.nrrd", "--min-feature-mm", "1.2"],
    ],
)
def test_removed_self_intersection_flag_is_a_usage_error(argv):
    result = runner.invoke(cli.app, argv, prog_name="medsurface")

    assert result.exit_code == 2
    assert "No such option" in result.stderr


@pytest.mark.parametrize(
    "command, argv",
    [
        ("list", ["list", "missing"]),
        ("convert", ["convert", "missing", "-o", "out.nrrd"]),
        ("extract", ["extract", "missing", "-o", "out.stl"]),
        ("fuse", ["fuse", "missing", "also-missing", "-o", "out.nrrd"]),
        ("validate", ["validate", "missing.stl"]),
        ("repair", ["repair", "missing.stl", "-o", "fixed.stl"]),
    ],
)
def test_required_inputs_use_typer_path_validation(command, argv):
    result = runner.invoke(cli.app, argv, prog_name="medsurface")

    assert result.exit_code == 2
    assert "does not exist" in result.stderr
    assert command in result.stderr


def test_list_json_writes_unique_ids_to_a_file(tmp_path, monkeypatch):
    found = [
        _candidate(
            description="[bold red]literal[/bold red]",
            kernel_values=("Hr68f", "1"),
        )
    ]
    monkeypatch.setattr(cli, "_discover", lambda _root: found)

    report = tmp_path / "volumes.json"
    result = runner.invoke(
        cli.app,
        ["list", str(tmp_path), "--json", str(report), "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0
    payload = json.loads(report.read_text())
    assert payload[0]["id"] == 1
    assert payload[0]["default"] is True
    assert payload[0]["format"] == "DICOM"
    assert payload[0]["dicom"]["series_number"] == 6
    assert payload[0]["dicom"]["part"] == 1
    assert payload[0]["dicom"]["n_parts"] == 1
    assert payload[0]["dicom"]["kernel"] == ["Hr68f", "1"]
    assert payload[0]["description"] == "[bold red]literal[/bold red]"
    assert "ident" not in payload[0]
    assert result.stdout == ""
    assert result.stderr == ""


def test_list_human_output_shows_discovery_progress(tmp_path, monkeypatch):
    candidate = _candidate(kernel_values=("Hr68f", "1"))
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: [candidate],
    )

    result = runner.invoke(cli.app, ["list", str(tmp_path)], prog_name="medsurface")

    assert result.exit_code == 0
    assert "Discovering volumes ..." in result.stdout
    assert "DICOM #" in result.stdout
    assert "axial" in result.stdout
    assert cli._volume_status(candidate, candidate) == "default; usable; sharp kernel"

    stream = io.StringIO()
    console = Console(file=stream, width=160, color_system=None, force_terminal=False)
    console.print(cli._volume_table([candidate], candidate))
    rendered = stream.getvalue()
    assert "Kernel: Hr68f, 1" in rendered
    assert "default; usable; sharp kernel" in rendered
    assert "sharp kernel Hr68f, 1" not in rendered


def test_list_command_hint_is_shell_quoted_without_hard_wrapping(tmp_path, monkeypatch):
    input_path = tmp_path / "scan folder;$(unsafe)"
    input_path.mkdir()
    monkeypatch.setattr(cli, "_discover", lambda _root: [_candidate()])

    result = runner.invoke(cli.app, ["list", str(input_path)], prog_name="medsurface")

    assert result.exit_code == 0
    expected = _expected_shell_command(
        ["medsurface", "extract", input_path, "-o", "out.stl"]
    )
    assert "Extract the default with:  " + expected in result.stdout


def test_list_file_volume_omits_unavailable_optional_details():
    candidate = VolumeCandidate(
        id=1,
        source=FileSource(Path("scan.nii.gz"), "NIfTI"),
        format="NIfTI",
        source_name="scan.nii.gz",
        modality=None,
        description=None,
        size=None,
        spacing=None,
        direction=None,
        origin=None,
        pixel_type="16-bit signed integer",
        components=1,
        plane="unknown",
    )
    stream = io.StringIO()
    console = Console(file=stream, width=160, color_system=None, force_terminal=False)

    console.print(cli._volume_table([candidate], candidate))

    rendered = stream.getvalue()
    assert "Format: NIfTI" in rendered
    assert "Source: scan.nii.gz" in rendered
    for unavailable in (
        "DICOM #:",
        "Modality:",
        "Description:",
        "Kernel:",
        "Slices:",
        "Voxel (mm):",
        "Plane:",
    ):
        assert unavailable not in rendered


def test_detail_block_omits_empty_values():
    assert (
        cli._detail_block(
            [
                ("None", None),
                ("Empty", ""),
                ("Whitespace", "  "),
                ("Value", "CT"),
            ]
        )
        == "Value: CT"
    )


def test_list_explains_why_mixed_dicom_modalities_have_no_default(
    tmp_path, monkeypatch
):
    ct = _candidate(row_id=1, uid="1.2.3", modality="CT")
    mr = _candidate(row_id=2, uid="1.2.4", modality="MR")
    monkeypatch.setattr(cli, "_discover", lambda _root: [ct, mr])

    result = runner.invoke(cli.app, ["list", str(tmp_path)], prog_name="medsurface")

    assert result.exit_code == 0
    assert "No automatic default" in result.stdout
    assert "multiple modalities: CT, MR" in result.stdout
    expected = _expected_shell_command(
        ["medsurface", "extract", tmp_path, "--volume", "ID", "-o", "out.stl"]
    )
    assert "Choose a volume with:  " + expected in result.stdout


def test_discovery_warnings_are_deduplicated_and_do_not_corrupt_json(
    tmp_path, monkeypatch
):
    from medsurface import catalog

    def fake_discover(_root):
        for _ in range(3):
            warnings.warn("Invalid value for VR UI: broken UID", UserWarning)
        warnings.warn("Invalid value for VR DS: broken spacing", UserWarning)
        return [_candidate()]

    monkeypatch.setattr(catalog, "discover", fake_discover)
    report = tmp_path / "volumes.json"
    result = runner.invoke(
        cli.app,
        ["list", str(tmp_path), "--json", str(report), "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0
    assert json.loads(report.read_text())[0]["id"] == 1
    assert result.stderr.count("Invalid value for VR UI") == 1
    assert "repeated 3 times" in result.stderr
    assert result.stderr.count("Invalid value for VR DS") == 1
    assert ".py:" not in result.stderr
    assert "\x1b" not in result.stdout + result.stderr


def test_series_table_is_responsive_safe_and_ansi_free():
    recommended = _candidate(
        row_id=1,
        description="[bold]recommended acquisition with a deliberately long description[/bold]",
    )
    dose = _candidate(
        row_id=2,
        uid="1.2.4",
        number=6,
        modality="RTDOSE",
        description="radiotherapy dose object with another long description",
    )
    dose.unusable_reason = "not an image series"
    plan = _candidate(
        row_id=3,
        uid="1.2.5",
        number=17,
        modality="RTPLAN",
        description="plan notes",
    )
    assert plan.dicom is not None
    plan.dicom.n_parts = 2
    plan.dicom.part = 2

    for width in (160, 84):
        stream = io.StringIO()
        console = Console(
            file=stream,
            width=width,
            color_system=None,
            force_terminal=False,
            highlight=False,
            markup=False,
        )
        console.print(cli._volume_table([recommended, dose, plan], recommended))
        rendered = stream.getvalue()
        normalized = " ".join(rendered.split())

        assert "ID" in rendered
        assert "DICOM #" in rendered
        assert "RTDOSE" in rendered
        assert "RTPLAN" in rendered
        assert "default" in normalized
        assert "not an" in normalized
        assert "image" in normalized
        assert "series" in normalized
        assert "orientation" in normalized
        assert "[bold]" in normalized
        assert "\x1b" not in rendered
        assert all(len(line) <= width for line in rendered.splitlines())


def test_presets_renders_tissue_table_without_ansi():
    result = runner.invoke(cli.app, ["presets"], prog_name="medsurface")

    assert result.exit_code == 0
    assert "Tissue presets (--preset)" in result.stdout
    assert "bone" in result.stdout
    assert "mask" in result.stdout
    assert "smooth off" in result.stdout
    assert "pre-mesh" in result.stdout
    assert "post-mesh" in result.stdout
    assert "smooth 20 iter" in result.stdout
    assert "smooth 40" in result.stdout
    assert "smooth 25 iter" in result.stdout
    assert "smooth 10 iter" in result.stdout
    assert "\x1b" not in result.stdout


def test_convert_accepts_storage_flags_and_writes_json_file(tmp_path, monkeypatch):
    chosen = _candidate()
    captured = {}
    output = tmp_path / "converted.nii.gz"
    json_file = tmp_path / "conversion.json"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from medsurface import volume

    def fake_convert(**kwargs):
        captured.update(kwargs)
        return _conversion_result(str(output))

    monkeypatch.setattr(volume, "convert", fake_convert)
    result = runner.invoke(
        cli.app,
        [
            "convert",
            str(tmp_path),
            "-o",
            str(output),
            "--volume",
            "1",
            "--strip-metadata",
            "--allow-large-volume",
            "--json",
            str(json_file),
            "-q",
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert result.stdout == ""
    assert result.stderr == ""
    assert captured["candidate"] is chosen
    assert captured["strip_metadata"] is True
    assert captured["allow_large_volume"] is True
    payload = json.loads(json_file.read_text())
    assert payload["result"] == {
        "output": str(output),
        "format": "NIfTI",
        "compression": "gzip",
        "dimensions": [8, 9, 10],
        "pixel_type": "16-bit signed integer",
        "components": 1,
        "spacing_mm": [0.7, 0.8, 0.9],
        "origin_mm": [10.0, 20.0, 30.0],
        "direction": [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        "seconds": 1.5,
        "warnings": [],
        "metadata_policy": "preserved-best-effort",
    }
    assert payload["provenance"]["input"]["format"] == "DICOM"


def test_extract_accepts_all_flags_and_writes_json_file(tmp_path, monkeypatch):
    chosen = _candidate()
    captured = {}
    output = tmp_path / "surface.stl"
    json_file = tmp_path / "result.json"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from medsurface import pipeline

    def fake_extract(**kwargs):
        captured.update(kwargs)
        return _extract_result(str(output))

    monkeypatch.setattr(pipeline, "extract", fake_extract)
    result = runner.invoke(
        cli.app,
        [
            "extract",
            str(tmp_path),
            "-o",
            str(output),
            "--volume",
            "1",
            "--preset",
            "skin",
            "--threshold",
            "auto",
            "--median-mm",
            "1.1",
            "--closing-mm",
            "2.2",
            "--opening-mm",
            "3.3",
            "--min-island-mm3",
            "4.4",
            "--all-islands",
            "--components",
            "all",
            "--resample-mm",
            "0.8",
            "--mask-smooth-mm",
            "0.9",
            "--mesh-smooth-iters",
            "17",
            "--simplify-error-mm",
            "0.18",
            "--post-mesh-smooth-iters",
            "19",
            "--no-cap",
            "--allow-large-volume",
            "--json",
            str(json_file),
            "-q",
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert result.stdout == ""
    assert result.stderr == ""
    assert captured["candidate"] is chosen
    assert captured["preset"].name == "skin"
    assert captured["preset"].threshold == pytest.approx(-300.0)
    assert captured["preset"].median_mm == pytest.approx(1.1)
    assert captured["preset"].opening_mm == pytest.approx(3.3)
    assert captured["preset"].mask_smooth_mm == pytest.approx(0.9)
    assert captured["preset"].surface_smooth_iters == 17
    assert captured["preset"].simplify_error_mm == pytest.approx(0.18)
    assert captured["preset"].post_surface_smooth_iters == 19
    assert captured["preset"].destep is None
    assert not captured["preset"].keep_largest_island
    assert not captured["preset"].keep_largest_component
    assert captured["threshold"] == "auto"
    assert not captured["cap_field_of_view"]
    assert captured["allow_large_volume"]
    payload = json.loads(json_file.read_text())
    assert payload["result"]["output"] == str(output)
    assert payload["quality"]["valid"]


def test_extract_forwards_destep_band_settings(tmp_path, monkeypatch):
    chosen = _candidate()
    captured = {}
    output = tmp_path / "surface.stl"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from medsurface import pipeline
    from medsurface.presets import DestepSettings

    def fake_extract(**kwargs):
        captured.update(kwargs)
        return _extract_result(str(output))

    monkeypatch.setattr(pipeline, "extract", fake_extract)
    result = runner.invoke(
        cli.app,
        [
            "extract",
            str(tmp_path),
            "-o",
            str(output),
            "--destep",
            "band",
            "--destep-axis",
            "y",
            "--destep-full-mm",
            "-555",
            "--destep-frozen-mm",
            "-600",
            "--destep-iters",
            "300",
            "--destep-max-mm",
            "0.6",
            "-q",
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["preset"].destep == DestepSettings(
        region="band",
        iterations=300,
        max_displacement_mm=0.6,
        axis="y",
        full_mm=-555.0,
        frozen_mm=-600.0,
    )


def test_labelmap_extract_forwards_destep_defaults(tmp_path, monkeypatch):
    source = tmp_path / "labels.nii.gz"
    source.write_bytes(b"labelmap")
    chosen = _file_candidate(source)
    output = tmp_path / "surface.stl"
    captured = {}
    monkeypatch.setattr(
        cli, "_select_labelmap", lambda _path, _role=None: (chosen, [chosen])
    )

    from medsurface import defaults
    from medsurface import labelmap as labelmap_mod
    from medsurface.presets import DestepSettings

    def fake_extract(**kwargs):
        captured.update(kwargs)
        return _extract_result(str(output))

    monkeypatch.setattr(labelmap_mod, "extract", fake_extract)
    result = runner.invoke(
        cli.app,
        ["labelmap", "extract", str(source), "-o", str(output), "--destep", "auto"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["destep"] == DestepSettings(
        region="auto",
        iterations=defaults.DEFAULT_DESTEP_ITERS,
        max_displacement_mm=defaults.DEFAULT_DESTEP_MAX_MM,
    )


@pytest.mark.parametrize("command", [["extract"], ["labelmap", "extract"]])
@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["--destep-iters", "100"], "--destep-iters requires --destep"),
        (["--destep-iters", "600"], "--destep-iters requires --destep"),
        (["--destep-full-mm", "1"], "--destep-full-mm requires --destep"),
        (
            ["--destep", "all", "--destep-axis", "z"],
            "--destep-axis requires --destep band",
        ),
        (
            ["--destep", "auto", "--destep-axis", "x"],
            "--destep-axis requires --destep band",
        ),
        (
            ["--destep", "band", "--destep-full-mm", "1"],
            "--destep band requires --destep-full-mm and --destep-frozen-mm",
        ),
        (
            [
                "--destep",
                "band",
                "--destep-full-mm",
                "1",
                "--destep-frozen-mm",
                "1",
            ],
            "must differ",
        ),
        (
            ["--destep-axis", "x", "--destep-max-mm", "1"],
            "--destep-axis, --destep-max-mm require --destep",
        ),
        (["--destep", "all", "--destep-iters", "0"], "--destep-iters must be at least 1"),
        (
            ["--destep", "all", "--destep-max-mm", "0"],
            "--destep-max-mm must be finite and greater than zero",
        ),
        (
            ["--destep", "all", "--destep-max-mm", "nan"],
            "--destep-max-mm must be finite and greater than zero",
        ),
    ],
)
def test_destep_option_errors_are_usage_errors(
    tmp_path, monkeypatch, command, arguments, message
):
    source = tmp_path / "labels.nii.gz"
    source.write_bytes(b"labelmap")
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail("destep options must be checked before discovery"),
    )
    monkeypatch.setattr(
        cli,
        "_select_labelmap",
        lambda *_args: pytest.fail("destep options must be checked before discovery"),
    )

    result = runner.invoke(
        cli.app,
        [*command, str(source), "-o", str(tmp_path / "out.stl"), *arguments],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert message in result.stderr


def test_extract_defaults_quality_output_and_invalid_exit(tmp_path, monkeypatch):
    chosen = _candidate()
    captured = {}
    output = tmp_path / "surface.stl"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from medsurface import pipeline

    def fake_extract(**kwargs):
        captured.update(kwargs)
        return _extract_result(str(output), quality=_quality(valid=False))

    monkeypatch.setattr(pipeline, "extract", fake_extract)
    result = runner.invoke(
        cli.app,
        ["extract", str(tmp_path), "-o", str(output)],
        prog_name="medsurface",
    )

    assert result.exit_code == 1
    assert captured["preset"].name == "bone"
    assert captured["preset"].keep_largest_component
    assert captured["threshold"] is None
    assert captured["cap_field_of_view"]
    assert "Success: wrote" in result.stdout
    assert "Mesh quality" in result.stdout
    assert "invalid" in result.stdout
    assert "Mesh problems" in result.stdout
    assert "failed validation" in result.stderr


def test_extract_omitted_components_preserves_teeth_preset_default(
    tmp_path, monkeypatch
):
    chosen = _candidate()
    captured = {}
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from medsurface import pipeline

    def fake_extract(**kwargs):
        captured.update(kwargs)
        return _extract_result("out.stl")

    monkeypatch.setattr(pipeline, "extract", fake_extract)

    result = runner.invoke(
        cli.app,
        ["extract", str(tmp_path), "-o", "out.stl", "--preset", "teeth", "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["preset"].name == "teeth"
    assert captured["preset"].keep_largest_component is False


def test_valid_quality_output_omits_empty_problems_panel():
    stream = io.StringIO()
    console = Console(
        file=stream,
        width=240,
        color_system=None,
        force_terminal=False,
        highlight=False,
        markup=False,
    )

    cli._print_quality(_quality(), console)

    rendered = stream.getvalue()
    assert "Mesh quality" in rendered
    assert "problems" not in rendered.casefold()


def test_extract_without_selector_keeps_dicom_best_stack_selection(
    tmp_path, monkeypatch
):
    coarse = _candidate(row_id=1, uid="1.2.3", description="coarse")
    fine = _candidate(row_id=2, uid="1.2.4", description="fine")
    assert coarse.dicom is not None and fine.dicom is not None
    coarse.dicom.pixel_spacing = (1.0, 1.0)
    coarse.dicom.slice_spacing = 2.0
    fine.dicom.pixel_spacing = (0.4, 0.4)
    fine.dicom.slice_spacing = 0.8
    monkeypatch.setattr(cli, "_discover", lambda _root: [coarse, fine])

    captured = {}
    from medsurface import pipeline

    def fake_extract(**kwargs):
        captured.update(kwargs)
        return _extract_result("out.stl")

    monkeypatch.setattr(pipeline, "extract", fake_extract)
    result = runner.invoke(
        cli.app,
        ["extract", str(tmp_path), "-o", "out.stl", "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["candidate"] is fine


def test_invalid_extract_threshold_is_exit_two_before_discovery(tmp_path, monkeypatch):
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail("discovery should not run for an invalid threshold"),
    )
    result = runner.invoke(
        cli.app,
        ["extract", str(tmp_path), "-o", "out.stl", "--threshold", "wat"],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert result.stdout == ""
    assert "number or 'auto'" in result.stderr


@pytest.mark.parametrize(
    "argv,option",
    [
        (["extract", "INPUT", "-o", "out.stl", "--median-mm", "-1"], "--median-mm"),
        (
            ["extract", "INPUT", "-o", "out.stl", "--resample-mm", "nan"],
            "--resample-mm",
        ),
        (
            ["extract", "INPUT", "-o", "out.stl", "--mask-smooth-mm", "-1"],
            "--mask-smooth-mm",
        ),
        (
            ["extract", "INPUT", "-o", "out.stl", "--mesh-smooth-iters", "-1"],
            "--mesh-smooth-iters",
        ),
        (
            [
                "extract",
                "INPUT",
                "-o",
                "out.stl",
                "--post-mesh-smooth-iters",
                "-1",
            ],
            "--post-mesh-smooth-iters",
        ),
        (["fuse", "INPUT", "MOVING", "-o", "out.nrrd", "--grid-mm", "0"], "--grid-mm"),
        (
            ["fuse", "INPUT", "MOVING", "-o", "out.nrrd", "--grid-mm", "inf"],
            "--grid-mm",
        ),
        (
            ["fuse", "INPUT", "MOVING", "-o", "out.nrrd", "--opening-mm", "-0.1"],
            "--opening-mm",
        ),
    ],
)
def test_invalid_processing_numbers_are_usage_errors_before_discovery(
    tmp_path, monkeypatch, argv, option
):
    moving = tmp_path / "moving"
    moving.mkdir()
    resolved = [
        str(tmp_path)
        if value == "INPUT"
        else str(moving)
        if value == "MOVING"
        else value
        for value in argv
    ]
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail(
            "discovery should not run for invalid processing options"
        ),
    )

    result = runner.invoke(cli.app, resolved, prog_name="medsurface")

    assert result.exit_code == 2
    assert option in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize("suffix", [".stl", ".dcm", ".nhdr", ".mhd", ".xyz"])
def test_convert_rejects_non_atomic_volume_output_before_discovery(
    tmp_path, monkeypatch, suffix
):
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail("discovery must not run for an unsupported output"),
    )

    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", str(tmp_path / ("out" + suffix))],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "unsupported volume output extension" in result.stderr


@pytest.mark.parametrize("suffix", [".stl", ".dcm", ".nhdr", ".mhd", ".xyz"])
def test_fuse_rejects_non_atomic_volume_output_before_discovery(
    tmp_path, monkeypatch, suffix
):
    moving = tmp_path / "moving"
    moving.mkdir()
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail("discovery must not run for an unsupported output"),
    )

    result = runner.invoke(
        cli.app,
        ["fuse", str(tmp_path), str(moving), "-o", str(tmp_path / ("out" + suffix))],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "unsupported volume output extension" in result.stderr


def test_repair_rejects_output_extension_before_loading(tmp_path, monkeypatch):
    mesh = tmp_path / "input.stl"
    mesh.write_text("placeholder")
    from medsurface import repair as repair_mod

    monkeypatch.setattr(
        repair_mod,
        "repair",
        lambda *_args, **_kwargs: pytest.fail(
            "repair must not start for an unsupported output"
        ),
    )

    result = runner.invoke(
        cli.app,
        ["repair", str(mesh), "-o", str(tmp_path / "out.xyz")],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "unsupported output extension" in result.stderr


def test_convert_rejects_report_aliases_before_processing(tmp_path, monkeypatch):
    source = tmp_path / "scan.mha"
    source.write_bytes(b"original medical image")
    alias = tmp_path / "report.json"
    alias.hardlink_to(source)
    candidate = _file_candidate(source)
    monkeypatch.setattr(cli, "_discover", lambda _root: [candidate])
    from medsurface import volume

    monkeypatch.setattr(
        volume,
        "convert",
        lambda **_kwargs: pytest.fail(
            "conversion must not start for a colliding report path"
        ),
    )

    result = runner.invoke(
        cli.app,
        ["convert", str(source), "-o", str(tmp_path / "out.nii.gz"), "--json", str(alias)],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "JSON report must not overwrite input" in result.stderr
    assert source.read_bytes() == b"original medical image"


def test_convert_rejects_report_overwriting_a_detached_payload(tmp_path, monkeypatch):
    import numpy as np
    import SimpleITK as sitk

    source = tmp_path / "scan.mhd"
    sitk.WriteImage(
        sitk.GetImageFromArray(np.ones((8, 8, 8), dtype=np.int16)), str(source)
    )
    candidate = cli._discover(source)[0]
    assert isinstance(candidate.source, FileSource)
    payload = candidate.source.payload_paths[0]
    original = payload.read_bytes()
    monkeypatch.setattr(cli, "_discover", lambda _root: [candidate])
    from medsurface import volume

    monkeypatch.setattr(
        volume,
        "convert",
        lambda **_kwargs: pytest.fail(
            "conversion must not start for a colliding payload"
        ),
    )

    result = runner.invoke(
        cli.app,
        [
            "convert",
            str(source),
            "-o",
            str(tmp_path / "out.nrrd"),
            "--json",
            str(payload),
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "JSON report must not overwrite input" in result.stderr
    assert payload.read_bytes() == original


def test_convert_rejects_one_path_for_volume_and_json_before_processing(
    tmp_path, monkeypatch
):
    chosen = _candidate()
    destination = tmp_path / "out.nii.gz"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])
    from medsurface import volume

    monkeypatch.setattr(
        volume,
        "convert",
        lambda **_kwargs: pytest.fail(
            "conversion must not start for colliding outputs"
        ),
    )

    result = runner.invoke(
        cli.app,
        [
            "convert",
            str(tmp_path),
            "-o",
            str(destination),
            "--json",
            str(destination),
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "output and JSON report must be different files" in result.stderr
    assert not destination.exists()


def test_selection_error_is_exit_two(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "_discover", lambda _root: [_candidate()])
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", "out.nii.gz", "--volume", "99"],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "no volume has ID" in result.stderr


def test_quiet_suppresses_progress_but_not_warnings(tmp_path, monkeypatch):
    chosen = _candidate()
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from medsurface import volume

    def fake_convert(**kwargs):
        kwargs["warn"]("sampling warning")
        return _conversion_result("out.nii.gz", warnings_=["sampling warning"])

    monkeypatch.setattr(volume, "convert", fake_convert)
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", "out.nii.gz", "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0
    assert result.stdout == ""
    assert "Warning: sampling warning" in result.stderr
    assert result.stderr.count("sampling warning") == 1


def test_labelmap_extract_accepts_surface_flags_and_writes_json(tmp_path, monkeypatch):
    source = tmp_path / "labels.nii.gz"
    source.write_bytes(b"labelmap")
    chosen = _file_candidate(source)
    output = tmp_path / "surface.stl"
    json_file = tmp_path / "result.json"
    captured = {}
    monkeypatch.setattr(
        cli, "_select_labelmap", lambda _path, _role=None: (chosen, [chosen])
    )

    from medsurface import labelmap as labelmap_mod

    def fake_extract(**kwargs):
        captured.update(kwargs)
        result = _extract_result(str(output))
        result.provenance = {
            "input": {"input_kind": "labelmap", "foreground": "all nonzero voxels"}
        }
        return result

    monkeypatch.setattr(labelmap_mod, "extract", fake_extract)
    result = runner.invoke(
        cli.app,
        [
            "labelmap",
            "extract",
            str(source),
            "-o",
            str(output),
            "--resample-mm",
            "0.8",
            "--mask-smooth-mm",
            "0.9",
            "--mesh-smooth-iters",
            "17",
            "--simplify-error-mm",
            "0.18",
            "--post-mesh-smooth-iters",
            "19",
            "--components",
            "largest",
            "--no-cap",
            "--allow-large-volume",
            "--json",
            str(json_file),
            "-q",
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert result.stdout == ""
    assert result.stderr == ""
    assert captured["candidate"] is chosen
    assert captured["resample_mm"] == pytest.approx(0.8)
    assert captured["mask_smooth_mm"] == pytest.approx(0.9)
    assert captured["surface_smooth_iters"] == 17
    assert captured["simplify_error_mm"] == pytest.approx(0.18)
    assert captured["post_surface_smooth_iters"] == 19
    assert captured["destep"] is None
    assert captured["keep_largest_component"] is True
    assert not captured["cap_field_of_view"]
    assert captured["allow_large_volume"]
    payload = json.loads(json_file.read_text())
    assert payload["provenance"]["input"]["input_kind"] == "labelmap"
    assert payload["quality"]["valid"]


def test_labelmap_fuse_accepts_flags_and_writes_generic_json(tmp_path, monkeypatch):
    fixed_path = tmp_path / "fixed.nii.gz"
    moving_path = tmp_path / "moving.nii.gz"
    fixed_path.write_bytes(b"fixed")
    moving_path.write_bytes(b"moving")
    fixed = _file_candidate(fixed_path, row_id=1)
    moving = _file_candidate(moving_path, row_id=1)
    output = tmp_path / "fused.nrrd"
    json_file = tmp_path / "result.json"
    captured = {}

    def select(path, role=None):
        candidate = fixed if path == fixed_path else moving
        return candidate, [candidate]

    monkeypatch.setattr(cli, "_select_labelmap", select)
    from medsurface import labelmap as labelmap_mod

    def fake_fuse(fixed_candidate, movings, output_path, **kwargs):
        captured.update(kwargs, fixed=fixed_candidate, movings=movings)
        kwargs["warn"]("rigid anatomy warning")
        return _label_fusion_result(
            output_path,
            warnings_=["rigid anatomy warning"],
        )

    monkeypatch.setattr(labelmap_mod, "fuse_labels", fake_fuse)
    result = runner.invoke(
        cli.app,
        [
            "labelmap",
            "fuse",
            str(fixed_path),
            str(moving_path),
            "-o",
            str(output),
            "--grid-mm",
            "0.7",
            "--force",
            "--allow-large-volume",
            "--json",
            str(json_file),
            "-q",
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert result.stdout == ""
    assert result.stderr.count("rigid anatomy warning") == 1
    assert captured["fixed"] is fixed
    assert captured["movings"] == [moving]
    assert not captured["preserve_labels"]
    assert captured["grid_mm"] == pytest.approx(0.7)
    assert captured["force"]
    assert captured["allow_large_volume"]
    payload = json.loads(json_file.read_text())
    assert payload["result"]["grid_mm"] == pytest.approx(0.8)
    assert payload["result"]["format"] == "NRRD"
    assert payload["result"]["compression"] == "gzip"
    assert payload["result"]["labels"] == {
        "1": {"name": None, "voxels": 293, "volume_mm3": 150.0}
    }
    assert payload["result"]["volume_fused_mm3"] == pytest.approx(150.0)
    assert payload["result"]["grid_origin_mm"] == [1.0, 2.0, 3.0]
    assert "quality" not in payload
    assert "triangles" not in payload["result"]
    assert "surface_components" not in payload["result"]


def test_labelmap_fuse_prints_shell_safe_extract_hint(tmp_path, monkeypatch):
    fixed_path = tmp_path / "fixed.nii.gz"
    moving_path = tmp_path / "moving.nii.gz"
    fixed_path.write_bytes(b"fixed")
    moving_path.write_bytes(b"moving")
    fixed = _file_candidate(fixed_path, row_id=1)
    moving = _file_candidate(moving_path, row_id=1)
    output = tmp_path / "fused labelmap;$(safe).mha"

    monkeypatch.setattr(
        cli,
        "_select_labelmap",
        lambda path, _role=None: (
            (fixed, [fixed]) if path == fixed_path else (moving, [moving])
        ),
    )
    from medsurface import labelmap as labelmap_mod

    def fake_fuse(_fixed, _movings, output_path, **_kwargs):
        return _label_fusion_result(output_path)

    monkeypatch.setattr(labelmap_mod, "fuse_labels", fake_fuse)
    result = runner.invoke(
        cli.app,
        [
            "labelmap",
            "fuse",
            str(fixed_path),
            str(moving_path),
            "-o",
            str(output),
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert "1 label(s)" in result.stdout
    assert "Mesh quality" not in result.stdout
    assert "triangles" not in result.stdout
    expected = _expected_shell_command(
        ["medsurface", "labelmap", "extract", output, "-o", "MODEL.stl"]
    )
    assert "Extract a surface with:  " + expected in result.stdout


@pytest.mark.parametrize(
    "argv,option",
    [
        (
            [
                "labelmap",
                "extract",
                "INPUT",
                "-o",
                "out.stl",
                "--mask-smooth-mm",
                "-1",
            ],
            "--mask-smooth-mm",
        ),
        (
            [
                "labelmap",
                "fuse",
                "INPUT",
                "MOVING",
                "-o",
                "out.nrrd",
                "--mask-smooth-mm",
                "nan",
            ],
            "--mask-smooth-mm",
        ),
        (
            [
                "labelmap",
                "extract",
                "INPUT",
                "-o",
                "out.stl",
                "--mesh-smooth-iters",
                "-1",
            ],
            "--mesh-smooth-iters",
        ),
        (
            [
                "labelmap",
                "extract",
                "INPUT",
                "-o",
                "out.stl",
                "--post-mesh-smooth-iters",
                "-1",
            ],
            "--post-mesh-smooth-iters",
        ),
        (
            ["labelmap", "fuse", "INPUT", "MOVING", "-o", "out.nrrd", "--grid-mm", "0"],
            "--grid-mm",
        ),
    ],
)
def test_invalid_labelmap_processing_options_fail_before_discovery(
    tmp_path, monkeypatch, argv, option
):
    source = tmp_path / "input.nii.gz"
    moving = tmp_path / "moving.nii.gz"
    source.write_bytes(b"input")
    moving.write_bytes(b"moving")
    resolved = [
        str(source) if value == "INPUT" else str(moving) if value == "MOVING" else value
        for value in argv
    ]
    monkeypatch.setattr(
        cli,
        "_select_labelmap",
        lambda *_args: pytest.fail("labelmap discovery must not start"),
    )

    result = runner.invoke(cli.app, resolved, prog_name="medsurface")

    assert result.exit_code == 2
    assert option in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize(
    "command,option",
    [
        ("fuse", "--mask-smooth-mm"),
        ("fuse", "--mesh-smooth-iters"),
        ("fuse", "--simplify-error-mm"),
        ("fuse", "--post-mesh-smooth-iters"),
        ("fuse", "--components"),
        ("fuse", "--no-cap"),
        ("labelmap fuse", "--mask-smooth-mm"),
        ("labelmap fuse", "--mesh-smooth-iters"),
        ("labelmap fuse", "--simplify-error-mm"),
        ("labelmap fuse", "--components"),
        ("labelmap fuse", "--no-cap"),
    ],
)
def test_fusion_commands_reject_surface_options_before_discovery(
    tmp_path, monkeypatch, command, option
):
    moving = tmp_path / "moving.nii.gz"
    fixed = tmp_path / "fixed.nii.gz"
    fixed.write_bytes(b"fixed")
    moving.write_bytes(b"moving")
    option_value = "all" if option == "--components" else "1"
    prefix = command.split()
    argv = [
        *prefix,
        str(fixed),
        str(moving),
        "-o",
        str(tmp_path / "merged.nii.gz"),
        option,
        option_value,
    ]
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda *_args: pytest.fail("volume discovery must not start"),
    )
    monkeypatch.setattr(
        cli,
        "_select_labelmap",
        lambda *_args: pytest.fail("labelmap discovery must not start"),
    )

    result = runner.invoke(cli.app, argv, prog_name="medsurface")

    assert result.exit_code == 2
    assert option in result.stderr
    assert "No such option" in result.stderr


@pytest.mark.parametrize("command", ["extract", "labelmap extract"])
def test_extract_commands_reject_volume_output_before_discovery(
    tmp_path, monkeypatch, command
):
    source = tmp_path / "input.nii.gz"
    source.write_bytes(b"input")
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda *_args: pytest.fail("volume discovery must not start"),
    )
    monkeypatch.setattr(
        cli,
        "_select_labelmap",
        lambda *_args: pytest.fail("labelmap discovery must not start"),
    )

    result = runner.invoke(
        cli.app,
        [*command.split(), str(source), "-o", str(tmp_path / "out.nii.gz")],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "unsupported output extension '.nii.gz'" in result.stderr


@pytest.mark.parametrize("command", ["validate", "repair"])
def test_mesh_tools_reject_nifti_input_before_loading(tmp_path, command):
    source = tmp_path / "input.nii.gz"
    source.write_bytes(b"input")
    argv = [command, str(source)]
    if command == "repair":
        argv.extend(["-o", str(tmp_path / "out.stl")])

    result = runner.invoke(cli.app, argv, prog_name="medsurface")

    assert result.exit_code == 2
    assert "unsupported mesh extension '.nii.gz'" in result.stderr


def test_convert_requires_id_for_multiple_dicom_modalities(tmp_path, monkeypatch):
    ct = _candidate(row_id=1, uid="1.2.3", modality="CT")
    mr = _candidate(row_id=2, uid="1.2.4", modality="MR")
    monkeypatch.setattr(cli, "_discover", lambda _root: [ct, mr])

    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", str(tmp_path / "out.nii.gz")],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "multiple modalities: CT, MR" in result.stderr
    assert "volume ID" in result.stderr


def test_fuse_identifies_the_input_with_ambiguous_dicom_modalities(
    tmp_path, monkeypatch
):
    moving_dir = tmp_path / "moving"
    moving_dir.mkdir()
    fixed = _candidate(row_id=1, uid="1.2.3", modality="CT")
    moving_ct = _candidate(row_id=1, uid="1.3.3", modality="CT")
    moving_mr = _candidate(row_id=2, uid="1.3.4", modality="MR")

    def fake_discover(path: Path):
        return [moving_ct, moving_mr] if path == moving_dir else [fixed]

    monkeypatch.setattr(cli, "_discover", fake_discover)
    result = runner.invoke(
        cli.app,
        ["fuse", str(tmp_path), str(moving_dir), "-o", str(tmp_path / "out.nrrd")],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "moving input" in result.stderr
    assert "multiple modalities: CT, MR" in result.stderr


def test_fuse_accepts_all_flags_and_safety_errors_exit_three(tmp_path, monkeypatch):
    directory_b = tmp_path / "moving"
    directory_b.mkdir()
    fixed = _candidate(uid="1.2.3", description="fixed")
    moving = _candidate(uid="1.2.4", description="moving")
    captured = {}
    output = tmp_path / "fused.mha"
    json_file = tmp_path / "fused.json"

    def fake_discover(path: Path):
        return [moving] if path == directory_b else [fixed]

    monkeypatch.setattr(cli, "_discover", fake_discover)
    from medsurface import fusion

    def fake_fuse(**kwargs):
        captured.update(kwargs)
        return _fusion_result(str(output))

    monkeypatch.setattr(fusion, "fuse", fake_fuse)
    result = runner.invoke(
        cli.app,
        [
            "fuse",
            str(tmp_path),
            str(directory_b),
            "-o",
            str(output),
            "--fixed-volume",
            "1",
            "--moving-volume",
            "1",
            "--preset",
            "teeth",
            "--fixed-threshold",
            "250",
            "--moving-threshold",
            "auto",
            "--median-mm",
            "1.1",
            "--closing-mm",
            "2.2",
            "--opening-mm",
            "1.7",
            "--min-island-mm3",
            "3.3",
            "--all-islands",
            "--grid-mm",
            "0.8",
            "--force",
            "--allow-large-volume",
            "--json",
            str(json_file),
            "-q",
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["fixed"] is fixed
    assert captured["moving"] is moving
    assert captured["preset"].name == "teeth"
    assert captured["preset"].opening_mm == pytest.approx(1.7)
    assert not captured["preset"].keep_largest_island
    assert captured["fixed_threshold"] == pytest.approx(250.0)
    assert captured["moving_threshold"] == "auto"
    assert captured["grid_mm"] == pytest.approx(0.8)
    assert captured["force"]
    assert captured["allow_large_volume"]
    payload = json.loads(json_file.read_text())
    assert payload["result"]["grid_size"] == [10, 20, 30]
    assert payload["result"]["foreground_fused_voxels"] == 293
    assert "quality" not in payload
    assert "surface_components" not in payload["result"]

    def refuse(**_kwargs):
        raise fusion.FusionError("registration gate refused the pair")

    monkeypatch.setattr(fusion, "fuse", refuse)
    refused = runner.invoke(
        cli.app,
        ["fuse", str(tmp_path), str(directory_b), "-o", str(output)],
        prog_name="medsurface",
    )
    assert refused.exit_code == 3
    assert "registration gate refused" in refused.stderr


def test_fuse_writes_generic_volume_result_without_mesh_fields(tmp_path, monkeypatch):
    moving_dir = tmp_path / "moving"
    moving_dir.mkdir()
    fixed = _candidate(uid="1.2.3", description="fixed")
    moving = _candidate(uid="1.2.4", description="moving")
    output = tmp_path / "fused.nii"
    json_file = tmp_path / "fused.json"
    captured = {}

    monkeypatch.setattr(
        cli,
        "_discover",
        lambda path: [moving] if path == moving_dir else [fixed],
    )
    from medsurface import fusion

    def fake_fuse(**kwargs):
        captured.update(kwargs)
        result = _fusion_result(str(output))
        result.output_format = "NIfTI"
        result.compression = "none"
        return result

    monkeypatch.setattr(fusion, "fuse", fake_fuse)
    result = runner.invoke(
        cli.app,
        [
            "fuse",
            str(tmp_path),
            str(moving_dir),
            "-o",
            str(output),
            "--preset",
            "bone",
            "--grid-mm",
            "0.8",
            "--json",
            str(json_file),
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["preset"].name == "bone"
    assert "foreground" in result.stdout
    assert "Mesh quality" not in result.stdout
    payload = json.loads(json_file.read_text())
    assert payload["result"]["format"] == "NIfTI"
    assert payload["result"]["compression"] == "none"
    assert payload["result"]["grid_size"] == [10, 20, 30]
    assert "quality" not in payload
    assert "surface_components" not in payload["result"]
    assert "foreground_voxels" not in payload["result"]
    assert "foreground_volume_mm3" not in payload["result"]


def test_fuse_does_not_override_surface_fields_on_the_segmentation_preset(
    tmp_path, monkeypatch
):
    moving_dir = tmp_path / "moving"
    moving_dir.mkdir()
    fixed = _candidate(uid="1.2.3", description="fixed")
    moving = _candidate(uid="1.2.4", description="moving")
    captured = {}

    def fake_discover(path: Path):
        return [moving] if path == moving_dir else [fixed]

    monkeypatch.setattr(cli, "_discover", fake_discover)
    from medsurface import fusion

    def fake_fuse(**kwargs):
        captured.update(kwargs)
        return _fusion_result("out.nrrd")

    monkeypatch.setattr(fusion, "fuse", fake_fuse)

    result = runner.invoke(
        cli.app,
        [
            "fuse",
            str(tmp_path),
            str(moving_dir),
            "-o",
            "out.nrrd",
            "--preset",
            "teeth",
            "-q",
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["preset"].name == "teeth"
    assert captured["preset"].keep_largest_component is False


def test_fuse_rejects_json_overwriting_a_dicom_instance(tmp_path, monkeypatch):
    moving_dir = tmp_path / "moving"
    moving_dir.mkdir()
    dicom_instance = tmp_path / "slice-001.dcm"
    dicom_instance.write_bytes(b"original dicom")
    fixed = _candidate(uid="1.2.3", description="fixed")
    moving = _candidate(uid="1.2.4", description="moving")
    assert fixed.dicom is not None
    fixed.dicom.files[0] = str(dicom_instance)

    monkeypatch.setattr(
        cli,
        "_discover",
        lambda path: [moving] if path == moving_dir else [fixed],
    )
    from medsurface import fusion

    monkeypatch.setattr(
        fusion,
        "fuse",
        lambda **_kwargs: pytest.fail(
            "fusion must not start for a colliding report path"
        ),
    )

    result = runner.invoke(
        cli.app,
        [
            "fuse",
            str(tmp_path),
            str(moving_dir),
            "-o",
            str(tmp_path / "out.nrrd"),
            "--json",
            str(dicom_instance),
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "JSON report must not overwrite input" in result.stderr
    assert dicom_instance.read_bytes() == b"original dicom"


def test_json_file_is_published_atomically(tmp_path, monkeypatch):
    destination = tmp_path / "report.json"
    destination.write_text("original", encoding="utf-8")

    def fail(*_args, **_kwargs):
        raise OSError("simulated write failure")

    monkeypatch.setattr(cli.json, "dump", fail)

    with pytest.raises(OSError, match="simulated"):
        cli._write_json_file(destination, {"result": "new"})

    assert destination.read_text(encoding="utf-8") == "original"
    assert not list(tmp_path.glob(".report.json.*"))


@pytest.mark.parametrize("hard_links_available", [True, False])
def test_report_failure_restores_the_existing_volume_output(
    tmp_path, monkeypatch, hard_links_available
):
    chosen = _candidate()
    output = tmp_path / "converted.nii.gz"
    report = tmp_path / "report.json"
    output.write_bytes(b"original volume")
    report.write_bytes(b"original report")
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])
    if not hard_links_available:
        monkeypatch.setattr(
            cli.os,
            "link",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                OSError("hard links unavailable")
            ),
        )

    from medsurface import volume

    def fake_convert(**kwargs):
        destination = Path(kwargs["output_path"])
        assert destination.read_bytes() == b"original volume"
        temporary = destination.with_name(".new-volume.nii.gz")
        temporary.write_bytes(b"new volume")
        os.replace(temporary, destination)
        return _conversion_result(kwargs["output_path"])

    monkeypatch.setattr(volume, "convert", fake_convert)
    monkeypatch.setattr(
        cli,
        "_write_json_file",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            OSError("simulated report failure")
        ),
    )

    result = runner.invoke(
        cli.app,
        [
            "convert",
            str(tmp_path),
            "-o",
            str(output),
            "--json",
            str(report),
            "-q",
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 1
    assert "simulated report failure" in result.stderr
    assert output.read_bytes() == b"original volume"
    assert report.read_bytes() == b"original report"
    assert not list(tmp_path.glob(".converted.nii.gz.rollback.*"))


def test_fuse_without_selectors_ranks_each_dicom_directory(tmp_path, monkeypatch):
    moving_dir = tmp_path / "moving"
    moving_dir.mkdir()
    fixed_coarse = _candidate(row_id=1, uid="1.2.1", description="fixed coarse")
    fixed_fine = _candidate(row_id=2, uid="1.2.2", description="fixed fine")
    moving_coarse = _candidate(row_id=1, uid="1.3.1", description="moving coarse")
    moving_fine = _candidate(row_id=2, uid="1.3.2", description="moving fine")
    for candidate, pixel_spacing, slice_spacing in (
        (fixed_coarse, (1.0, 1.0), 2.0),
        (fixed_fine, (0.4, 0.4), 0.8),
        (moving_coarse, (1.2, 1.2), 2.0),
        (moving_fine, (0.5, 0.5), 0.7),
    ):
        assert candidate.dicom is not None
        candidate.dicom.pixel_spacing = pixel_spacing
        candidate.dicom.slice_spacing = slice_spacing

    def fake_discover(path: Path):
        if path == moving_dir:
            return [moving_coarse, moving_fine]
        return [fixed_coarse, fixed_fine]

    monkeypatch.setattr(cli, "_discover", fake_discover)
    captured = {}
    from medsurface import fusion

    def fake_fuse(**kwargs):
        captured.update(kwargs)
        return _fusion_result("fused.nrrd")

    monkeypatch.setattr(fusion, "fuse", fake_fuse)
    result = runner.invoke(
        cli.app,
        ["fuse", str(tmp_path), str(moving_dir), "-o", "fused.nrrd", "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["fixed"] is fixed_fine
    assert captured["moving"] is moving_fine


def test_validate_json_writes_a_file_and_invalid_quality_exits_one(
    tmp_path, monkeypatch
):
    mesh = tmp_path / "mesh.stl"
    mesh.write_text("placeholder")
    from medsurface import validate as validate_mod

    monkeypatch.setattr(validate_mod, "validate", lambda _path: _quality(valid=False))
    report = tmp_path / "quality.json"
    result = runner.invoke(
        cli.app,
        ["validate", str(mesh), "--json", str(report), "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 1
    assert json.loads(report.read_text())["valid"] is False
    assert result.stdout == ""
    assert result.stderr == ""

    human = runner.invoke(
        cli.app,
        ["validate", str(mesh)],
        prog_name="medsurface",
    )
    assert human.exit_code == 1
    assert "Loading validation engine ..." in human.stdout
    assert "Validating mesh structure and self-intersections ..." in human.stdout
    assert "Mesh quality" in human.stdout


def test_repair_json_writes_a_file_and_errors_are_concise(tmp_path, monkeypatch):
    mesh = tmp_path / "mesh.stl"
    output = tmp_path / "fixed.stl"
    mesh.write_text("placeholder")
    from medsurface import repair as repair_mod

    monkeypatch.setattr(
        repair_mod,
        "repair",
        lambda *_args, **_kwargs: repair_mod.RepairResult(
            stats={"holes_filled": 1},
            quality=_quality(),
        ),
    )
    report = tmp_path / "repair.json"
    result = runner.invoke(
        cli.app,
        ["repair", str(mesh), "-o", str(output), "--json", str(report), "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0
    assert json.loads(report.read_text())["repair"]["holes_filled"] == 1
    assert result.stdout == ""
    assert result.stderr == ""

    def fail(*_args, **_kwargs):
        raise ValueError("bad mesh syntax")

    monkeypatch.setattr(repair_mod, "repair", fail)
    failed = runner.invoke(
        cli.app,
        ["repair", str(mesh), "-o", str(output)],
        prog_name="medsurface",
    )
    assert failed.exit_code == 1
    assert "Loading repair engine ..." in failed.stdout
    assert "Error: cannot repair" in failed.stderr
    assert "Traceback" not in failed.stderr
