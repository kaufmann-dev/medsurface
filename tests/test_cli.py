"""Cross-command Typer and Rich CLI contracts."""

from __future__ import annotations

import io
import json
import os
import pty
import select
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

runner = CliRunner()


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
) -> Series:
    value = Series(
        uid=uid,
        modality=modality,
        description=description,
        series_number=number,
        id=row_id,
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


def _convert_result(output: str, *, warnings_: list[str] | None = None, quality: dict | None = None):
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


def _merge_result(output: str, *, warnings_: list[str] | None = None, quality: dict | None = None):
    return SimpleNamespace(
        output_path=output,
        triangles=1400,
        vertices=702,
        bounds_mm=(0.0, 1.0, 0.0, 2.0, 0.0, 3.0),
        grid_mm=0.8,
        grid_size=(10, 20, 30),
        volume_fixed_mm3=100.0,
        volume_moving_mm3=110.0,
        volume_union_mm3=150.0,
        surface_components=2,
        seconds=2.5,
        warnings=warnings_ or [],
        provenance={"fixed": {"format": "DICOM"}, "moving": {"format": "DICOM"}},
        quality=_quality() if quality is None else quality,
    )


def _run_cli_in_clean_interpreter(argv):
    code = """
import json
import sys
from typer.testing import CliRunner
from medsurface import cli

heavy_modules = {
    "medsurface.merge",
    "medsurface.pipeline",
    "medsurface.registration",
    "medsurface.repair",
    "medsurface.surface",
    "medsurface.validate",
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
    for command in ("list", "presets", "convert", "merge", "validate", "repair"):
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
        ["merge", "--help"],
        ["validate", "--help"],
        ["repair", "--help"],
    ],
)
def test_every_help_surface_is_lightweight(argv):
    result = _run_cli_in_clean_interpreter(argv)

    assert result.returncode == 0
    assert "Usage: medsurface" in result.stdout
    assert result.stderr == ""
    assert "\x1b" not in result.stdout
    assert "Traceback" not in result.stdout


def test_merge_help_describes_shared_processing_options():
    result = _run_cli_in_clean_interpreter(["merge", "--help"])

    assert result.returncode == 0
    for option, description_start in (
        ("--median-mm", "Despeckle kernel"),
        ("--closing-mm", "Pore-sealing kernel"),
        ("--opening-mm", "Bridge-breaking kernel"),
        ("--min-island-mm3", "Drop blobs smaller"),
        ("--smooth-iters", "MeshLib relaxation"),
        ("--smooth-force", "MeshLib relaxation"),
        ("--post-smooth-iters", "Smoothing after"),
    ):
        assert option in result.stdout
        assert description_start in result.stdout


@pytest.mark.parametrize(
    "argv",
    [
        ["unknown"],
        ["list"],
        ["presets", "unexpected"],
        ["convert", "scan"],
        ["convert", "scan", "-o", "out.stl", "--preset", "unknown"],
        ["merge"],
        ["validate"],
        ["repair", "broken.stl"],
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
        ["convert", str(tmp_path), "-o", "out.stl", "--preset", "bone-detail"],
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
    with cli._ProgressDisplay(True, "Starting [literal] ...", plain_console) as progress:
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
        assert any(isinstance(column, cli.TimeElapsedColumn) for column in progress.progress.columns)
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


@pytest.mark.parametrize(
    "argv",
    [
        ["convert", ".", "-o", "out.stl", "--self-intersections"],
        ["merge", ".", ".", "-o", "out.stl", "--self-intersections"],
        ["validate", __file__, "--self-intersections"],
        ["repair", __file__, "-o", "fixed.stl", "--self-intersections"],
        ["convert", ".", "-o", "out.stl", "--target-faces", "100"],
        ["merge", ".", ".", "-o", "out.stl", "--target-faces", "100"],
        ["convert", ".", "-o", "out.stl", "--no-validate"],
        ["merge", ".", ".", "-o", "out.stl", "--no-validate"],
        ["convert", ".", "-o", "out.stl", "--passband", "0.1"],
        ["merge", ".", ".", "-o", "out.stl", "--passband", "0.1"],
        ["convert", ".", "-o", "out.stl", "--print-profile", "resin"],
        ["merge", ".", ".", "-o", "out.stl", "--print-profile", "fdm"],
        ["convert", ".", "-o", "out.stl", "--min-feature-mm", "0.6"],
        ["merge", ".", ".", "-o", "out.stl", "--min-feature-mm", "1.2"],
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
        ("convert", ["convert", "missing", "-o", "out.stl"]),
        ("merge", ["merge", "missing", "also-missing", "-o", "out.stl"]),
        ("validate", ["validate", "missing.stl"]),
        ("repair", ["repair", "missing.stl", "-o", "fixed.stl"]),
    ],
)
def test_required_inputs_use_typer_path_validation(command, argv):
    result = runner.invoke(cli.app, argv, prog_name="medsurface")

    assert result.exit_code == 2
    assert "does not exist" in result.stderr
    assert command in result.stderr


def test_list_json_uses_unique_id_and_plain_stdout(tmp_path, monkeypatch):
    found = [_candidate(description="[bold red]literal[/bold red]")]
    monkeypatch.setattr(cli, "_discover", lambda _root: found)

    result = runner.invoke(cli.app, ["list", str(tmp_path), "--json"], prog_name="medsurface")

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload[0]["id"] == 1
    assert payload[0]["default"] is True
    assert payload[0]["format"] == "DICOM"
    assert payload[0]["dicom"]["series_number"] == 6
    assert payload[0]["dicom"]["part"] == 1
    assert payload[0]["dicom"]["n_parts"] == 1
    assert payload[0]["description"] == "[bold red]literal[/bold red]"
    assert "ident" not in payload[0]
    assert "\x1b" not in result.stdout
    assert result.stderr == ""


def test_list_human_output_shows_discovery_progress(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "_discover", lambda _root: [_candidate()])

    result = runner.invoke(cli.app, ["list", str(tmp_path)], prog_name="medsurface")

    assert result.exit_code == 0
    assert "Discovering volumes ..." in result.stdout
    assert "DICOM #" in result.stdout
    assert "axial" in result.stdout


def test_list_file_volume_uses_dash_for_missing_metadata():
    candidate = VolumeCandidate(
        id=1,
        source=FileSource(Path("scan.nii.gz"), "NIfTI"),
        format="NIfTI",
        source_name="scan.nii.gz",
        modality=None,
        description=None,
        size=(10, 20, 30),
        spacing=(0.5, 0.5, 1.0),
        direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        origin=(0.0, 0.0, 0.0),
        pixel_type="16-bit signed integer",
        components=1,
        plane="axial",
    )
    stream = io.StringIO()
    console = Console(file=stream, width=160, color_system=None, force_terminal=False)

    console.print(cli._volume_table([candidate], candidate))

    rendered = stream.getvalue()
    assert "Format: NIfTI" in rendered
    assert "DICOM #: -" in rendered
    assert "Modality: -" in rendered
    assert "Description: -" in rendered
    assert "Plane: axial" in rendered


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
    assert "Choose a volume with" in result.stdout


def test_discovery_warnings_are_deduplicated_and_do_not_corrupt_json(tmp_path, monkeypatch):
    from medsurface import catalog

    def fake_discover(_root):
        for _ in range(3):
            warnings.warn("Invalid value for VR UI: broken UID", UserWarning)
        warnings.warn("Invalid value for VR DS: broken spacing", UserWarning)
        return [_candidate()]

    monkeypatch.setattr(catalog, "discover", fake_discover)
    result = runner.invoke(cli.app, ["list", str(tmp_path), "--json"], prog_name="medsurface")

    assert result.exit_code == 0
    assert json.loads(result.stdout)[0]["id"] == 1
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
    assert "\x1b" not in result.stdout


def test_convert_accepts_all_flags_and_writes_json_file(tmp_path, monkeypatch):
    chosen = _candidate()
    captured = {}
    output = tmp_path / "surface.stl"
    json_file = tmp_path / "result.json"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from medsurface import pipeline

    def fake_convert(**kwargs):
        captured.update(kwargs)
        return _convert_result(str(output))

    monkeypatch.setattr(pipeline, "convert", fake_convert)
    result = runner.invoke(
        cli.app,
        [
            "convert",
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
            "--all-components",
            "--resample-mm",
            "0.8",
            "--smooth-iters",
            "11",
            "--smooth-force",
            "0.12",
            "--simplify-error-mm",
            "0.18",
            "--post-smooth-iters",
            "9",
            "--no-cap",
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
    assert captured["preset"].simplify_error_mm == pytest.approx(0.18)
    assert not captured["preset"].keep_largest_island
    assert not captured["preset"].keep_largest_component
    assert captured["threshold"] == "auto"
    assert not captured["cap_field_of_view"]
    payload = json.loads(json_file.read_text())
    assert payload["result"]["output"] == str(output)
    assert payload["quality"]["valid"]


def test_convert_defaults_quality_output_and_invalid_exit(tmp_path, monkeypatch):
    chosen = _candidate()
    captured = {}
    output = tmp_path / "surface.stl"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from medsurface import pipeline

    def fake_convert(**kwargs):
        captured.update(kwargs)
        return _convert_result(str(output), quality=_quality(valid=False))

    monkeypatch.setattr(pipeline, "convert", fake_convert)
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", str(output)],
        prog_name="medsurface",
    )

    assert result.exit_code == 1
    assert captured["preset"].name == "bone"
    assert captured["threshold"] is None
    assert captured["cap_field_of_view"]
    assert "Success: wrote" in result.stdout
    assert "Mesh quality" in result.stdout
    assert "invalid" in result.stdout
    assert "Mesh problems" in result.stdout
    assert "failed validation" in result.stderr


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


def test_convert_without_selector_keeps_dicom_best_stack_selection(tmp_path, monkeypatch):
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

    def fake_convert(**kwargs):
        captured.update(kwargs)
        return _convert_result("out.stl")

    monkeypatch.setattr(pipeline, "convert", fake_convert)
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", "out.stl", "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["candidate"] is fine


def test_invalid_threshold_is_exit_two_before_discovery(tmp_path, monkeypatch):
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail("discovery should not run for an invalid threshold"),
    )
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", "out.stl", "--threshold", "wat"],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert result.stdout == ""
    assert "number or 'auto'" in result.stderr


@pytest.mark.parametrize(
    "argv,option",
    [
        (["convert", "INPUT", "-o", "out.stl", "--median-mm", "-1"], "--median-mm"),
        (["convert", "INPUT", "-o", "out.stl", "--resample-mm", "nan"], "--resample-mm"),
        (["convert", "INPUT", "-o", "out.stl", "--smooth-iters", "-1"], "--smooth-iters"),
        (["convert", "INPUT", "-o", "out.stl", "--smooth-force", "2"], "--smooth-force"),
        (["merge", "INPUT", "MOVING", "-o", "out.stl", "--grid-mm", "0"], "--grid-mm"),
        (["merge", "INPUT", "MOVING", "-o", "out.stl", "--grid-mm", "inf"], "--grid-mm"),
        (["merge", "INPUT", "MOVING", "-o", "out.stl", "--opening-mm", "-0.1"], "--opening-mm"),
    ],
)
def test_invalid_processing_numbers_are_usage_errors_before_discovery(
    tmp_path, monkeypatch, argv, option
):
    moving = tmp_path / "moving"
    moving.mkdir()
    resolved = [
        str(tmp_path) if value == "INPUT" else str(moving) if value == "MOVING" else value
        for value in argv
    ]
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail("discovery should not run for invalid processing options"),
    )

    result = runner.invoke(cli.app, resolved, prog_name="medsurface")

    assert result.exit_code == 2
    assert option in result.stderr
    assert "Traceback" not in result.stderr


def test_convert_rejects_output_extension_before_discovery(tmp_path, monkeypatch):
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail("discovery must not run for an unsupported output"),
    )

    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", str(tmp_path / "out.xyz")],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "unsupported output extension" in result.stderr


def test_merge_rejects_output_extension_before_discovery(tmp_path, monkeypatch):
    moving = tmp_path / "moving"
    moving.mkdir()
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail("discovery must not run for an unsupported output"),
    )

    result = runner.invoke(
        cli.app,
        ["merge", str(tmp_path), str(moving), "-o", str(tmp_path / "out.xyz")],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "unsupported output extension" in result.stderr


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
    from medsurface import pipeline

    monkeypatch.setattr(
        pipeline,
        "convert",
        lambda **_kwargs: pytest.fail("conversion must not start for a colliding report path"),
    )

    result = runner.invoke(
        cli.app,
        ["convert", str(source), "-o", str(tmp_path / "out.stl"), "--json", str(alias)],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "JSON report must not overwrite input" in result.stderr
    assert source.read_bytes() == b"original medical image"


def test_convert_rejects_report_overwriting_a_detached_payload(tmp_path, monkeypatch):
    import numpy as np
    import SimpleITK as sitk

    source = tmp_path / "scan.mhd"
    sitk.WriteImage(sitk.GetImageFromArray(np.ones((8, 8, 8), dtype=np.int16)), str(source))
    candidate = cli._discover(source)[0]
    assert isinstance(candidate.source, FileSource)
    payload = candidate.source.payload_paths[0]
    original = payload.read_bytes()
    monkeypatch.setattr(cli, "_discover", lambda _root: [candidate])
    from medsurface import pipeline

    monkeypatch.setattr(
        pipeline,
        "convert",
        lambda **_kwargs: pytest.fail("conversion must not start for a colliding payload"),
    )

    result = runner.invoke(
        cli.app,
        ["convert", str(source), "-o", str(tmp_path / "out.stl"), "--json", str(payload)],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "JSON report must not overwrite input" in result.stderr
    assert payload.read_bytes() == original


def test_convert_rejects_one_path_for_mesh_and_json_before_processing(tmp_path, monkeypatch):
    chosen = _candidate()
    destination = tmp_path / "out.stl"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])
    from medsurface import pipeline

    monkeypatch.setattr(
        pipeline,
        "convert",
        lambda **_kwargs: pytest.fail("conversion must not start for colliding outputs"),
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
    assert "mesh output and JSON report must be different files" in result.stderr
    assert not destination.exists()


def test_selection_error_is_exit_two(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "_discover", lambda _root: [_candidate()])
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", "out.stl", "--volume", "99"],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "no volume has ID" in result.stderr


def test_quiet_suppresses_progress_but_not_warnings(tmp_path, monkeypatch):
    chosen = _candidate()
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from medsurface import pipeline

    def fake_convert(**kwargs):
        kwargs["warn"]("sampling warning")
        return _convert_result("out.stl", warnings_=["sampling warning"])

    monkeypatch.setattr(pipeline, "convert", fake_convert)
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", "out.stl", "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0
    assert result.stdout == ""
    assert "Warning: sampling warning" in result.stderr
    assert result.stderr.count("sampling warning") == 1


def test_convert_requires_id_for_multiple_dicom_modalities(tmp_path, monkeypatch):
    ct = _candidate(row_id=1, uid="1.2.3", modality="CT")
    mr = _candidate(row_id=2, uid="1.2.4", modality="MR")
    monkeypatch.setattr(cli, "_discover", lambda _root: [ct, mr])

    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", str(tmp_path / "out.stl")],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "multiple modalities: CT, MR" in result.stderr
    assert "volume ID" in result.stderr


def test_merge_identifies_the_input_with_ambiguous_dicom_modalities(
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
        ["merge", str(tmp_path), str(moving_dir), "-o", str(tmp_path / "out.stl")],
        prog_name="medsurface",
    )

    assert result.exit_code == 2
    assert "moving input" in result.stderr
    assert "multiple modalities: CT, MR" in result.stderr


def test_merge_accepts_all_flags_and_safety_errors_exit_three(tmp_path, monkeypatch):
    directory_b = tmp_path / "moving"
    directory_b.mkdir()
    fixed = _candidate(uid="1.2.3", description="fixed")
    moving = _candidate(uid="1.2.4", description="moving")
    captured = {}
    output = tmp_path / "merged.stl"
    json_file = tmp_path / "merged.json"

    def fake_discover(path: Path):
        return [moving] if path == directory_b else [fixed]

    monkeypatch.setattr(cli, "_discover", fake_discover)
    from medsurface import merge as merge_mod

    def fake_merge(**kwargs):
        captured.update(kwargs)
        return _merge_result(str(output))

    monkeypatch.setattr(merge_mod, "merge", fake_merge)
    result = runner.invoke(
        cli.app,
        [
            "merge",
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
            "--all-components",
            "--grid-mm",
            "0.8",
            "--smooth-iters",
            "12",
            "--smooth-force",
            "0.1",
            "--simplify-error-mm",
            "0.2",
            "--post-smooth-iters",
            "8",
            "--force",
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
    assert not captured["preset"].keep_largest_component
    assert captured["fixed_threshold"] == pytest.approx(250.0)
    assert captured["moving_threshold"] == "auto"
    assert captured["grid_mm"] == pytest.approx(0.8)
    assert captured["smooth_iters"] == 12
    assert captured["smooth_force"] == pytest.approx(0.1)
    assert captured["simplify_error_mm"] == pytest.approx(0.2)
    assert captured["post_smooth_iters"] == 8
    assert captured["force"]
    payload = json.loads(json_file.read_text())
    assert payload["result"]["grid_size"] == [10, 20, 30]
    assert payload["result"]["surface_components"] == 2

    def refuse(**_kwargs):
        raise merge_mod.MergeError("registration gate refused the pair")

    monkeypatch.setattr(merge_mod, "merge", refuse)
    refused = runner.invoke(
        cli.app,
        ["merge", str(tmp_path), str(directory_b), "-o", str(output)],
        prog_name="medsurface",
    )
    assert refused.exit_code == 3
    assert "registration gate refused" in refused.stderr


def test_merge_rejects_json_overwriting_a_dicom_instance(tmp_path, monkeypatch):
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
    from medsurface import merge as merge_mod

    monkeypatch.setattr(
        merge_mod,
        "merge",
        lambda **_kwargs: pytest.fail("merge must not start for a colliding report path"),
    )

    result = runner.invoke(
        cli.app,
        [
            "merge",
            str(tmp_path),
            str(moving_dir),
            "-o",
            str(tmp_path / "out.stl"),
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


def test_merge_without_selectors_ranks_each_dicom_directory(tmp_path, monkeypatch):
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
    from medsurface import merge as merge_mod

    def fake_merge(**kwargs):
        captured.update(kwargs)
        return _merge_result("merged.stl")

    monkeypatch.setattr(merge_mod, "merge", fake_merge)
    result = runner.invoke(
        cli.app,
        ["merge", str(tmp_path), str(moving_dir), "-o", "merged.stl", "-q"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert captured["fixed"] is fixed_fine
    assert captured["moving"] is moving_fine


def test_validate_json_is_plain_and_invalid_quality_exits_one(tmp_path, monkeypatch):
    mesh = tmp_path / "mesh.stl"
    mesh.write_text("placeholder")
    from medsurface import validate as validate_mod

    monkeypatch.setattr(validate_mod, "validate", lambda _path: _quality(valid=False))
    result = runner.invoke(
        cli.app,
        ["validate", str(mesh), "--json"],
        prog_name="medsurface",
    )

    assert result.exit_code == 1
    assert json.loads(result.stdout)["valid"] is False
    assert result.stderr == ""
    assert "\x1b" not in result.stdout

    human = runner.invoke(
        cli.app,
        ["validate", str(mesh)],
        prog_name="medsurface",
    )
    assert human.exit_code == 1
    assert "Loading validation engine ..." in human.stdout
    assert "Validating mesh structure and self-intersections ..." in human.stdout
    assert "Mesh quality" in human.stdout


def test_repair_json_is_plain_and_errors_are_concise(tmp_path, monkeypatch):
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
    result = runner.invoke(
        cli.app,
        ["repair", str(mesh), "-o", str(output), "--json"],
        prog_name="medsurface",
    )

    assert result.exit_code == 0
    assert json.loads(result.stdout)["repair"]["holes_filled"] == 1
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
    assert "Loading mesh for repair ..." in failed.stdout
    assert "Error: cannot repair" in failed.stderr
    assert "Traceback" not in failed.stderr
