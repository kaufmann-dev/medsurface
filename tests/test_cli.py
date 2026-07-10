"""Cross-command Typer and Rich CLI contracts."""

from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest
from rich.console import Console
from typer.testing import CliRunner

from dicom_surface import cli
from dicom_surface.series import Series


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
        "nonmanifold_edge_uses": 0,
        "degenerate_faces": 0,
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
        provenance={"series_uid": "1.2.3"},
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
        volume_a_mm3=100.0,
        volume_b_mm3=110.0,
        volume_union_mm3=150.0,
        seconds=2.5,
        warnings=warnings_ or [],
        provenance={"fixed": {"uid": "1.2.3"}, "moving": {"uid": "1.2.4"}},
        quality=_quality() if quality is None else quality,
    )


def _run_cli_in_clean_interpreter(argv):
    code = """
import json
import sys
from typer.testing import CliRunner
from dicom_surface import cli

heavy_modules = {
    "dicom_surface.merge",
    "dicom_surface.pipeline",
    "dicom_surface.registration",
    "dicom_surface.repair",
    "dicom_surface.surface",
    "dicom_surface.validate",
    "meshlib",
    "vtk",
}
assert heavy_modules.isdisjoint(sys.modules)
result = CliRunner().invoke(
    cli.app,
    json.loads(sys.argv[1]),
    prog_name="dicom-surface",
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
    assert "Usage: dicom-surface [OPTIONS] [COMMAND]" in result.stdout
    for command in ("list", "presets", "convert", "merge", "validate", "repair"):
        assert command in result.stdout


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
    assert "Usage: dicom-surface" in result.stdout
    assert result.stderr == ""
    assert "\x1b" not in result.stdout
    assert "Traceback" not in result.stdout


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
    assert result.stderr.startswith("Usage: dicom-surface")
    assert "\x1b" not in result.stderr
    assert "Error" in result.stderr
    assert "Traceback" not in result.stderr


def test_python_module_uses_public_program_name():
    result = subprocess.run(
        [sys.executable, "-m", "dicom_surface", "--help"],
        capture_output=True,
        text=True,
        check=False,
        env=_plain_subprocess_env(),
    )

    assert result.returncode == 0
    assert "Usage: dicom-surface" in result.stdout
    assert "\x1b" not in result.stdout


def test_registry_enums_match_presets_exactly():
    assert {choice.value for choice in cli.PresetChoice} == set(cli.PRESETS)
    assert {choice.value for choice in cli.PrintProfileChoice} == set(cli.PRINT_PROFILES)


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
    ],
)
def test_removed_self_intersection_flag_is_a_usage_error(argv):
    result = runner.invoke(cli.app, argv, prog_name="dicom-surface")

    assert result.exit_code == 2
    assert "No such option" in result.stderr


@pytest.mark.parametrize(
    "command, argv",
    [
        ("list", ["list", "missing"]),
        ("convert", ["convert", "missing", "-o", "out.stl"]),
        ("merge", ["merge", "missing", "-o", "out.stl"]),
        ("validate", ["validate", "missing.stl"]),
        ("repair", ["repair", "missing.stl", "-o", "fixed.stl"]),
    ],
)
def test_required_inputs_use_typer_path_validation(command, argv):
    result = runner.invoke(cli.app, argv, prog_name="dicom-surface")

    assert result.exit_code == 2
    assert "does not exist" in result.stderr
    assert command in result.stderr


def test_list_json_uses_unique_id_and_plain_stdout(tmp_path, monkeypatch):
    found = [_series(description="[bold red]literal[/bold red]")]
    monkeypatch.setattr(cli, "_discover", lambda _root: found)

    result = runner.invoke(cli.app, ["list", str(tmp_path), "--json"], prog_name="dicom-surface")

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload[0]["id"] == 1
    assert payload[0]["series_number"] == 6
    assert payload[0]["part"] == 1
    assert payload[0]["n_parts"] == 1
    assert payload[0]["description"] == "[bold red]literal[/bold red]"
    assert "ident" not in payload[0]
    assert "\x1b" not in result.stdout
    assert result.stderr == ""


def test_list_human_output_shows_discovery_progress(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "_discover", lambda _root: [_series()])

    result = runner.invoke(cli.app, ["list", str(tmp_path)], prog_name="dicom-surface")

    assert result.exit_code == 0
    assert "Discovering DICOM series ..." in result.stdout
    assert "DICOM #" in result.stdout


def test_discovery_warnings_are_deduplicated_and_do_not_corrupt_json(tmp_path, monkeypatch):
    from dicom_surface import series as series_mod

    def fake_discover(_root):
        for _ in range(3):
            warnings.warn("Invalid value for VR UI: broken UID", UserWarning)
        warnings.warn("Invalid value for VR DS: broken spacing", UserWarning)
        return [_series()]

    monkeypatch.setattr(series_mod, "discover", fake_discover)
    result = runner.invoke(cli.app, ["list", str(tmp_path), "--json"], prog_name="dicom-surface")

    assert result.exit_code == 0
    assert json.loads(result.stdout)[0]["id"] == 1
    assert result.stderr.count("Invalid value for VR UI") == 1
    assert "repeated 3 times" in result.stderr
    assert result.stderr.count("Invalid value for VR DS") == 1
    assert ".py:" not in result.stderr
    assert "\x1b" not in result.stdout + result.stderr


def test_series_table_is_responsive_safe_and_ansi_free():
    recommended = _series(
        row_id=1,
        description="[bold]recommended acquisition with a deliberately long description[/bold]",
    )
    dose = _series(
        row_id=2,
        uid="1.2.4",
        number=6,
        modality="RTDOSE",
        description="radiotherapy dose object with another long description",
    )
    dose.files = ["one"]
    plan = _series(
        row_id=3,
        uid="1.2.5",
        number=17,
        modality="RTPLAN",
        description="plan notes",
    )
    plan.files = ["one"]
    plan.n_parts = 2
    plan.part = 2

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
        console.print(cli._series_table([recommended, dose, plan], recommended))
        rendered = stream.getvalue()
        normalized = " ".join(rendered.split())

        assert "ID" in rendered
        assert "DICOM #" in rendered
        assert "RTDOSE" in rendered
        assert "RTPLAN" in rendered
        assert "default" in normalized
        assert "not an" in normalized
        assert "image series" in normalized
        assert "orientation" in normalized
        assert "[bold]" in normalized
        assert "\x1b" not in rendered
        assert all(len(line) <= width for line in rendered.splitlines())


def test_presets_renders_two_rich_tables_without_ansi():
    result = runner.invoke(cli.app, ["presets"], prog_name="dicom-surface")

    assert result.exit_code == 0
    assert "Tissue presets (--preset)" in result.stdout
    assert "Print profiles (--print-profile)" in result.stdout
    assert "bone-detail" in result.stdout
    assert "anatomical" in result.stdout
    assert "\x1b" not in result.stdout


def test_convert_accepts_all_flags_and_writes_json_file(tmp_path, monkeypatch):
    chosen = _series()
    captured = {}
    output = tmp_path / "surface.stl"
    json_file = tmp_path / "result.json"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from dicom_surface import pipeline

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
            "--series",
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
            "--print-profile",
            "resin",
            "--min-feature-mm",
            "0.7",
            "--all-islands",
            "--all-components",
            "--resample-mm",
            "0.8",
            "--smooth-iters",
            "11",
            "--passband",
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
        prog_name="dicom-surface",
    )

    assert result.exit_code == 0, result.output
    assert result.stdout == ""
    assert result.stderr == ""
    assert captured["series"] is chosen
    assert captured["preset"].name == "skin"
    assert captured["preset"].threshold == "auto"
    assert captured["preset"].median_mm == pytest.approx(1.1)
    assert captured["preset"].opening_mm == pytest.approx(3.3)
    assert captured["preset"].simplify_error_mm == pytest.approx(0.18)
    assert not captured["preset"].keep_largest_island
    assert not captured["preset"].keep_largest_component
    assert captured["threshold"] is None
    assert not captured["cap_field_of_view"]
    assert captured["print_profile"].name == "resin+flags"
    assert captured["print_profile"].min_feature_mm == pytest.approx(0.7)
    payload = json.loads(json_file.read_text())
    assert payload["result"]["output"] == str(output)
    assert payload["quality"]["valid"]


def test_convert_defaults_quality_output_and_invalid_exit(tmp_path, monkeypatch):
    chosen = _series()
    captured = {}
    output = tmp_path / "surface.stl"
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from dicom_surface import pipeline

    def fake_convert(**kwargs):
        captured.update(kwargs)
        return _convert_result(str(output), quality=_quality(valid=False))

    monkeypatch.setattr(pipeline, "convert", fake_convert)
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", str(output)],
        prog_name="dicom-surface",
    )

    assert result.exit_code == 1
    assert captured["preset"].name == "bone"
    assert captured["threshold"] is None
    assert captured["cap_field_of_view"]
    assert captured["print_profile"] == cli.presets_mod.ANATOMICAL
    assert "Success: wrote" in result.stdout
    assert "Mesh quality" in result.stdout
    assert "invalid" in result.stdout
    assert "Problems" in result.stdout
    assert "failed validation" in result.stderr


def test_invalid_threshold_is_exit_two_before_discovery(tmp_path, monkeypatch):
    monkeypatch.setattr(
        cli,
        "_discover",
        lambda _root: pytest.fail("discovery should not run for an invalid threshold"),
    )
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", "out.stl", "--threshold", "wat"],
        prog_name="dicom-surface",
    )

    assert result.exit_code == 2
    assert result.stdout == ""
    assert "number or 'auto'" in result.stderr


def test_selection_error_is_exit_two(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "_discover", lambda _root: [_series()])
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", "out.stl", "--series", "99"],
        prog_name="dicom-surface",
    )

    assert result.exit_code == 2
    assert "no series matches" in result.stderr


def test_quiet_suppresses_progress_but_not_warnings(tmp_path, monkeypatch):
    chosen = _series()
    monkeypatch.setattr(cli, "_discover", lambda _root: [chosen])

    from dicom_surface import pipeline

    monkeypatch.setattr(
        pipeline,
        "convert",
        lambda **_kwargs: _convert_result("out.stl", warnings_=["sampling warning"]),
    )
    result = runner.invoke(
        cli.app,
        ["convert", str(tmp_path), "-o", "out.stl", "-q"],
        prog_name="dicom-surface",
    )

    assert result.exit_code == 0
    assert result.stdout == ""
    assert "Warning: sampling warning" in result.stderr


def test_merge_accepts_all_flags_and_safety_errors_exit_three(tmp_path, monkeypatch):
    directory_b = tmp_path / "moving"
    directory_b.mkdir()
    fixed = _series(uid="1.2.3", description="fixed")
    moving = _series(uid="1.2.4", description="moving")
    captured = {}
    output = tmp_path / "merged.stl"
    json_file = tmp_path / "merged.json"

    def fake_discover(path: Path):
        return [moving] if path == directory_b else [fixed]

    monkeypatch.setattr(cli, "_discover", fake_discover)
    from dicom_surface import merge as merge_mod

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
            "--series-a",
            "1",
            "--series-b",
            "1",
            "--preset",
            "bone-detail",
            "--threshold",
            "250",
            "--median-mm",
            "1.1",
            "--closing-mm",
            "2.2",
            "--min-island-mm3",
            "3.3",
            "--print-profile",
            "fdm",
            "--min-feature-mm",
            "1.4",
            "--grid-mm",
            "0.8",
            "--smooth-iters",
            "12",
            "--passband",
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
        prog_name="dicom-surface",
    )

    assert result.exit_code == 0, result.output
    assert captured["series_a"] is fixed
    assert captured["series_b"] is moving
    assert captured["preset"].name == "bone-detail"
    assert captured["threshold"] == pytest.approx(250.0)
    assert captured["grid_mm"] == pytest.approx(0.8)
    assert captured["smooth_iters"] == 12
    assert captured["passband"] == pytest.approx(0.1)
    assert captured["simplify_error_mm"] == pytest.approx(0.2)
    assert captured["post_smooth_iters"] == 8
    assert captured["force"]
    assert captured["print_profile"].name == "fdm+flags"
    assert json.loads(json_file.read_text())["result"]["grid_size"] == [10, 20, 30]

    def refuse(**_kwargs):
        raise merge_mod.MergeError("registration gate refused the pair")

    monkeypatch.setattr(merge_mod, "merge", refuse)
    refused = runner.invoke(
        cli.app,
        ["merge", str(tmp_path), str(directory_b), "-o", str(output)],
        prog_name="dicom-surface",
    )
    assert refused.exit_code == 3
    assert "registration gate refused" in refused.stderr


def test_validate_json_is_plain_and_invalid_quality_exits_one(tmp_path, monkeypatch):
    mesh = tmp_path / "mesh.stl"
    mesh.write_text("placeholder")
    from dicom_surface import validate as validate_mod

    monkeypatch.setattr(validate_mod, "validate", lambda _path: _quality(valid=False))
    result = runner.invoke(
        cli.app,
        ["validate", str(mesh), "--json"],
        prog_name="dicom-surface",
    )

    assert result.exit_code == 1
    assert json.loads(result.stdout)["valid"] is False
    assert result.stderr == ""
    assert "\x1b" not in result.stdout

    human = runner.invoke(
        cli.app,
        ["validate", str(mesh)],
        prog_name="dicom-surface",
    )
    assert human.exit_code == 1
    assert "Loading validation engine ..." in human.stdout
    assert "Validating mesh structure and self-intersections ..." in human.stdout
    assert "Mesh quality" in human.stdout


def test_repair_json_is_plain_and_errors_are_concise(tmp_path, monkeypatch):
    mesh = tmp_path / "mesh.stl"
    output = tmp_path / "fixed.stl"
    mesh.write_text("placeholder")
    from dicom_surface import repair as repair_mod, validate as validate_mod

    monkeypatch.setattr(repair_mod, "repair", lambda *_args, **_kwargs: {"holes_filled": 1})
    monkeypatch.setattr(validate_mod, "validate", lambda _path: _quality())
    result = runner.invoke(
        cli.app,
        ["repair", str(mesh), "-o", str(output), "--json"],
        prog_name="dicom-surface",
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
        prog_name="dicom-surface",
    )
    assert failed.exit_code == 1
    assert "Loading repair engine ..." in failed.stdout
    assert "Loading mesh for repair ..." in failed.stdout
    assert "Error: cannot repair" in failed.stderr
    assert "Traceback" not in failed.stderr
