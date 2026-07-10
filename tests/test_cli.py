"""Cross-command CLI contracts."""

import json
import subprocess
import sys

import pytest

from dicom_surface import cli


def _run_cli_in_clean_interpreter(argv):
    code = """
import json
import sys
from dicom_surface import cli

heavy_modules = {
    "dicom_surface.merge",
    "dicom_surface.pipeline",
    "dicom_surface.surface",
    "dicom_surface.validate",
    "vtk",
}
assert heavy_modules.isdisjoint(sys.modules)
try:
    status = cli.main(json.loads(sys.argv[1]))
except SystemExit as exc:
    status = exc.code
assert heavy_modules.isdisjoint(sys.modules)
raise SystemExit(status)
"""
    return subprocess.run(
        [sys.executable, "-c", code, json.dumps(argv)],
        capture_output=True,
        text=True,
        check=False,
    )


def test_root_command_without_arguments_is_lightweight_help():
    result = _run_cli_in_clean_interpreter([])

    assert result.returncode == 0
    assert result.stderr == ""
    assert result.stdout.startswith("usage: dicom-surface")
    assert "{list,presets,convert,merge,validate,repair}" in result.stdout


@pytest.mark.parametrize(
    "argv",
    [
        ["unknown"],
        ["--unknown"],
        ["list"],
        ["presets", "unexpected"],
        ["convert", "scan"],
        ["convert", "scan", "-o", "out.stl", "--preset", "unknown"],
        ["merge"],
        ["validate"],
        ["repair", "broken.stl"],
    ],
)
def test_invalid_command_lines_fail_before_heavy_imports(argv):
    result = _run_cli_in_clean_interpreter(argv)

    assert result.returncode == 2
    assert result.stdout == ""
    assert result.stderr.startswith("usage: dicom-surface")
    assert ": error:" in result.stderr
    assert "Traceback" not in result.stderr


def test_every_validating_command_uses_the_same_quality_status():
    assert cli._quality_status(None) == 0
    assert cli._quality_status({"valid": True}) == 0
    assert cli._quality_status({"valid": False}) == 1


@pytest.mark.parametrize(
    "argv",
    [
        ["convert", "scan", "-o", "out.stl", "--self-intersections"],
        ["merge", "a", "b", "-o", "out.stl", "--self-intersections"],
        ["validate", "mesh.stl", "--self-intersections"],
        ["repair", "mesh.stl", "-o", "fixed.stl", "--self-intersections"],
    ],
)
def test_self_intersection_flag_is_not_a_second_validation_mode(argv):
    with pytest.raises(SystemExit, match="2"):
        cli.build_parser().parse_args(argv)


@pytest.mark.parametrize(
    "argv, action",
    [
        (["validate", "missing.stl"], "cannot validate"),
        (["repair", "missing.stl", "-o", "fixed.stl"], "cannot repair"),
    ],
)
def test_mesh_file_errors_are_concise(argv, action, capsys):
    assert cli.main(argv) == 1
    captured = capsys.readouterr()

    assert captured.out == ""
    assert captured.err.startswith("error: %s" % action)
    assert "Traceback" not in captured.err
