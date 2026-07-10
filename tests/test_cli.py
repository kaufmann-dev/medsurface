"""Cross-command CLI contracts."""

import pytest

from dicom_surface import cli


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
