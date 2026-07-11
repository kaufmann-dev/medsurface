"""Mesh-repair behavior."""

import json
import os

import pytest
from typer.testing import CliRunner

from dicom_surface import cli, repair, validate
from tests.mesh_helpers import box, write


def test_repair_closes_an_open_mesh(tmp_path):
    vertices, faces = box()

    source = str(tmp_path / "open.stl")
    output = str(tmp_path / "repaired.stl")
    write(source, (vertices, faces[2:]))

    before = validate.validate(source)
    assert not before["watertight"]

    stats = repair.repair(source, output)
    after = validate.validate(output)

    assert stats["holes_in"] > 0
    assert stats["holes_out"] == 0
    assert stats["holes_filled"] > 0
    assert after["watertight"]


def test_repair_refuses_to_overwrite_its_input_or_a_hard_link(tmp_path):
    source = tmp_path / "source.stl"
    alias = tmp_path / "alias.stl"
    write(source, box())
    original = source.read_bytes()

    with pytest.raises(ValueError, match="not in-place"):
        repair.repair(str(source), str(source))
    assert source.read_bytes() == original

    os.link(source, alias)
    with pytest.raises(ValueError, match="not in-place"):
        repair.repair(str(source), str(alias))
    assert source.read_bytes() == original


def test_repair_announces_mesh_loading_before_it_starts(monkeypatch):
    messages = []

    class StopLoading(Exception):
        pass

    def stop(_path):
        assert messages[-1] == "load mesh ..."
        raise StopLoading

    monkeypatch.setattr(repair.mm, "loadMesh", stop)
    with pytest.raises(StopLoading):
        repair.repair("input.stl", "output.stl", log=messages.append)


def test_repair_json_mode_writes_only_json(tmp_path):
    vertices, faces = box()
    source = tmp_path / "open.stl"
    output = tmp_path / "fixed.stl"
    write(source, (vertices, faces[2:]))

    result = CliRunner().invoke(
        cli.app,
        ["repair", str(source), "-o", str(output), "--json"],
        prog_name="dicom-surface",
    )
    assert result.exit_code == 0
    payload = json.loads(result.stdout)

    assert result.stderr == ""
    assert payload["output"] == str(output)
    assert payload["repair"]["holes_filled"] > 0
    assert payload["quality"]["valid"]
