"""Every command streams machine-readable progress with ``--progress json``."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import SimpleITK as sitk
from typer.testing import CliRunner

from medsurface import cli

runner = CliRunner()


def _events(stderr: str) -> list[dict]:
    return [json.loads(line) for line in stderr.splitlines() if line.strip()]


def _run(argv: list[str]) -> list[dict]:
    result = runner.invoke(cli.app, [*argv, "--progress", "json"], prog_name="medsurface")
    assert result.exit_code == 0, result.output
    assert result.stdout == ""
    events = _events(result.stderr)
    assert events[0]["event"] == "status"
    final = events[-1]
    assert final["event"] == "result" and final["exit_code"] == 0
    assert final["result"] is not None
    assert all(event["event"] != "result" for event in events[:-1])
    return events


def _stages(events: list[dict]) -> set[str]:
    return {event["stage"] for event in events if event["event"] == "stage_end"}


def _volume(path: Path) -> Path:
    pixels = np.zeros((32, 32, 32), dtype=np.int16)
    pixels[4:12, 4:12, 4:12] = 2000
    pixels[18:28, 16:28, 20:28] = 2000
    sitk.WriteImage(sitk.GetImageFromArray(pixels), str(path))
    return path


def test_every_command_streams_json_progress(tmp_path):
    source = _volume(tmp_path / "scan.nii.gz")
    moving = _volume(tmp_path / "moving.nrrd")
    mesh = tmp_path / "model.stl"

    assert _run(["list", str(source)])

    assert {"load volume", "write and verify volume"} <= _stages(
        _run(["convert", str(source), "-o", str(tmp_path / "copy.mha")])
    )

    assert {"segment", "marching cubes", "validate and publish mesh"} <= _stages(
        _run(["extract", str(source), "-o", str(mesh), "--threshold", "1000"])
    )

    fuse_events = _run(
        [
            "fuse",
            str(source),
            str(moving),
            "-o",
            str(tmp_path / "fused.nrrd"),
            "--fixed-threshold",
            "1000",
            "--moving-threshold",
            "1000",
            "--grid-mm",
            "1",
        ]
    )
    assert "register moving to fixed" in _stages(fuse_events)
    assert any(event["event"] == "warning" for event in fuse_events)

    validate_events = _run(["validate", str(mesh)])
    assert any(
        event["event"] == "status" and "Validating" in event["message"]
        for event in validate_events
    )

    repair_stages = _stages(
        _run(["repair", str(mesh), "-o", str(tmp_path / "repaired.stl")])
    )
    assert {"load mesh", "fill boundary holes"} <= repair_stages


def test_usage_errors_before_processing_are_json_events(tmp_path):
    source = _volume(tmp_path / "labels.nii.gz")
    for argv, message in (
        (["labelmap", "extract", str(source)], "pass -o/--output, --split, or both"),
        (["labelmap", "extract", str(source), "-o", "x.stl", "--labels", "a"], "--labels"),
        (["extract", str(source), "-o", "model.txt"], "unsupported output extension"),
    ):
        result = runner.invoke(
            cli.app, [*argv, "--progress", "json"], prog_name="medsurface"
        )
        assert result.exit_code == 2
        assert result.stdout == ""
        events = _events(result.stderr)
        assert events[0]["event"] == "error" and message in events[0]["message"]
        assert events[-1] == {"event": "result", "exit_code": 2, "result": None}


def test_failed_quality_still_reports_the_result(tmp_path, monkeypatch):
    from medsurface import validate as validate_mod

    mesh = tmp_path / "mesh.stl"
    mesh.write_text("placeholder")
    report = {"valid": False, "problems": ["not watertight"]}
    monkeypatch.setattr(validate_mod, "validate", lambda _path: report)

    result = runner.invoke(
        cli.app, ["validate", str(mesh), "--progress", "json"], prog_name="medsurface"
    )

    assert result.exit_code == 1
    assert result.stdout == ""
    final = _events(result.stderr)[-1]
    assert final == {"event": "result", "exit_code": 1, "result": report}
