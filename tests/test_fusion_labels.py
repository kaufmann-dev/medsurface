"""Label-preserving N-way fusion of partially overlapping labelmaps."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk
from typer.testing import CliRunner

from medsurface import catalog, cli, fusion, labelmap, labelnames

runner = CliRunner()


def _ellipsoid(shape, center, radii) -> np.ndarray:
    zz, yy, xx = np.indices(shape)
    return (
        ((zz - center[0]) / radii[0]) ** 2
        + ((yy - center[1]) / radii[1]) ** 2
        + ((xx - center[2]) / radii[2]) ** 2
    ) <= 1.0


def _body() -> np.ndarray:
    """Labelled structures along z with no translational self-similarity (z, y, x)."""
    values = np.zeros((96, 64, 64), dtype=np.uint8)
    # Irregular spacing and sizes: a periodic stack would let rigid registration
    # lock on one element off, which real anatomy avoids through its variation.
    centers = (8, 19, 33, 44, 60, 71, 86)
    radii = ((3.5, 8, 10), (4.5, 9, 12), (4, 10, 11), (5, 8, 13), (4, 11, 9), (3.5, 9, 12), (4.5, 7, 10))
    for index, (z, radius) in enumerate(zip(centers, radii)):
        values[_ellipsoid(values.shape, (z, 30, 32 + index), radius)] = 10 + index
    values[_ellipsoid(values.shape, (30, 48, 20), (14, 6, 7))] = 2
    values[_ellipsoid(values.shape, (70, 14, 44), (12, 5, 8))] = 3
    # A curved vessel whose position changes along z breaks residual ambiguity.
    zz, yy, xx = np.indices(values.shape)
    center_y = 50 + 6 * np.sin(zz / 9.0)
    center_x = 8 + 0.35 * zz
    vessel = ((yy - center_y) ** 2 + (xx - center_x) ** 2 <= 9) & (zz >= 4) & (zz < 92)
    values[vessel] = 4
    return values


def _write(path: Path, values: np.ndarray, z0: int, *, shift=(0.0, 0.0, 0.0)):
    image = sitk.GetImageFromArray(values)
    image.SetSpacing((1.0, 1.0, 1.0))
    image.SetOrigin((shift[0], shift[1], z0 + shift[2]))
    sitk.WriteImage(image, str(path))
    return catalog.discover(path)[0]


def _dice(a: np.ndarray, b: np.ndarray) -> float:
    total = int(a.sum() + b.sum())
    return 1.0 if total == 0 else 2.0 * float((a & b).sum()) / total


def _truth_on_grid(result, body: np.ndarray) -> np.ndarray:
    fused = sitk.ReadImage(result.output_path)
    truth = np.zeros(sitk.GetArrayViewFromImage(fused).shape, dtype=np.uint8)
    offset = [int(round(-value)) for value in fused.GetOrigin()]  # x, y, z
    zs, ys, xs = body.shape
    truth[
        offset[2] : offset[2] + zs,
        offset[1] : offset[1] + ys,
        offset[0] : offset[0] + xs,
    ] = body
    return truth


def test_three_partial_scans_fuse_into_the_full_labelled_body(tmp_path):
    body = _body()
    fixed = _write(tmp_path / "fixed.nii.gz", body[0:64], 0)
    moving_a = _write(tmp_path / "moving-a.nii.gz", body[24:96], 24)
    # Same anatomy described in a translated physical frame: registration must undo it.
    moving_b = _write(tmp_path / "moving-b.nii.gz", body[16:80], 16, shift=(4.0, -3.0, 2.0))
    output = tmp_path / "fused.nii.gz"
    events: list[dict] = []

    result = labelmap.fuse_labels(
        fixed,
        [moving_a, moving_b],
        str(output),
        grid_mm=1.0,
        names={2: "liver", 3: "spleen"},
        progress=events.append,
    )

    assert result.preserve_labels
    assert result.pixel_type == "uint8"
    assert len(result.registrations) == 2
    assert np.allclose(result.registrations[0].transform[:3, 3], [0.0, 0.0, 0.0], atol=0.5)
    # The moving frame was shifted by +shift, so moving->fixed undoes it.
    assert np.allclose(result.registrations[1].transform[:3, 3], [-4.0, 3.0, -2.0], atol=0.5)
    fused = sitk.GetArrayFromImage(sitk.ReadImage(str(output)))
    truth = _truth_on_grid(result, body)
    for label in (2, 3, 4, *range(10, 17)):
        assert _dice(fused == label, truth == label) > 0.9, label
    assert set(result.labels) == {2, 3, 4, *range(10, 17)}
    assert labelnames.read_label_names(output) == {2: "liver", 3: "spleen"}
    assert result.provenance["output"]["label_names_embedded"]
    assert any(event["event"] == "label_end" for event in events)


def test_binary_n_way_fusion_is_a_union(tmp_path):
    body = _body()
    fixed = _write(tmp_path / "fixed.nii.gz", body[0:64], 0)
    moving = _write(tmp_path / "moving.nii.gz", body[24:96], 24)
    other = _write(tmp_path / "other.nii.gz", body[16:80], 16)
    output = tmp_path / "union.nrrd"

    result = labelmap.fuse_labels(
        fixed, [moving, other], str(output), grid_mm=1.0, preserve_labels=False
    )

    fused = sitk.GetArrayFromImage(sitk.ReadImage(str(output)))
    assert set(np.unique(fused)) == {0, 1}
    truth = _truth_on_grid(result, body) > 0
    assert _dice(fused == 1, truth) > 0.95
    assert result.pixel_type == "uint8"


def test_label_selection_limits_fused_labels(tmp_path):
    body = _body()
    fixed = _write(tmp_path / "fixed.nii.gz", body[0:64], 0)
    moving = _write(tmp_path / "moving.nii.gz", body[24:96], 24)
    result = labelmap.fuse_labels(
        fixed, [moving], str(tmp_path / "subset.nii.gz"), grid_mm=1.0, labels=[2, 3, 4]
    )
    assert set(result.labels) == {2, 3, 4}


def test_wide_labels_use_uint16_output(tmp_path):
    body = _body().astype(np.uint16)
    body[body == 4] = 1000
    fixed = _write(tmp_path / "fixed.nii.gz", body[0:64], 0)
    moving = _write(tmp_path / "moving.nii.gz", body[24:96], 24)
    output = tmp_path / "wide.nii.gz"
    result = labelmap.fuse_labels(fixed, [moving], str(output), grid_mm=1.0)
    assert result.pixel_type == "uint16"
    assert 1000 in result.labels
    assert sitk.ReadImage(str(output)).GetPixelID() == sitk.sitkUInt16


def test_unrelated_anatomy_fails_the_quality_gates(tmp_path):
    body = _body()
    fixed = _write(tmp_path / "fixed.nii.gz", body[0:50], 0)
    unrelated = np.zeros((50, 64, 64), dtype=np.uint8)
    unrelated[_ellipsoid(unrelated.shape, (25, 32, 32), (20, 25, 5))] = 9
    moving = _write(tmp_path / "moving.nii.gz", unrelated, 0)
    with pytest.raises(fusion.FusionError, match="quality gates"):
        labelmap.fuse_labels(fixed, [moving], str(tmp_path / "bad.nii.gz"), grid_mm=1.0)


def test_duplicate_inputs_are_rejected(tmp_path):
    body = _body()
    fixed = _write(tmp_path / "fixed.nii.gz", body[0:50], 0)
    with pytest.raises(fusion.FusionError, match="same labelmap"):
        labelmap.fuse_labels(fixed, [fixed], str(tmp_path / "dup.nii.gz"))


def test_cli_fuses_three_inputs_with_preserved_labels(tmp_path):
    body = _body()
    _write(tmp_path / "a.nii.gz", body[0:64], 0)
    _write(tmp_path / "b.nii.gz", body[24:96], 24)
    _write(tmp_path / "c.nii.gz", body[16:80], 16)
    report = tmp_path / "fused.json"

    result = runner.invoke(
        cli.app,
        [
            "labelmap",
            "fuse",
            str(tmp_path / "a.nii.gz"),
            str(tmp_path / "b.nii.gz"),
            str(tmp_path / "c.nii.gz"),
            "-o",
            str(tmp_path / "fused.nii.gz"),
            "--preserve-labels",
            "--grid-mm",
            "1",
            "--json",
            str(report),
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert "--split" in result.stdout
    payload = json.loads(report.read_text())
    assert payload["result"]["labels_preserved"] is True
    assert sorted(int(label) for label in payload["result"]["labels"]) == [
        2,
        3,
        4,
        *range(10, 17),
    ]
    assert len(payload["provenance"]["registrations"]) == 2


def test_cli_two_inputs_without_flags_write_a_binary_union(tmp_path):
    body = _body()
    _write(tmp_path / "a.nii.gz", body[0:64], 0)
    _write(tmp_path / "b.nii.gz", body[24:96], 24)
    result = runner.invoke(
        cli.app,
        [
            "labelmap",
            "fuse",
            str(tmp_path / "a.nii.gz"),
            str(tmp_path / "b.nii.gz"),
            "-o",
            str(tmp_path / "binary.nrrd"),
            "--grid-mm",
            "1",
        ],
        prog_name="medsurface",
    )
    assert result.exit_code == 0, result.output
    stored = sitk.GetArrayFromImage(sitk.ReadImage(str(tmp_path / "binary.nrrd")))
    assert set(np.unique(stored)) == {0, 1}


def test_cli_json_progress_streams_fusion_events(tmp_path):
    body = _body()
    _write(tmp_path / "a.nii.gz", body[0:64], 0)
    _write(tmp_path / "b.nii.gz", body[24:96], 24)
    result = runner.invoke(
        cli.app,
        [
            "labelmap",
            "fuse",
            str(tmp_path / "a.nii.gz"),
            str(tmp_path / "b.nii.gz"),
            "-o",
            str(tmp_path / "binary.nrrd"),
            "--grid-mm",
            "1",
            "--progress",
            "json",
            "-q",
        ],
        prog_name="medsurface",
    )
    assert result.exit_code == 0, result.output
    events = [json.loads(line) for line in result.stderr.splitlines() if line.strip()]
    stages = {event.get("stage") for event in events if event["event"] == "stage_end"}
    assert "register moving 1 to fixed" in stages
    assert any(event["event"] == "warning" for event in events)


@pytest.mark.parametrize("grid_mm", [0.7, 0.45, 0.4, 0.33])
def test_touching_labels_keep_the_whole_union(grid_mm):
    values = np.zeros((40, 40, 40), dtype=np.uint8)
    values[8:32, 8:32, 8:20] = 1
    values[8:32, 8:20, 20:32] = 2
    values[8:32, 20:32, 20:32] = 3
    image = sitk.GetImageFromArray(values)

    labelled = fusion.fuse_label_fields([image], [np.eye(4)], grid_mm).labels
    union = fusion.fuse_label_fields([image], [np.eye(4)], grid_mm, binary=True).labels

    np.testing.assert_array_equal(labelled > 0, union > 0)
    assert set(np.unique(labelled)) == {0, 1, 2, 3}
