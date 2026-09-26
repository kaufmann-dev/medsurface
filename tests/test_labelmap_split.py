"""Label selection, per-label split extraction, and label-name tables."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk
from typer.testing import CliRunner

from medsurface import catalog, cli, labelmap, labelnames

runner = CliRunner()


def _sphere(shape, center, radius) -> np.ndarray:
    zz, yy, xx = np.indices(shape)
    distance = (zz - center[0]) ** 2 + (yy - center[1]) ** 2 + (xx - center[2]) ** 2
    return distance <= radius**2


def _phantom() -> np.ndarray:
    """Three spheres, two touching cubes, and one label cut by the boundary."""
    values = np.zeros((48, 48, 48), dtype=np.uint16)
    values[_sphere(values.shape, (12, 12, 12), 6)] = 3
    values[_sphere(values.shape, (12, 34, 12), 5)] = 7
    values[_sphere(values.shape, (34, 12, 30), 6)] = 12
    values[28:36, 28:36, 28:34] = 20
    values[28:36, 28:36, 34:40] = 21
    values[40:48, 36:44, 6:14] = 300
    return values


def _write(path: Path, values: np.ndarray, spacing=(1.0, 1.0, 1.0)):
    image = sitk.GetImageFromArray(values)
    image.SetSpacing(spacing)
    image.SetOrigin((-10.0, 5.0, 2.5))
    sitk.WriteImage(image, str(path))
    return catalog.discover(path)[0]


def _fast_settings(**overrides):
    return labelmap.resolve_surface_settings(
        mask_smooth_mm=overrides.get("mask_smooth_mm", 0.0),
        surface_smooth_iters=0,
        simplify_error_mm=0,
    )


def test_normalize_labels_sorts_and_rejects_invalid_ids():
    assert labelmap.normalize_labels([7, 3, 7]) == (3, 7)
    assert labelmap.normalize_labels(None) is None
    for invalid in ([0], [-1], [1.5], []):
        with pytest.raises(ValueError):
            labelmap.normalize_labels(invalid)


def test_binary_mask_can_select_labels():
    image = sitk.GetImageFromArray(_phantom())
    mask = labelmap._binary_mask(image, [3, 300])
    values = sitk.GetArrayViewFromImage(mask)
    expected = np.isin(_phantom(), [3, 300])
    assert np.array_equal(values.astype(bool), expected)
    with pytest.raises(ValueError, match="none of the selected labels"):
        labelmap._binary_mask(image, [99])


def test_inspect_labels_reports_boxes_volume_and_boundary():
    values = _phantom()
    image = sitk.GetImageFromArray(values)
    image.SetSpacing((0.5, 1.0, 2.0))
    stats = labelmap.inspect_labels(image)

    assert sorted(stats) == [3, 7, 12, 20, 21, 300]
    assert stats[300].touches_boundary
    assert not stats[3].touches_boundary
    assert stats[20].voxels == 8 * 8 * 6
    assert stats[20].volume_mm3 == pytest.approx(8 * 8 * 6 * 1.0)
    # SimpleITK index order is x, y, z.
    assert stats[20].bbox_index == (28, 28, 28)
    assert stats[20].bbox_size == (6, 8, 8)


def test_extract_with_selected_labels_matches_subset(tmp_path):
    candidate = _write(tmp_path / "labels.nii.gz", _phantom())
    output = tmp_path / "selected.stl"

    result = labelmap.extract(
        candidate,
        str(output),
        mask_smooth_mm=0,
        surface_smooth_iters=0,
        simplify_error_mm=0,
        labels=[3, 7],
    )

    assert result.quality["valid"]
    assert result.quality["components"] == 2
    assert result.provenance["input"]["foreground"] == "labels 3, 7"


def test_all_labels_selection_equals_default_extraction(tmp_path):
    candidate = _write(tmp_path / "labels.nii.gz", _phantom())
    default = labelmap.extract(
        candidate, str(tmp_path / "all.stl"), mask_smooth_mm=0.8
    )
    selected = labelmap.extract(
        candidate,
        str(tmp_path / "selected.stl"),
        mask_smooth_mm=0.8,
        labels=[3, 7, 12, 20, 21, 300],
    )
    assert (tmp_path / "all.stl").read_bytes() == (tmp_path / "selected.stl").read_bytes()
    assert default.triangles == selected.triangles


@pytest.mark.parametrize("mask_smooth_mm", [0.0, 0.8])
def test_split_crop_matches_uncropped_extraction(tmp_path, mask_smooth_mm):
    candidate = _write(tmp_path / "labels.nii.gz", _phantom())
    settings = labelmap.resolve_surface_settings(
        mask_smooth_mm=mask_smooth_mm, surface_smooth_iters=0, simplify_error_mm=0
    )
    split = labelmap.extract_labels(
        candidate, str(tmp_path / "parts"), labels=[3, 12], settings=settings
    )
    assert not split.failed
    for mesh in split.meshes:
        whole = labelmap.extract(
            candidate,
            str(tmp_path / ("whole-%d.stl" % mesh.label)),
            mask_smooth_mm=mask_smooth_mm,
            surface_smooth_iters=0,
            simplify_error_mm=0,
            labels=[mesh.label],
        )
        assert mesh.quality["valid"]
        assert not mesh.capped_field_of_view
        assert np.allclose(mesh.bounds_mm, whole.bounds_mm, atol=0.1)
        assert mesh.quality["volume_mm3"] == pytest.approx(
            whole.quality["volume_mm3"], rel=1e-3
        )


def test_split_writes_named_files_and_caps_only_the_boundary_label(tmp_path):
    candidate = _write(tmp_path / "labels.nii.gz", _phantom())
    events: list[dict] = []
    seen: list[int] = []

    result = labelmap.extract_labels(
        candidate,
        str(tmp_path / "parts"),
        names={3: "Left Femur", 300: "skull"},
        settings=_fast_settings(),
        combined_path=str(tmp_path / "combined.stl"),
        progress=events.append,
        on_label=lambda mesh: seen.append(mesh.label),
    )

    names = sorted(path.name for path in (tmp_path / "parts").iterdir())
    assert names == [
        "003_left-femur.stl",
        "007_label-7.stl",
        "012_label-12.stl",
        "020_label-20.stl",
        "021_label-21.stl",
        "300_skull.stl",
    ]
    assert seen == [3, 7, 12, 20, 21, 300]
    assert result.valid
    capped = {mesh.label for mesh in result.meshes if mesh.capped_field_of_view}
    assert capped == {300}
    assert all(mesh.quality["watertight"] for mesh in result.meshes)
    assert result.combined is not None
    assert result.combined.quality["valid"]
    # Touching cubes 20 and 21 merge into one shell in the union.
    assert result.combined.quality["components"] == 5
    kinds = {event["event"] for event in events}
    assert {"stage_start", "stage_end", "label_start", "label_end"} <= kinds
    assert result.provenance["labels"]["300"]["touches_boundary"] is True


def test_split_skips_missing_labels_and_reports_them(tmp_path):
    candidate = _write(tmp_path / "labels.nii.gz", _phantom())
    result = labelmap.extract_labels(
        candidate, str(tmp_path / "parts"), labels=[3, 99], settings=_fast_settings()
    )
    assert [mesh.label for mesh in result.meshes] == [3]
    assert result.skipped == [{"label": 99, "reason": "no voxels"}]
    assert any("label 99" in warning for warning in result.warnings)
    with pytest.raises(ValueError, match="none of the selected labels"):
        labelmap.extract_labels(
            candidate, str(tmp_path / "none"), labels=[98], settings=_fast_settings()
        )


def test_split_accepts_per_label_settings(tmp_path):
    candidate = _write(tmp_path / "labels.nii.gz", _phantom())
    requested: list[int] = []

    def per_label(label: int):
        requested.append(label)
        return labelmap.resolve_surface_settings(
            mask_smooth_mm=0.0 if label == 3 else 0.8,
            surface_smooth_iters=0,
            simplify_error_mm=0,
        )

    result = labelmap.extract_labels(
        candidate, str(tmp_path / "parts"), labels=[3, 7], settings=per_label
    )
    assert requested == [3, 7]
    smoothing = {mesh.label: mesh.surface["mask_smooth_mm"] for mesh in result.meshes}
    assert smoothing == {3: 0.0, 7: 0.8}


def test_label_table_round_trips_through_nifti(tmp_path):
    path = tmp_path / "labels.nii.gz"
    _write(path, _phantom())
    before = sitk.GetArrayFromImage(sitk.ReadImage(str(path)))

    assert labelnames.embed_label_names(path, {3: "femur_left", 300: "skull <main>"})
    assert labelnames.read_label_names(path) == {3: "femur_left", 300: "skull <main>"}
    after = sitk.ReadImage(str(path))
    assert np.array_equal(before, sitk.GetArrayFromImage(after))
    with gzip.open(path, "rb") as handle:
        assert b"<![CDATA[femur_left]]>" in handle.read(4096)


def test_totalsegmentator_style_table_names_split_files(tmp_path):
    path = tmp_path / "labels.nii"
    _write(path, _phantom())
    labelnames.embed_label_names(path, {12: "vertebrae_L1"})
    candidate = catalog.discover(path)[0]
    result = labelmap.extract_labels(
        candidate, str(tmp_path / "parts"), labels=[12], settings=_fast_settings()
    )
    assert Path(result.meshes[0].output_path).name == "012_vertebrae-l1.stl"
    assert result.meshes[0].name == "vertebrae_L1"


def test_read_label_names_tolerates_other_formats(tmp_path):
    assert labelnames.read_label_names(tmp_path / "missing.nii.gz") == {}
    nrrd = tmp_path / "labels.nrrd"
    _write(nrrd, _phantom())
    assert labelnames.read_label_names(nrrd) == {}
    assert labelnames.embed_label_names(nrrd, {1: "x"}) is False


def test_cli_split_writes_json_and_json_progress(tmp_path):
    source = tmp_path / "labels.nii.gz"
    _write(source, _phantom())
    names = tmp_path / "names.json"
    names.write_text(json.dumps({"3": "femur", "7": "hip"}))
    report = tmp_path / "report.json"

    result = runner.invoke(
        cli.app,
        [
            "labelmap",
            "extract",
            str(source),
            "--split",
            str(tmp_path / "parts"),
            "-o",
            str(tmp_path / "combined.stl"),
            "--labels",
            "3,7",
            "--label-names",
            str(names),
            "--mask-smooth-mm",
            "0",
            "--mesh-smooth-iters",
            "0",
            "--simplify-error-mm",
            "0",
            "--json",
            str(report),
            "--progress",
            "json",
            "-q",
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert result.stdout == ""
    events = [json.loads(line) for line in result.stderr.splitlines() if line.strip()]
    assert any(event["event"] == "label_end" and event["label"] == 7 for event in events)
    payload = json.loads(report.read_text())
    assert [Path(mesh["output"]).name for mesh in payload["result"]["meshes"]] == [
        "003_femur.stl",
        "007_hip.stl",
    ]
    assert payload["result"]["combined"]["quality"]["valid"]
    assert payload["provenance"]["input"]["foreground"] == "labels 3, 7"


def test_cli_requires_output_or_split(tmp_path):
    source = tmp_path / "labels.nii.gz"
    _write(source, _phantom())
    result = runner.invoke(
        cli.app, ["labelmap", "extract", str(source)], prog_name="medsurface"
    )
    assert result.exit_code == 2
    assert "pass -o/--output, --split, or both" in result.stderr


@pytest.mark.parametrize("value", ["0", "a,b", "5-3", ","])
def test_cli_rejects_invalid_label_selection(tmp_path, value):
    source = tmp_path / "labels.nii.gz"
    _write(source, _phantom())
    result = runner.invoke(
        cli.app,
        ["labelmap", "extract", str(source), "-o", str(tmp_path / "x.stl"), "--labels", value],
        prog_name="medsurface",
    )
    assert result.exit_code == 2
    assert "--labels expects positive label IDs" in result.stderr


def test_cli_label_ranges_expand(tmp_path):
    assert cli._parse_label_selection("10-12, 3") == [3, 10, 11, 12]


def test_cli_extract_maps_runtime_errors_to_clean_exit(tmp_path, monkeypatch):
    source = tmp_path / "labels.nii.gz"
    _write(source, _phantom())

    def explode(**_kwargs):
        raise RuntimeError("meshlib exploded")

    monkeypatch.setattr(labelmap, "extract", explode)
    result = runner.invoke(
        cli.app,
        ["labelmap", "extract", str(source), "-o", str(tmp_path / "x.stl")],
        prog_name="medsurface",
    )
    assert result.exit_code == 1
    assert "meshlib exploded" in result.stderr
    assert "Traceback" not in result.output
