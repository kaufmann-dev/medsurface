"""Preset resolution and format-neutral threshold behavior."""

from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk

from medsurface import pipeline as pipeline_mod, presets, validate
from medsurface.catalog import DicomSource, FileSource, VolumeCandidate
from medsurface.pipeline import resolve_threshold
from medsurface.series import Series
from tests.mesh_helpers import write


def _series(modality="CT"):
    return Series(uid="1.2.3", modality=modality, description="test", series_number=1)


def _candidate(*, modality="CT", file=False):
    series = _series(modality)
    source = FileSource(Path("scan.nii"), "NIfTI") if file else DicomSource(Path("scans"), series)
    return VolumeCandidate(
        id=1,
        source=source,
        format="NIfTI" if file else "DICOM",
        source_name="scan.nii" if file else "scans",
        modality=None if file else modality,
        description=None,
        size=(32, 32, 32),
        spacing=(1.0, 1.0, 1.0),
        direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        origin=(0.0, 0.0, 0.0),
        pixel_type="32-bit float",
        components=1,
        plane="axial",
    )


def _ct_image():
    rng = np.random.default_rng(3)
    arr = np.concatenate(
        [np.full(20_000, -1000.0), rng.normal(40, 20, 8_000), rng.normal(900, 80, 3_000)]
    ).astype(np.float32)
    side = 32
    vol = np.full(side**3, -1000.0, dtype=np.float32)
    vol[: arr.size] = arr[: side**3]
    return sitk.GetImageFromArray(vol.reshape(side, side, side))


def test_conversion_announces_volume_loading_before_it_starts(monkeypatch):
    messages = []

    class StopLoading(Exception):
        pass

    def stop(_candidate):
        assert messages[-1] == "load volume ..."
        raise StopLoading

    monkeypatch.setattr(pipeline_mod.volume_mod, "load", stop)
    with pytest.raises(StopLoading):
        pipeline_mod.convert(
            _candidate(),
            presets.get("bone"),
            "unused.stl",
            log=messages.append,
        )


def test_explicit_threshold_wins_over_preset():
    v, src = resolve_threshold(_ct_image(), presets.get("bone"), 123.0)
    assert v == 123.0 and src == "explicit"


def test_preset_threshold_used_for_ct():
    v, src = resolve_threshold(_ct_image(), presets.get("bone"), None)
    assert v == 300.0 and src == "preset:bone"


def test_explicit_auto_uses_otsu():
    value, source = resolve_threshold(_ct_image(), presets.get("bone"), "auto")
    assert source == "otsu"
    assert np.isfinite(value)


def test_auto_preset_works_on_any_modality():
    v, src = resolve_threshold(_ct_image(), presets.get("auto"), None)
    assert src == "otsu"
    assert np.isfinite(v)


def test_explicit_threshold_bypasses_the_modality_guard():
    v, src = resolve_threshold(_ct_image(), presets.get("bone"), 500.0)
    assert v == 500.0 and src == "explicit"


def test_hu_warning_depends_on_verified_calibration_not_processing():
    bone = presets.get("bone")
    assert pipeline_mod.threshold_warnings(_candidate(), bone, None) == []
    warnings = pipeline_mod.threshold_warnings(_candidate(file=True), bone, None)
    assert len(warnings) == 1 and "cannot be verified" in warnings[0]
    assert "--threshold" in warnings[0]
    merge_warning = pipeline_mod.threshold_warnings(
        _candidate(file=True), bone, None, "--moving-threshold"
    )
    assert "--moving-threshold" in merge_warning[0]
    assert pipeline_mod.threshold_warnings(_candidate(file=True), bone, 300.0) == []


def test_every_preset_is_self_consistent():
    for name, p in presets.PRESETS.items():
        assert p.name == name
        assert p.description
        assert p.median_mm >= 0 and p.closing_mm >= 0 and p.opening_mm >= 0
        assert 0.0 < p.smooth_force <= 0.5
        assert p.simplify_error_mm >= 0
        if isinstance(p.threshold, str):
            assert p.threshold == "auto"
            assert p.threshold_unit == "auto"
        else:
            assert p.threshold_unit == "HU"


def test_override_ignores_none_and_applies_values():
    base = presets.get("bone")
    same = presets.override(base, median_mm=None, closing_mm=None)
    assert same == base
    changed = presets.override(base, median_mm=2.0)
    assert changed.median_mm == 2.0
    assert changed.closing_mm == base.closing_mm


def test_unknown_preset_lists_alternatives():
    with pytest.raises(KeyError, match="available:"):
        presets.get("nope")


def test_validate_does_not_repair_the_mesh_it_measures(tmp_path):
    """A validator must not silently fill an open input mesh."""

    # An open box: five faces of a cube, one side missing.
    v = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
         [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    f = [[0, 2, 1], [0, 3, 2],  # bottom
         [4, 5, 6], [4, 6, 7],  # top
         [0, 1, 5], [0, 5, 4],  # front
         [1, 2, 6], [1, 6, 5],  # right
         [2, 3, 7], [2, 7, 6]]  # back  (left face omitted)
    p = str(tmp_path / "openbox.stl")
    write(p, (np.asarray(v, dtype=float), np.asarray(f, dtype=np.int32)))

    report = validate.validate(p)
    assert report["watertight"] is False
    assert report["boundary_edges"] > 0
    assert report["components"] == 1
    assert report["volume_mm3"] is None


def test_validate_reports_open_mesh_volume_as_none(tmp_path):
    """An open mesh has no well-defined volume."""

    # a single triangle: maximally open
    p = str(tmp_path / "tri.stl")
    write(
        p,
        (
            np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float),
            np.asarray([[0, 1, 2]], dtype=np.int32),
        ),
    )
    report = validate.validate(p)
    assert report["watertight"] is False
    assert report["volume_mm3"] is None
    assert report["boundary_edges"] == 3
    assert "genus" not in report
