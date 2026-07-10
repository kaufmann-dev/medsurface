"""Preset resolution and the modality guard."""

import numpy as np
import pytest
import SimpleITK as sitk

from dicom_surface import presets, validate
from dicom_surface.pipeline import ModalityMismatch, resolve_threshold
from dicom_surface.series import Series


def _series(modality="CT"):
    return Series(uid="1.2.3", modality=modality, description="test", series_number=1)


def _ct_image():
    rng = np.random.default_rng(3)
    arr = np.concatenate(
        [np.full(20_000, -1000.0), rng.normal(40, 20, 8_000), rng.normal(900, 80, 3_000)]
    ).astype(np.float32)
    side = 32
    vol = np.full(side**3, -1000.0, dtype=np.float32)
    vol[: arr.size] = arr[: side**3]
    return sitk.GetImageFromArray(vol.reshape(side, side, side))


def test_explicit_threshold_wins_over_preset():
    v, src = resolve_threshold(_ct_image(), _series(), presets.get("bone"), 123.0)
    assert v == 123.0 and src == "explicit"


def test_preset_threshold_used_for_ct():
    v, src = resolve_threshold(_ct_image(), _series("CT"), presets.get("bone"), None)
    assert v == 300.0 and src == "preset:bone"


def test_hu_preset_rejected_for_non_ct():
    """A Hounsfield threshold is meaningless on MR. Refuse rather than emit
    a confidently wrong mesh."""
    for modality in ("MR", "US", "PT", "XA"):
        with pytest.raises(ModalityMismatch, match="only.*meaningful for CT"):
            resolve_threshold(_ct_image(), _series(modality), presets.get("bone"), None)


def test_auto_preset_works_on_any_modality():
    for modality in ("CT", "MR", "US"):
        v, src = resolve_threshold(_ct_image(), _series(modality), presets.get("auto"), None)
        assert src == "otsu"
        assert np.isfinite(v)


def test_explicit_threshold_bypasses_the_modality_guard():
    v, src = resolve_threshold(_ct_image(), _series("MR"), presets.get("bone"), 500.0)
    assert v == 500.0 and src == "explicit"


def test_every_preset_is_self_consistent():
    for name, p in presets.PRESETS.items():
        assert p.name == name
        assert p.description
        assert p.median_mm >= 0 and p.closing_mm >= 0 and p.opening_mm >= 0
        assert 0.0 < p.passband <= 2.0
        assert p.target_faces >= 0
        if isinstance(p.threshold, str):
            assert p.threshold == "auto"
            assert p.modalities == ()
        else:
            assert p.modalities, "an HU threshold must declare its modalities"


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


def test_accepts_modality():
    assert presets.accepts_modality(presets.get("auto"), "MR")
    assert presets.accepts_modality(presets.get("bone"), "CT")
    assert not presets.accepts_modality(presets.get("bone"), "MR")


def test_validate_does_not_repair_the_mesh_it_measures(tmp_path):
    """Regression: trimesh's ``split()`` builds submeshes with ``repair=True``,
    which fills holes. A validator that silently repairs would report an open
    mesh as closed."""
    import trimesh

    # An open box: five faces of a cube, one side missing.
    v = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
         [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    f = [[0, 2, 1], [0, 3, 2],  # bottom
         [4, 5, 6], [4, 6, 7],  # top
         [0, 1, 5], [0, 5, 4],  # front
         [1, 2, 6], [1, 6, 5],  # right
         [2, 3, 7], [2, 7, 6]]  # back  (left face omitted)
    p = str(tmp_path / "openbox.stl")
    trimesh.Trimesh(vertices=v, faces=f, process=False).export(p)

    report = validate.validate(p)
    assert report["watertight"] is False
    assert report["boundary_edges"] > 0
    assert report["components"] == 1
    assert report["volume_mm3"] is None


def test_validate_reports_open_mesh_volume_as_none(tmp_path):
    """An open mesh has no well-defined volume; the report must say so rather
    than emit trimesh's meaningless divergence-theorem number."""
    import trimesh

    # a single triangle: maximally open
    m = trimesh.Trimesh(vertices=[[0, 0, 0], [1, 0, 0], [0, 1, 0]], faces=[[0, 1, 2]])
    p = str(tmp_path / "tri.stl")
    m.export(p)
    report = validate.validate(p)
    assert report["watertight"] is False
    assert report["volume_mm3"] is None
    assert report["boundary_edges"] == 3
    assert "genus" not in report
