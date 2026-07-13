"""Thresholding, morphology and island handling."""

import numpy as np
import pytest
import SimpleITK as sitk

from medsurface import segment


def test_otsu_separates_two_modes():
    rng = np.random.default_rng(0)
    low = rng.normal(-100, 15, 50_000)
    high = rng.normal(600, 40, 20_000)
    t = segment.otsu_threshold(np.concatenate([low, high]))
    assert -50 < t < 550


def test_otsu_on_constant_input_does_not_crash():
    assert segment.otsu_threshold(np.full(100, 42.0)) == pytest.approx(42.0)


def test_otsu_rejects_empty():
    with pytest.raises(ValueError):
        segment.otsu_threshold(np.array([]))


def test_auto_threshold_ignores_ct_air():
    """Whole-volume Otsu on a CT just finds the air/tissue edge. Masking air out
    makes it find the tissue/bone edge instead."""
    rng = np.random.default_rng(1)
    air = np.full(200_000, -1000.0)
    soft = rng.normal(40, 20, 40_000)
    bone = rng.normal(800, 100, 10_000)
    arr = np.concatenate([air, soft, bone]).astype(np.float32)
    # pad to a 3D shape
    side = int(np.ceil(arr.size ** (1 / 3))) + 1
    vol = np.zeros(side**3, dtype=np.float32)
    vol[: arr.size] = arr
    vol[arr.size :] = -1000.0
    img = sitk.GetImageFromArray(vol.reshape(side, side, side))

    t = segment.auto_threshold(img)
    assert 100 < t < 700, t


def test_binarize_bounds():
    arr = np.array([[[-1000, 0, 300, 900, 3000]]], dtype=np.int16)
    img = sitk.GetImageFromArray(arr)
    out = sitk.GetArrayFromImage(segment.binarize(img, 300, None))
    assert out.ravel().tolist() == [0, 0, 1, 1, 1]

    capped = sitk.GetArrayFromImage(segment.binarize(img, 300, 1000))
    assert capped.ravel().tolist() == [0, 0, 1, 1, 0]


def _two_blobs():
    arr = np.zeros((10, 10, 10), dtype=np.uint8)
    arr[1:6, 1:6, 1:6] = 1  # 125 voxels
    arr[8:10, 8:10, 8:10] = 1  # 8 voxels
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((1.0, 1.0, 1.0))
    return img


def test_keep_largest_island():
    out = sitk.GetArrayFromImage(segment.islands(_two_blobs(), keep_largest=True, min_mm3=0))
    assert out.sum() == 125
    assert out[9, 9, 9] == 0


def test_min_island_volume_drops_small_blobs():
    out = sitk.GetArrayFromImage(
        segment.islands(_two_blobs(), keep_largest=False, min_mm3=50.0)
    )
    assert out.sum() == 125  # the 8-voxel (8 mm3) blob is gone


def test_min_island_volume_keeps_both_when_permissive():
    out = sitk.GetArrayFromImage(
        segment.islands(_two_blobs(), keep_largest=False, min_mm3=1.0)
    )
    assert out.sum() == 133


def test_islands_raises_when_nothing_survives():
    empty = sitk.GetImageFromArray(np.zeros((5, 5, 5), dtype=np.uint8))
    empty.SetSpacing((1.0, 1.0, 1.0))
    with pytest.raises(ValueError, match="no foreground"):
        segment.islands(empty, keep_largest=True, min_mm3=0)


def test_closing_seals_a_small_pore():
    arr = np.zeros((12, 12, 12), dtype=np.uint8)
    arr[2:10, 2:10, 2:10] = 1
    arr[6, 6, 6] = 0  # one-voxel pore
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((1.0, 1.0, 1.0))

    before = sitk.GetArrayFromImage(img)[6, 6, 6]
    after = sitk.GetArrayFromImage(segment.closing(img, 3.0))[6, 6, 6]
    assert before == 0 and after == 1


def test_median_removes_an_isolated_speck():
    arr = np.zeros((12, 12, 12), dtype=np.uint8)
    arr[2:10, 2:10, 2:10] = 1
    arr[0, 0, 0] = 1  # lone noise voxel
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((1.0, 1.0, 1.0))
    out = sitk.GetArrayFromImage(segment.median(img, 3.0))
    assert out[0, 0, 0] == 0
    assert out[5, 5, 5] == 1


def test_filters_are_noop_when_kernel_collapses():
    """A 1 mm kernel on 3 mm slices must not silently become a 3 mm kernel."""
    arr = np.zeros((6, 12, 12), dtype=np.uint8)
    arr[2:4, 2:10, 2:10] = 1
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((0.4, 0.4, 3.0))
    same = sitk.GetArrayFromImage(segment.median(img, 1.0))
    assert np.array_equal(same, arr)


def test_pad_adds_background_border():
    arr = np.ones((3, 3, 3), dtype=np.uint8)
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((1.0, 1.0, 1.0))
    padded = sitk.GetArrayFromImage(segment.pad(img, 1))
    assert padded.shape == (5, 5, 5)
    assert padded[0, 0, 0] == 0
    assert padded[2, 2, 2] == 1
    assert padded.sum() == 27


@pytest.mark.parametrize("spacing", [0.0, -0.1, float("nan"), float("inf")])
def test_resample_rejects_invalid_target_spacing(spacing):
    with pytest.raises(ValueError, match="finite and greater than zero"):
        segment.resample_isotropic(_two_blobs(), spacing)


def test_resample_rejects_an_excessive_target_grid_before_allocating():
    with pytest.raises(ValueError, match="resampled grid would hold"):
        segment.resample_isotropic(_two_blobs(), 0.001)
