"""Thresholding, morphology and island handling."""

import numpy as np
import pytest
import SimpleITK as sitk

from medsurface import segment, surface

from .mesh_helpers import mean_adjacency_angle, volume


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


def test_physical_occupancy_smoothing_removes_anisotropic_terracing():
    spacing = (0.4, 0.4, 0.8)
    shape = (40, 80, 80)
    zz, yy, xx = np.meshgrid(
        np.arange(shape[0]) * spacing[2],
        np.arange(shape[1]) * spacing[1],
        np.arange(shape[2]) * spacing[0],
        indexing="ij",
    )
    values = (
        (zz - 16.0) ** 2 + (yy - 16.0) ** 2 + (xx - 16.0) ** 2 <= 12.0**2
    ).astype(np.uint8)
    image = sitk.GetImageFromArray(values)
    image.SetSpacing(spacing)
    image = segment.pad(image, 4)
    affine = surface.index_to_physical(image)

    raw = surface.transform(surface.marching_cubes(image), affine)
    relaxed = surface.smooth(raw, 20, 0.1)
    field = segment.smooth_occupancy(image, 0.8)
    physically_smoothed = surface.transform(
        surface.marching_cubes(field, segment.ISO_OCCUPANCY),
        affine,
    )
    physically_smoothed = surface.smooth(physically_smoothed, 20, 0.1)

    relaxed_arrays = surface.to_arrays(relaxed)
    smoothed_arrays = surface.to_arrays(physically_smoothed)
    assert mean_adjacency_angle(smoothed_arrays) < 0.6 * mean_adjacency_angle(
        relaxed_arrays
    )
    assert volume(smoothed_arrays) == pytest.approx(volume(relaxed_arrays), rel=0.02)


@pytest.mark.parametrize("sigma", [-0.1, float("nan"), float("inf")])
def test_occupancy_smoothing_rejects_invalid_sigma(sigma):
    with pytest.raises(ValueError, match="finite and non-negative"):
        segment.smooth_occupancy(_two_blobs(), sigma)


def test_zero_occupancy_smoothing_is_a_noop():
    image = _two_blobs()
    assert segment.smooth_occupancy(image, 0) is image


def test_occupancy_smoothing_supports_minimum_volume_dimensions():
    image = sitk.Image((2, 2, 2), sitk.sitkUInt8)
    assert segment.smooth_occupancy(image, 0.8).GetSize() == (2, 2, 2)


@pytest.mark.parametrize("spacing", [0.0, -0.1, float("nan"), float("inf")])
def test_resample_rejects_invalid_target_spacing(spacing):
    with pytest.raises(ValueError, match="finite and greater than zero"):
        segment.resample_isotropic(_two_blobs(), spacing)


def test_resample_rejects_an_excessive_target_grid_before_allocating():
    with pytest.raises(ValueError, match="resampled grid would hold"):
        segment.resample_isotropic(_two_blobs(), 0.001)


def test_resample_accepts_exactly_the_default_voxel_limit(monkeypatch):
    image = sitk.Image((2, 2, 2), sitk.sitkUInt8)
    image.SetSpacing((250.0, 500.0, 500.0))

    class AllocationReached(Exception):
        pass

    monkeypatch.setattr(
        segment,
        "antialias_for_grid",
        lambda _image, _spacing: (_ for _ in ()).throw(AllocationReached),
    )

    with pytest.raises(AllocationReached):
        segment.resample_isotropic(image, 1.0)


def test_resample_override_reaches_the_large_grid_allocation(monkeypatch):
    class AllocationReached(Exception):
        pass

    def stop_before_allocation(_image, _spacing):
        raise AllocationReached

    monkeypatch.setattr(segment, "antialias_for_grid", stop_before_allocation)

    with pytest.raises(AllocationReached):
        segment.resample_isotropic(
            _two_blobs(),
            0.001,
            allow_large_volume=True,
        )
