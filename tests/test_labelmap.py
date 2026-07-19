"""External labelmap validation, extraction, and registered binary fusion."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk

from medsurface import catalog, fusion, labelmap, registration


def _image(values: np.ndarray, *, origin=(0.0, 0.0, 0.0)) -> sitk.Image:
    image = sitk.GetImageFromArray(values)
    image.SetSpacing((1.0, 1.0, 1.0))
    image.SetOrigin(origin)
    return image


def _write(path: Path, values: np.ndarray, *, origin=(0.0, 0.0, 0.0)):
    sitk.WriteImage(_image(values, origin=origin), str(path))
    return catalog.discover(path)[0]


def _two_labels(dtype=np.uint8) -> np.ndarray:
    values = np.zeros((28, 28, 28), dtype=dtype)
    values[3:11, 3:11, 3:11] = 5
    values[17:25, 17:25, 17:25] = 91
    return values


def test_binary_mask_unions_every_nonzero_label():
    mask = labelmap._binary_mask(_image(_two_labels()))
    values = sitk.GetArrayViewFromImage(mask)

    assert set(np.unique(values)) == {0, 1}
    assert int(values.sum()) == 2 * 8**3


def test_integer_valued_float_labelmap_is_supported():
    mask = labelmap._binary_mask(_image(_two_labels(np.float32)))
    assert int(sitk.GetArrayViewFromImage(mask).sum()) == 2 * 8**3


@pytest.mark.parametrize(
    "values,message",
    [
        (np.zeros((4, 4, 4), dtype=np.uint8), "no nonzero foreground"),
        (np.full((4, 4, 4), -1, dtype=np.int16), "non-negative"),
        (np.full((4, 4, 4), 0.5, dtype=np.float32), "fractional"),
        (np.full((4, 4, 4), np.nan, dtype=np.float32), "non-finite"),
        (np.full((4, 4, 4), np.inf, dtype=np.float32), "non-finite"),
    ],
)
def test_invalid_labelmap_values_are_rejected(values, message):
    with pytest.raises(ValueError, match=message):
        labelmap._binary_mask(_image(values))


@pytest.mark.parametrize(
    "shape,dimensions",
    [
        ((2, 2, 2), "2x2x2"),
        ((4, 4, 3), "3x4x4"),
        ((4, 3, 4), "4x3x4"),
        ((3, 4, 4), "4x4x3"),
    ],
)
def test_labelmap_extract_rejects_short_dimensions_before_meshing(
    tmp_path, monkeypatch, shape, dimensions
):
    candidate = _write(tmp_path / "short.nii.gz", np.ones(shape, dtype=np.uint8))
    monkeypatch.setattr(
        labelmap.pipeline,
        "mesh_binary_mask",
        lambda *_args, **_kwargs: pytest.fail("meshing must not start"),
    )

    with pytest.raises(
        ValueError,
        match=r"each volume axis must contain at least 4 voxels; got %s$" % dimensions,
    ):
        labelmap.extract(
            candidate,
            str(tmp_path / "out.stl"),
            mask_smooth_mm=0,
        )


@pytest.mark.parametrize("invalid_role", ["fixed", "moving"])
def test_labelmap_fuse_rejects_short_dimensions_before_registration(
    tmp_path, monkeypatch, invalid_role
):
    valid = _write(tmp_path / "valid.nii.gz", np.ones((4, 4, 4), dtype=np.uint8))
    short = _write(tmp_path / "short.nii.gz", np.ones((4, 4, 2), dtype=np.uint8))
    fixed, moving = (short, valid) if invalid_role == "fixed" else (valid, short)
    monkeypatch.setattr(
        labelmap.fusion,
        "fuse_masks",
        lambda *_args, **_kwargs: pytest.fail("registration must not start"),
    )

    with pytest.raises(
        ValueError,
        match=(r"each volume axis must contain at least 4 voxels; got 2x4x4$"),
    ):
        labelmap.fuse(
            fixed,
            moving,
            str(tmp_path / "out.nii.gz"),
        )


def test_four_voxel_labelmap_is_accepted_with_recursive_smoothing(tmp_path):
    candidate = _write(
        tmp_path / "minimum.nii.gz",
        np.ones((4, 4, 4), dtype=np.uint8),
    )
    output = tmp_path / "minimum.stl"

    result = labelmap.extract(
        candidate,
        str(output),
        simplify_error_mm=0,
    )

    assert output.exists()
    assert result.quality["valid"]


def test_default_finishing_is_independent_and_keeps_all_shells():
    settings = labelmap.default_surface_settings()

    assert settings.resample_mm == 0
    assert settings.mask_smooth_mm == pytest.approx(0.8)
    assert settings.surface_smooth_iters == 20
    assert settings.simplify_error_mm == pytest.approx(0.25)
    assert settings.post_surface_smooth_iters == 0
    assert settings.keep_largest_component is False


def test_disabling_mask_smoothing_keeps_surface_smoothing():
    settings = labelmap.resolve_surface_settings(mask_smooth_mm=0)

    assert settings.mask_smooth_mm == 0
    assert settings.surface_smooth_iters == 20
    assert settings.simplify_error_mm == pytest.approx(0.25)
    assert settings.post_surface_smooth_iters == 0


def test_component_selection_can_override_the_labelmap_default():
    settings = labelmap.resolve_surface_settings(keep_largest_component=True)

    assert settings.keep_largest_component is True


def test_default_labelmap_finishing_uses_no_post_simplification_relaxation(tmp_path):
    source = tmp_path / "labels.nii.gz"
    output = tmp_path / "labels.stl"
    candidate = _write(source, _two_labels())

    result = labelmap.extract(candidate, str(output))
    finishing = result.provenance["surface_finishing"]

    assert result.quality["valid"]
    assert result.provenance["surface"]["mask_smooth_mm"] == pytest.approx(0.8)
    assert result.provenance["surface"]["surface_smooth_iters"] == 20
    assert result.provenance["surface"]["post_surface_smooth_iters"] == 0
    assert finishing["pre_smoothing"]["requested_iterations"] == 20
    assert finishing["decimation"]["simplify_error_mm"] == pytest.approx(0.25)
    assert finishing["decimation"]["error_introduced_mm"] <= 0.25
    assert finishing["post_smoothing"]["requested_iterations"] == 0
    assert set(finishing) == {"pre_smoothing", "decimation", "post_smoothing"}
    assert any("Gaussian sigma of 0.80 mm" in warning for warning in result.warnings)


@pytest.mark.parametrize("extension", [".nii.gz", ".nrrd", ".mha"])
def test_real_labelmap_formats_extract_all_components(tmp_path, extension):
    source = tmp_path / ("labels" + extension)
    output = tmp_path / ("labels-%s.stl" % extension.replace(".", ""))
    candidate = _write(source, _two_labels())

    result = labelmap.extract(
        candidate,
        str(output),
        mask_smooth_mm=0,
        surface_smooth_iters=0,
        simplify_error_mm=0,
    )

    assert output.exists()
    assert result.quality["valid"]
    assert result.labelmap_components == 2
    assert result.surface_components == 2
    assert result.quality["components"] == 2
    assert result.provenance["input"]["foreground"] == "all nonzero voxels"
    assert result.provenance["surface"]["keep_largest_component"] is False
    assert (
        result.provenance["surface_finishing"]["pre_smoothing"]["requested_iterations"]
        == 0
    )


def test_labelmap_extract_can_keep_only_the_largest_component(tmp_path):
    source = tmp_path / "labels.nii.gz"
    output = tmp_path / "largest.stl"
    candidate = _write(source, _two_labels())

    result = labelmap.extract(
        candidate,
        str(output),
        mask_smooth_mm=0,
        surface_smooth_iters=0,
        simplify_error_mm=0,
        keep_largest_component=True,
    )

    assert result.quality["valid"]
    assert result.surface_components == 2
    assert result.quality["components"] == 1
    assert result.provenance["surface"]["keep_largest_component"] is True


def test_labelmap_extraction_preserves_left_handed_physical_geometry(tmp_path):
    source = tmp_path / "left-handed.nii.gz"
    output = tmp_path / "left-handed.stl"
    image = _image(_two_labels())
    image.SetDirection((-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    sitk.WriteImage(image, str(source))

    result = labelmap.extract(
        catalog.discover(source)[0],
        str(output),
        mask_smooth_mm=0,
        surface_smooth_iters=0,
        simplify_error_mm=0,
    )

    assert result.quality["valid"]
    assert result.quality["winding_consistent"]
    assert result.quality["disoriented_faces"] == 0


def test_labelmap_boundary_is_capped_and_reported(tmp_path):
    values = np.zeros((16, 16, 16), dtype=np.uint8)
    values[:8, 4:12, 4:12] = 1
    source = tmp_path / "clipped.nii.gz"
    output = tmp_path / "clipped.stl"
    candidate = _write(source, values)

    result = labelmap.extract(
        candidate,
        str(output),
        mask_smooth_mm=0,
        surface_smooth_iters=0,
        simplify_error_mm=0,
    )

    assert result.capped_field_of_view
    assert result.quality["watertight"]
    assert any("capped flat" in warning for warning in result.warnings)


def test_labelmap_fuse_registers_rigid_masks_and_uses_fixed_frame(tmp_path):
    zz, yy, xx = np.indices((64, 64, 64))
    values = (
        ((xx - 32) / 18) ** 2 + ((yy - 31) / 15) ** 2 + ((zz - 30) / 12) ** 2 <= 1
    ).astype(np.uint8)
    values[25:37, 28:35, 45:50] = 1
    fixed_path = tmp_path / "fixed.nii.gz"
    moving_path = tmp_path / "moving.nii.gz"
    fixed = _write(fixed_path, values * 5)
    moving = _write(moving_path, values * 91, origin=(4.0, -3.0, 2.0))
    output = tmp_path / "fused.nrrd"

    result = labelmap.fuse(
        fixed,
        moving,
        str(output),
        grid_mm=1.0,
    )

    assert output.exists()
    stored = sitk.ReadImage(str(output))
    assert stored.GetPixelID() == sitk.sitkUInt8
    assert set(np.unique(sitk.GetArrayViewFromImage(stored))) == {0, 1}
    assert int(sitk.GetArrayViewFromImage(stored).sum()) == result.foreground_fused_voxels
    assert result.registration.surface_overlap > 0.9
    assert result.registration.shared_fov_dice > 0.9
    assert "fixed labelmap" in result.provenance["coordinate_system"]
    assert result.provenance["fixed"]["foreground"] == "all nonzero voxels"
    assert result.provenance["moving"]["foreground"] == "all nonzero voxels"
    assert result.provenance["segmentation"] == {
        "input_kind": "labelmap",
        "validation": "finite, discrete, and non-negative",
        "foreground": "all nonzero source values normalized to 1",
    }
    assert "surface" not in result.provenance
    assert "surface_finishing" not in result.provenance
    assert any("rigid registration" in warning for warning in result.warnings)
    assert any("same rigid structures" in warning for warning in result.warnings)


def _registration_result() -> registration.RegistrationResult:
    return registration.RegistrationResult(
        transform=np.eye(4),
        fft_translation_mm=np.zeros(3),
        rotation_deg=0.0,
        translation_mm=np.zeros(3),
        inlier_rms_mm=0.1,
        inlier_median_mm=0.1,
        overlap_moving_in_fixed=0.95,
        overlap_fixed_in_moving=0.95,
        shared_fov_dice=0.9,
        shared_fov_mm3=100_000.0,
    )


@pytest.mark.parametrize(
    "extension,expected_format,expected_compression",
    [
        (".nii", "NIfTI", "none"),
        (".nii.gz", "NIfTI", "gzip"),
        (".nrrd", "NRRD", "gzip"),
        (".mha", "MetaImage", "zlib"),
    ],
)
def test_labelmap_fuse_normalizes_multilabel_values_in_every_output_format(
    tmp_path,
    monkeypatch,
    extension,
    expected_format,
    expected_compression,
):
    fixed_values = _two_labels()
    moving_values = np.zeros_like(fixed_values)
    moving_values[fixed_values == 5] = 17
    moving_values[fixed_values == 91] = 203
    fixed = _write(tmp_path / "fixed.nii.gz", fixed_values)
    moving = _write(tmp_path / "moving.mha", moving_values)
    output = tmp_path / ("fused" + extension)
    monkeypatch.setattr(
        fusion.registration,
        "rigid_register",
        lambda *_args, **_kwargs: _registration_result(),
    )

    result = labelmap.fuse(fixed, moving, str(output), grid_mm=1.0)
    stored = sitk.ReadImage(str(output))
    stored_values = sitk.GetArrayViewFromImage(stored)

    assert isinstance(result, fusion.FusionResult)
    assert result.output_format == expected_format
    assert result.compression == expected_compression
    assert set(np.unique(stored_values)) == {0, 1}
    assert 5 not in stored_values
    assert 17 not in stored_values
    assert 91 not in stored_values
    assert 203 not in stored_values
    assert result.foreground_fixed_voxels == 2 * 8**3
    assert result.foreground_moving_voxels == 2 * 8**3
    assert result.foreground_fused_voxels == 2 * 8**3
    assert result.provenance["output"]["kind"] == "binary labelmap"
    assert "surface" not in result.provenance
    assert "surface_finishing" not in result.provenance
    assert any("fused labelmap" in warning for warning in result.warnings)


@pytest.mark.parametrize(
    "invalid_values,message",
    [
        (np.zeros((4, 4, 4), dtype=np.uint8), "no nonzero foreground"),
        (np.full((4, 4, 4), 0.25, dtype=np.float32), "fractional"),
    ],
)
def test_labelmap_fuse_rejects_invalid_moving_values_before_registration(
    tmp_path, monkeypatch, invalid_values, message
):
    valid = _write(tmp_path / "valid.nii.gz", np.ones((4, 4, 4), dtype=np.uint8))
    invalid = _write(tmp_path / "invalid.nrrd", invalid_values)
    monkeypatch.setattr(
        fusion,
        "fuse_masks",
        lambda *_args, **_kwargs: pytest.fail("registration must not start"),
    )

    with pytest.raises(ValueError, match=message):
        labelmap.fuse(valid, invalid, str(tmp_path / "out.nii"))


def test_labelmap_fuse_rejects_the_same_input_before_loading(tmp_path, monkeypatch):
    source = tmp_path / "mask.nii.gz"
    candidate = _write(source, _two_labels())
    monkeypatch.setattr(
        labelmap,
        "load",
        lambda *_args, **_kwargs: pytest.fail(
            "duplicate inputs must fail before loading"
        ),
    )

    with pytest.raises(fusion.FusionError, match="same labelmap"):
        labelmap.fuse(candidate, candidate, str(tmp_path / "out.nii"))


def test_labelmap_extract_does_not_call_segmentation(monkeypatch, tmp_path):
    source = tmp_path / "mask.nii.gz"
    candidate = _write(source, _two_labels())

    from medsurface import segment

    for name in ("binarize", "islands", "median", "opening", "closing"):
        monkeypatch.setattr(
            segment,
            name,
            lambda *_args, _name=name, **_kwargs: pytest.fail(
                "labelmap extraction called segmentation step %s" % _name
            ),
        )

    result = labelmap.extract(
        candidate,
        str(tmp_path / "out.stl"),
        mask_smooth_mm=0,
        surface_smooth_iters=0,
        simplify_error_mm=0,
    )
    assert result.quality["valid"]
