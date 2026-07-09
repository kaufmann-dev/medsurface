"""Print profiles: printability morphology applied where the mask lives.

Doing this on the mask, rather than on the finished STL, is the whole point.
A mesh -> voxel -> mesh round trip costs fidelity before any morphology runs at
all: re-rasterising a good surface onto a hard binary grid, then running marching
cubes, windowed-sinc smoothing and quadric decimation a second time. Measured on a
600k-triangle skull, with the morphology switched off entirely, the round trip
still moved the surface by up to 0.40 mm and dropped the mean dihedral angle from
11.90 to 10.22 degrees.
"""

import math

import numpy as np
import pytest
import SimpleITK as sitk

from dicom_surface import pipeline, presets, segment
from dicom_surface.geometry import dilation_extent_mm, dilation_radius_voxels
from dicom_surface.presets import ANATOMICAL, PRINT_PROFILES


HEAD_CT = (0.315, 0.315, 0.8)


def _slab(thickness_voxels, shape=(20, 40, 40), spacing=(1.0, 1.0, 1.0)):
    arr = np.zeros(shape, dtype=np.uint8)
    mid = shape[1] // 2
    arr[3:shape[0] - 3, mid:mid + thickness_voxels, 3:shape[2] - 3] = 1
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(spacing)
    return img


def _box(shape=(24, 24, 24), spacing=(1.0, 1.0, 1.0)):
    arr = np.zeros(shape, dtype=np.uint8)
    arr[8:16, 8:16, 8:16] = 1
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(spacing)
    return img


# ------------------------------------------------------- the units contract
def test_thickening_never_undershoots():
    """A dilation is a distance the surface must move, not a kernel extent.

    ``kernel_radius_voxels`` floors, because overshooting a filter kernel erases
    anatomy. ``dilation_radius_voxels`` ceils, because undershooting a thickening
    leaves a wall too thin to print. Opposite rounding, on purpose.
    """
    for mm in (0.2, 0.4, 0.8, 1.2, 2.0):
        for spacing in (HEAD_CT, (1.0, 1.0, 1.0), (0.4, 0.4, 0.4)):
            radii = dilation_radius_voxels(mm, spacing)
            grown = dilation_extent_mm(radii, spacing)
            assert all(r >= 1 for r in radii), (mm, spacing, radii)
            for g in grown:
                assert g >= mm - 1e-9, (mm, spacing, radii, g)


def test_thickening_of_zero_is_a_no_op():
    assert dilation_radius_voxels(0.0, HEAD_CT) == [0, 0, 0]
    assert dilation_radius_voxels(-1.0, HEAD_CT) == [0, 0, 0]


def test_a_sub_voxel_request_still_grows_one_voxel():
    """0.4 mm cannot be realised on 0.8 mm slices; it must round up, not vanish."""
    radii = dilation_radius_voxels(0.4, HEAD_CT)
    assert radii[2] == 1
    assert dilation_extent_mm(radii, HEAD_CT)[2] == pytest.approx(0.8)


def test_dilate_grows_the_object_by_the_requested_radius():
    image = _box(spacing=(1.0, 1.0, 1.0))
    before = np.count_nonzero(sitk.GetArrayViewFromImage(image))
    grown = segment.dilate(image, 2.0)
    after = np.count_nonzero(sitk.GetArrayViewFromImage(grown))
    assert after > before

    arr = sitk.GetArrayViewFromImage(grown)
    idx = np.argwhere(arr)
    # an 8-voxel box dilated by radius 2 spans 8 + 2*2 = 12 on every axis
    assert (idx.max(0) - idx.min(0) + 1).tolist() == [12, 12, 12]


def test_dilate_of_zero_returns_the_input_untouched():
    image = _box()
    assert segment.dilate(image, 0.0) is image


# ---------------------------------------------------- the thin-feature metric
def test_thin_fraction_flags_a_thin_slab_and_clears_a_thick_one():
    """Thin material is what a ball of the minimum feature size cannot reach."""
    thin = _slab(thickness_voxels=1)     # 1 mm wall
    thick = _slab(thickness_voxels=8)    # 8 mm wall

    assert segment.thin_fraction(thin, 3.0) == pytest.approx(1.0)
    assert segment.thin_fraction(thick, 3.0) < 0.5


def test_thin_fraction_is_zero_when_not_asked_for():
    assert segment.thin_fraction(_slab(1), 0.0) == 0.0


def test_thickening_reduces_the_thin_fraction():
    """The point of the whole exercise."""
    wall = _slab(thickness_voxels=1)
    before = segment.thin_fraction(wall, 3.0)
    after = segment.thin_fraction(segment.dilate(wall, 2.0), 3.0)
    assert before > after
    assert after < 0.5


# ------------------------------------------------------------- the profiles
def test_anatomical_is_a_true_identity():
    """Every field is a no-op at zero, so there is no 'printing disabled' branch
    that could drift out of sync with the enabled one."""
    assert ANATOMICAL.closing_mm == 0.0
    assert ANATOMICAL.thicken_mm == 0.0
    assert ANATOMICAL.min_island_mm3 == 0.0
    assert ANATOMICAL.min_feature_mm == 0.0
    assert PRINT_PROFILES["anatomical"] == ANATOMICAL


def test_anatomical_mask_matches_the_profile_free_mask():
    """The hard requirement: the default output must not move."""
    rng = np.random.default_rng(0)
    arr = (rng.random((24, 24, 24)) * 1000 - 200).astype(np.float32)
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((0.5, 0.5, 0.5))
    preset = presets.get("bone")

    default = pipeline.build_mask(image, preset, 300.0)
    explicit = pipeline.build_mask(image, preset, 300.0, ANATOMICAL)
    assert np.array_equal(sitk.GetArrayViewFromImage(default),
                          sitk.GetArrayViewFromImage(explicit))


@pytest.mark.parametrize("name", ["fdm", "resin"])
def test_profiles_only_add_printability(name):
    """Composition is max(), never assignment: a profile may raise a preset's
    closing but never relax it. `skin` already closes 3.2 mm."""
    profile = PRINT_PROFILES[name]
    skin = presets.get("skin")
    assert max(skin.closing_mm, profile.closing_mm) >= skin.closing_mm
    assert profile.thicken_mm > 0
    assert profile.min_feature_mm > 0


def test_fdm_is_coarser_than_resin():
    fdm, resin = PRINT_PROFILES["fdm"], PRINT_PROFILES["resin"]
    assert fdm.closing_mm > resin.closing_mm
    assert fdm.thicken_mm > resin.thicken_mm
    assert fdm.min_island_mm3 >= resin.min_island_mm3
    assert fdm.min_feature_mm > resin.min_feature_mm


def test_unknown_profile_lists_alternatives():
    with pytest.raises(KeyError, match="available:"):
        presets.get_print_profile("sla")


def test_build_mask_seals_and_thickens_under_a_profile():
    """A pitted wall: the pore is sealed by the closing and the wall grown by the
    dilation.

    Note the phantom must be bigger than the profile's island floor (100 mm3 for
    resin), or the profile correctly deletes it as an unprintable fragment.
    """
    arr = np.zeros((28, 28, 28), dtype=np.int16) - 1000
    arr[8:18, 10:16, 8:18] = 1000     # 10 x 6 x 10 mm slab of bone = 600 mm3
    arr[12:14, 10:16, 12:14] = -1000  # punched clean through: a pore
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))

    preset = presets.override(presets.get("bone"), median_mm=0.0)
    plain = pipeline.build_mask(image, preset, 300.0)
    printed = pipeline.build_mask(image, preset, 300.0, PRINT_PROFILES["resin"])

    plain_n = np.count_nonzero(sitk.GetArrayViewFromImage(plain))
    printed_n = np.count_nonzero(sitk.GetArrayViewFromImage(printed))
    assert plain_n > 0 and printed_n > plain_n, "the resin profile must add material"

    # the pore is gone
    assert sitk.GetArrayViewFromImage(printed)[13, 13, 13] == 1


def test_a_profile_deletes_fragments_it_cannot_print():
    """A 30 mm3 speck is below the resin profile's 100 mm3 island floor."""
    arr = np.zeros((28, 28, 28), dtype=np.int16) - 1000
    arr[8:18, 8:18, 8:18] = 1000   # 1000 mm3 body
    arr[2:5, 2:5, 2:5] = 1000      # 27 mm3 speck, far away
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))

    preset = presets.override(presets.get("bone"), median_mm=0.0, min_island_mm3=1.0,
                              keep_largest_island=False)
    plain = pipeline.build_mask(image, preset, 300.0)
    printed = pipeline.build_mask(image, preset, 300.0, PRINT_PROFILES["resin"])

    assert sitk.GetArrayViewFromImage(plain)[3, 3, 3] == 1
    assert sitk.GetArrayViewFromImage(printed)[3, 3, 3] == 0


def test_printability_warnings_are_silent_for_anatomical():
    assert pipeline.printability_warnings(_box(), ANATOMICAL) == []


def test_printability_warnings_report_thickening_and_thin_material():
    wall = _slab(thickness_voxels=1)
    profile = presets.PrintProfile(name="t", description="d",
                                   thicken_mm=1.0, min_feature_mm=3.0)
    messages = " ".join(pipeline.printability_warnings(wall, profile))
    assert "thickened" in messages or "thickening" in messages
    assert "thinner than" in messages


def test_anisotropic_thickening_is_reported_not_hidden():
    """On 0.8 mm slices a 0.4 mm request realises as 0.8 mm along z. Say so."""
    image = _box(spacing=HEAD_CT)
    profile = presets.PrintProfile(name="t", description="d", thicken_mm=0.4)
    messages = " ".join(pipeline.printability_warnings(image, profile))
    assert "realised as" in messages


# ------------------------------------------------------------------- merge
class _Stop(Exception):
    """Halts `merge` at a chosen stage so the masks built so far can be inspected."""


def test_merge_registers_on_anatomy_not_on_the_print_mask(monkeypatch):
    """The safety invariant.

    Thickening both scans inflates Dice and surface overlap -- the very numbers
    the impostor gates are calibrated against -- so registration must never see a
    print mask. Confirmed end to end too: the bar phantom is refused with exactly
    the same 24.5% overlap and 0.314 Dice whether or not `--print-profile fdm` is
    given.
    """
    from dicom_surface import merge as merge_mod
    from dicom_surface.series import Series
    from dicom_surface.volume import Volume

    arr = np.zeros((16, 16, 16), dtype=np.int16) - 1000
    arr[4:12, 4:12, 4:12] = 1000
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))

    seen: list[str] = []
    real_build_mask = pipeline.build_mask

    def spy(img, preset, threshold, profile=ANATOMICAL, log=None):
        seen.append(profile.name)
        return real_build_mask(img, preset, threshold, profile, log)

    monkeypatch.setattr(merge_mod, "check_compatible", lambda a, b, force=False: [])
    monkeypatch.setattr(merge_mod.volume_mod, "load", lambda s: Volume(image=image, series=s))
    monkeypatch.setattr(merge_mod.volume_mod, "warnings_for", lambda v: [])
    monkeypatch.setattr(merge_mod.pipeline, "build_mask", spy)

    # Registration always succeeds and always passes its gates, so the only thing
    # under test is which masks reach which stage.
    def fake_register(mask_a, mask_b, **kw):
        return registration_result()

    def registration_result():
        from dicom_surface.registration import RegistrationResult
        return RegistrationResult(
            transform=np.eye(4), fft_translation_mm=np.zeros(3), rotation_deg=0.0,
            translation_mm=np.zeros(3), inlier_rms_mm=0.1, inlier_median_mm=0.1,
            overlap_moving_in_fixed=0.99, overlap_fixed_in_moving=0.99,
            shared_fov_dice=0.9, shared_fov_mm3=1e6,
        )

    def stop(*_a, **_k):
        raise _Stop

    monkeypatch.setattr(merge_mod.registration, "rigid_register", fake_register)
    monkeypatch.setattr(merge_mod, "_common_grid", stop)

    def run(profile):
        seen.clear()
        a = Series(uid="a", modality="CT", description="fixed", series_number=1)
        b = Series(uid="b", modality="CT", description="moving", series_number=2)
        with pytest.raises(_Stop):
            merge_mod.merge(series_a=a, series_b=b, preset=presets.get("bone"),
                            output_path="/dev/null", print_profile=profile)
        return list(seen)

    # A print profile: register on anatomy, then re-segment for the union.
    assert run(PRINT_PROFILES["fdm"]) == ["anatomical", "anatomical", "fdm", "fdm"]

    # No profile: the registration masks are the union masks. No second pass.
    assert run(ANATOMICAL) == ["anatomical", "anatomical"]
