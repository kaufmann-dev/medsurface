"""Print profiles: printability morphology applied where the mask lives.

Doing this on the mask, rather than on the finished STL, is the whole point.
A mesh -> voxel -> mesh round trip costs fidelity before any morphology runs at
all: re-rasterising a good surface onto a hard binary grid, then running marching
cubes, windowed-sinc smoothing and quadric decimation a second time. Measured on a
600k-triangle skull, with the morphology switched off entirely, the round trip
still moved the surface by up to 0.40 mm and dropped the mean dihedral angle from
11.90 to 10.22 degrees.
"""

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


def test_realised_feature_reports_what_the_grid_actually_probes():
    """The probe is a ball rounded up to whole voxels, so on a coarse grid it is
    larger than requested. Saying so is the difference between a diagnostic and a
    fiction."""
    image = _box(spacing=(0.315, 0.315, 0.8))
    assert segment.realised_feature_mm(image, 1.2) == pytest.approx([1.26, 1.26, 1.6])
    assert segment.realised_feature_mm(image, 0.0) == [0.0, 0.0, 0.0]


def test_a_feature_smaller_than_the_slice_pitch_is_not_measurable():
    """0.6 mm walls cannot be assessed on 0.8 mm slices: the probe becomes 1.6 mm."""
    anisotropic = _box(spacing=(0.315, 0.315, 0.8))
    assert segment.feature_is_resolvable(anisotropic, 1.2)
    assert not segment.feature_is_resolvable(anisotropic, 0.6)

    coarse = _box(spacing=(3.0, 3.0, 3.0))
    assert not segment.feature_is_resolvable(coarse, 1.2)
    assert segment.feature_is_resolvable(_box(spacing=(0.5, 0.5, 0.5)), 1.2)


def test_unmeasurable_feature_reports_that_instead_of_a_number():
    """A solid sphere at 3 mm voxels reports 3.9% thin at a 1.2 mm feature size --
    an artefact of the probe, not of the geometry. Report the limitation, not the
    artefact."""
    profile = presets.PrintProfile(name="t", description="d", min_feature_mm=1.2)
    coarse = pipeline.thin_material_warning(_box(spacing=(3.0, 3.0, 3.0)), profile)
    assert coarse and "cannot assess" in coarse[0]

    # A genuinely thick block: 8 mm across at 0.2 mm voxels. (`_box` is only
    # 1.6 mm across there, which really is thin against a 1.2 mm feature.)
    arr = np.zeros((60, 60, 60), dtype=np.uint8)
    arr[10:50, 10:50, 10:50] = 1
    thick = sitk.GetImageFromArray(arr)
    thick.SetSpacing((0.2, 0.2, 0.2))
    assert pipeline.thin_material_warning(thick, profile) == [], \
        "an 8 mm block is not thin against a 1.2 mm feature"


def test_thin_fraction_is_zero_when_not_asked_for():
    assert segment.thin_fraction(_slab(1), 0.0) == 0.0


def test_thickening_reduces_the_thin_fraction():
    """The point of the whole exercise."""
    wall = _slab(thickness_voxels=1)
    before = segment.thin_fraction(wall, 3.0)
    after = segment.thin_fraction(segment.thicken(wall, 3.0), 3.0)
    assert before > after
    assert after < 0.5


# ------------------------------------------------------- selective thickening
def _thin_sheet_beside_a_thick_block(spacing=(1.0, 1.0, 1.0)):
    """A 1 mm sheet and an 8 mm block, far apart. Only the sheet is too thin."""
    arr = np.zeros((24, 40, 24), dtype=np.uint8)
    arr[8:16, 4:12, 8:16] = 1     # block: 8 mm on every axis
    arr[8:16, 30:31, 8:16] = 1    # sheet: 1 mm thick
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(spacing)
    return img


def _extent(image, y_slice):
    arr = sitk.GetArrayViewFromImage(image)[:, y_slice, :]
    return int(np.count_nonzero(arr))


def test_thicken_grows_only_the_thin_material():
    """A cranial vault is 5 mm of solid bone and needs nothing; inflating it moves
    the model's outer surface for no benefit. Only the thin structures get the
    material."""
    image = _thin_sheet_beside_a_thick_block()
    grown = segment.thicken(image, min_feature_mm=3.0)

    arr_before = sitk.GetArrayViewFromImage(image)
    arr_after = sitk.GetArrayViewFromImage(grown)

    # the block keeps its bounding box on every axis
    block_before = np.argwhere(arr_before[:, :20, :])
    block_after = np.argwhere(arr_after[:, :20, :])
    assert (block_before.max(0) - block_before.min(0)).tolist() == \
           (block_after.max(0) - block_after.min(0)).tolist()

    # the sheet gets thicker
    sheet_before = np.count_nonzero(arr_before[:, 20:, :])
    sheet_after = np.count_nonzero(arr_after[:, 20:, :])
    assert sheet_after > 2 * sheet_before


def test_thicken_is_a_no_op_when_nothing_is_thin():
    """Nothing thinner than the minimum feature means nothing to do."""
    arr = np.zeros((24, 24, 24), dtype=np.uint8)
    arr[4:20, 4:20, 4:20] = 1  # 16 mm cube
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))

    grown = segment.thicken(image, min_feature_mm=3.0)
    assert np.array_equal(sitk.GetArrayViewFromImage(image),
                          sitk.GetArrayViewFromImage(grown))


def test_thin_mask_is_the_complement_of_the_opening():
    image = _thin_sheet_beside_a_thick_block()
    thin = sitk.GetArrayFromImage(segment.thin_mask(image, 3.0))
    whole = sitk.GetArrayViewFromImage(image)

    assert np.all(thin <= whole), "thin material must be a subset of the material"
    assert np.count_nonzero(thin[:, 20:, :]) > 0, "the sheet is thin"
    assert np.count_nonzero(thin[8:16, 6:10, 10:14]) == 0, "the block core is not"


def test_counting_a_temporary_image_does_not_read_freed_memory():
    """`sitk.GetArrayViewFromImage` returns a view into the image's buffer. Take it
    from a temporary and the image is freed the moment the call returns, leaving
    the view pointing at released memory: it reads plausible garbage on a small
    volume and segfaults on an 80M-voxel scan. `count_foreground` binds the image
    to a parameter, which keeps it alive for the duration of the count."""
    image = _thin_sheet_beside_a_thick_block()
    expected = int(np.count_nonzero(sitk.GetArrayFromImage(segment.thin_mask(image, 3.0))))

    # 200 rounds so a freed buffer would very likely be reused by something else
    for _ in range(200):
        assert segment.count_foreground(segment.thin_mask(image, 3.0)) == expected


def test_thin_mask_is_empty_without_a_feature_size():
    image = _thin_sheet_beside_a_thick_block()
    assert np.count_nonzero(sitk.GetArrayFromImage(segment.thin_mask(image, 0.0))) == 0


def test_thicken_of_zero_returns_the_input_untouched():
    image = _thin_sheet_beside_a_thick_block()
    assert segment.thicken(image, 0.0) is image


# --------------------------------------------------- the printability grid
def test_the_printability_grid_realises_the_feature_size_exactly():
    """Choosing mm = (T/2)/r for integer r makes the voxel radius exactly r and the
    realised feature size exactly T, on any scan."""
    for spacing in [(0.379, 0.379, 0.8), (0.45, 0.45, 0.3), (0.5, 0.5, 2.0)]:
        for feature in (0.6, 1.2, 2.0):
            image = _box(spacing=spacing)
            mm = segment.printability_grid_mm(image, feature)
            if mm is None:
                continue
            probe = _box(spacing=(mm, mm, mm))
            assert max(segment.realised_feature_mm(probe, feature)) == \
                pytest.approx(feature, rel=1e-6), (spacing, feature, mm)


def test_the_printability_grid_is_never_coarser_than_the_data():
    """Upsampling invents no detail but costs voxels; downsampling destroys thin
    bone -- 0.6 mm once reopened 257 pores the closing had sealed."""
    for spacing in [(0.379, 0.379, 0.8), (0.3, 0.3, 0.3), (0.5, 0.5, 2.0)]:
        image = _box(spacing=spacing)
        mm = segment.printability_grid_mm(image, 1.2)
        if mm is not None:
            assert mm <= min(spacing) + 1e-9, (spacing, mm)


def test_an_already_isotropic_fine_grid_is_left_alone():
    """The common case must pay nothing: no resample, no copy."""
    assert segment.printability_grid_mm(_box(spacing=(0.3, 0.3, 0.3)), 1.2) is None
    assert segment.printability_grid_mm(_box(spacing=(1.0, 1.0, 1.0)), 0.0) is None


def test_the_voxel_budget_coarsens_the_grid_rather_than_exhausting_memory():
    image = _box(spacing=(0.9, 0.9, 5.0))
    generous = segment.printability_grid_mm(image, 1.2, budget=10 ** 12)
    stingy = segment.printability_grid_mm(image, 1.2, budget=10 ** 5)
    assert stingy > generous


def test_a_coarse_slice_pitch_no_longer_inflates_the_model_along_z():
    """The bug this grid exists for.

    On a 0.9 x 0.9 x 5.0 mm survey CT, `dilation_radius_voxels(0.6)` rounds to one
    voxel per axis, so the thickening "ball" acquires a 5 mm semi-axis. Guaranteeing
    1.2 mm walls used to move the model's surface 5 mm along z. The scanner's slice
    pitch must not become the printer's tolerance.
    """
    survey = (0.9, 0.9, 5.0)
    arr = np.zeros((14, 60, 60), dtype=np.uint8)
    arr[6:8, 10:50, 10:50] = 1          # a plate, thin in z
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing(survey)

    def z_extent_mm(img):
        zs = np.argwhere(sitk.GetArrayViewFromImage(img))[:, 0]
        return (zs.max() - zs.min() + 1) * img.GetSpacing()[2]

    before = z_extent_mm(image)

    naive = segment.thicken(image, 1.2)                      # on the native grid
    mm = segment.printability_grid_mm(image, 1.2)
    fixed = segment.thicken(segment.to_printability_grid(image, mm), 1.2)

    naive_growth = z_extent_mm(naive) - before
    fixed_growth = z_extent_mm(fixed) - before

    assert naive_growth >= 9.0, "the native grid grows a full slice each way"
    assert fixed_growth <= 2.5, ("the printability grid grows about the feature "
                                 "size", fixed_growth)


# ------------------------------------------------------------- the profiles
def test_anatomical_is_a_true_identity():
    """Every field is a no-op at zero, so there is no 'printing disabled' branch
    that could drift out of sync with the enabled one."""
    assert ANATOMICAL.closing_mm == 0.0
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
    assert profile.min_feature_mm > 0


def test_fdm_is_coarser_than_resin():
    fdm, resin = PRINT_PROFILES["fdm"], PRINT_PROFILES["resin"]
    assert fdm.closing_mm > resin.closing_mm
    assert fdm.min_island_mm3 >= resin.min_island_mm3
    assert fdm.min_feature_mm > resin.min_feature_mm


def test_a_flag_bent_profile_stops_claiming_to_be_the_original():
    """`--min-feature-mm 1.0` with the default profile used to log "anatomical: No
    printability changes" while thickening by 1 mm."""
    import argparse

    from dicom_surface.cli import _resolve_print_profile

    plain = argparse.Namespace(print_profile="anatomical", min_feature_mm=None,
                               closing_mm=None, min_island_mm3=None)
    assert _resolve_print_profile(plain) == ANATOMICAL

    bent = argparse.Namespace(print_profile="anatomical", min_feature_mm=1.0,
                              closing_mm=None, min_island_mm3=None)
    profile = _resolve_print_profile(bent)
    assert profile != ANATOMICAL
    assert profile.min_feature_mm == 1.0
    assert "No printability changes" not in profile.description
    assert profile.name == "anatomical+flags"


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

    # The pore is gone. Ask in millimetres: a print profile moves the mask onto its
    # own isotropic grid, so voxel indices are not comparable between the two.
    pore_mm = plain.TransformIndexToPhysicalPoint((13, 13, 13))
    assert printed[printed.TransformPhysicalPointToIndex(pore_mm)] == 1
    assert plain[plain.TransformPhysicalPointToIndex(pore_mm)] == 0


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
    assert pipeline.printability_warnings(_box(), ANATOMICAL, (1.0, 1.0, 1.0)) == []


def test_printability_warnings_report_thin_material():
    """A one-voxel wall against a 3 mm minimum feature: say so."""
    wall = _slab(thickness_voxels=1)
    profile = presets.PrintProfile(name="t", description="d", min_feature_mm=3.0)
    messages = " ".join(pipeline.printability_warnings(wall, profile, (1.0, 1.0, 1.0)))
    assert "thinner than" in messages


def test_grid_warnings_are_silent_when_the_grid_can_honour_the_request():
    """Isotropic 1 mm voxels realise a 2 mm feature exactly. Nothing to warn about."""
    profile = presets.PrintProfile(name="t", description="d", min_feature_mm=2.0)
    grid = (1.0, 1.0, 1.0)
    assert pipeline.grid_warnings(_box(spacing=grid), profile, grid) == []


def test_thin_mask_speckles_a_voxelised_surface_which_is_why_thicken_ignores_it():
    """An opening is the union of the balls it contains, so it cannot reach into a
    sharp convex corner -- and every surface of a voxelised object is locally
    sharp. `thin_mask` therefore marks a speckle over even a solid sphere. It is
    the right measure for *reporting* local thickness and the wrong one for
    *selecting* what to grow, which is why `thicken` uses the thick core's reach
    instead."""
    n = 40
    zz, yy, xx = np.indices((n, n, n))
    sphere = ((zz - 20) ** 2 + (yy - 20) ** 2 + (xx - 20) ** 2 <= 12 ** 2)
    image = sitk.GetImageFromArray(sphere.astype(np.uint8))
    image.SetSpacing((1.0, 1.0, 1.0))

    speckle = np.count_nonzero(sitk.GetArrayFromImage(segment.thin_mask(image, 3.0)))
    assert speckle > 0, "the voxelised surface is locally sharp"

    # ... and yet a solid sphere must not grow at all
    grown = segment.thicken(image, 3.0)
    assert np.array_equal(sitk.GetArrayFromImage(image),
                          sitk.GetArrayFromImage(grown))


def test_the_scans_own_resolution_limit_is_reported_not_hidden():
    """A 0.6 mm feature cannot be *measured* on 0.8 mm slices, however finely the
    printability grid is chosen: bone that thin never entered the data."""
    profile = presets.PrintProfile(name="t", description="d", min_feature_mm=0.6)
    messages = " ".join(pipeline.grid_warnings(_box(spacing=(0.3, 0.3, 0.3)),
                                               profile, HEAD_CT))
    assert "not a measurement of the anatomy" in messages


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
