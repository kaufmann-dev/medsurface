"""Registration and fusion, including the refusals.

Registration cannot fail on its own -- FFT always has a peak, ICP always
converges somewhere. The tests that matter here are the ones asserting that a
confident-looking transform over unrelated anatomy is *rejected*.
"""

import math
from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk

from medsurface import merge as merge_mod
from medsurface import presets, registration, volume
from medsurface.catalog import DicomSource, FileSource, VolumeCandidate
from medsurface.merge import MergeError
from medsurface.series import Series


def _image(arr, spacing=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    img = sitk.GetImageFromArray(arr.astype(np.uint8))
    img.SetSpacing(spacing)
    img.SetOrigin(origin)
    return img


def _ball(shape, center, radius):
    zz, yy, xx = np.indices(shape)
    d2 = (zz - center[0]) ** 2 + (yy - center[1]) ** 2 + (xx - center[2]) ** 2
    return (d2 <= radius ** 2).astype(np.uint8)


def _ellipsoid(shape, center, radii):
    zz, yy, xx = np.indices(shape)
    d = (((zz - center[0]) / radii[0]) ** 2
         + ((yy - center[1]) / radii[1]) ** 2
         + ((xx - center[2]) / radii[2]) ** 2)
    return (d <= 1.0).astype(np.uint8)


def _lumpy_shell(shape=(70, 70, 70)):
    """A hollow, chiral blob whose pose is unambiguous.

    Every symmetry must be broken, or the tests below measure nothing. A sphere
    with a bump on the z axis is rotationally symmetric about z: ICP settles 10
    degrees off with a 0.18 mm residual and 100% surface overlap, and it is
    right to -- the pose is genuinely ambiguous. Unequal semi-axes plus two
    off-axis lobes fix the pose.
    """
    c = (35, 35, 35)
    body = _ellipsoid(shape, c, (26, 20, 16)) & ~_ellipsoid(shape, c, (21, 15, 11))
    body |= _ball(shape, (20, 24, 46), 8)   # off every axis
    body |= _ball(shape, (48, 44, 26), 6)   # and a different one, elsewhere
    body |= _ellipsoid(shape, (35, 50, 20), (4, 9, 4))  # an elongated wing
    return body.astype(np.uint8)


def _rigid(rotation_deg, axis, translation):
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    theta = math.radians(rotation_deg)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rot = np.eye(3) + math.sin(theta) * k + (1 - math.cos(theta)) * (k @ k)
    t = np.eye(4)
    t[:3, :3] = rot
    t[:3, 3] = translation
    return t


def _apply_to_image(mask, transform):
    """Resample a mask so its anatomy moves by `transform` in world coordinates.

    Occupancy is resampled linearly and re-thresholded, not nearest-neighbour.
    Nearest-neighbour quantises the ground truth to the voxel grid, which put a
    floor of ~1 degree on how accurately any registration could be *seen* to
    recover a small rotation -- the fixture, not the algorithm.
    """
    inv = registration.inverse_transform(transform)
    r = sitk.ResampleImageFilter()
    r.SetReferenceImage(mask)
    size = np.asarray(mask.GetSize()) + 40
    r.SetSize([int(v) for v in size])
    r.SetOutputOrigin([o - 20.0 * s for o, s in zip(mask.GetOrigin(), mask.GetSpacing())])
    r.SetInterpolator(sitk.sitkLinear)
    r.SetDefaultPixelValue(0.0)
    r.SetTransform(inv)
    resampled = r.Execute(sitk.Cast(mask, sitk.sitkFloat32))
    return sitk.BinaryThreshold(resampled, 0.5, 1e9, 1, 0)


# ------------------------------------------------------------- registration
def test_recovers_a_known_rigid_transform():
    fixed = _image(_lumpy_shell())
    truth = _rigid(12.0, (0.2, 0.3, 1.0), (7.0, -5.0, 11.0))
    moving = _apply_to_image(fixed, truth)

    result = registration.rigid_register(fixed, moving, samples=20000)

    # rigid_register maps moving -> fixed, i.e. it should invert `truth`
    composed = result.transform @ truth
    assert composed[:3, 3] == pytest.approx([0, 0, 0], abs=1.0)
    angle = math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(composed[:3, :3]) - 1) / 2))))
    assert angle < 2.0, angle

    assert result.inlier_rms_mm < 0.6
    assert result.surface_overlap > 0.9
    assert result.shared_fov_dice > 0.85


@pytest.mark.parametrize("degrees", [3.0, 8.0, 25.0, 40.0])
def test_recovers_rotations_across_the_capture_range(degrees):
    """Both ends of the range are load-bearing.

    Small rotations fail if the trim is fixed: the correspondences a hard trim
    discards are the ones furthest from the rotation axis, which is where the
    rotational signal lives. At a fixed 0.45 trim, a 5 degree misalignment of two
    identical volumes converges 4.49 degrees off.

    Large rotations fail if the fit is bootstrapped with point-to-point ICP,
    which at 40 degrees settles into a wrong basin 33 degrees away.
    Point-to-plane from the FFT translation handles the whole range.
    """
    fixed = _image(_lumpy_shell())
    truth = _rigid(degrees, (0.3, 0.4, 1.0), (6.0, -4.0, 8.0))
    moving = _apply_to_image(fixed, truth)

    result = registration.rigid_register(fixed, moving, samples=20000)

    composed = result.transform @ truth
    angle = math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(composed[:3, :3]) - 1) / 2))))
    assert angle < 1.0, "rotation error %.2f deg at %.0f deg" % (angle, degrees)
    assert np.linalg.norm(composed[:3, 3]) < 1.0


def test_identical_input_registers_to_the_identity():
    mask = _image(_lumpy_shell())
    result = registration.rigid_register(mask, mask, samples=20000)
    assert result.rotation_deg < 0.5
    assert np.allclose(result.translation_mm, 0.0, atol=0.5)
    assert result.shared_fov_dice > 0.98
    assert result.surface_overlap > 0.98


def test_empty_mask_is_reported_not_silently_registered():
    empty = _image(np.zeros((30, 30, 30), dtype=np.uint8))
    solid = _image(_ball((30, 30, 30), (15, 15, 15), 8))
    with pytest.raises(registration.RegistrationError, match="empty"):
        registration.rigid_register(solid, empty)


def test_merge_translates_registration_failure_to_its_public_error(monkeypatch):
    image = _image(_ball((30, 30, 30), (15, 15, 15), 8))

    def load(candidate, *, allow_large_volume):
        assert not allow_large_volume
        return volume.Volume(image=image, candidate=candidate)

    def fail_registration(*_args, **_kwargs):
        raise registration.RegistrationError("one mask vanished on the registration grid")

    monkeypatch.setattr(merge_mod.volume_mod, "load", load)
    monkeypatch.setattr(merge_mod.pipeline, "build_mask", lambda *_args, **_kwargs: image)
    monkeypatch.setattr(merge_mod.registration, "rigid_register", fail_registration)

    with pytest.raises(MergeError, match="one mask vanished"):
        merge_mod.merge(
            _candidate(uid="a"),
            _candidate(uid="b"),
            presets.get("bone"),
            "unused.stl",
            fixed_threshold=1.0,
            moving_threshold=1.0,
        )


# --------------------------------------------------------------- the refusals
def _result(**kw):
    base = dict(
        transform=np.eye(4), fft_translation_mm=np.zeros(3), rotation_deg=0.0,
        translation_mm=np.zeros(3), inlier_rms_mm=0.2, inlier_median_mm=0.1,
        overlap_moving_in_fixed=0.95, overlap_fixed_in_moving=0.92,
        shared_fov_dice=0.85, shared_fov_mm3=1e6,
    )
    base.update(kw)
    return registration.RegistrationResult(**base)


def test_good_registration_passes_the_gates():
    merge_mod.check_registration(_result())  # must not raise


def test_surface_overlap_is_the_weaker_direction():
    """A small object buried in a large one matches the large one almost
    everywhere in one direction and hardly at all in the other. The impostor
    scored 0.964 forward and 0.270 backward."""
    r = _result(overlap_moving_in_fixed=0.964, overlap_fixed_in_moving=0.270)
    assert r.surface_overlap == pytest.approx(0.270)
    with pytest.raises(MergeError, match="fixed->moving"):
        merge_mod.check_registration(r)


@pytest.mark.parametrize(
    "kw,expected",
    [
        ({"overlap_fixed_in_moving": 0.02}, "surfaces agree over only"),
        ({"overlap_moving_in_fixed": 0.02}, "surfaces agree over only"),
        ({"shared_fov_dice": 0.05}, "bone agreement"),
        ({"shared_fov_mm3": 100.0}, "share only"),
    ],
)
def test_bad_registration_is_refused_with_a_reason(kw, expected):
    with pytest.raises(MergeError, match=expected):
        merge_mod.check_registration(_result(**kw))


def test_residual_is_reported_but_not_gated():
    """The residual does not discriminate: an impostor scores 0.43 mm, a
    5%-oversized skull 0.41 mm, true pairs 0.09-0.21 mm. ICP drives some residual
    down whatever it is fitting, so it informs but never decides."""
    merge_mod.check_registration(_result(inlier_rms_mm=9.9))  # must not raise
    assert "rms" in _result().summary()


def test_geometry_alone_cannot_reject_a_similar_body():
    """Measured on a real skull scaled by 3%, within person-to-person variation:
    overlap 0.989, dice 0.675 -- both gates pass. Subject identity therefore
    remains an explicit user responsibility."""
    similar_body = _result(overlap_moving_in_fixed=0.989, overlap_fixed_in_moving=0.989,
                           shared_fov_dice=0.675)
    merge_mod.check_registration(similar_body)  # geometry waves it through


def test_gates_have_margin_against_real_measurements():
    """Thresholds must sit between the worst true pair and the best impostor,
    not hug either. Numbers measured on real studies; see merge.py."""
    worst_true_overlap, worst_true_dice = 0.922, 0.835
    impostor_overlap, impostor_dice = 0.270, 0.323

    assert impostor_overlap < merge_mod.MIN_SURFACE_OVERLAP < worst_true_overlap
    assert impostor_dice < merge_mod.MIN_SHARED_FOV_DICE < worst_true_dice
    # and not by a hair
    assert merge_mod.MIN_SURFACE_OVERLAP - impostor_overlap > 0.2
    assert worst_true_overlap - merge_mod.MIN_SURFACE_OVERLAP > 0.2
    assert merge_mod.MIN_SHARED_FOV_DICE - impostor_dice > 0.15
    assert worst_true_dice - merge_mod.MIN_SHARED_FOV_DICE > 0.15


def test_surface_stage_comes_from_the_preset():
    """A fused mesh must be finished exactly as a single-scan one is.

    Slice terracing is baked into each scan's mask by its own slice pitch, so an
    isotropic fused grid does not remove it -- it merely samples it more finely.
    Smoothing the fused surface any more lightly than `convert` does leaves it
    visibly rougher than the scans it was built from.
    """
    import inspect

    from medsurface import presets

    signature = inspect.signature(merge_mod.merge)
    for name in ("smooth_iters", "smooth_force", "simplify_error_mm", "post_smooth_iters"):
        assert signature.parameters[name].default is None, (
            "%s must default to the preset, not to a merge-specific constant" % name
        )

    bone = presets.get("bone")
    assert (bone.smooth_iters, bone.smooth_force, bone.post_smooth_iters) == (20, 0.1, 40)


@pytest.mark.parametrize("grid_mm", [0.0, -0.4, float("nan"), float("inf")])
def test_merge_rejects_invalid_grid_before_loading(grid_mm, monkeypatch):
    monkeypatch.setattr(
        merge_mod.volume_mod,
        "load",
        lambda _candidate: pytest.fail("volume loading must not start for an invalid grid"),
    )

    with pytest.raises(ValueError, match="grid_mm must be finite and greater than zero"):
        merge_mod.merge(
            _candidate(uid="a"),
            _candidate(uid="b"),
            presets.get("bone"),
            "unused.stl",
            grid_mm=grid_mm,
        )


def test_common_grid_rejects_tiny_spacing_without_integer_overflow():
    mask = _image(_lumpy_shell(), spacing=(1.0, 1.0, 1.0))

    with pytest.raises(MergeError, match="raise --grid-mm"):
        merge_mod._common_grid(mask, mask, np.eye(4), 1e-12)


def test_common_grid_accepts_exactly_the_default_voxel_limit():
    mask = sitk.Image((2, 2, 2), sitk.sitkUInt8)
    mask.SetSpacing((495.0, 995.0, 995.0))

    size, _origin = merge_mod._common_grid(mask, mask, np.eye(4), 1.0)

    assert size == (500, 1_000, 1_000)
    assert math.prod(size) == 500_000_000


def test_common_grid_override_allows_a_grid_above_the_default_limit():
    mask = _image(_lumpy_shell(), spacing=(1.0, 1.0, 1.0))

    size, _origin = merge_mod._common_grid(
        mask,
        mask,
        np.eye(4),
        1e-12,
        allow_large_volume=True,
    )

    assert math.prod(size) > 500_000_000


def test_force_overrides_the_gates():
    merge_mod.check_registration(
        _result(overlap_moving_in_fixed=0.0, overlap_fixed_in_moving=0.0,
                shared_fov_dice=0.0), force=True)


def test_unrelated_anatomy_fails_the_gates():
    """The headline safety property: two scans that share no anatomy must not
    silently fuse. Registration still returns a transform -- it always does."""
    shell = _image(_lumpy_shell())
    # a solid bar, far away, nothing like the shell
    bar = np.zeros((70, 70, 70), dtype=np.uint8)
    bar[10:60, 30:38, 30:38] = 1
    other = _image(bar, origin=(400.0, -250.0, 900.0))

    result = registration.rigid_register(shell, other, samples=20000)
    with pytest.raises(MergeError):
        merge_mod.check_registration(result)


# ------------------------------------------------------------- compatibility
def _candidate(uid="1.2.3", modality="CT", part=1, *, path=None):
    series = Series(uid=uid, modality=modality, description="d", series_number=1)
    series.files = ["slice"]
    series.part = part
    source = FileSource(Path(path), "NIfTI") if path else DicomSource(Path("scans"), series)
    return VolumeCandidate(
        id=1,
        source=source,
        format="NIfTI" if path else "DICOM",
        source_name=path or "scans",
        modality=None if path else modality,
        description="d",
        size=(10, 10, 10),
        spacing=(1.0, 1.0, 1.0),
        direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        origin=(0.0, 0.0, 0.0),
        pixel_type="16-bit signed integer",
        components=1,
        plane="axial",
    )


def test_merge_announces_volume_loading_before_it_starts(monkeypatch):
    messages = []
    emitted_warnings = []

    class StopLoading(Exception):
        pass

    def stop(_candidate, *, allow_large_volume):
        assert messages[-1] == "load fixed volume ..."
        assert allow_large_volume
        assert emitted_warnings[0] == (
            "subject identity is not verified; confirm that fixed and moving volumes "
            "show the same subject before using the fused surface"
        )
        assert any("--fixed-threshold" in message for message in emitted_warnings)
        assert any("--moving-threshold" in message for message in emitted_warnings)
        raise StopLoading

    monkeypatch.setattr(merge_mod.volume_mod, "load", stop)
    with pytest.raises(StopLoading):
        merge_mod.merge(
            _candidate(uid="a"),
            _candidate(uid="b"),
            presets.get("bone"),
            "unused.stl",
            allow_large_volume=True,
            log=messages.append,
            warn=emitted_warnings.append,
        )


def test_merge_emits_fixed_volume_warnings_before_loading_moving(monkeypatch):
    fixed = _candidate(uid="a")
    moving = _candidate(uid="b")
    emitted_warnings = []

    class StopMovingLoad(Exception):
        pass

    fixed_loaded = object()

    def load(candidate, *, allow_large_volume):
        assert not allow_large_volume
        if candidate is fixed:
            return fixed_loaded
        assert "fixed geometry warning" in emitted_warnings
        raise StopMovingLoad

    monkeypatch.setattr(merge_mod.volume_mod, "load", load)
    monkeypatch.setattr(
        merge_mod.volume_mod,
        "warnings_for",
        lambda volume: ["fixed geometry warning"] if volume is fixed_loaded else [],
    )

    with pytest.raises(StopMovingLoad):
        merge_mod.merge(
            fixed,
            moving,
            presets.get("bone"),
            "unused.stl",
            warn=emitted_warnings.append,
        )


def test_same_volume_twice_is_refused():
    candidate = _candidate()
    with pytest.raises(MergeError, match="same volume"):
        merge_mod.check_compatible(candidate, candidate)


def test_modality_differences_do_not_change_merge_compatibility():
    fixed = _candidate(uid="a", modality="CT")
    moving = _candidate(uid="b", modality="MR")
    warnings = merge_mod.check_compatible(fixed, moving)
    assert any("identity is not verified" in warning for warning in warnings)


def test_same_file_via_hard_link_is_refused(tmp_path):
    source = tmp_path / "source.nii"
    alias = tmp_path / "alias.nii"
    source.write_bytes(b"data")
    alias.hardlink_to(source)
    fixed = _candidate(path=str(source))
    moving = _candidate(path=str(alias))
    with pytest.raises(MergeError, match="same volume"):
        merge_mod.check_compatible(fixed, moving)
