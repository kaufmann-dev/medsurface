"""Registration and fusion, including the refusals.

Registration cannot fail on its own -- FFT always has a peak, ICP always
converges somewhere. The tests that matter here are the ones asserting that a
confident-looking transform over unrelated anatomy is *rejected*.
"""

import math
import os

import numpy as np
import pytest
import SimpleITK as sitk

from dicom_surface import merge as merge_mod
from dicom_surface import registration
from dicom_surface.merge import MergeError
from dicom_surface.series import Series


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
    overlap 0.989, dice 0.675 -- both gates pass. This is why the demographic
    check exists, and why removing it would not be a simplification."""
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

    from dicom_surface import presets

    signature = inspect.signature(merge_mod.merge)
    for name in ("smooth_iters", "passband", "target_faces", "post_smooth_iters"):
        assert signature.parameters[name].default is None, (
            "%s must default to the preset, not to a merge-specific constant" % name
        )

    bone = presets.get("bone")
    assert (bone.smooth_iters, bone.passband, bone.post_smooth_iters) == (20, 0.1, 25)


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
def _series(uid="1.2.3", modality="CT", files=(), part=1):
    s = Series(uid=uid, modality=modality, description="d", series_number=1)
    s.files = list(files)
    s.part = part
    return s


def test_same_series_twice_is_refused():
    s = _series()
    with pytest.raises(MergeError, match="same series"):
        merge_mod.check_compatible(s, s)


def test_modality_mismatch_is_refused_but_forceable():
    a, b = _series(uid="a", modality="CT"), _series(uid="b", modality="MR")
    with pytest.raises(MergeError, match="modalities differ"):
        merge_mod.check_compatible(a, b)
    warnings = merge_mod.check_compatible(a, b, force=True)
    assert any("modalities differ" in w for w in warnings)


def _dicom_with_patient(tmp_path, name, uid, **patient):
    """A minimal 6-slice CT series carrying the given patient attributes."""
    from pydicom.dataset import Dataset, FileMetaDataset
    from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, generate_uid

    d = str(tmp_path / name)
    os.makedirs(d, exist_ok=True)
    paths = []
    for i in range(6):
        fm = FileMetaDataset()
        fm.MediaStorageSOPClassUID = CTImageStorage
        fm.MediaStorageSOPInstanceUID = generate_uid()
        fm.TransferSyntaxUID = ExplicitVRLittleEndian
        ds = Dataset()
        ds.file_meta = fm
        ds.SOPClassUID = CTImageStorage
        ds.SOPInstanceUID = fm.MediaStorageSOPInstanceUID
        ds.SeriesInstanceUID = uid
        ds.StudyInstanceUID = uid
        ds.Modality = "CT"
        for key, value in patient.items():
            setattr(ds, key, value)
        ds.Rows = ds.Columns = 4
        ds.PixelSpacing = [0.5, 0.5]
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        ds.ImagePositionPatient = [0.0, 0.0, float(i)]
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.BitsAllocated = ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 0
        ds.PixelData = np.zeros((4, 4), dtype=np.uint16).tobytes()
        p = os.path.join(d, "s%d" % i)
        ds.save_as(p, enforce_file_format=True)
        paths.append(p)
    return _series(uid=uid, files=paths)


def test_same_person_across_institutions_is_accepted(tmp_path):
    """The real case that broke a naive check: one skull, two hospitals. Medical
    record numbers are institution-scoped, so PatientID differs; the name is
    formatted differently; birth date and sex agree exactly."""
    a = _dicom_with_patient(tmp_path, "a", "1.2.3", PatientID="0012345",
                            PatientName="Doe^Jane", PatientBirthDate="19800101",
                            PatientSex="F")
    b = _dicom_with_patient(tmp_path, "b", "1.2.4", PatientID="X-99887766-A",
                            PatientName="DOE  JANE", PatientBirthDate="19800101",
                            PatientSex="F", IssuerOfPatientID="OTHER-HOSPITAL")

    match = merge_mod.compare_patients(a.files[0], b.files[0])
    assert match.verdict == "same", match
    assert "name" in match.matched and "birth_date" in match.matched
    assert merge_mod.check_compatible(a, b) == []


def test_conflicting_birth_date_is_refused(tmp_path):
    a = _dicom_with_patient(tmp_path, "a", "1.2.3", PatientID="SAME",
                            PatientName="Doe^Jane", PatientBirthDate="19800101")
    b = _dicom_with_patient(tmp_path, "b", "1.2.4", PatientID="SAME",
                            PatientName="Doe^Jane", PatientBirthDate="19731224")
    assert merge_mod.compare_patients(a.files[0], b.files[0]).verdict == "different"
    with pytest.raises(MergeError, match="birth_date"):
        merge_mod.check_compatible(a, b)


def test_conflicting_name_is_refused(tmp_path):
    a = _dicom_with_patient(tmp_path, "a", "1.2.3", PatientName="Doe^Jane")
    b = _dicom_with_patient(tmp_path, "b", "1.2.4", PatientName="Roe^Richard")
    with pytest.raises(MergeError, match="name"):
        merge_mod.check_compatible(a, b)


def test_differing_ids_with_no_corroboration_are_refused(tmp_path):
    a = _dicom_with_patient(tmp_path, "a", "1.2.3", PatientID="PAT-AAA")
    b = _dicom_with_patient(tmp_path, "b", "1.2.4", PatientID="PAT-BBB")
    assert merge_mod.compare_patients(a.files[0], b.files[0]).verdict == "different"
    with pytest.raises(MergeError, match="different patients"):
        merge_mod.check_compatible(a, b)
    assert any("different patients" in w
               for w in merge_mod.check_compatible(a, b, force=True))


def test_matching_id_is_enough(tmp_path):
    a = _dicom_with_patient(tmp_path, "a", "1.2.3", PatientID="PAT-SAME")
    b = _dicom_with_patient(tmp_path, "b", "1.2.4", PatientID="PAT-SAME")
    assert merge_mod.compare_patients(a.files[0], b.files[0]).verdict == "same"
    assert merge_mod.check_compatible(a, b) == []


def test_matching_id_from_different_issuers_is_not_enough(tmp_path):
    """Two hospitals can both issue MRN '12345' to different people."""
    a = _dicom_with_patient(tmp_path, "a", "1.2.3", PatientID="12345",
                            IssuerOfPatientID="HOSPITAL-A")
    b = _dicom_with_patient(tmp_path, "b", "1.2.4", PatientID="12345",
                            IssuerOfPatientID="HOSPITAL-B")
    assert merge_mod.compare_patients(a.files[0], b.files[0]).verdict == "unknown"
    assert any("cannot verify" in w for w in merge_mod.check_compatible(a, b))


def test_patient_values_never_appear_in_output(tmp_path):
    """Identifiers are compared, never surfaced -- not in the error, not in the
    fingerprint."""
    a = _dicom_with_patient(tmp_path, "a", "1.2.3", PatientID="PAT-SECRET",
                            PatientName="Secretname^Alice")
    b = _dicom_with_patient(tmp_path, "b", "1.2.4", PatientID="PAT-OTHER",
                            PatientName="Othername^Bob")

    match = merge_mod.compare_patients(a.files[0], b.files[0])
    assert match.verdict == "different"
    assert match.conflicts == ["name"]          # field names only, never values

    with pytest.raises(MergeError) as exc:
        merge_mod.check_compatible(a, b)
    text = str(exc.value).upper()
    for secret in ("SECRET", "OTHER", "ALICE", "BOB"):
        assert secret not in text


def test_deidentified_scans_warn_rather_than_refuse(tmp_path):
    a = _dicom_with_patient(tmp_path, "a", "1.2.3")
    b = _dicom_with_patient(tmp_path, "b", "1.2.4")
    warnings = merge_mod.check_compatible(a, b)
    assert any("cannot verify" in w for w in warnings)
