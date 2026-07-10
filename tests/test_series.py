"""Series discovery, slice ordering and selection.

The headline test is :func:`test_slices_are_ordered_by_position_not_filename`.
In the real head CT this tool was built against, file ``1`` held instance 222 and
file ``231`` held instance 226. Sorting by filename yields a scrambled volume that
still renders as a plausible-looking blob.
"""

import os

import numpy as np
import pydicom
import pytest
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, generate_uid

from dicom_surface import series as series_mod

SERIES_UID = generate_uid()


AXIAL = [1, 0, 0, 0, 1, 0]
CORONAL = [1, 0, 0, 0, 0, -1]
SAGITTAL = [0, 1, 0, 0, 0, -1]


def _write_slice(path, z, series_uid=SERIES_UID, *, modality="CT", localizer=False,
                 description="TEST SERIES", kernel="Hr68", series_number=6, rows=8,
                 orientation=None, pixel_spacing=(0.5, 0.5), patient_id=None):
    fm = FileMetaDataset()
    fm.MediaStorageSOPClassUID = CTImageStorage
    fm.MediaStorageSOPInstanceUID = generate_uid()
    fm.TransferSyntaxUID = ExplicitVRLittleEndian

    ds = Dataset()
    ds.file_meta = fm
    ds.SOPClassUID = CTImageStorage
    ds.SOPInstanceUID = fm.MediaStorageSOPInstanceUID
    ds.SeriesInstanceUID = series_uid
    ds.StudyInstanceUID = generate_uid()
    ds.Modality = modality
    if patient_id is not None:
        ds.PatientID = patient_id
    ds.SeriesDescription = description
    ds.SeriesNumber = series_number
    ds.ConvolutionKernel = kernel
    ds.ImageType = ["DERIVED", "SECONDARY", "LOCALIZER"] if localizer else ["ORIGINAL", "PRIMARY", "AXIAL"]

    ds.Rows = rows
    ds.Columns = rows
    ds.PixelSpacing = list(pixel_spacing)
    ds.SliceThickness = 1.0
    ds.ImageOrientationPatient = list(orientation or AXIAL)
    # Offset the slice along its own normal so ordering works in any plane.
    o = np.asarray(ds.ImageOrientationPatient, dtype=float)
    normal = np.cross(o[0:3], o[3:6])
    ds.ImagePositionPatient = [float(v) for v in (normal * float(z))]
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0
    ds.RescaleIntercept = -1024
    ds.RescaleSlope = 1
    ds.PixelData = np.full((rows, rows), int(1024 + z), dtype=np.uint16).tobytes()

    ds.save_as(path, enforce_file_format=True)


def test_slices_are_ordered_by_position_not_filename(tmp_path):
    d = str(tmp_path)
    z_values = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    # Deliberately adversarial: filename order is the reverse of position order.
    names = ["1", "2", "3", "4", "5", "6"]
    for name, z in zip(names, reversed(z_values)):
        _write_slice(os.path.join(d, name), z)

    found = series_mod.discover(d)
    assert len(found) == 1
    s = found[0]
    assert s.n_slices == 6

    got_z = []
    for f in s.files:
        ds = pydicom.dcmread(f, stop_before_pixels=True)
        got_z.append(float(ds.ImagePositionPatient[2]))
    assert got_z == sorted(z_values)

    # And filename order would have been wrong.
    assert [os.path.basename(f) for f in s.files] == list(reversed(names))


def test_slice_spacing_is_measured_from_positions(tmp_path):
    d = str(tmp_path)
    for i, z in enumerate([0.0, 0.8, 1.6, 2.4]):
        _write_slice(os.path.join(d, "s%d" % i), z)
    s = series_mod.discover(d)[0]
    assert s.slice_spacing == pytest.approx(0.8)
    assert s.spacing_uniform
    # SliceThickness (1.0) disagrees with the true spacing (0.8) -- overlapping
    # reconstruction. The measured value must win.
    assert s.slice_thickness == pytest.approx(1.0)


def test_irregular_spacing_is_flagged(tmp_path):
    d = str(tmp_path)
    for i, z in enumerate([0.0, 1.0, 2.0, 5.0]):
        _write_slice(os.path.join(d, "s%d" % i), z)
    s = series_mod.discover(d)[0]
    assert not s.spacing_uniform
    assert s.spacing_spread_mm == pytest.approx(2.0)


def test_series_are_grouped_and_localizers_excluded(tmp_path):
    d = str(tmp_path)
    main = generate_uid()
    scout = generate_uid()
    for i in range(6):
        _write_slice(os.path.join(d, "m%d" % i), i, series_uid=main, series_number=6)
    _write_slice(os.path.join(d, "scout0"), 0, series_uid=scout, series_number=1,
                 localizer=True, description="Topogramm")

    found = series_mod.discover(d)
    assert len(found) == 2
    by_uid = {s.uid: s for s in found}
    assert by_uid[main].usable
    assert not by_uid[scout].usable  # localizer AND too few slices

    assert series_mod.select(found, None).uid == main


def test_nested_directories_are_walked(tmp_path):
    d = str(tmp_path)
    for i in range(6):
        sub = os.path.join(d, "1", "2", "6", str(i))
        os.makedirs(os.path.dirname(sub), exist_ok=True)
        _write_slice(sub, i)
    s = series_mod.discover(d)
    assert len(s) == 1 and s[0].n_slices == 6


def test_sharp_kernel_detection(tmp_path):
    d = str(tmp_path)
    for i in range(6):
        _write_slice(os.path.join(d, "a%d" % i), i, kernel="Hr68")
    assert series_mod.discover(d)[0].sharp_kernel

    d2 = str(tmp_path / "soft")
    os.makedirs(d2)
    for i in range(6):
        _write_slice(os.path.join(d2, "b%d" % i), i, kernel="Hr40")
    assert not series_mod.discover(d2)[0].sharp_kernel


@pytest.mark.parametrize(
    "kernel,sharp",
    [("Hr68", True), ("Hr40", False), ("B70f", True), ("B30f", False),
     ("BONE", False), ("STANDARD", False), ("U90u", True), ("I70f", False)],
)
def test_sharp_kernel_patterns(kernel, sharp):
    assert bool(series_mod.SHARP_KERNEL_RE.search(kernel)) is sharp


def test_select_by_row_id_uid_and_description(tmp_path):
    d = str(tmp_path)
    for i in range(6):
        _write_slice(os.path.join(d, "x%d" % i), i, description="GS nativ 1,00 ax")
    found = series_mod.discover(d)
    uid = found[0].uid

    assert series_mod.select(found, uid).uid == uid
    assert series_mod.select(found, "1").uid == uid
    assert series_mod.select(found, "nativ").uid == uid
    with pytest.raises(ValueError, match="no series matches"):
        series_mod.select(found, "6")
    with pytest.raises(ValueError, match="no series matches"):
        series_mod.select(found, "does-not-exist")


def test_uid_holding_two_orientations_is_split(tmp_path):
    """Regression, from a real study: SeriesInstanceUID 1021 held 64 axial frames
    and one perpendicular frame. Nothing in DICOM forbids that."""
    d = str(tmp_path)
    uid = generate_uid()
    # The odd frame sorts FIRST by filename, so a naive walk picks its normal.
    _write_slice(os.path.join(d, "aaa_odd"), 0, series_uid=uid, series_number=1021,
                 orientation=SAGITTAL)
    for i in range(20):
        _write_slice(os.path.join(d, "zzz_%02d" % i), i * 2.0, series_uid=uid,
                     series_number=1021, orientation=AXIAL)

    found = series_mod.discover(d)
    assert len(found) == 2, [(s.id, s.uid, s.part) for s in found]
    assert {s.n_parts for s in found} == {2}

    big = [s for s in found if s.n_slices == 20][0]
    small = [s for s in found if s.n_slices == 1][0]

    assert big.id == 1
    assert small.id == 2
    assert big.part == 1     # largest stack is orientation part 1
    assert small.part == 2
    assert big.plane == "axial"
    assert big.slice_spacing == pytest.approx(2.0)
    assert big.usable
    assert not small.usable
    assert "only 1 slice" in small.unusable_reason


def test_wrong_normal_would_collapse_the_spacing(tmp_path):
    """The failure this prevents: projecting axial positions onto a sagittal
    normal collapses every depth to ~0, giving a sub-micron slice spacing that
    then wins the ranking and scrambles the slice order."""
    d = str(tmp_path)
    uid = generate_uid()
    _write_slice(os.path.join(d, "aaa_odd"), 0, series_uid=uid, series_number=1021,
                 orientation=SAGITTAL)
    for i in range(20):
        _write_slice(os.path.join(d, "zzz_%02d" % i), i * 2.0, series_uid=uid,
                     series_number=1021, orientation=AXIAL)
    # a genuine, finer axial acquisition elsewhere in the study
    real = generate_uid()
    for i in range(30):
        _write_slice(os.path.join(d, "real_%02d" % i), i * 0.5, series_uid=real,
                     series_number=2, orientation=AXIAL, pixel_spacing=(0.4, 0.4))

    found = series_mod.discover(d)
    for s in found:
        if s.usable:
            assert s.slice_spacing > series_mod.MIN_SLICE_SPACING_MM
            assert s.voxel_volume_mm3 > 1e-4

    # the real acquisition wins; the reformat does not sneak in on a bogus voxel size
    chosen = series_mod.select(found, None)
    assert chosen.uid == real
    assert chosen.n_slices == 30


def test_discovery_is_independent_of_filesystem_order(tmp_path):
    """The slice normal must not depend on which file os.walk yields first."""
    import os as _os

    orders = []
    for run, names in enumerate((["a", "b", "c", "d", "e", "f"],
                                 ["f", "e", "d", "c", "b", "a"])):
        d = str(tmp_path / ("run%d" % run))
        _os.makedirs(d)
        for name, z in zip(names, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]):
            _write_slice(_os.path.join(d, name), z)
        s = series_mod.discover(d)[0]
        orders.append([_os.path.basename(f) for f in s.files])
        assert s.normal == pytest.approx((0.0, 0.0, 1.0))
    # both runs order slices by position, giving mirrored filename sequences
    assert orders[0] == list(reversed(orders[1]))


def test_geometry_guards_reject_degenerate_stacks(tmp_path):
    d = str(tmp_path)
    for i in range(8):
        _write_slice(os.path.join(d, "s%d" % i), i * 1.0)
    s = series_mod.discover(d)[0]
    assert s.usable and s.unusable_reason is None

    s.slice_spacing = 3.9e-07  # only a mis-derived slice normal yields this
    assert not s.usable
    assert "implausible slice spacing" in s.unusable_reason

    s.slice_spacing = 1.0
    s.spacing_spread_mm = 73.77  # the real reformat's spread
    assert not s.usable
    assert "irregular spacing" in s.unusable_reason


def test_rank_sends_degenerate_voxel_volume_last(tmp_path):
    d = str(tmp_path)
    good = generate_uid()
    for i in range(10):
        _write_slice(os.path.join(d, "g%02d" % i), i * 1.0, series_uid=good, series_number=2)
    found = series_mod.discover(d)
    bad = found[0]
    ranked = series_mod.rank([bad])
    assert ranked  # sanity

    # a stack whose geometry did not survive: it must never outrank a real one
    broken = series_mod.Series(uid="x", modality="CT", description="broken",
                               series_number=9, normal=(0.0, 0.0, 1.0))
    broken.files = list(bad.files)
    broken.pixel_spacing = (1.0, 1.0)
    broken.slice_spacing = 3.9e-07
    assert series_mod.rank([broken, bad])[0] is bad


def test_select_each_orientation_split_by_row_id(tmp_path):
    d = str(tmp_path)
    uid = generate_uid()
    _write_slice(os.path.join(d, "aaa_odd"), 0, series_uid=uid, series_number=1021,
                 orientation=SAGITTAL)
    for i in range(20):
        _write_slice(os.path.join(d, "zzz_%02d" % i), i * 2.0, series_uid=uid,
                     series_number=1021, orientation=AXIAL)
    found = series_mod.discover(d)

    assert series_mod.select(found, "1").n_slices == 20
    assert series_mod.select(found, "2").n_slices == 1
    with pytest.raises(ValueError, match="no series matches"):
        series_mod.select(found, "1021.1")
    with pytest.raises(ValueError, match=r"ambiguous \(2 rows: 1, 2\)"):
        series_mod.select(found, uid)


def test_ambiguous_uid_lists_the_row_ids(tmp_path):
    d = str(tmp_path)
    uid = generate_uid()
    for i in range(8):
        _write_slice(os.path.join(d, "ax%02d" % i), i * 2.0, series_uid=uid,
                     series_number=7, orientation=AXIAL)
    for i in range(8):
        _write_slice(os.path.join(d, "sg%02d" % i), i * 2.0, series_uid=uid,
                     series_number=7, orientation=SAGITTAL)
    found = series_mod.discover(d)
    assert len(found) == 2 and all(s.usable for s in found)
    with pytest.raises(ValueError, match=r"2 rows: 1, 2"):
        series_mod.select(found, uid)


def test_duplicate_dicom_series_numbers_get_unique_row_ids(tmp_path):
    d = str(tmp_path)
    left, right = generate_uid(), generate_uid()
    for i in range(6):
        _write_slice(os.path.join(d, "left%d" % i), i, series_uid=left,
                     description="shared left", series_number=7)
        _write_slice(os.path.join(d, "right%d" % i), i, series_uid=right,
                     description="shared right", series_number=7)

    found = series_mod.discover(d)

    assert [series.id for series in found] == [1, 2]
    assert [series.series_number for series in found] == [7, 7]
    assert {series_mod.select(found, "1").uid, series_mod.select(found, "2").uid} == {
        left,
        right,
    }
    with pytest.raises(ValueError, match="no series matches"):
        series_mod.select(found, "7")
    with pytest.raises(ValueError, match=r"description 'shared' is ambiguous \(2 rows: 1, 2\)"):
        series_mod.select(found, "shared")
    assert series_mod.select(found, left).uid == left


def test_orientation_jitter_does_not_split_a_series(tmp_path):
    """Oblique reformats carry float noise in ImageOrientationPatient; a series
    must not shatter into one group per slice."""
    d = str(tmp_path)
    uid = generate_uid()
    for i in range(10):
        eps = 1e-4 * (i - 5)
        _write_slice(os.path.join(d, "s%02d" % i), i * 1.0, series_uid=uid,
                     orientation=[1, eps, 0, -eps, 1, 0])
    found = series_mod.discover(d)
    assert len(found) == 1
    assert found[0].n_slices == 10
    assert found[0].n_parts == 1


def test_select_rejects_ambiguous_description(tmp_path):
    d = str(tmp_path)
    bone, soft = generate_uid(), generate_uid()
    for i in range(6):
        _write_slice(os.path.join(d, "a%d" % i), i, series_uid=bone,
                     description="axial bone", series_number=1)
    # second series, same description substring
    for i in range(6):
        _write_slice(os.path.join(d, "b%d" % i), i, series_uid=soft,
                     description="axial soft", series_number=2)
    found = series_mod.discover(d)
    assert len(found) == 2
    with pytest.raises(ValueError, match="is ambiguous"):
        series_mod.select(found, "axial")


def test_rank_prefers_smaller_voxels(tmp_path):
    d = str(tmp_path)
    thin = generate_uid()
    thick = generate_uid()
    for i in range(20):
        _write_slice(os.path.join(d, "t%d" % i), i * 0.5, series_uid=thin, series_number=2)
    for i in range(8):
        _write_slice(os.path.join(d, "k%d" % i), i * 3.0, series_uid=thick, series_number=3)
    found = series_mod.discover(d)
    assert series_mod.rank(found)[0].uid == thin


@pytest.mark.parametrize(
    "orientation,expected",
    [(AXIAL, "axial"), (CORONAL, "coronal"), (SAGITTAL, "sagittal")],
)
def test_plane_detection(tmp_path, orientation, expected):
    d = str(tmp_path / expected)
    os.makedirs(d)
    for i in range(6):
        _write_slice(os.path.join(d, "s%d" % i), i, orientation=orientation)
    assert series_mod.discover(d)[0].plane == expected


def test_oblique_plane_detected(tmp_path):
    d = str(tmp_path)
    # 30 degrees off axial
    c, s = 0.866, 0.5
    for i in range(6):
        _write_slice(os.path.join(d, "o%d" % i), i, orientation=[1, 0, 0, 0, c, s])
    assert series_mod.discover(d)[0].plane == "oblique"


def test_significant_figure_bucketing_of_voxel_volume():
    """Real numbers from the study this tool was built on: the sagittal reformat
    measured 0.7997661776617022 mm spacing against the axial's 0.8. Comparing raw
    floats let 0.03% of rounding noise decide the default series."""
    sig = series_mod._significant
    axial = 0.315 * 0.315 * 0.8
    sagittal = 0.315 * 0.315 * 0.7997661776617022
    assert axial != sagittal
    assert sig(axial) == sig(sagittal)
    # but a genuinely coarser series stays distinct
    coronal = 0.3374322917 * 0.3374317343 * 0.8
    assert sig(coronal) != sig(axial)
    assert sig(0.0) == 0.0


def test_rank_ignores_subpercent_spacing_noise(tmp_path):
    d = str(tmp_path)
    axial, sagittal = generate_uid(), generate_uid()
    for i in range(10):
        _write_slice(os.path.join(d, "ax%02d" % i), i * 0.8, series_uid=axial,
                     series_number=6, orientation=AXIAL, pixel_spacing=(0.315, 0.315))
    for i in range(12):  # more slices AND a hair smaller voxel, from float noise
        _write_slice(os.path.join(d, "sg%02d" % i), i * 0.7997661776617022,
                     series_uid=sagittal, series_number=8, orientation=SAGITTAL,
                     pixel_spacing=(0.315, 0.315))
    found = series_mod.discover(d)
    assert series_mod.rank(found)[0].uid == axial


def test_rank_prefers_axial_over_a_reformat_with_more_slices(tmp_path):
    """Regression: a sagittal reformat of the same acquisition can carry more
    slices while holding no more information. Slice count is the wrong tiebreak;
    the acquired (axial) plane wins."""
    d = str(tmp_path)
    axial, sagittal = generate_uid(), generate_uid()
    for i in range(231):
        _write_slice(os.path.join(d, "ax%03d" % i), i * 0.8, series_uid=axial,
                     series_number=6, orientation=AXIAL, rows=4)
    for i in range(242):  # MORE slices, same voxel size
        _write_slice(os.path.join(d, "sg%03d" % i), i * 0.8, series_uid=sagittal,
                     series_number=8, orientation=SAGITTAL, rows=4)

    found = series_mod.discover(d)
    best = series_mod.rank(found)[0]
    assert best.uid == axial
    assert best.plane == "axial"
    assert series_mod.select(found, None).uid == axial


def test_no_usable_series_raises(tmp_path):
    d = str(tmp_path)
    _write_slice(os.path.join(d, "only"), 0, localizer=True)
    with pytest.raises(ValueError, match="no usable image series"):
        series_mod.select(series_mod.discover(d), None)


def test_non_dicom_files_are_ignored(tmp_path):
    d = str(tmp_path)
    with open(os.path.join(d, "README.txt"), "w") as fh:
        fh.write("not dicom")
    for i in range(6):
        _write_slice(os.path.join(d, "s%d" % i), i)
    assert len(series_mod.discover(d)) == 1
