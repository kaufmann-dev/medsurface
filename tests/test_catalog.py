from __future__ import annotations

from pathlib import Path

import numpy as np
import pydicom
import pytest
import SimpleITK as sitk
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

from medsurface import catalog, volume
from medsurface.defaults import MAX_VOXELS


def _image() -> sitk.Image:
    image = sitk.GetImageFromArray(np.arange(4 * 5 * 6, dtype=np.int16).reshape(4, 5, 6))
    image.SetSpacing((0.7, 0.8, 1.2))
    image.SetOrigin((10.0, -4.0, 2.5))
    image.SetMetaData("modality", "MR")
    image.SetMetaData("description", "synthetic volume")
    return image


def _write_dicom_series(
    root: Path,
    *,
    uid: str | None = None,
    modality: str | None = "CT",
    series_number: int = 1,
    slices: int = 6,
    slice_spacing: float = 1.0,
    pixel_spacing: tuple[float, float] = (0.5, 0.5),
) -> str:
    uid = uid or generate_uid()
    root.mkdir(parents=True, exist_ok=True)
    for index in range(slices):
        path = root / ("slice-%03d.dcm" % index)
        meta = FileMetaDataset()
        meta.TransferSyntaxUID = ExplicitVRLittleEndian
        meta.MediaStorageSOPClassUID = generate_uid()
        meta.MediaStorageSOPInstanceUID = generate_uid()
        ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
        ds.SOPClassUID = meta.MediaStorageSOPClassUID
        ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
        ds.SeriesInstanceUID = uid
        if modality is not None:
            ds.Modality = modality
        ds.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
        ds.SeriesDescription = "axial source"
        ds.SeriesNumber = series_number
        ds.InstanceNumber = index + 1
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        ds.ImagePositionPatient = [0, 0, index * slice_spacing]
        ds.PixelSpacing = list(pixel_spacing)
        ds.SliceThickness = slice_spacing
        ds.Rows = 4
        ds.Columns = 5
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 1
        ds.RescaleSlope = 1
        ds.RescaleIntercept = -1024
        ds.PixelData = np.full((4, 5), index, dtype=np.int16).tobytes()
        pydicom.dcmwrite(path, ds, enforce_file_format=True)
    return uid


@pytest.mark.parametrize(
    "suffix,format_name",
    [
        (".nii", "NIfTI"),
        (".nii.gz", "NIfTI"),
        (".nrrd", "NRRD"),
        (".nhdr", "NRRD"),
        (".mha", "MetaImage"),
        (".mhd", "MetaImage"),
    ],
)
def test_direct_file_is_one_physical_volume(tmp_path, suffix, format_name):
    path = tmp_path / ("scan" + suffix)
    sitk.WriteImage(_image(), str(path), True)

    found = catalog.discover(path)

    assert len(found) == 1
    candidate = found[0]
    assert candidate.id == 1
    assert candidate.format == format_name
    assert candidate.source_name == path.name
    assert candidate.size == (6, 5, 4)
    assert candidate.spacing == pytest.approx((0.7, 0.8, 1.2))
    assert candidate.slices == 4
    assert candidate.plane == "axial"
    assert candidate.usable
    assert catalog.select(found, None) is candidate
    assert volume.load(candidate).image.GetSize() == (6, 5, 4)


@pytest.mark.parametrize(
    "size",
    [
        (1_000, 1_000, 499),
        (1_000, 1_000, 500),
    ],
)
def test_native_volume_limit_accepts_counts_at_or_below_the_boundary(
    tmp_path, monkeypatch, size
):
    path = tmp_path / "scan.mha"
    sitk.WriteImage(_image(), str(path))
    candidate = catalog.discover(path)[0]
    candidate.size = size
    loaded = _image()
    monkeypatch.setattr(sitk, "ReadImage", lambda _path: loaded)

    assert volume.load(candidate).image is loaded
    assert np.prod(size, dtype=object) <= MAX_VOXELS


def test_native_volume_limit_refuses_before_reading_pixels(tmp_path, monkeypatch):
    path = tmp_path / "scan.mha"
    sitk.WriteImage(_image(), str(path))
    candidate = catalog.discover(path)[0]
    candidate.size = (1_000, 1_000, 501)
    monkeypatch.setattr(
        sitk,
        "ReadImage",
        lambda _path: pytest.fail("pixel reading must not start above the voxel limit"),
    )

    with pytest.raises(ValueError, match="500,000,000.*--allow-large-volume"):
        volume.load(candidate)


def test_native_volume_override_allows_pixel_loading(tmp_path, monkeypatch):
    path = tmp_path / "scan.mha"
    sitk.WriteImage(_image(), str(path))
    candidate = catalog.discover(path)[0]
    candidate.size = (1_000, 1_000, 501)
    loaded = _image()
    monkeypatch.setattr(sitk, "ReadImage", lambda _path: loaded)

    assert volume.load(candidate, allow_large_volume=True).image is loaded


def test_file_metadata_is_shown_when_the_format_preserves_it(tmp_path):
    path = tmp_path / "scan.mha"
    sitk.WriteImage(_image(), str(path))

    candidate = catalog.discover(path)[0]

    assert candidate.modality == "MR"
    assert candidate.description == "synthetic volume"


def test_file_plane_is_derived_from_stored_direction(tmp_path):
    image = _image()
    image.SetDirection((1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0))
    path = tmp_path / "coronal.nrrd"
    sitk.WriteImage(image, str(path))

    candidate = catalog.discover(path)[0]

    assert candidate.plane == "coronal"


def test_mixed_directory_has_one_id_space_and_no_implicit_choice(tmp_path):
    _write_dicom_series(tmp_path / "dicom")
    sitk.WriteImage(_image(), str(tmp_path / "other.nrrd"))

    found = catalog.discover(tmp_path)

    assert [candidate.id for candidate in found] == [1, 2]
    assert {candidate.format for candidate in found} == {"DICOM", "NRRD"}
    assert catalog.recommended(found) is None
    with pytest.raises(ValueError, match="multiple usable volumes"):
        catalog.select(found, None)
    assert catalog.select(found, 1).id == 1
    assert catalog.select(found, 2).id == 2


def test_dicom_only_catalog_still_automatically_ranks_the_best_stack(tmp_path):
    coarse = _write_dicom_series(
        tmp_path / "coarse",
        series_number=1,
        slices=8,
        slice_spacing=2.0,
        pixel_spacing=(1.0, 1.0),
    )
    fine = _write_dicom_series(
        tmp_path / "fine",
        series_number=2,
        slices=6,
        slice_spacing=0.8,
        pixel_spacing=(0.4, 0.4),
    )

    found = catalog.discover(tmp_path)
    chosen = catalog.select(found, None)

    assert len(found) == 2
    assert chosen.dicom is not None
    assert chosen.dicom.uid == fine
    assert chosen.dicom.uid != coarse
    assert chosen.dicom.has_calibrated_hu


@pytest.mark.parametrize(
    "other_modality,expected",
    [("MR", "CT, MR"), (None, "CT, unknown")],
)
def test_dicom_catalog_spanning_modalities_requires_an_id(
    tmp_path, other_modality, expected
):
    _write_dicom_series(tmp_path / "ct", modality="CT", series_number=1)
    _write_dicom_series(
        tmp_path / "other",
        modality=other_modality,
        series_number=2,
    )

    found = catalog.discover(tmp_path)

    assert len(found) == 2
    assert catalog.recommended(found) is None
    with pytest.raises(ValueError, match="multiple modalities") as caught:
        catalog.select(found, None)
    assert expected in str(caught.value)
    assert catalog.select(found, found[0].id) is found[0]
    assert catalog.select(found, found[1].id) is found[1]


def test_only_integer_catalog_ids_are_accepted(tmp_path):
    _write_dicom_series(tmp_path)
    found = catalog.discover(tmp_path)

    assert catalog.select(found, 1) is found[0]
    with pytest.raises(ValueError, match="no volume has ID"):
        catalog.select(found, 99)


@pytest.mark.parametrize(
    "image,reason",
    [
        (sitk.Image([5, 5], sitk.sitkInt16), "expected a 3D volume"),
        (sitk.Image([5, 5, 5], sitk.sitkVectorUInt8, 3), "one scalar component"),
    ],
)
def test_non_scalar_or_non_3d_files_are_listed_as_unusable(tmp_path, image, reason):
    path = tmp_path / "invalid.mha"
    sitk.WriteImage(image, str(path))

    candidate = catalog.discover(path)[0]

    assert not candidate.usable
    assert reason in candidate.unusable_reason
    with pytest.raises(ValueError, match="no usable volumes"):
        catalog.select([candidate], None)


def test_file_source_identity_detects_the_same_file(tmp_path):
    path = tmp_path / "scan.mha"
    sitk.WriteImage(_image(), str(path))
    alias = tmp_path / "alias.mha"
    alias.hardlink_to(path)

    original = catalog.discover(path)[0]
    linked = catalog.discover(alias)[0]

    assert catalog.same_source(original, linked)


@pytest.mark.parametrize("suffix", [".mhd", ".nhdr"])
def test_detached_header_with_missing_payload_is_listed_as_unusable(tmp_path, suffix):
    path = tmp_path / ("scan" + suffix)
    sitk.WriteImage(_image(), str(path))
    complete = catalog.discover(path)[0]
    assert isinstance(complete.source, catalog.FileSource)
    assert len(complete.source.payload_paths) == 1
    complete.source.payload_paths[0].unlink()

    candidate = catalog.discover(path)[0]

    assert not candidate.usable
    assert "referenced payload is missing" in candidate.unusable_reason
    with pytest.raises(ValueError, match="no usable volumes"):
        catalog.select([candidate], None)
