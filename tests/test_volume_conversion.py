"""Strict one-volume storage conversion contracts."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pydicom
import pytest
import SimpleITK as sitk
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, generate_uid

from medsurface import catalog, volume

INPUT_EXTENSIONS = (".nii", ".nii.gz", ".nrrd", ".nhdr", ".mha", ".mhd")
OUTPUT_EXTENSIONS = (".nii", ".nii.gz", ".nrrd", ".mha")


def _image() -> sitk.Image:
    values = (np.arange(4 * 5 * 6, dtype=np.int16) - 50).reshape(4, 5, 6)
    image = sitk.GetImageFromArray(values)
    image.SetSpacing((0.7, 0.8, 1.2))
    image.SetOrigin((10.0, -4.0, 2.5))
    image.SetDirection((0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0))
    return image


def _read(path: Path) -> sitk.Image:
    name = path.name.casefold()
    if name.endswith((".nii", ".nii.gz")):
        image_io = "NiftiImageIO"
    elif name.endswith((".nrrd", ".nhdr")):
        image_io = "NrrdImageIO"
    else:
        image_io = "MetaImageIO"
    return sitk.ReadImage(str(path), imageIO=image_io)


def _assert_same_image(actual: sitk.Image, expected: sitk.Image) -> None:
    assert actual.GetDimension() == expected.GetDimension()
    assert actual.GetSize() == expected.GetSize()
    assert actual.GetPixelID() == expected.GetPixelID()
    assert (
        actual.GetNumberOfComponentsPerPixel()
        == expected.GetNumberOfComponentsPerPixel()
    )
    assert actual.GetSpacing() == pytest.approx(expected.GetSpacing(), abs=1e-5)
    assert actual.GetOrigin() == pytest.approx(expected.GetOrigin(), abs=1e-5)
    assert actual.GetDirection() == pytest.approx(expected.GetDirection(), abs=1e-5)
    np.testing.assert_array_equal(
        sitk.GetArrayViewFromImage(actual),
        sitk.GetArrayViewFromImage(expected),
    )


@pytest.mark.parametrize("input_extension", INPUT_EXTENSIONS)
@pytest.mark.parametrize("output_extension", OUTPUT_EXTENSIONS)
def test_every_file_input_converts_to_every_atomic_volume_output(
    tmp_path, input_extension, output_extension
):
    source = tmp_path / ("source" + input_extension)
    sitk.WriteImage(_image(), str(source), useCompression=True)
    candidate = catalog.select(catalog.discover(source), None)
    expected = volume.load(candidate).image
    destination = tmp_path / ("converted" + output_extension)

    result = volume.convert(candidate, str(destination))
    stored = _read(destination)

    _assert_same_image(stored, expected)
    assert result.output_path == str(destination)
    assert result.dimensions == expected.GetSize()
    assert result.pixel_type == expected.GetPixelIDTypeAsString()
    assert result.components == expected.GetNumberOfComponentsPerPixel()
    assert result.spacing == pytest.approx(expected.GetSpacing())
    assert result.origin == pytest.approx(expected.GetOrigin())
    assert result.direction == pytest.approx(expected.GetDirection())
    assert result.metadata_policy == "preserved-best-effort"
    assert result.provenance["path"] == str(source)
    assert result.provenance["components"] == 1


@pytest.mark.parametrize("extension", OUTPUT_EXTENSIONS)
def test_same_format_conversion_uses_a_distinct_atomic_destination(tmp_path, extension):
    source = tmp_path / ("source" + extension)
    sitk.WriteImage(_image(), str(source), useCompression=True)
    candidate = catalog.select(catalog.discover(source), None)
    expected = volume.load(candidate).image
    destination = tmp_path / ("copy" + extension)

    volume.convert(candidate, str(destination))

    _assert_same_image(_read(destination), expected)
    assert source.exists()
    assert destination.exists()


@pytest.mark.parametrize("output_extension", (".NII", ".NII.GZ", ".NRRD", ".MHA"))
def test_conversion_accepts_case_insensitive_output_suffixes(
    tmp_path, output_extension
):
    source = tmp_path / "source.nrrd"
    sitk.WriteImage(_image(), str(source), useCompression=True)
    candidate = catalog.select(catalog.discover(source), None)
    expected = volume.load(candidate).image
    destination = tmp_path / ("UPPER" + output_extension)

    volume.convert(candidate, str(destination))

    _assert_same_image(_read(destination), expected)


@pytest.mark.parametrize("output_extension", (".nrrd", ".mha"))
def test_conversion_preserves_representable_source_metadata(
    tmp_path, output_extension
):
    image = _image()
    image.SetMetaData("custom_note", "keep this value")
    image.SetMetaData("modality", "MR")
    source = tmp_path / "source.nrrd"
    sitk.WriteImage(image, str(source), useCompression=True)
    candidate = catalog.select(catalog.discover(source), None)
    destination = tmp_path / ("preserved" + output_extension)

    result = volume.convert(candidate, str(destination))
    stored = _read(destination)

    assert stored.GetMetaData("custom_note") == "keep this value"
    assert stored.GetMetaData("modality") == "MR"
    assert result.metadata_policy == "preserved-best-effort"


@pytest.mark.parametrize("output_extension", OUTPUT_EXTENSIONS)
def test_conversion_strips_source_metadata_but_keeps_required_headers(
    tmp_path, output_extension
):
    image = _image()
    image.SetMetaData("custom_note", "remove this value")
    source = tmp_path / "source.nrrd"
    sitk.WriteImage(image, str(source), useCompression=True)
    candidate = catalog.select(catalog.discover(source), None)
    destination = tmp_path / ("stripped" + output_extension)

    result = volume.convert(candidate, str(destination), strip_metadata=True)
    stored = _read(destination)

    assert not stored.HasMetaDataKey("custom_note")
    assert result.metadata_policy == "stripped"
    assert stored.GetSize() == image.GetSize()
    assert stored.GetSpacing() == pytest.approx(image.GetSpacing(), abs=1e-5)
    assert stored.GetMetaDataKeys()


def test_conversion_warns_when_the_destination_cannot_represent_all_metadata(
    tmp_path,
):
    image = _image()
    image.SetMetaData("custom_note", "NIfTI cannot represent this custom field")
    source = tmp_path / "source.nrrd"
    sitk.WriteImage(image, str(source), useCompression=True)
    candidate = catalog.select(catalog.discover(source), None)
    emitted = []

    result = volume.convert(
        candidate,
        str(tmp_path / "converted.nii.gz"),
        warn=emitted.append,
    )

    assert any("metadata" in warning and "omitted" in warning for warning in emitted)
    assert result.warnings == emitted


def _write_ct_stack(root: Path) -> None:
    root.mkdir()
    series_uid = generate_uid()
    for index in range(6):
        path = root / ("slice-%03d.dcm" % index)
        sop_uid = generate_uid()
        file_meta = FileMetaDataset()
        file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
        file_meta.MediaStorageSOPClassUID = CTImageStorage
        file_meta.MediaStorageSOPInstanceUID = sop_uid
        dataset = FileDataset(
            str(path), {}, file_meta=file_meta, preamble=b"\0" * 128
        )
        dataset.SOPClassUID = CTImageStorage
        dataset.SOPInstanceUID = sop_uid
        dataset.SeriesInstanceUID = series_uid
        dataset.Modality = "CT"
        dataset.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
        dataset.SeriesDescription = "conversion CT"
        dataset.SeriesNumber = 1
        dataset.InstanceNumber = index + 1
        dataset.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        dataset.ImagePositionPatient = [0, 0, index]
        dataset.PixelSpacing = [1, 1]
        dataset.SliceThickness = 1
        dataset.Rows = 8
        dataset.Columns = 8
        dataset.SamplesPerPixel = 1
        dataset.PhotometricInterpretation = "MONOCHROME2"
        dataset.BitsAllocated = 16
        dataset.BitsStored = 16
        dataset.HighBit = 15
        dataset.PixelRepresentation = 1
        dataset.RescaleSlope = 1
        dataset.RescaleIntercept = -1024
        dataset.PixelData = np.full((8, 8), index, dtype=np.int16).tobytes()
        pydicom.dcmwrite(path, dataset, enforce_file_format=True)


@pytest.mark.parametrize("output_extension", OUTPUT_EXTENSIONS)
def test_synthetic_dicom_converts_without_changing_the_loaded_volume(
    tmp_path, output_extension
):
    source = tmp_path / "dicom"
    _write_ct_stack(source)
    candidate = catalog.select(catalog.discover(source), None)
    expected = volume.load(candidate).image
    destination = tmp_path / ("dicom-converted" + output_extension)

    result = volume.convert(candidate, str(destination))

    _assert_same_image(_read(destination), expected)
    assert result.provenance["format"] == "DICOM"
    assert result.provenance["dicom"]["series_uid"] == candidate.dicom.uid


def test_conversion_rejects_an_unsupported_output_before_loading(tmp_path, monkeypatch):
    source = tmp_path / "source.nrrd"
    sitk.WriteImage(_image(), str(source))
    candidate = catalog.select(catalog.discover(source), None)
    monkeypatch.setattr(
        volume,
        "load",
        lambda *_args, **_kwargs: pytest.fail(
            "pixel loading must not start for an unsupported output"
        ),
    )

    with pytest.raises(ValueError, match=r"unsupported.*\.mhd"):
        volume.convert(candidate, str(tmp_path / "detached.mhd"))


def test_conversion_rejects_a_core_contract_the_target_cannot_preserve(tmp_path):
    image = _image()
    image.SetDirection((1.0, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    source = tmp_path / "sheared.nrrd"
    sitk.WriteImage(image, str(source))
    candidate = catalog.select(catalog.discover(source), None)
    destination = tmp_path / "converted.nii"
    destination.write_bytes(b"original")

    with pytest.raises(ValueError, match="direction changed"):
        volume.convert(candidate, str(destination))

    assert destination.read_bytes() == b"original"
    assert not list(tmp_path.glob(".converted.nii.*"))
