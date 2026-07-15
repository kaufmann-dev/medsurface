"""Real-reader end-to-end workflows across supported input families."""

from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pydicom
import SimpleITK as sitk
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, generate_uid

from medsurface import catalog, merge, pipeline, presets


def _minimal_preset(name: str):
    return replace(
        presets.get(name),
        median_mm=0.0,
        closing_mm=0.0,
        opening_mm=0.0,
        min_island_mm3=0.0,
        resample_mm=0.0,
        smooth_mm=0,
        simplify_error_mm=0.0,
    )


def _write_ct_stack(root: Path) -> None:
    root.mkdir()
    series_uid = generate_uid()
    for index in range(6):
        path = root / ("slice-%03d.dcm" % index)
        sop_uid = generate_uid()
        meta = FileMetaDataset()
        meta.TransferSyntaxUID = ExplicitVRLittleEndian
        meta.MediaStorageSOPClassUID = CTImageStorage
        meta.MediaStorageSOPInstanceUID = sop_uid
        ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
        ds.SOPClassUID = CTImageStorage
        ds.SOPInstanceUID = sop_uid
        ds.SeriesInstanceUID = series_uid
        ds.Modality = "CT"
        ds.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
        ds.SeriesDescription = "workflow CT"
        ds.SeriesNumber = 1
        ds.InstanceNumber = index + 1
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        ds.ImagePositionPatient = [0, 0, index]
        ds.PixelSpacing = [1, 1]
        ds.SliceThickness = 1
        ds.Rows = 8
        ds.Columns = 8
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 1
        ds.RescaleSlope = 1
        ds.RescaleIntercept = -1024
        pixels = np.full((8, 8), index, dtype=np.int16)
        ds.PixelData = pixels.tobytes()
        pydicom.dcmwrite(path, ds, enforce_file_format=True)


def test_real_dicom_stack_converts_to_a_valid_mesh(tmp_path):
    source = tmp_path / "dicom"
    _write_ct_stack(source)
    candidate = catalog.discover(source)[0]
    output = tmp_path / "dicom.stl"

    result = pipeline.convert(
        candidate,
        _minimal_preset("bone"),
        str(output),
        threshold=-1022.0,
    )

    assert output.exists()
    assert result.quality["valid"]
    assert result.quality["triangles"] > 0


def test_left_handed_nifti_converts_to_a_valid_outward_wound_mesh(tmp_path):
    zz, yy, xx = np.indices((24, 24, 24))
    pixels = np.zeros((24, 24, 24), dtype=np.int16)
    pixels[(xx - 12) ** 2 + (yy - 12) ** 2 + (zz - 12) ** 2 <= 7**2] = 2_000
    image = sitk.GetImageFromArray(pixels)
    image.SetDirection((-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    source = tmp_path / "left-handed.nii.gz"
    output = tmp_path / "left-handed.stl"
    sitk.WriteImage(image, str(source))

    candidate = catalog.discover(source)[0]
    result = pipeline.convert(
        candidate,
        _minimal_preset("bone"),
        str(output),
        threshold=1_000.0,
        allow_large_volume=True,
    )

    assert output.exists()
    assert result.quality["valid"]
    assert result.quality["winding_consistent"]
    assert result.quality["disoriented_faces"] == 0
    assert result.provenance["allow_large_volume"] is True


def test_real_cross_format_merge_preserves_teeth_components(tmp_path):
    pixels = np.zeros((32, 32, 32), dtype=np.int16)
    pixels[4:12, 4:12, 4:12] = 2000
    pixels[20:28, 20:28, 20:28] = 2000
    image = sitk.GetImageFromArray(pixels)
    image.SetSpacing((1.0, 1.0, 1.0))
    fixed_path = tmp_path / "fixed.mha"
    moving_path = tmp_path / "moving.nrrd"
    sitk.WriteImage(image, str(fixed_path))
    sitk.WriteImage(image, str(moving_path))
    fixed = catalog.discover(fixed_path)[0]
    moving = catalog.discover(moving_path)[0]
    output = tmp_path / "merged.stl"

    preset = _minimal_preset("teeth")
    effective_preset = replace(
        preset,
        smooth_mm=0.2,
        simplify_error_mm=0.01,
    )
    result = merge.merge(
        fixed,
        moving,
        preset,
        str(output),
        fixed_threshold=1000.0,
        moving_threshold=1000.0,
        grid_mm=1.0,
        smooth_mm=effective_preset.smooth_mm,
        simplify_error_mm=effective_preset.simplify_error_mm,
        allow_large_volume=True,
    )

    assert output.exists()
    assert result.quality["valid"]
    assert result.quality["components"] == 2
    assert result.surface_components == 2
    assert result.provenance["preset"] == asdict(effective_preset)
    assert result.provenance["allow_large_volume"] is True
    finishing = result.provenance["surface_finishing"]
    assert finishing["smoothing"]["requested_iterations"] == 20
    assert finishing["decimation"]["simplify_error_mm"] == 0.01
