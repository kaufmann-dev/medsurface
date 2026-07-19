"""Real-reader end-to-end workflows across supported input families."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pydicom
import SimpleITK as sitk
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, generate_uid

from medsurface import catalog, fusion, labelmap, pipeline, presets


def _minimal_preset(name: str):
    return replace(
        presets.get(name),
        median_mm=0.0,
        closing_mm=0.0,
        opening_mm=0.0,
        min_island_mm3=0.0,
        resample_mm=0.0,
        mask_smooth_mm=0,
        surface_smooth_iters=0,
        simplify_error_mm=0.0,
        post_surface_smooth_iters=0,
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
        ds.add_new((0x0008, 0x0080), "LO", b"Clinic \xfc")
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


def test_real_dicom_stack_extracts_to_a_valid_mesh(tmp_path):
    source = tmp_path / "dicom"
    _write_ct_stack(source)
    candidate = catalog.discover(source)[0]
    output = tmp_path / "dicom.stl"

    result = pipeline.extract(
        candidate,
        _minimal_preset("bone"),
        str(output),
        threshold=-1022.0,
    )

    assert output.exists()
    assert result.quality["valid"]
    assert result.quality["triangles"] > 0


def test_left_handed_nifti_extracts_to_a_valid_outward_wound_mesh(tmp_path):
    zz, yy, xx = np.indices((24, 24, 24))
    pixels = np.zeros((24, 24, 24), dtype=np.int16)
    pixels[(xx - 12) ** 2 + (yy - 12) ** 2 + (zz - 12) ** 2 <= 7**2] = 2_000
    image = sitk.GetImageFromArray(pixels)
    image.SetDirection((-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    source = tmp_path / "left-handed.nii.gz"
    output = tmp_path / "left-handed.stl"
    sitk.WriteImage(image, str(source))

    candidate = catalog.discover(source)[0]
    result = pipeline.extract(
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


def test_intensity_fuse_then_labelmap_extract_produces_a_valid_mesh(tmp_path):
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
    fused_output = tmp_path / "fused.nrrd"
    mesh_output = tmp_path / "fused.stl"

    fusion_preset = replace(
        _minimal_preset("teeth"),
        mask_smooth_mm=9.0,
        surface_smooth_iters=77,
        simplify_error_mm=4.0,
        post_surface_smooth_iters=88,
    )
    fused = fusion.fuse(
        fixed,
        moving,
        fusion_preset,
        str(fused_output),
        fixed_threshold=1000.0,
        moving_threshold=1000.0,
        grid_mm=1.0,
        allow_large_volume=True,
    )
    fused_image = sitk.ReadImage(str(fused_output))
    extracted = labelmap.extract(
        catalog.discover(fused_output)[0],
        str(mesh_output),
    )

    assert fused_output.exists()
    assert set(np.unique(sitk.GetArrayViewFromImage(fused_image))) == {0, 1}
    assert fused.foreground_fused_voxels == int(
        sitk.GetArrayViewFromImage(fused_image).sum()
    )
    assert set(fused.provenance["segmentation"]) == {
        "name",
        "description",
        "threshold",
        "threshold_unit",
        "threshold_max",
        "median_mm",
        "closing_mm",
        "opening_mm",
        "min_island_mm3",
        "keep_largest_island",
    }
    assert "surface_finishing" not in fused.provenance
    assert fused.provenance["allow_large_volume"] is True
    assert mesh_output.exists()
    assert extracted.quality["valid"]
    assert extracted.quality["components"] == 2
    assert extracted.provenance["surface"]["mask_smooth_mm"] == 0.8
    assert extracted.provenance["surface"]["surface_smooth_iters"] == 20
    assert extracted.provenance["surface"]["post_surface_smooth_iters"] == 0


def test_labelmap_fuse_then_labelmap_extract_produces_a_valid_mesh(tmp_path):
    zz, yy, xx = np.indices((48, 48, 48))
    foreground = (
        ((xx - 24) / 14) ** 2
        + ((yy - 23) / 11) ** 2
        + ((zz - 22) / 9) ** 2
        <= 1
    )
    foreground[16:24, 29:36, 33:40] = True
    fixed_image = sitk.GetImageFromArray((foreground * 5).astype(np.uint16))
    moving_image = sitk.GetImageFromArray((foreground * 203).astype(np.uint16))
    moving_image.SetOrigin((3.0, -2.0, 1.0))
    fixed_path = tmp_path / "fixed-labels.nii.gz"
    moving_path = tmp_path / "moving-labels.mha"
    fused_path = tmp_path / "fused-labels.mha"
    mesh_path = tmp_path / "fused-labels.stl"
    sitk.WriteImage(fixed_image, str(fixed_path))
    sitk.WriteImage(moving_image, str(moving_path))

    fused = labelmap.fuse(
        catalog.discover(fixed_path)[0],
        catalog.discover(moving_path)[0],
        str(fused_path),
        grid_mm=1.0,
    )
    fused_image = sitk.ReadImage(str(fused_path))
    extracted = labelmap.extract(
        catalog.discover(fused_path)[0],
        str(mesh_path),
    )

    assert set(np.unique(sitk.GetArrayViewFromImage(fused_image))) == {0, 1}
    assert fused.registration.surface_overlap > 0.9
    assert fused.registration.shared_fov_dice > 0.9
    assert fused.provenance["segmentation"]["foreground"].endswith(
        "normalized to 1"
    )
    assert mesh_path.exists()
    assert extracted.quality["valid"]
    assert extracted.quality["components"] == 1
