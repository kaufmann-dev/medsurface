"""Output classification and atomic NIfTI publishing."""

import numpy as np
import pytest
import SimpleITK as sitk

from medsurface import outputs, volume


@pytest.mark.parametrize("name", ["mask.nii", "mask.NII", "mask.nii.gz", "mask.NII.GZ"])
def test_merge_output_recognizes_nifti_extensions(name):
    assert outputs.merge_output_kind(name) is outputs.OutputKind.NIFTI


@pytest.mark.parametrize("name", ["mesh.stl", "mesh.PLY", "mesh.obj"])
def test_merge_output_recognizes_mesh_extensions(name):
    assert outputs.merge_output_kind(name) is outputs.OutputKind.MESH


def test_merge_output_reports_the_full_compound_extension():
    with pytest.raises(ValueError, match=r"'\.nii\.gz'"):
        outputs.validate_mesh_output("mask.nii.gz")


@pytest.mark.parametrize("extension", [".nii", ".nii.gz"])
def test_binary_nifti_is_atomic_and_preserves_geometry(tmp_path, extension):
    values = np.zeros((6, 7, 8), dtype=np.float32)
    values[1:5, 2:6, 3:7] = 0.75
    image = sitk.GetImageFromArray(values)
    image.SetSpacing((0.7, 0.8, 0.9))
    image.SetOrigin((10.0, 20.0, 30.0))
    image.SetDirection((0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0))
    destination = tmp_path / ("merged" + extension)

    written = volume.write_binary_nifti(image, str(destination))
    stored = sitk.ReadImage(str(destination))

    assert written.GetPixelID() == sitk.sitkUInt8
    assert stored.GetPixelID() == sitk.sitkUInt8
    assert stored.GetSize() == image.GetSize()
    assert stored.GetSpacing() == pytest.approx(image.GetSpacing(), abs=1e-5)
    assert stored.GetOrigin() == pytest.approx(image.GetOrigin(), abs=1e-5)
    assert stored.GetDirection() == pytest.approx(image.GetDirection(), abs=1e-5)
    assert set(np.unique(sitk.GetArrayViewFromImage(stored))) == {0, 1}
    assert not list(tmp_path.glob(".merged*"))


def test_binary_nifti_preserves_existing_destination_after_validation_failure(
    tmp_path, monkeypatch
):
    destination = tmp_path / "merged.nii.gz"
    destination.write_bytes(b"original")
    image = sitk.Image((4, 4, 4), sitk.sitkUInt8)
    monkeypatch.setattr(
        volume.sitk,
        "ReadImage",
        lambda _path: (_ for _ in ()).throw(RuntimeError("simulated read failure")),
    )

    with pytest.raises(RuntimeError, match="simulated read failure"):
        volume.write_binary_nifti(image, str(destination))

    assert destination.read_bytes() == b"original"
    assert not list(tmp_path.glob(".merged.nii.gz.*"))
