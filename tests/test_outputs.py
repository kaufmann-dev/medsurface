"""Volume-output classification and verified atomic publication."""

from __future__ import annotations

import os
import weakref
from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk

from medsurface import outputs, volume


def _image() -> sitk.Image:
    values = np.arange(4 * 5 * 6, dtype=np.int16).reshape(4, 5, 6)
    image = sitk.GetImageFromArray(values)
    image.SetSpacing((0.7, 0.8, 1.2))
    image.SetOrigin((10.0, -4.0, 2.5))
    image.SetDirection((0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0))
    image.SetMetaData("custom_note", "preserve me")
    return image


def _read(path: Path) -> sitk.Image:
    name = path.name.casefold()
    if name.endswith((".nii", ".nii.gz")):
        image_io = "NiftiImageIO"
    elif name.endswith(".nrrd"):
        image_io = "NrrdImageIO"
    else:
        image_io = "MetaImageIO"
    return sitk.ReadImage(str(path), imageIO=image_io)


@pytest.mark.parametrize(
    "name,extension,format_name,compression",
    [
        ("volume.nii", ".nii", "NIfTI", "none"),
        ("volume.NII", ".nii", "NIfTI", "none"),
        ("volume.nii.gz", ".nii.gz", "NIfTI", "gzip"),
        ("volume.NII.GZ", ".nii.gz", "NIfTI", "gzip"),
        ("volume.nrrd", ".nrrd", "NRRD", "gzip"),
        ("volume.NRRD", ".nrrd", "NRRD", "gzip"),
        ("volume.mha", ".mha", "MetaImage", "zlib"),
        ("volume.MHA", ".mha", "MetaImage", "zlib"),
    ],
)
def test_volume_output_contract_is_suffix_driven(
    name, extension, format_name, compression
):
    contract = outputs.volume_output(name)

    assert contract.extension == extension
    assert contract.format == format_name
    assert contract.compression == compression
    assert contract.compressed is (compression != "none")


@pytest.mark.parametrize(
    "name",
    [
        "mesh.stl",
        "dicom.dcm",
        "detached.nhdr",
        "detached.mhd",
        "volume.raw",
        "volume",
    ],
)
def test_volume_output_rejects_non_atomic_or_unsupported_formats(name):
    with pytest.raises(ValueError, match="unsupported volume output extension"):
        outputs.volume_output(name)


def test_mesh_validation_reports_the_full_compound_extension():
    with pytest.raises(ValueError, match=r"'\.nii\.gz'"):
        outputs.validate_mesh_output("mask.nii.gz")


@pytest.mark.parametrize(
    "extension,format_name,compression",
    [
        (".nii", "NIfTI", "none"),
        (".nii.gz", "NIfTI", "gzip"),
        (".nrrd", "NRRD", "gzip"),
        (".mha", "MetaImage", "zlib"),
    ],
)
def test_verified_writer_preserves_the_exact_core_image_contract(
    tmp_path, extension, format_name, compression
):
    owner = [_image()]
    expected_dimension = owner[0].GetDimension()
    expected_size = owner[0].GetSize()
    expected_pixel_id = owner[0].GetPixelID()
    expected_components = owner[0].GetNumberOfComponentsPerPixel()
    expected_spacing = owner[0].GetSpacing()
    expected_origin = owner[0].GetOrigin()
    expected_direction = owner[0].GetDirection()
    expected_values = np.array(sitk.GetArrayViewFromImage(owner[0]), copy=True)
    destination = tmp_path / ("converted" + extension)

    result = volume.write_verified_volume(
        owner[0], str(destination), release_source=owner.clear
    )
    stored = _read(destination)

    assert result.output.format == format_name
    assert result.output.compression == compression
    assert stored.GetDimension() == expected_dimension
    assert stored.GetSize() == expected_size
    assert stored.GetPixelID() == expected_pixel_id
    assert stored.GetNumberOfComponentsPerPixel() == expected_components
    assert stored.GetSpacing() == pytest.approx(expected_spacing, abs=1e-5)
    assert stored.GetOrigin() == pytest.approx(expected_origin, abs=1e-5)
    assert stored.GetDirection() == pytest.approx(expected_direction, abs=1e-5)
    np.testing.assert_array_equal(
        sitk.GetArrayViewFromImage(stored),
        expected_values,
    )
    assert owner == []
    assert not list(tmp_path.glob(".converted*"))


@pytest.mark.parametrize("extension", [".NII", ".NII.GZ", ".NRRD", ".MHA"])
def test_verified_writer_accepts_case_insensitive_destinations(tmp_path, extension):
    owner = [_image()]
    expected_values = np.array(sitk.GetArrayViewFromImage(owner[0]), copy=True)
    destination = tmp_path / ("UPPER" + extension)

    volume.write_verified_volume(
        owner[0], str(destination), release_source=owner.clear
    )

    assert destination.is_file()
    stored = _read(destination)
    np.testing.assert_array_equal(
        sitk.GetArrayViewFromImage(stored),
        expected_values,
    )
    assert not list(tmp_path.glob(".UPPER*"))


@pytest.mark.parametrize(
    "extension,assert_compressed",
    [
        (".nii", lambda header: not header.startswith(b"\x1f\x8b")),
        (".nii.gz", lambda header: header.startswith(b"\x1f\x8b")),
        (".nrrd", lambda header: b"encoding: gzip" in header.lower()),
        (".mha", lambda header: b"compresseddata = true" in header.lower()),
    ],
)
def test_verified_writer_uses_suffix_selected_compression(
    tmp_path, extension, assert_compressed
):
    destination = tmp_path / ("compressed" + extension)
    owner = [_image()]

    volume.write_verified_volume(
        owner[0], str(destination), release_source=owner.clear
    )

    assert assert_compressed(destination.read_bytes()[:4096])


def _assert_failed_publication_is_clean(destination: Path) -> None:
    assert destination.read_bytes() == b"original"
    assert not list(destination.parent.glob(".%s.*" % destination.name))


@pytest.mark.parametrize(
    "failure",
    ["write", "read", "digest", "geometry", "replace", "cancel"],
)
def test_verified_writer_preserves_existing_destination_and_cleans_temp_files(
    tmp_path, monkeypatch, failure
):
    destination = tmp_path / "converted.nrrd"
    destination.write_bytes(b"original")
    owner = [_image()]

    if failure == "write":
        monkeypatch.setattr(
            volume.sitk,
            "WriteImage",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                RuntimeError("simulated write failure")
            ),
        )
        expected = RuntimeError
    elif failure == "read":
        monkeypatch.setattr(
            volume.sitk,
            "ReadImage",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                RuntimeError("simulated read failure")
            ),
        )
        expected = RuntimeError
    elif failure == "digest":
        real_digest = volume._voxel_digest
        calls = 0

        def mismatched_digest(stored):
            nonlocal calls
            calls += 1
            return real_digest(stored) if calls == 1 else "mismatch"

        monkeypatch.setattr(volume, "_voxel_digest", mismatched_digest)
        expected = ValueError
    elif failure == "geometry":
        real_read = volume.sitk.ReadImage

        def changed_geometry(path):
            stored = real_read(path)
            stored.SetOrigin((100.0, 200.0, 300.0))
            return stored

        monkeypatch.setattr(volume.sitk, "ReadImage", changed_geometry)
        expected = ValueError
    elif failure == "replace":
        monkeypatch.setattr(
            volume.os,
            "replace",
            lambda *_args: (_ for _ in ()).throw(OSError("simulated replace failure")),
        )
        expected = OSError
    else:
        monkeypatch.setattr(
            volume.sitk,
            "ReadImage",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(KeyboardInterrupt()),
        )
        expected = KeyboardInterrupt

    with pytest.raises(expected):
        volume.write_verified_volume(
            owner[0], str(destination), release_source=owner.clear
        )

    _assert_failed_publication_is_clean(destination)


def test_verified_writer_does_not_retry_or_mutate_metadata_after_write_failure(
    tmp_path, monkeypatch
):
    destination = tmp_path / "converted.nrrd"
    destination.write_bytes(b"original")
    owner = [_image()]
    calls = []

    def fail_once(*_args, **_kwargs):
        calls.append(True)
        raise RuntimeError("simulated one-shot write failure")

    monkeypatch.setattr(volume.sitk, "WriteImage", fail_once)

    with pytest.raises(RuntimeError, match="one-shot"):
        volume.write_verified_volume(
            owner[0], str(destination), release_source=owner.clear
        )

    assert calls == [True]
    assert owner[0].GetMetaData("custom_note") == "preserve me"
    _assert_failed_publication_is_clean(destination)


def test_verified_writer_releases_the_source_before_readback(tmp_path, monkeypatch):
    owner = [_image()]
    source_reference = weakref.ref(owner[0])
    real_read = volume.sitk.ReadImage
    observed = []

    def read_after_release(path):
        observed.append(True)
        assert owner == []
        assert source_reference() is None
        return real_read(path)

    monkeypatch.setattr(volume.sitk, "ReadImage", read_after_release)

    volume.write_verified_volume(
        owner[0],
        str(tmp_path / "released.nrrd"),
        release_source=owner.clear,
    )

    assert observed == [True]


def test_verified_writer_uses_a_same_filesystem_lowercase_suffix_temp(
    tmp_path, monkeypatch
):
    destination = tmp_path / "FINAL.NII.GZ"
    real_write = volume.sitk.WriteImage
    temporary_paths = []

    def capture_path(image, path, **kwargs):
        temporary_paths.append(Path(path))
        return real_write(image, path, **kwargs)

    monkeypatch.setattr(volume.sitk, "WriteImage", capture_path)

    owner = [_image()]
    volume.write_verified_volume(
        owner[0], str(destination), release_source=owner.clear
    )

    assert len(temporary_paths) == 1
    temporary = temporary_paths[0]
    assert temporary.parent == destination.parent
    assert temporary.name.endswith(".nii.gz")
    assert not temporary.exists()
    assert destination.exists()
    assert os.stat(destination).st_dev == os.stat(destination.parent).st_dev
