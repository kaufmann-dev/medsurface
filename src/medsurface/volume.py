"""Load a selected catalog volume into a SimpleITK image."""

from __future__ import annotations

import hashlib
import math
import os
import tempfile
import time
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
import SimpleITK as sitk

from .catalog import DicomSource, FileSource, VolumeCandidate, validate_image
from .defaults import MAX_VOXELS, MIN_VOLUME_AXIS_VOXELS
from .outputs import VolumeOutput, volume_output


@dataclass
class Volume:
    image: sitk.Image
    candidate: VolumeCandidate
    metadata_omissions: int = 0

    @property
    def spacing(self) -> tuple[float, float, float]:
        return tuple(self.image.GetSpacing())  # type: ignore[return-value]

    @property
    def size(self) -> tuple[int, int, int]:
        return tuple(self.image.GetSize())  # type: ignore[return-value]

    def intensity_range(self) -> tuple[float, float]:
        f = sitk.MinimumMaximumImageFilter()
        f.Execute(self.image)
        return float(f.GetMinimum()), float(f.GetMaximum())

    def array(self) -> np.ndarray:
        """Voxels as ``(z, y, x)``."""
        return sitk.GetArrayViewFromImage(self.image)


@dataclass(frozen=True)
class VolumeWriteResult:
    """Facts observed while serializing and reading back one volume."""

    output: VolumeOutput
    source_metadata: int
    preserved_metadata: int


@dataclass
class ConversionResult:
    """Strict one-volume storage-format conversion result."""

    output_path: str
    format: str
    compression: str
    dimensions: tuple[int, ...]
    pixel_type: str
    components: int
    spacing: tuple[float, ...]
    origin: tuple[float, ...]
    direction: tuple[float, ...]
    seconds: float
    warnings: list[str]
    metadata_policy: str
    provenance: dict[str, Any]


@dataclass(frozen=True)
class _ImageContract:
    digest: str
    dimension: int
    size: tuple[int, ...]
    pixel_id: int
    pixel_type: str
    components: int
    spacing: tuple[float, ...]
    origin: tuple[float, ...]
    direction: tuple[float, ...]


def _voxel_count(candidate: VolumeCandidate) -> int:
    size = candidate.size
    if size is None or len(size) != 3:
        raise ValueError(
            "volume ID %d has no complete positive 3D dimensions; "
            "pixel data will not be loaded without a usable header" % candidate.id
        )
    if any(value < MIN_VOLUME_AXIS_VOXELS for value in size):
        raise ValueError(
            "each volume axis must contain at least %d voxels; got %s"
            % (
                MIN_VOLUME_AXIS_VOXELS,
                "x".join(str(value) for value in size),
            )
        )
    return math.prod(int(value) for value in size)


def _copy_dicom_metadata(
    reader: sitk.ImageSeriesReader,
    image: sitk.Image,
) -> int:
    """Promote valid first-slice metadata without guessing malformed text."""
    omissions = 0
    for key in reader.GetMetaDataKeys(0):
        try:
            value = reader.GetMetaData(0, key)
            value.encode("utf-8", errors="strict")
        except (RuntimeError, UnicodeEncodeError):
            omissions += 1
            continue
        image.SetMetaData(key, value)
    return omissions


def _erase_metadata(image: sitk.Image) -> None:
    for key in image.GetMetaDataKeys():
        image.EraseMetaData(key)


def load(
    candidate: VolumeCandidate,
    *,
    allow_large_volume: bool = False,
    preserve_metadata: bool = False,
) -> Volume:
    """Load pixels for one already-selected candidate."""
    if not candidate.usable:
        raise ValueError("volume ID %d is unusable: %s" % (candidate.id, candidate.unusable_reason))
    voxels = _voxel_count(candidate)
    if voxels > MAX_VOXELS and not allow_large_volume:
        raise ValueError(
            "volume contains %s voxels, above the default limit of %s; resample it "
            "to a coarser spacing before processing, or pass --allow-large-volume "
            "to attempt it (this may exhaust memory)"
            % (f"{voxels:,}", f"{MAX_VOXELS:,}")
        )
    metadata_omissions = 0
    try:
        if isinstance(candidate.source, DicomSource):
            series = candidate.source.series
            if series.n_slices < 2:
                raise ValueError(
                    "series %r has %d slice(s)" % (series.description, series.n_slices)
                )
            reader = sitk.ImageSeriesReader()
            reader.SetFileNames(series.files)  # already ordered by physical position
            if preserve_metadata:
                reader.MetaDataDictionaryArrayUpdateOn()
                reader.LoadPrivateTagsOn()
            image = reader.Execute()
            if preserve_metadata and series.files:
                metadata_omissions = _copy_dicom_metadata(reader, image)
        else:
            source_path = str(candidate.source.path)
            if candidate.source.path.name == candidate.source.path.name.casefold():
                image = sitk.ReadImage(source_path)
            else:
                image = sitk.ReadImage(
                    source_path,
                    imageIO=_image_io(candidate.source.path),
                )
            if not preserve_metadata:
                _erase_metadata(image)
    except RuntimeError as exc:
        detail = next(
            (line.strip() for line in reversed(str(exc).splitlines()) if line.strip()),
            "the image reader rejected the input",
        )
        companion = ""
        if isinstance(candidate.source, FileSource) and candidate.source.path.name.casefold().endswith(
            (".mhd", ".nhdr")
        ):
            companion = "; ensure the referenced payload exists and is readable"
        raise ValueError("cannot load %s: %s%s" % (candidate.source_name, detail, companion)) from None

    reason = validate_image(image)
    if reason:
        raise ValueError("volume ID %d is unusable: %s" % (candidate.id, reason))
    return Volume(
        image=image,
        candidate=candidate,
        metadata_omissions=metadata_omissions,
    )


def has_calibrated_hu(candidate: VolumeCandidate) -> bool:
    """Whether the loader can verify that values are Hounsfield units."""
    return bool(candidate.dicom and candidate.dicom.has_calibrated_hu)


def loading_warnings_for(volume: Volume) -> list[str]:
    """Non-fatal observations relevant to loading without image processing."""
    out: list[str] = []
    series = volume.candidate.dicom

    if volume.metadata_omissions:
        noun = "value" if volume.metadata_omissions == 1 else "values"
        out.append(
            "omitted %d DICOM metadata %s with invalid text encoding; voxel data "
            "and geometry are unchanged" % (volume.metadata_omissions, noun)
        )

    if series is not None and not series.spacing_uniform:
        out.append(
            "slice spacing is not uniform (spread %.3f mm); the volume will be "
            "resampled onto a regular grid and geometry may shift slightly"
            % series.spacing_spread_mm
        )
    return out


def warnings_for(volume: Volume) -> list[str]:
    """Non-fatal data-quality observations relevant to surface workflows."""
    out = loading_warnings_for(volume)
    candidate = volume.candidate
    series = candidate.dicom

    sx, sy, sz = volume.spacing
    aniso = max(sx, sy, sz) / min(sx, sy, sz)
    if aniso > 2.0:
        out.append(
            "voxels are strongly anisotropic (%.2f x %.2f x %.2f mm, ratio %.1f:1); "
            "expect stair-stepping along the thick axis" % (sx, sy, sz, aniso)
        )

    if series is not None and series.sharp_kernel:
        out.append(
            "reconstruction kernel %s is a sharp/edge-enhancing kernel: it amplifies "
            "noise, so a low bone threshold will produce a spiky surface. Prefer a "
            "higher threshold (~300 HU rather than ~200 HU), or a smoother kernel "
            "series if the study has one." % series.kernel_display
        )

    return out


def _image_io(path: os.PathLike[str] | str) -> str:
    name = os.fspath(path).casefold()
    if name.endswith((".nii", ".nii.gz")):
        return "NiftiImageIO"
    if name.endswith((".nrrd", ".nhdr")):
        return "NrrdImageIO"
    if name.endswith((".mha", ".mhd")):
        return "MetaImageIO"
    return ""


def provenance_for(volume: Volume) -> dict[str, Any]:
    """Format-neutral provenance shared by every volume workflow."""
    candidate = volume.candidate
    if isinstance(candidate.source, DicomSource):
        path = str(candidate.source.catalog_path)
    else:
        path = str(candidate.source.path)
    record: dict[str, Any] = {
        "id": candidate.id,
        "format": candidate.format,
        "path": path,
        "source": candidate.source_name,
        "size": list(volume.size),
        "pixel_type": volume.image.GetPixelIDTypeAsString(),
        "components": volume.image.GetNumberOfComponentsPerPixel(),
        "spacing_mm": list(volume.spacing),
        "origin_mm": list(volume.image.GetOrigin()),
        "direction": list(volume.image.GetDirection()),
        "modality": candidate.modality,
        "description": candidate.description,
        "plane": candidate.plane,
    }
    if candidate.dicom is not None:
        series = candidate.dicom
        record["dicom"] = {
            "series_uid": series.uid,
            "series_orientation_part": [series.part, series.n_parts],
            "series_number": series.series_number,
            "description": series.description,
            "convolution_kernel": list(series.kernel_values),
            "slices": series.n_slices,
            "image_type": list(series.image_type),
            "rescale_type": series.rescale_type,
            "multi_energy_ct_acquisition": series.multi_energy_ct_acquisition,
            "hu_calibration_consistent": series.hu_calibration_consistent,
        }
    return record


def touches_boundary(image: sitk.Image, value: int = 1) -> bool:
    """True if any foreground voxel lies on the outer face of the volume.

    When this holds, marching cubes leaves the surface open there unless the
    volume is padded -- the classic cause of a "watertight" pipeline emitting an
    open mesh. Anatomy cut off by the scanner's field of view triggers it.
    """
    a = sitk.GetArrayViewFromImage(image)
    if a.size == 0:
        return False
    faces = (
        a[0, :, :], a[-1, :, :],
        a[:, 0, :], a[:, -1, :],
        a[:, :, 0], a[:, :, -1],
    )
    return any(bool(np.any(f == value)) for f in faces)


def _voxel_digest(image: sitk.Image) -> str:
    """Hash voxel bytes without materializing another complete volume."""
    values = sitk.GetArrayViewFromImage(image)
    digest = hashlib.sha256()
    if values.ndim <= 1:
        digest.update(values.tobytes(order="C"))
    else:
        for plane in values:
            digest.update(np.ascontiguousarray(plane).tobytes(order="C"))
    return digest.hexdigest()


def _image_contract(image: sitk.Image) -> _ImageContract:
    return _ImageContract(
        digest=_voxel_digest(image),
        dimension=image.GetDimension(),
        size=tuple(int(value) for value in image.GetSize()),
        pixel_id=image.GetPixelID(),
        pixel_type=image.GetPixelIDTypeAsString(),
        components=image.GetNumberOfComponentsPerPixel(),
        spacing=tuple(float(value) for value in image.GetSpacing()),
        origin=tuple(float(value) for value in image.GetOrigin()),
        direction=tuple(float(value) for value in image.GetDirection()),
    )


def _metadata(image: sitk.Image) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for key in image.GetMetaDataKeys():
        try:
            metadata[key] = image.GetMetaData(key)
        except RuntimeError:
            continue
    return metadata


def _geometry_tolerance(
    expected: tuple[float, ...],
    output: VolumeOutput,
) -> np.ndarray:
    tolerance = np.full(len(expected), 1e-5, dtype=np.float64)
    expected_values = np.asarray(expected, dtype=np.float64)
    if (
        output.format == "NIfTI"
        and np.all(np.abs(expected_values) <= np.finfo(np.float32).max)
    ):
        encoded = expected_values.astype(np.float32)
        ulps = np.abs(np.spacing(encoded).astype(np.float64))
        tolerance = np.maximum(tolerance, ulps)
    return tolerance


def _verify_contract(
    stored: sitk.Image,
    expected: _ImageContract,
    output: VolumeOutput,
) -> None:
    checks = (
        ("dimension", stored.GetDimension(), expected.dimension),
        ("dimensions", tuple(stored.GetSize()), expected.size),
        ("pixel type", stored.GetPixelID(), expected.pixel_id),
        (
            "component count",
            stored.GetNumberOfComponentsPerPixel(),
            expected.components,
        ),
    )
    for label, actual, wanted in checks:
        if actual != wanted:
            raise ValueError(
                "serialized volume %s changed during writing (%r != %r)"
                % (label, actual, wanted)
            )
    for geometry_label, geometry_actual, geometry_wanted in (
        ("spacing", stored.GetSpacing(), expected.spacing),
        ("origin", stored.GetOrigin(), expected.origin),
        ("direction", stored.GetDirection(), expected.direction),
    ):
        delta = np.abs(
            np.asarray(geometry_actual, dtype=np.float64)
            - np.asarray(geometry_wanted, dtype=np.float64)
        )
        tolerance = _geometry_tolerance(geometry_wanted, output)
        if np.any(delta > tolerance):
            raise ValueError(
                "serialized volume %s changed during writing (maximum delta %.6g, "
                "allowed %.6g)"
                % (geometry_label, float(delta.max()), float(tolerance.max()))
            )
    if _voxel_digest(stored) != expected.digest:
        raise ValueError("serialized volume voxel values changed during writing")


def _verify_compression(path: str, output: VolumeOutput) -> None:
    with open(path, "rb") as handle:
        header = handle.read(4096)
    lowered = header.lower()
    if output.extension == ".nii" and header.startswith(b"\x1f\x8b"):
        raise ValueError("serialized .nii volume was unexpectedly compressed")
    if output.extension == ".nii.gz" and not header.startswith(b"\x1f\x8b"):
        raise ValueError("serialized .nii.gz volume is not gzip-compressed")
    if output.extension == ".nrrd" and b"encoding: gzip" not in lowered:
        raise ValueError("serialized .nrrd volume is not losslessly compressed")
    if output.extension == ".mha" and b"compresseddata = true" not in lowered:
        raise ValueError("serialized .mha volume is not losslessly compressed")


def write_verified_volume(
    image: sitk.Image,
    path: str,
    *,
    release_source: Callable[[], None],
) -> VolumeWriteResult:
    """Verify and atomically publish one exact single-file volume.

    ``release_source`` clears the caller's owning reference after serialization
    and before read-back, so verification never requires two complete pixel
    buffers.
    """
    output = volume_output(path)
    expected = _image_contract(image)
    source_metadata = _metadata(image)
    destination = os.path.abspath(path)
    parent = os.path.dirname(destination)
    os.makedirs(parent, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=".%s." % os.path.basename(destination),
        # ITK dispatch is case-sensitive for NIfTI and treats uppercase .MHA as
        # detached MetaImage. A normalized final compound suffix guarantees the
        # requested single-file writer while the final rename preserves the
        # caller's spelling.
        suffix=output.extension,
        dir=parent,
    )
    os.close(descriptor)
    try:
        sitk.WriteImage(
            image,
            temporary,
            useCompression=output.compressed,
            compressionLevel=9 if output.compressed else -1,
        )
        release_source()
        image = sitk.Image()

        stored = sitk.ReadImage(temporary)
        _verify_contract(stored, expected, output)
        _verify_compression(temporary, output)
        preserved_metadata = sum(
            stored.HasMetaDataKey(key) and stored.GetMetaData(key) == value
            for key, value in source_metadata.items()
        )
        stored = sitk.Image()
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return VolumeWriteResult(
        output=output,
        source_metadata=len(source_metadata),
        preserved_metadata=preserved_metadata,
    )


def convert(
    candidate: VolumeCandidate,
    output_path: str,
    *,
    strip_metadata: bool = False,
    allow_large_volume: bool = False,
    log: Callable[[str], None] | None = None,
    warn: Callable[[str], None] | None = None,
) -> ConversionResult:
    """Preserve one selected image while changing only its storage format."""
    output = volume_output(output_path)
    if os.path.isdir(output_path):
        raise ValueError("volume output path is a directory: %s" % output_path)
    started = time.time()

    def say(message: str) -> None:
        if log:
            log(message)

    def step(message: str, function):
        say("%s ..." % message)
        before = time.time()
        value = function()
        say("  %-34s %6.1fs" % (message, time.time() - before))
        return value

    warnings: list[str] = []

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    loaded = step(
        "load volume",
        lambda: load(
            candidate,
            allow_large_volume=allow_large_volume,
            preserve_metadata=not strip_metadata,
        ),
    )
    provenance = provenance_for(loaded)
    for message in loading_warnings_for(loaded):
        add_warning(message)

    dimensions = tuple(int(value) for value in loaded.image.GetSize())
    pixel_type = loaded.image.GetPixelIDTypeAsString()
    components = loaded.image.GetNumberOfComponentsPerPixel()
    spacing = tuple(float(value) for value in loaded.image.GetSpacing())
    origin = tuple(float(value) for value in loaded.image.GetOrigin())
    direction = tuple(float(value) for value in loaded.image.GetDirection())

    written = step(
        "write and verify volume",
        lambda: write_verified_volume(
            loaded.image,
            output_path,
            release_source=lambda: setattr(loaded, "image", sitk.Image()),
        ),
    )
    if (
        not strip_metadata
        and written.preserved_metadata < written.source_metadata
    ):
        add_warning(
            "%s preserved %d of %d source metadata entries; metadata not "
            "representable by the destination format was omitted"
            % (
                output.format,
                written.preserved_metadata,
                written.source_metadata,
            )
        )

    return ConversionResult(
        output_path=output_path,
        format=output.format,
        compression=output.compression,
        dimensions=dimensions,
        pixel_type=pixel_type,
        components=components,
        spacing=spacing,
        origin=origin,
        direction=direction,
        seconds=time.time() - started,
        warnings=warnings,
        metadata_policy="stripped" if strip_metadata else "preserved-best-effort",
        provenance=provenance,
    )
