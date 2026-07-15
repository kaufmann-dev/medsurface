"""Load a selected catalog volume into a SimpleITK image."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import SimpleITK as sitk

from .catalog import DicomSource, FileSource, VolumeCandidate, validate_image
from .defaults import MAX_VOXELS, MIN_VOLUME_AXIS_VOXELS


@dataclass
class Volume:
    image: sitk.Image
    candidate: VolumeCandidate

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


def load(candidate: VolumeCandidate, *, allow_large_volume: bool = False) -> Volume:
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
    try:
        if isinstance(candidate.source, DicomSource):
            series = candidate.source.series
            if series.n_slices < 2:
                raise ValueError(
                    "series %r has %d slice(s)" % (series.description, series.n_slices)
                )
            reader = sitk.ImageSeriesReader()
            reader.SetFileNames(series.files)  # already ordered by physical position
            image = reader.Execute()
        else:
            image = sitk.ReadImage(str(candidate.source.path))
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
    return Volume(image=image, candidate=candidate)


def has_calibrated_hu(candidate: VolumeCandidate) -> bool:
    """Whether the loader can verify that values are Hounsfield units."""
    return bool(candidate.dicom and candidate.dicom.has_calibrated_hu)


def warnings_for(volume: Volume) -> list[str]:
    """Non-fatal data-quality observations worth surfacing to the user."""
    out: list[str] = []
    candidate = volume.candidate
    series = candidate.dicom

    if series is not None and not series.spacing_uniform:
        out.append(
            "slice spacing is not uniform (spread %.3f mm); the volume will be "
            "resampled onto a regular grid and geometry may shift slightly"
            % series.spacing_spread_mm
        )

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
