"""Load an ordered DICOM series into a SimpleITK image.

SimpleITK works in DICOM patient space (LPS), which is also the convention STL
readers expect, so no coordinate flip is needed anywhere downstream. The reader
applies RescaleSlope/RescaleIntercept, meaning CT pixel values arrive as
Hounsfield units.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import SimpleITK as sitk

from .series import Series


@dataclass
class Volume:
    image: sitk.Image
    series: Series

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


def load(series: Series) -> Volume:
    if series.n_slices < 2:
        raise ValueError("series %r has %d slice(s)" % (series.description, series.n_slices))
    reader = sitk.ImageSeriesReader()
    reader.SetFileNames(series.files)  # already ordered by physical position
    image = reader.Execute()
    if image.GetDimension() != 3:
        raise ValueError("expected a 3D volume, got %dD" % image.GetDimension())
    return Volume(image=image, series=series)


def has_calibrated_hu(series: Series) -> bool:
    """Whether pixel values can be interpreted as Hounsfield units."""
    return series.modality == "CT"


def warnings_for(volume: Volume) -> list[str]:
    """Non-fatal data-quality observations worth surfacing to the user."""
    out: list[str] = []
    s = volume.series

    if not s.spacing_uniform:
        out.append(
            "slice spacing is not uniform (spread %.3f mm); the volume will be "
            "resampled onto a regular grid and geometry may shift slightly"
            % s.spacing_spread_mm
        )

    sx, sy, sz = volume.spacing
    aniso = max(sx, sy, sz) / min(sx, sy, sz)
    if aniso > 2.0:
        out.append(
            "voxels are strongly anisotropic (%.2f x %.2f x %.2f mm, ratio %.1f:1); "
            "expect stair-stepping along the thick axis" % (sx, sy, sz, aniso)
        )

    if s.sharp_kernel:
        out.append(
            "reconstruction kernel %r is a sharp/edge-enhancing kernel: it amplifies "
            "noise, so a low bone threshold will produce a spiky surface. Prefer a "
            "higher threshold (~300 HU rather than ~200 HU), or a smoother kernel "
            "series if the study has one." % s.kernel
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
