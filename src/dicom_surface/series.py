"""DICOM series discovery, ordering and ranking.

Two things routinely go wrong when reading a DICOM directory:

1. **Filenames carry no ordering guarantee.** Instances are frequently written in
   acquisition or arbitrary order. Sorting by filename produces a volume with
   shuffled slices that still looks plausible in a thumbnail. We sort by
   ``ImagePositionPatient`` projected onto the slice normal, which is the only
   ordering the standard actually defines.

2. **One directory holds many series.** A single study typically contains scout
   images, several reconstruction kernels, and multi-planar reformats. Picking
   "the DICOM files in this folder" silently mixes them.

This module never reads patient identifiers. Only geometry, modality and
acquisition parameters are extracted.
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np
import pydicom

#: Siemens/GE/Philips sharp ("bone", "edge-enhancing") reconstruction kernels.
#: These amplify high-frequency noise by design; thresholding them at the usual
#: bone value yields a surface covered in spurious spikes.
SHARP_KERNEL_RE = re.compile(r"(?:^|[^0-9])(?:[BHUY]r?|BONE|EDGE|LUNG)\s*_?([6-9]\d)", re.I)

#: Series we should never try to reconstruct.
_NON_IMAGE_MODALITIES = {"SR", "PR", "KO", "SEG", "RTSTRUCT", "RTPLAN", "RTDOSE", "DOC"}


@dataclass
class Series:
    uid: str
    modality: str
    description: str
    series_number: int | None
    files: list[str] = field(default_factory=list)
    rows: int | None = None
    columns: int | None = None
    pixel_spacing: tuple[float, float] | None = None
    slice_thickness: float | None = None
    kernel: str | None = None
    is_localizer: bool = False
    #: Unit normal of the slice plane, in patient coordinates.
    normal: tuple[float, float, float] | None = None
    #: Populated by :meth:`finalise`.
    slice_spacing: float | None = None
    spacing_uniform: bool = True
    spacing_spread_mm: float = 0.0
    gantry_tilted: bool = False

    # ---------------------------------------------------------------- helpers
    @property
    def n_slices(self) -> int:
        return len(self.files)

    @property
    def voxel_volume_mm3(self) -> float | None:
        if not self.pixel_spacing or not self.slice_spacing:
            return None
        return self.pixel_spacing[0] * self.pixel_spacing[1] * self.slice_spacing

    @property
    def sharp_kernel(self) -> bool:
        return bool(self.kernel and SHARP_KERNEL_RE.search(self.kernel))

    @property
    def plane(self) -> str:
        """Slice plane in patient coordinates: axial, coronal, sagittal, oblique."""
        if self.normal is None:
            return "unknown"
        n = np.abs(np.asarray(self.normal, dtype=float))
        axis = int(np.argmax(n))
        if n[axis] < 0.95:  # more than ~18 degrees off a cardinal axis
            return "oblique"
        return ("sagittal", "coronal", "axial")[axis]

    @property
    def usable(self) -> bool:
        return (
            self.n_slices >= 5
            and not self.is_localizer
            and self.modality not in _NON_IMAGE_MODALITIES
        )

    def label(self) -> str:
        parts = [self.modality, self.description or "(no description)"]
        return " ".join(p for p in parts if p)


def _slice_normal(orientation: Iterable[float]) -> np.ndarray:
    o = np.asarray(list(orientation), dtype=float)
    return np.cross(o[0:3], o[3:6])


def _read_header(path: str):
    try:
        return pydicom.dcmread(path, stop_before_pixels=True, force=False)
    except Exception:
        return None


def discover(root: str) -> list[Series]:
    """Walk ``root`` and group every readable DICOM instance by SeriesInstanceUID."""
    by_uid: dict[str, Series] = {}
    positions: dict[str, list[tuple[float, str]]] = {}
    normals: dict[str, np.ndarray] = {}

    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            path = os.path.join(dirpath, fn)
            ds = _read_header(path)
            if ds is None or not hasattr(ds, "SeriesInstanceUID"):
                continue
            uid = str(ds.SeriesInstanceUID)

            if uid not in by_uid:
                ps = getattr(ds, "PixelSpacing", None)
                image_type = [str(x).upper() for x in getattr(ds, "ImageType", [])]
                by_uid[uid] = Series(
                    uid=uid,
                    modality=str(getattr(ds, "Modality", "?")),
                    description=str(getattr(ds, "SeriesDescription", "")).strip(),
                    series_number=_as_int(getattr(ds, "SeriesNumber", None)),
                    rows=_as_int(getattr(ds, "Rows", None)),
                    columns=_as_int(getattr(ds, "Columns", None)),
                    pixel_spacing=(float(ps[0]), float(ps[1])) if ps else None,
                    slice_thickness=_as_float(getattr(ds, "SliceThickness", None)),
                    kernel=_kernel_of(ds),
                    is_localizer="LOCALIZER" in image_type,
                )
                normals[uid] = _slice_normal(
                    getattr(ds, "ImageOrientationPatient", [1, 0, 0, 0, 1, 0])
                )
                by_uid[uid].normal = tuple(float(v) for v in normals[uid])  # type: ignore[assignment]
                positions[uid] = []

            by_uid[uid].files.append(path)

            ipp = getattr(ds, "ImagePositionPatient", None)
            if ipp is not None:
                depth = float(np.dot(np.asarray(ipp, dtype=float), normals[uid]))
            else:
                # No position: fall back to InstanceNumber, which at least beats
                # filename order.
                depth = float(_as_int(getattr(ds, "InstanceNumber", 0)) or 0)
            positions[uid].append((depth, path))

    out = []
    for uid, series in by_uid.items():
        series.finalise = None  # type: ignore[attr-defined]
        _finalise(series, positions[uid])
        out.append(series)

    out.sort(key=lambda s: (s.series_number if s.series_number is not None else 1 << 30, s.uid))
    return out


def _finalise(series: Series, positions: list[tuple[float, str]]) -> None:
    """Order slices by physical position and characterise the spacing."""
    positions.sort(key=lambda t: t[0])
    series.files = [p for _d, p in positions]

    depths = np.asarray([d for d, _p in positions], dtype=float)
    if depths.size >= 2:
        diffs = np.diff(depths)
        diffs = diffs[np.abs(diffs) > 1e-9]
        if diffs.size:
            series.slice_spacing = float(np.median(np.abs(diffs)))
            series.spacing_spread_mm = float(np.max(np.abs(diffs)) - np.min(np.abs(diffs)))
            # 1 micron of jitter is float noise; 10 microns is a real gap.
            series.spacing_uniform = series.spacing_spread_mm < 0.01
    if series.slice_spacing is None:
        series.slice_spacing = series.slice_thickness


def _kernel_of(ds) -> str | None:
    k = getattr(ds, "ConvolutionKernel", None)
    if k is None:
        return None
    if isinstance(k, (list, tuple, pydicom.multival.MultiValue)):
        return "\\".join(str(x) for x in k)
    return str(k)


def _as_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _significant(x: float, digits: int = 3) -> float:
    """Round to ``digits`` significant figures.

    Used to bucket voxel volumes before comparing them. Reformats of one
    acquisition differ in measured slice spacing by rounding noise (0.80000 mm
    vs 0.79977 mm in a real study), and comparing raw floats lets that noise
    decide which series a user gets.
    """
    if x == 0 or not math.isfinite(x):
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (digits - 1))


def rank(series: Iterable[Series]) -> list[Series]:
    """Best-first ordering for automatic selection.

    Prefers, in order:

    1. usable image series;
    2. smaller voxels, compared to 3 significant figures;
    3. the axial plane. A study often carries one acquisition reformatted into
       three planes. The axial stack is normally the acquired one; coronal and
       sagittal reformats have been interpolated a second time.
    4. more slices. This is deliberately the *last* tiebreak: a sagittal reformat
       of the same volume can have more slices while carrying no more information.

    Sharp kernels are *not* penalised: they are the higher-resolution data. The
    right response is a higher threshold, which the CLI warns about.
    """

    def key(s: Series):
        vv = s.voxel_volume_mm3 or 1e9
        return (
            0 if s.usable else 1,
            _significant(vv),
            0 if s.plane == "axial" else 1,
            -s.n_slices,
        )

    return sorted(series, key=key)


def select(series: list[Series], wanted: str | None) -> Series:
    """Resolve ``wanted`` (a UID, a series number, or a description substring)."""
    usable = [s for s in series if s.usable]
    if not usable:
        raise ValueError("no usable image series found")

    if wanted is None:
        return rank(usable)[0]

    for s in series:
        if s.uid == wanted:
            return s
    if wanted.isdigit():
        n = int(wanted)
        hits = [s for s in series if s.series_number == n]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            raise ValueError("series number %d is ambiguous (%d matches)" % (n, len(hits)))
    lowered = wanted.lower()
    hits = [s for s in series if lowered in s.description.lower()]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise ValueError(
            "description %r matches %d series: %s"
            % (wanted, len(hits), ", ".join(s.description for s in hits))
        )
    raise ValueError("no series matches %r" % wanted)
