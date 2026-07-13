"""DICOM series discovery, ordering and ranking.

Three things routinely go wrong when reading a DICOM directory:

1. **Filenames carry no ordering guarantee.** Instances are frequently written in
   acquisition or arbitrary order. Sorting by filename produces a volume with
   shuffled slices that still looks plausible in a thumbnail. We sort by
   ``ImagePositionPatient`` projected onto the slice normal, which is the only
   ordering the standard actually defines.

2. **One directory holds many series.** A single study typically contains scout
   images, several reconstruction kernels, and multi-planar reformats. Picking
   "the DICOM files in this folder" silently mixes them.

3. **One SeriesInstanceUID may hold many orientations.** Nothing in DICOM forbids
   it, and reformat/secondary-capture series do it routinely. A real study seen
   here had 64 axial frames and one perpendicular frame under a single UID.
   Taking the slice normal from whichever instance the filesystem yields first
   then projects every position onto the wrong axis: slice spacing collapses from
   2.0 mm to 3.9e-07 mm, ordering is scrambled, and nothing raises. Instances are
   therefore grouped by ``(SeriesInstanceUID, orientation)``.

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

#: Minimum slices for a stack to be worth reconstructing.
MIN_SLICES = 5

#: No clinical scanner reconstructs slices thinner than this. A value below it
#: means the geometry was computed wrong, not that the scan is very fine.
MIN_SLICE_SPACING_MM = 0.01

#: Instances whose orientations agree to within this angle belong to the same
#: stack. Exact equality is too strict: oblique reformats carry float jitter in
#: ImageOrientationPatient.
ORIENTATION_TOLERANCE_DEG = 2.0


@dataclass
class Series:
    uid: str
    modality: str
    description: str
    series_number: int | None
    #: Unique 1-based row identifier within one deterministic discovery result.
    #: It is intentionally local to that result, unlike ``uid``.
    id: int = 0
    files: list[str] = field(default_factory=list)
    rows: int | None = None
    columns: int | None = None
    pixel_spacing: tuple[float, float] | None = None
    slice_thickness: float | None = None
    kernel: str | None = None
    image_type: tuple[str, ...] = ()
    rescale_type: str | None = None
    rescale_slope: float | None = None
    rescale_intercept: float | None = None
    multi_energy_ct_acquisition: str | None = None
    hu_calibration_consistent: bool = True
    is_localizer: bool = False
    #: Unit normal of the slice plane, in patient coordinates.
    normal: tuple[float, float, float] | None = None
    #: Which orientation group of its SeriesInstanceUID this is, 1-based, and how
    #: many groups that UID was split into. ``n_parts > 1`` means the UID mixed
    #: orientations.
    part: int = 1
    n_parts: int = 1
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
    def has_calibrated_hu(self) -> bool:
        """Whether every discovered CT instance provides sufficient HU evidence."""
        if self.modality.upper() != "CT" or not self.hu_calibration_consistent:
            return False
        if self.rescale_slope is None or self.rescale_intercept is None:
            return False
        if not math.isfinite(self.rescale_slope) or not math.isfinite(self.rescale_intercept):
            return False
        if self.rescale_type is not None:
            return self.rescale_type.upper() == "HU"
        if self.multi_energy_ct_acquisition == "YES":
            return False
        return bool(self.image_type and self.image_type[0] == "ORIGINAL" and not self.is_localizer)

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
    def unusable_reason(self) -> str | None:
        """Why this stack cannot be reconstructed, or None if it can.

        The geometry checks are not paranoia. A stack whose slice spacing is a
        fraction of a micron, or whose spacing varies by more than half its own
        median, cannot be placed on a regular grid; SimpleITK will resample it
        anyway and hand back a confidently wrong volume.
        """
        if self.is_localizer:
            return "localizer"
        if self.modality in _NON_IMAGE_MODALITIES:
            return "not an image series"
        if self.n_slices < MIN_SLICES:
            return "only %d slice(s)" % self.n_slices

        spacing = self.slice_spacing
        if spacing is None or not math.isfinite(spacing):
            return "no slice spacing"
        if spacing < MIN_SLICE_SPACING_MM:
            return "implausible slice spacing (%.2g mm)" % spacing
        if self.spacing_spread_mm > max(0.1, 0.5 * spacing):
            return "irregular spacing (spread %.2f mm)" % self.spacing_spread_mm
        return None

    @property
    def usable(self) -> bool:
        return self.unusable_reason is None

    def label(self) -> str:
        parts = [self.modality, self.description or "(no description)"]
        return " ".join(p for p in parts if p)


def _slice_normal(orientation: Iterable[float]) -> np.ndarray:
    o = np.asarray(list(orientation), dtype=float)
    return np.cross(o[0:3], o[3:6])


_DISCOVERY_TAGS = [
    "SeriesInstanceUID",
    "ImageOrientationPatient",
    "PixelSpacing",
    "ImageType",
    "RescaleType",
    "RescaleSlope",
    "RescaleIntercept",
    "MultienergyCTAcquisition",
    "Modality",
    "SeriesDescription",
    "SeriesNumber",
    "Rows",
    "Columns",
    "SliceThickness",
    "ConvolutionKernel",
    "ImagePositionPatient",
    "InstanceNumber",
]


def _read_header(path: str):
    try:
        return pydicom.dcmread(
            path,
            stop_before_pixels=True,
            force=False,
            specific_tags=_DISCOVERY_TAGS,
        )
    except Exception:
        return None


def _orientation_matches(a: np.ndarray, b: np.ndarray, tol_deg: float) -> bool:
    """Same slice plane *and* same in-plane axes, within an angular tolerance."""
    cos_lim = math.cos(math.radians(tol_deg))
    for i in (0, 1):  # row direction, column direction
        u, v = a[3 * i:3 * i + 3], b[3 * i:3 * i + 3]
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        if nu == 0 or nv == 0:
            return False
        if float(np.dot(u, v) / (nu * nv)) < cos_lim:
            return False
    return True


def discover_files(paths: Iterable[str]) -> list[Series]:
    """Group candidate files by (SeriesInstanceUID, orientation).

    Grouping on the UID alone is not enough: a single UID may carry several
    orientations, and the slice normal taken from an arbitrary member then
    scrambles the ordering of all the others. Paths are sorted so discovery does
    not depend on filesystem iteration order.
    """
    groups: dict[tuple[str, int], Series] = {}
    positions: dict[tuple[str, int], list[tuple[float, str]]] = {}
    normals: dict[tuple[str, int], np.ndarray] = {}
    #: Representative orientation vectors per UID, in first-seen order.
    orientations: dict[str, list[np.ndarray]] = {}

    for path in sorted(paths):
        ds = _read_header(path)
        if ds is None or not hasattr(ds, "SeriesInstanceUID"):
            continue
        uid = str(ds.SeriesInstanceUID)

        iop = np.asarray(
            getattr(ds, "ImageOrientationPatient", [1, 0, 0, 0, 1, 0]), dtype=float
        )
        reps = orientations.setdefault(uid, [])
        for index, rep in enumerate(reps):
            if _orientation_matches(iop, rep, ORIENTATION_TOLERANCE_DEG):
                break
        else:
            reps.append(iop)
            index = len(reps) - 1
        key = (uid, index)

        if key not in groups:
            ps = getattr(ds, "PixelSpacing", None)
            image_type = [str(x).upper() for x in getattr(ds, "ImageType", [])]
            groups[key] = Series(
                uid=uid,
                modality=str(getattr(ds, "Modality", "?")),
                description=str(getattr(ds, "SeriesDescription", "")).strip(),
                series_number=_as_int(getattr(ds, "SeriesNumber", None)),
                rows=_as_int(getattr(ds, "Rows", None)),
                columns=_as_int(getattr(ds, "Columns", None)),
                pixel_spacing=(float(ps[0]), float(ps[1])) if ps else None,
                slice_thickness=_as_float(getattr(ds, "SliceThickness", None)),
                kernel=_kernel_of(ds),
                image_type=tuple(image_type),
                rescale_type=_normalised_text(getattr(ds, "RescaleType", None)),
                rescale_slope=_as_float(getattr(ds, "RescaleSlope", None)),
                rescale_intercept=_as_float(getattr(ds, "RescaleIntercept", None)),
                multi_energy_ct_acquisition=_normalised_text(
                    getattr(ds, "MultienergyCTAcquisition", None)
                ),
                is_localizer="LOCALIZER" in image_type,
            )
            normals[key] = _slice_normal(iop)
            groups[key].normal = tuple(float(v) for v in normals[key])  # type: ignore[assignment]
            positions[key] = []

        series = groups[key]
        if series.files:
            current = (
                tuple(str(x).upper() for x in getattr(ds, "ImageType", [])),
                _normalised_text(getattr(ds, "RescaleType", None)),
                _as_float(getattr(ds, "RescaleSlope", None)),
                _as_float(getattr(ds, "RescaleIntercept", None)),
                _normalised_text(getattr(ds, "MultienergyCTAcquisition", None)),
            )
            expected = (
                series.image_type,
                series.rescale_type,
                series.rescale_slope,
                series.rescale_intercept,
                series.multi_energy_ct_acquisition,
            )
            if current != expected:
                series.hu_calibration_consistent = False
        series.files.append(path)

        ipp = getattr(ds, "ImagePositionPatient", None)
        if ipp is not None:
            depth = float(np.dot(np.asarray(ipp, dtype=float), normals[key]))
        else:
            # No position: fall back to InstanceNumber, which at least beats
            # filename order.
            depth = float(_as_int(getattr(ds, "InstanceNumber", 0)) or 0)
        positions[key].append((depth, path))

    for key, series in groups.items():
        _finalise(series, positions[key])

    # Number the orientation groups of a split UID, biggest stack first.
    by_uid: dict[str, list[Series]] = {}
    for (uid, _index), series in groups.items():
        by_uid.setdefault(uid, []).append(series)
    for uid, members in by_uid.items():
        if len(members) == 1:
            continue
        members.sort(key=lambda s: (-s.n_slices, s.files[0] if s.files else ""))
        for part, series in enumerate(members, start=1):
            series.part = part
            series.n_parts = len(members)

    out = list(groups.values())
    out.sort(
        key=lambda s: (
            s.series_number if s.series_number is not None else 1 << 30,
            s.uid,
            s.part,
        )
    )
    for row_id, series in enumerate(out, start=1):
        series.id = row_id
    return out


def discover(root: str) -> list[Series]:
    """Walk ``root`` and discover DICOM stacks."""
    paths: list[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        paths.extend(os.path.join(dirpath, name) for name in sorted(filenames))
    return discover_files(paths)


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


def _normalised_text(value) -> str | None:
    if value is None:
        return None
    text = str(value).strip().upper()
    return text or None


def _significant(x: float, digits: int = 3) -> float:
    """Round to ``digits`` significant figures.

    Buckets voxel volumes before they are compared. Reformats of one acquisition
    differ in measured slice spacing by rounding noise (0.80000 mm vs 0.79977 mm
    in a real study), and comparing raw floats lets that noise decide which
    series a user gets.
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
        vv = s.voxel_volume_mm3
        # A degenerate voxel volume must sort last, never first. Before geometry
        # was validated, a mis-derived 3.9e-07 mm slice spacing made a reformat
        # look like the finest series in the study and win this comparison.
        if vv is None or not math.isfinite(vv) or vv <= 0:
            vv = 1e9
        return (
            0 if s.usable else 1,
            _significant(vv),
            0 if s.plane == "axial" else 1,
            -s.n_slices,
        )

    return sorted(series, key=key)
