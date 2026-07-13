"""Discover format-neutral 3D volume candidates.

A directory is a catalog: DICOM instances are grouped into physical stacks and
every supported self-describing image file contributes one candidate.  A direct
file is simply a one-item catalog.  Candidate IDs are deliberately local to one
discovery result and are the only public selectors.
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias

import numpy as np
import SimpleITK as sitk

from . import series as series_mod
from .series import Series

SUPPORTED_EXTENSIONS = (".nii.gz", ".nii", ".nrrd", ".nhdr", ".mha", ".mhd")


@dataclass(frozen=True)
class DicomSource:
    catalog_path: Path
    series: Series


@dataclass(frozen=True)
class FileSource:
    path: Path
    format: str
    payload_paths: tuple[Path, ...] = ()


VolumeSource: TypeAlias = DicomSource | FileSource


@dataclass
class VolumeCandidate:
    id: int
    source: VolumeSource
    format: str
    source_name: str
    modality: str | None
    description: str | None
    size: tuple[int, ...] | None
    spacing: tuple[float, ...] | None
    direction: tuple[float, ...] | None
    origin: tuple[float, ...] | None
    pixel_type: str | None
    components: int | None
    plane: str | None
    unusable_reason: str | None = None

    @property
    def usable(self) -> bool:
        return self.unusable_reason is None

    @property
    def dicom(self) -> Series | None:
        return self.source.series if isinstance(self.source, DicomSource) else None

    @property
    def slices(self) -> int | None:
        if self.dicom is not None:
            return self.dicom.n_slices
        if self.size is not None and len(self.size) == 3:
            return self.size[2]
        return None

    @property
    def series_number(self) -> int | None:
        return self.dicom.series_number if self.dicom is not None else None

    @property
    def kernel_values(self) -> tuple[str, ...]:
        return self.dicom.kernel_values if self.dicom is not None else ()

    @property
    def sharp_kernel(self) -> bool:
        return bool(self.dicom and self.dicom.sharp_kernel)


def _extension(path: Path) -> str | None:
    name = path.name.casefold()
    return next((ext for ext in SUPPORTED_EXTENSIONS if name.endswith(ext)), None)


def _format_for(extension: str) -> str:
    if extension in (".nii", ".nii.gz"):
        return "NIfTI"
    if extension in (".nrrd", ".nhdr"):
        return "NRRD"
    return "MetaImage"


def _header_value(lines: list[str], pattern: re.Pattern[str]) -> tuple[int, str] | None:
    for index, line in enumerate(lines):
        match = pattern.match(line.strip())
        if match:
            return index, match.group(1).strip()
    return None


def _listed_payloads(lines: list[str], start: int) -> list[str]:
    names = []
    for line in lines[start:]:
        value = line.strip()
        if not value:
            break
        if not value.startswith("#"):
            names.append(value.strip('"'))
    return names


def _expanded_payloads(value: str) -> list[str]:
    """Expand the integer filename sequence supported by detached image headers."""
    parts = value.split()
    if len(parts) < 4 or "%" not in parts[0]:
        return [value.strip('"')]
    try:
        first, last, step = map(int, parts[1:4])
    except ValueError:
        return [value.strip('"')]
    if step == 0 or (last - first) * step < 0:
        return [value.strip('"')]
    count = abs((last - first) // step) + 1
    if count > 100_000:
        return [value.strip('"')]
    try:
        return [parts[0] % index for index in range(first, last + (1 if step > 0 else -1), step)]
    except (TypeError, ValueError):
        return [value.strip('"')]


def _detached_payloads(path: Path, extension: str) -> tuple[Path, ...]:
    if extension not in (".mhd", ".nhdr"):
        return ()
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if extension == ".mhd":
        found = _header_value(lines, re.compile(r"elementdatafile\s*=\s*(.+)", re.I))
    else:
        found = _header_value(lines, re.compile(r"data\s*file\s*:\s*(.+)", re.I))
    if found is None:
        return ()
    index, value = found
    upper = value.upper()
    if upper == "LOCAL":
        return ()
    if upper.startswith("LIST"):
        names = _listed_payloads(lines, index + 1)
    else:
        names = _expanded_payloads(value)
    return tuple(
        item if item.is_absolute() else path.parent / item
        for item in (Path(name) for name in names)
    )


def _payload_reason(payloads: tuple[Path, ...]) -> str | None:
    for payload in payloads:
        if not payload.is_file():
            return "referenced payload is missing: %s" % payload.name
        if not os.access(payload, os.R_OK):
            return "referenced payload is not readable: %s" % payload.name
    return None


def _metadata(reader: sitk.ImageFileReader, *names: str) -> str | None:
    wanted = {name.casefold() for name in names}
    for key in reader.GetMetaDataKeys():
        if key.casefold() not in wanted:
            continue
        value = reader.GetMetaData(key).strip()
        if value and value.casefold() not in {"unknown", "met_mod_unknown", "none"}:
            return value
    return None


def _plane(direction: tuple[float, ...] | None) -> str | None:
    if direction is None or len(direction) != 9:
        return None
    matrix = np.asarray(direction, dtype=float).reshape(3, 3)
    normal = np.abs(matrix[:, 2])
    axis = int(np.argmax(normal))
    if not math.isfinite(float(normal[axis])):
        return None
    if normal[axis] < 0.95:
        return "oblique"
    return ("sagittal", "coronal", "axial")[axis]


def _geometry_reason(
    dimension: int,
    size: tuple[int, ...],
    spacing: tuple[float, ...],
    direction: tuple[float, ...],
    origin: tuple[float, ...],
    pixel_type: str,
    components: int,
) -> str | None:
    if dimension != 3:
        return "expected a 3D volume, got %dD" % dimension
    if components != 1:
        return "expected one scalar component, got %d" % components
    if "complex" in pixel_type.casefold():
        return "complex-valued pixels are unsupported"
    if len(size) != 3 or any(value < 2 for value in size):
        return "each volume axis must contain at least 2 voxels"
    if len(spacing) != 3 or any(not math.isfinite(value) or value <= 0 for value in spacing):
        return "spacing must contain three finite positive values"
    if len(origin) != 3 or any(not math.isfinite(value) for value in origin):
        return "origin must contain three finite values"
    if len(direction) != 9 or not np.all(np.isfinite(direction)):
        return "direction must contain nine finite values"
    if abs(float(np.linalg.det(np.asarray(direction).reshape(3, 3)))) < 1e-8:
        return "direction matrix is singular"
    return None


def validate_image(image: sitk.Image) -> str | None:
    """Return why a loaded image is unusable, or ``None``."""
    return _geometry_reason(
        image.GetDimension(),
        tuple(int(value) for value in image.GetSize()),
        tuple(float(value) for value in image.GetSpacing()),
        tuple(float(value) for value in image.GetDirection()),
        tuple(float(value) for value in image.GetOrigin()),
        image.GetPixelIDTypeAsString(),
        image.GetNumberOfComponentsPerPixel(),
    )


def _file_candidate(path: Path, root: Path) -> VolumeCandidate:
    extension = _extension(path)
    if extension is None:
        raise ValueError("unsupported input extension %r" % path.name)
    format_name = _format_for(extension)
    payload_problem: str | None
    try:
        payloads = _detached_payloads(path, extension)
    except OSError as exc:
        payloads = ()
        payload_problem = "cannot inspect detached header: %s" % exc
    else:
        payload_problem = _payload_reason(payloads)
    source = FileSource(path=path, format=format_name, payload_paths=payloads)
    if payload_problem is not None:
        return VolumeCandidate(
            id=0,
            source=source,
            format=format_name,
            source_name=_relative_name(path, root),
            modality=None,
            description=None,
            size=None,
            spacing=None,
            direction=None,
            origin=None,
            pixel_type=None,
            components=None,
            plane=None,
            unusable_reason=payload_problem,
        )
    try:
        reader = sitk.ImageFileReader()
        reader.SetFileName(str(path))
        reader.ReadImageInformation()
        dimension = reader.GetDimension()
        size = tuple(int(value) for value in reader.GetSize())
        spacing = tuple(float(value) for value in reader.GetSpacing())
        direction = tuple(float(value) for value in reader.GetDirection())
        origin = tuple(float(value) for value in reader.GetOrigin())
        pixel_type = sitk.GetPixelIDValueAsString(reader.GetPixelID())
        components = reader.GetNumberOfComponents()
        reason = _geometry_reason(
            dimension, size, spacing, direction, origin, pixel_type, components
        )
        modality = _metadata(reader, "modality")
        description = _metadata(reader, "description", "descrip", "content", "intent_name")
    except RuntimeError as exc:
        return VolumeCandidate(
            id=0,
            source=source,
            format=format_name,
            source_name=_relative_name(path, root),
            modality=None,
            description=None,
            size=None,
            spacing=None,
            direction=None,
            origin=None,
            pixel_type=None,
            components=None,
            plane=None,
            unusable_reason=_reader_error(exc),
        )

    return VolumeCandidate(
        id=0,
        source=source,
        format=format_name,
        source_name=_relative_name(path, root),
        modality=modality,
        description=description,
        size=size,
        spacing=spacing,
        direction=direction,
        origin=origin,
        pixel_type=pixel_type,
        components=components,
        plane=_plane(direction),
        unusable_reason=reason,
    )


def _reader_error(exc: RuntimeError) -> str:
    lines = [line.strip() for line in str(exc).splitlines() if line.strip()]
    detail = lines[-1] if lines else "reader rejected the file"
    return "cannot read image header: %s" % detail


def _relative_name(path: Path, root: Path) -> str:
    if root.is_file():
        return path.name
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def _dicom_candidate(series: Series, root: Path) -> VolumeCandidate:
    size = None
    if series.columns is not None and series.rows is not None:
        size = (series.columns, series.rows, series.n_slices)
    spacing = None
    if series.pixel_spacing and series.slice_spacing:
        spacing = (
            float(series.pixel_spacing[1]),
            float(series.pixel_spacing[0]),
            float(series.slice_spacing),
        )
    representative = Path(min(series.files)) if series.files else root
    return VolumeCandidate(
        id=0,
        source=DicomSource(catalog_path=root, series=series),
        format="DICOM",
        source_name=_relative_name(representative.parent, root),
        modality=series.modality if series.modality != "?" else None,
        description=series.description or None,
        size=size,
        spacing=spacing,
        direction=None,
        origin=None,
        pixel_type=None,
        components=1,
        plane=series.plane if series.normal is not None else None,
        unusable_reason=series.unusable_reason,
    )


def _paths(root: Path) -> list[Path]:
    if root.is_file():
        return [root]
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        found.extend(Path(dirpath) / filename for filename in sorted(filenames))
    return found


def discover(input_path: str | Path) -> list[VolumeCandidate]:
    """Discover every supported volume at or below ``input_path``."""
    root = Path(input_path)
    paths = _paths(root)
    dicom_series = series_mod.discover_files([str(path) for path in paths])
    dicom_paths = {Path(path) for series in dicom_series for path in series.files}

    candidates = [_dicom_candidate(series, root) for series in dicom_series]
    for path in paths:
        if path in dicom_paths or _extension(path) is None:
            continue
        candidates.append(_file_candidate(path, root))

    def sort_key(candidate: VolumeCandidate):
        native = ""
        if candidate.dicom is not None:
            native = "%s:%04d" % (candidate.dicom.uid, candidate.dicom.part)
        return (candidate.source_name.casefold(), candidate.format, native)

    candidates.sort(key=sort_key)
    for candidate_id, candidate in enumerate(candidates, start=1):
        candidate.id = candidate_id
    return candidates


def recommended(candidates: list[VolumeCandidate]) -> VolumeCandidate | None:
    usable = [candidate for candidate in candidates if candidate.usable]
    if len(usable) == 1:
        return usable[0]
    if usable and all(candidate.dicom is not None for candidate in usable):
        if selection_ambiguity(usable) is not None:
            return None
        ranked = series_mod.rank(candidate.dicom for candidate in usable if candidate.dicom is not None)
        winner = ranked[0]
        return next(candidate for candidate in usable if candidate.dicom is winner)
    return None


def selection_ambiguity(candidates: list[VolumeCandidate]) -> str | None:
    """Explain a DICOM-only catalog that has no safe automatic winner."""
    usable = [candidate for candidate in candidates if candidate.usable]
    if len(usable) < 2 or not all(candidate.dicom is not None for candidate in usable):
        return None
    modalities = sorted(
        {
            candidate.modality.strip().upper()
            if candidate.modality and candidate.modality.strip()
            else "unknown"
            for candidate in usable
        }
    )
    if len(modalities) < 2:
        return None
    return "usable DICOM volumes span multiple modalities: %s" % ", ".join(modalities)


def select(candidates: list[VolumeCandidate], wanted: int | None) -> VolumeCandidate:
    usable = [candidate for candidate in candidates if candidate.usable]
    if not usable:
        if not candidates:
            raise ValueError("no supported volumes found")
        reasons = sorted({candidate.unusable_reason or "unusable" for candidate in candidates})
        raise ValueError("no usable volumes found: %s" % "; ".join(reasons))
    if wanted is not None:
        hits = [candidate for candidate in candidates if candidate.id == wanted]
        if not hits:
            raise ValueError(
                "no volume has ID %d; available IDs: %s"
                % (wanted, ", ".join(str(candidate.id) for candidate in candidates))
            )
        chosen = hits[0]
        if not chosen.usable:
            raise ValueError("volume ID %d is unusable: %s" % (wanted, chosen.unusable_reason))
        return chosen
    default = recommended(candidates)
    if default is not None:
        return default
    ambiguity = selection_ambiguity(candidates)
    if ambiguity is not None:
        raise ValueError(
            "%s; run 'medsurface list INPUT' and pass the displayed volume ID"
            % ambiguity
        )
    raise ValueError(
        "input contains multiple usable volumes; run 'medsurface list INPUT' and pass a volume ID"
    )


def same_source(a: VolumeCandidate, b: VolumeCandidate) -> bool:
    if isinstance(a.source, DicomSource) and isinstance(b.source, DicomSource):
        return a.source.series.uid == b.source.series.uid and a.source.series.part == b.source.series.part
    if isinstance(a.source, FileSource) and isinstance(b.source, FileSource):
        try:
            return os.path.samefile(a.source.path, b.source.path)
        except OSError:
            return a.source.path.resolve() == b.source.path.resolve()
    return False


def source_paths(candidate: VolumeCandidate) -> tuple[Path, ...]:
    """Every file required to load a selected volume."""
    if isinstance(candidate.source, DicomSource):
        return tuple(Path(path) for path in candidate.source.series.files)
    return (candidate.source.path, *candidate.source.payload_paths)
