"""Extract or fuse externally produced discrete segmentation labelmaps."""

from __future__ import annotations

import math
import os
import time
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import SimpleITK as sitk

from . import fusion, labelnames, pipeline, registration, surface
from . import volume as volume_mod
from .catalog import DicomSource, FileSource, VolumeCandidate, same_source
from .defaults import (
    DEFAULT_FUSION_GRID_MM,
    DEFAULT_LABELMAP_MASK_SMOOTH_MM,
    DEFAULT_LABELMAP_POST_SURFACE_SMOOTH_ITERS,
    DEFAULT_LABELMAP_SURFACE_SMOOTH_ITERS,
    SUPPORTED_MESH_EXTENSIONS,
)
from .outputs import volume_output
from .presets import DestepSettings

Logger = Callable[[str], None]
#: Receives machine-readable progress events such as
#: ``{"event": "stage_start", "stage": "marching cubes"}``.
ProgressSink = Callable[[dict[str, Any]], None]
SettingsForLabel = Callable[[int], pipeline.SurfaceSettings]


@dataclass
class LoadedLabelmap:
    volume: volume_mod.Volume
    mask: sitk.Image
    provenance: dict[str, Any]


@dataclass
class Result:
    output_path: str
    triangles: int
    vertices: int
    bounds_mm: tuple[float, ...]
    labelmap_components: int
    surface_components: int
    capped_field_of_view: bool
    seconds: float
    warnings: list[str]
    provenance: dict[str, Any]
    quality: dict[str, Any]


def default_surface_settings() -> pipeline.SurfaceSettings:
    """Return independent finishing defaults for an external segmentation."""
    return pipeline.SurfaceSettings(
        resample_mm=0.0,
        mask_smooth_mm=DEFAULT_LABELMAP_MASK_SMOOTH_MM,
        surface_smooth_iters=DEFAULT_LABELMAP_SURFACE_SMOOTH_ITERS,
        simplify_error_mm=0.25,
        post_surface_smooth_iters=DEFAULT_LABELMAP_POST_SURFACE_SMOOTH_ITERS,
        keep_largest_component=False,
    )


def resolve_surface_settings(
    *,
    resample_mm: float | None = None,
    mask_smooth_mm: float | None = None,
    surface_smooth_iters: int | None = None,
    simplify_error_mm: float | None = None,
    post_surface_smooth_iters: int | None = None,
    keep_largest_component: bool | None = None,
    destep: DestepSettings | None = None,
) -> pipeline.SurfaceSettings:
    base = default_surface_settings()
    settings = replace(
        base,
        resample_mm=base.resample_mm if resample_mm is None else resample_mm,
        mask_smooth_mm=(
            base.mask_smooth_mm if mask_smooth_mm is None else float(mask_smooth_mm)
        ),
        surface_smooth_iters=(
            base.surface_smooth_iters
            if surface_smooth_iters is None
            else surface_smooth_iters
        ),
        simplify_error_mm=(
            base.simplify_error_mm if simplify_error_mm is None else simplify_error_mm
        ),
        post_surface_smooth_iters=(
            base.post_surface_smooth_iters
            if post_surface_smooth_iters is None
            else post_surface_smooth_iters
        ),
        keep_largest_component=(
            base.keep_largest_component
            if keep_largest_component is None
            else keep_largest_component
        ),
        destep=destep,
    )
    pipeline.validate_surface_settings(settings)
    return settings


def surface_provenance(settings: pipeline.SurfaceSettings) -> dict[str, Any]:
    """Record the external labelmap surface controls."""
    return asdict(settings)


def normalize_labels(labels: Sequence[int] | None) -> tuple[int, ...] | None:
    """Validate a label selection: sorted, unique, positive integer IDs."""
    if labels is None:
        return None
    normalized: set[int] = set()
    for value in labels:
        if isinstance(value, bool) or int(value) != value:
            raise ValueError("label IDs must be integers; got %r" % (value,))
        if int(value) <= 0:
            raise ValueError("label IDs must be positive; got %d" % int(value))
        normalized.add(int(value))
    if not normalized:
        raise ValueError("the label selection is empty")
    return tuple(sorted(normalized))


def _validate_labelmap(image: sitk.Image) -> bool:
    """Reject non-discrete label volumes; return whether any voxel is nonzero."""
    values = sitk.GetArrayViewFromImage(image)
    is_float = np.issubdtype(values.dtype, np.floating)
    is_signed = np.issubdtype(values.dtype, np.signedinteger)
    has_foreground = False

    for plane in values:
        if is_float:
            if not np.all(np.isfinite(plane)):
                raise ValueError("labelmap contains non-finite voxel values")
            if np.any(plane != np.trunc(plane)):
                raise ValueError(
                    "labelmap contains fractional voxel values; probability maps are unsupported"
                )
        if (is_float or is_signed) and np.any(plane < 0):
            raise ValueError("labelmap values must be non-negative")
        has_foreground = has_foreground or bool(np.any(plane != 0))
    return has_foreground


def _selected_mask(image: sitk.Image, labels: tuple[int, ...]) -> sitk.Image:
    values = sitk.GetArrayViewFromImage(image)
    selected = np.zeros(values.shape, dtype=np.uint8)
    for index, plane in enumerate(values):
        selected[index] = np.isin(plane, labels)
    mask = sitk.GetImageFromArray(selected)
    mask.CopyInformation(image)
    return mask


def _binary_mask(image: sitk.Image, labels: Sequence[int] | None = None) -> sitk.Image:
    """Validate a discrete label volume and union every nonzero or selected label."""
    selection = normalize_labels(labels)
    if not _validate_labelmap(image):
        raise ValueError("labelmap contains no nonzero foreground voxels")
    if selection is None:
        return sitk.Cast(sitk.NotEqual(image, 0), sitk.sitkUInt8)
    mask = _selected_mask(image, selection)
    if not np.any(sitk.GetArrayViewFromImage(mask)):
        raise ValueError(
            "labelmap contains none of the selected labels: %s"
            % ", ".join(str(label) for label in selection)
        )
    return mask


def _foreground_description(labels: tuple[int, ...] | None) -> str:
    if labels is None:
        return "all nonzero voxels"
    return "labels %s" % ", ".join(str(label) for label in labels)


@dataclass(frozen=True)
class LabelStats:
    """Geometry of one label in its source grid."""

    label: int
    voxels: int
    volume_mm3: float
    bbox_index: tuple[int, int, int]
    bbox_size: tuple[int, int, int]
    touches_boundary: bool


def inspect_labels(image: sitk.Image) -> dict[int, LabelStats]:
    """Per-label voxel count, physical volume, bounding box, and boundary contact.

    The image must already be a valid discrete labelmap; label 0 is background.
    """
    if not _validate_labelmap(image):
        return {}
    discrete = sitk.Cast(image, sitk.sitkUInt32)
    shape = sitk.LabelShapeStatisticsImageFilter()
    shape.SetBackgroundValue(0)
    shape.ComputePerimeterOff()
    shape.ComputeFeretDiameterOff()
    shape.ComputeOrientedBoundingBoxOff()
    shape.Execute(discrete)
    size = image.GetSize()
    voxel_mm3 = float(np.prod(image.GetSpacing()))
    stats: dict[int, LabelStats] = {}
    for label in shape.GetLabels():
        box = shape.GetBoundingBox(label)
        index = (int(box[0]), int(box[1]), int(box[2]))
        extent = (int(box[3]), int(box[4]), int(box[5]))
        touches = any(
            index[axis] == 0 or index[axis] + extent[axis] >= size[axis]
            for axis in range(3)
        )
        voxels = int(shape.GetNumberOfPixels(label))
        stats[int(label)] = LabelStats(
            label=int(label),
            voxels=voxels,
            volume_mm3=voxels * voxel_mm3,
            bbox_index=index,
            bbox_size=extent,
            touches_boundary=touches,
        )
    return stats


def _load_volume(
    candidate: VolumeCandidate,
    *,
    allow_large_volume: bool = False,
) -> volume_mod.Volume:
    if isinstance(candidate.source, DicomSource):
        raise ValueError(
            "labelmap input must be a NIfTI, NRRD, or MetaImage file, not DICOM"
        )
    return volume_mod.load(candidate, allow_large_volume=allow_large_volume)


def load(
    candidate: VolumeCandidate,
    *,
    allow_large_volume: bool = False,
    labels: Sequence[int] | None = None,
) -> LoadedLabelmap:
    selection = normalize_labels(labels)
    volume = _load_volume(candidate, allow_large_volume=allow_large_volume)
    mask = _binary_mask(volume.image, selection)
    provenance = volume_mod.provenance_for(volume)
    provenance.update(
        {
            "input_kind": "labelmap",
            "foreground": _foreground_description(selection),
        }
    )
    return LoadedLabelmap(volume=volume, mask=mask, provenance=provenance)


def embedded_label_names(candidate: VolumeCandidate) -> dict[int, str]:
    """Label names stored in a NIfTI label-table extension, or ``{}``."""
    if isinstance(candidate.source, FileSource):
        return labelnames.read_label_names(candidate.source.path)
    return {}


def _stage_runner(
    say: Logger,
    progress: ProgressSink | None,
    width: int,
    **context: Any,
):
    def step(message: str, function):
        say("%s ..." % message)
        if progress is not None:
            progress({"event": "stage_start", "stage": message, **context})
        before = time.time()
        value = function()
        seconds = time.time() - before
        say("  %-*s %6.1fs" % (width, message, seconds))
        if progress is not None:
            progress(
                {"event": "stage_end", "stage": message, "seconds": seconds, **context}
            )
        return value

    return step


def extract(
    candidate: VolumeCandidate,
    output_path: str,
    *,
    resample_mm: float | None = None,
    mask_smooth_mm: float | None = None,
    surface_smooth_iters: int | None = None,
    simplify_error_mm: float | None = None,
    post_surface_smooth_iters: int | None = None,
    keep_largest_component: bool | None = None,
    destep: DestepSettings | None = None,
    cap_field_of_view: bool = True,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
    labels: Sequence[int] | None = None,
    progress: ProgressSink | None = None,
) -> Result:
    selection = normalize_labels(labels)
    settings = resolve_surface_settings(
        resample_mm=resample_mm,
        mask_smooth_mm=mask_smooth_mm,
        surface_smooth_iters=surface_smooth_iters,
        simplify_error_mm=simplify_error_mm,
        post_surface_smooth_iters=post_surface_smooth_iters,
        keep_largest_component=keep_largest_component,
        destep=destep,
    )
    surface.validate_output_path(output_path)
    started = time.time()

    def say(message: str) -> None:
        if log:
            log(message)

    step = _stage_runner(say, progress, 34)

    warnings: list[str] = []

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    for message in (
        pipeline.mask_smoothing_warning(settings.mask_smooth_mm),
        pipeline.destep_warning(settings.destep),
    ):
        if message:
            add_warning(message)

    loaded = step(
        "load labelmap",
        lambda: load(
            candidate, allow_large_volume=allow_large_volume, labels=selection
        ),
    )
    for message in volume_mod.warnings_for(loaded.volume):
        add_warning(message)
    say(
        "labelmap %s  spacing %s mm  foreground %s"
        % (
            "x".join(str(value) for value in loaded.volume.size),
            " x ".join("%.3f" % value for value in loaded.volume.spacing),
            "nonzero" if selection is None else _foreground_description(selection),
        )
    )

    meshed = pipeline.mesh_binary_mask(
        loaded.mask,
        output_path,
        settings=settings,
        cap_field_of_view=cap_field_of_view,
        allow_large_volume=allow_large_volume,
        step=step,
        log=say,
        warn=add_warning,
    )
    provenance = {
        "input": loaded.provenance,
        "surface": surface_provenance(settings),
        "surface_finishing": meshed.surface_finishing,
        "capped_field_of_view": meshed.capped_field_of_view,
        "allow_large_volume": bool(allow_large_volume),
        "coordinate_system": "SimpleITK physical space of the input labelmap",
    }
    return Result(
        output_path=output_path,
        triangles=meshed.triangles,
        vertices=meshed.vertices,
        bounds_mm=meshed.bounds_mm,
        labelmap_components=meshed.labelmap_components,
        surface_components=meshed.surface_components,
        capped_field_of_view=meshed.capped_field_of_view,
        seconds=time.time() - started,
        warnings=warnings,
        provenance=provenance,
        quality=meshed.quality,
    )


def fuse(
    fixed: VolumeCandidate,
    moving: VolumeCandidate,
    output_path: str,
    *,
    grid_mm: float = DEFAULT_FUSION_GRID_MM,
    force: bool = False,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
    progress: ProgressSink | None = None,
) -> fusion.FusionResult:
    volume_output(output_path)
    if os.path.isdir(output_path):
        raise ValueError("volume output path is a directory: %s" % output_path)
    if not np.isfinite(grid_mm) or grid_mm <= 0:
        raise ValueError("grid_mm must be finite and greater than zero")
    if same_source(fixed, moving):
        raise fusion.FusionError(
            "fixed and moving inputs resolve to the same labelmap"
        )
    started = time.time()

    def say(message: str) -> None:
        if log:
            log(message)

    step = _stage_runner(say, progress, 36)

    warnings: list[str] = []

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    for message in fusion.check_compatible(fixed, moving):
        add_warning(message)
    add_warning(
        "labelmap contents are not verified; confirm that fixed and moving masks "
        "represent the same rigid structures before using the fused labelmap"
    )

    say("fixed  ID %d  %s  %s" % (fixed.id, fixed.format, fixed.source_name))
    say("moving ID %d  %s  %s" % (moving.id, moving.format, moving.source_name))
    fixed_loaded = step(
        "load fixed labelmap",
        lambda: load(fixed, allow_large_volume=allow_large_volume),
    )
    for message in volume_mod.warnings_for(fixed_loaded.volume):
        add_warning(message)
    moving_loaded = step(
        "load moving labelmap",
        lambda: load(moving, allow_large_volume=allow_large_volume),
    )
    for message in volume_mod.warnings_for(moving_loaded.volume):
        add_warning(message)

    fixed_provenance = fixed_loaded.provenance
    moving_provenance = moving_loaded.provenance
    fixed_loaded.volume.image = sitk.Image()
    moving_loaded.volume.image = sitk.Image()
    result = fusion.fuse_masks(
        fixed_loaded.mask,
        moving_loaded.mask,
        output_path,
        grid_mm=grid_mm,
        force=force,
        allow_large_volume=allow_large_volume,
        log=say,
        warn=warn,
        progress=progress,
    )
    result.seconds = time.time() - started
    result.warnings = warnings + result.warnings
    result.provenance.update(
        {
            "fixed": fixed_provenance,
            "moving": moving_provenance,
            "segmentation": {
                "input_kind": "labelmap",
                "validation": "finite, discrete, and non-negative",
                "foreground": "all nonzero source values normalized to 1",
            },
            "forced": bool(force),
            "allow_large_volume": bool(allow_large_volume),
            "coordinate_system": (
                "axis-aligned isotropic lattice in the fixed labelmap's "
                "SimpleITK physical coordinate system"
            ),
        }
    )
    return result


# --------------------------------------------------------------- split output
@dataclass
class LabelMesh:
    """One published mesh from a per-label or combined extraction."""

    label: int | None
    name: str | None
    output_path: str
    triangles: int
    vertices: int
    bounds_mm: tuple[float, ...]
    labelmap_components: int
    surface_components: int
    capped_field_of_view: bool
    seconds: float
    quality: dict[str, Any]
    surface: dict[str, Any]

    def payload(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "name": self.name,
            "output": self.output_path,
            "triangles": self.triangles,
            "vertices": self.vertices,
            "bounds_mm": list(self.bounds_mm),
            "labelmap_components": self.labelmap_components,
            "surface_components": self.surface_components,
            "capped_field_of_view": self.capped_field_of_view,
            "seconds": self.seconds,
            "surface": self.surface,
            "quality": self.quality,
        }


@dataclass
class SplitResult:
    """Every mesh written by :func:`extract_labels`."""

    output_dir: str
    meshes: list[LabelMesh]
    combined: LabelMesh | None
    skipped: list[dict[str, Any]]
    failed: list[dict[str, Any]]
    seconds: float
    warnings: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def valid(self) -> bool:
        meshes = [*self.meshes, *([self.combined] if self.combined else [])]
        return not self.failed and all(mesh.quality.get("valid") for mesh in meshes)


def crop_margin_voxels(
    settings: pipeline.SurfaceSettings, spacing: Sequence[float]
) -> tuple[int, int, int]:
    """Background voxels kept around a label so cropping never changes its surface.

    The margin covers three Gaussian sigmas of mask smoothing, two voxels of
    resampling support, and one voxel so a label that does not reach the image
    boundary never reaches the crop boundary either (which would cap it).
    """
    reach_mm = 3.0 * max(settings.mask_smooth_mm, 0.0) + 2.0 * max(
        settings.resample_mm, 0.0
    )
    margin = []
    for value in spacing:
        voxel = float(value)
        margin.append(int(math.ceil((reach_mm + 2.0 * voxel) / voxel)) + 1)
    return (margin[0], margin[1], margin[2])


def _crop_region(
    size: Sequence[int],
    index: Sequence[int],
    extent: Sequence[int],
    margin: Sequence[int],
) -> tuple[list[int], list[int]]:
    start = [max(0, int(index[axis]) - int(margin[axis])) for axis in range(3)]
    stop = [
        min(int(size[axis]), int(index[axis]) + int(extent[axis]) + int(margin[axis]))
        for axis in range(3)
    ]
    return start, [stop[axis] - start[axis] for axis in range(3)]


def _union_region(stats: Sequence[LabelStats]) -> tuple[list[int], list[int]]:
    low = [min(item.bbox_index[axis] for item in stats) for axis in range(3)]
    high = [
        max(item.bbox_index[axis] + item.bbox_size[axis] for item in stats)
        for axis in range(3)
    ]
    return low, [high[axis] - low[axis] for axis in range(3)]


def split_file_name(label: int, name: str | None, mesh_format: str) -> str:
    stem = labelnames.slug(name) if name else "label-%d" % label
    return "%03d_%s.%s" % (label, stem, mesh_format)


def extract_labels(
    candidate: VolumeCandidate,
    output_dir: str,
    *,
    labels: Sequence[int] | None = None,
    names: Mapping[int, str] | None = None,
    settings: pipeline.SurfaceSettings | SettingsForLabel | None = None,
    combined_path: str | None = None,
    combined_settings: pipeline.SurfaceSettings | None = None,
    mesh_format: str = "stl",
    cap_field_of_view: bool = True,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
    progress: ProgressSink | None = None,
    on_label: Callable[[LabelMesh], None] | None = None,
) -> SplitResult:
    """Extract one validated mesh per label, plus an optional combined mesh.

    Each label is cropped to its bounding box plus :func:`crop_margin_voxels`,
    so results match an uncropped extraction while memory and time scale with
    the structure instead of the whole scan. A label that reaches the image
    boundary is capped exactly as :func:`extract` would cap it. The combined
    mesh is extracted from the union of the selected labels; concatenating the
    per-label meshes would not form one valid closed surface.

    ``settings`` may be one :class:`pipeline.SurfaceSettings` for every label or
    a callable returning the settings for a label ID. Labels that fail are
    reported in ``failed`` and do not stop the remaining labels.
    """
    mesh_format = mesh_format.lower().lstrip(".")
    if "." + mesh_format not in SUPPORTED_MESH_EXTENSIONS:
        raise ValueError(
            "unsupported mesh format %r; supported: %s"
            % (mesh_format, ", ".join(SUPPORTED_MESH_EXTENSIONS))
        )
    selection = normalize_labels(labels)
    if combined_path is not None:
        surface.validate_output_path(combined_path)
    directory = Path(output_dir)
    if directory.exists() and not directory.is_dir():
        raise ValueError("split output path is not a directory: %s" % directory)

    def settings_for(label: int) -> pipeline.SurfaceSettings:
        if settings is None:
            return default_surface_settings()
        if isinstance(settings, pipeline.SurfaceSettings):
            chosen = settings
        else:
            chosen = settings(label)
        pipeline.validate_surface_settings(chosen)
        return chosen

    started = time.time()
    warnings: list[str] = []

    def say(message: str) -> None:
        if log:
            log(message)

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    step = _stage_runner(say, progress, 34)
    volume = step(
        "load labelmap",
        lambda: _load_volume(candidate, allow_large_volume=allow_large_volume),
    )
    for message in volume_mod.warnings_for(volume):
        add_warning(message)
    image = volume.image
    stats = step("measure labels", lambda: inspect_labels(image))
    if not stats:
        raise ValueError("labelmap contains no nonzero foreground voxels")

    label_names: dict[int, str] = dict(embedded_label_names(candidate))
    if names:
        label_names.update({int(key): str(value) for key, value in names.items()})

    skipped: list[dict[str, Any]] = []
    if selection is None:
        chosen_labels = sorted(stats)
    else:
        chosen_labels = [label for label in selection if label in stats]
        for label in selection:
            if label not in stats:
                skipped.append({"label": label, "reason": "no voxels"})
                add_warning("label %d has no voxels and was skipped" % label)
        if not chosen_labels:
            raise ValueError(
                "labelmap contains none of the selected labels: %s"
                % ", ".join(str(label) for label in selection)
            )

    directory.mkdir(parents=True, exist_ok=True)
    say(
        "labelmap %s  spacing %s mm  %d label(s)"
        % (
            "x".join(str(value) for value in image.GetSize()),
            " x ".join("%.3f" % value for value in image.GetSpacing()),
            len(chosen_labels),
        )
    )

    def mesh_region(
        label: int | None,
        start: list[int],
        extent: list[int],
        path: str,
        label_settings: pipeline.SurfaceSettings,
    ) -> LabelMesh:
        before = time.time()
        cropped = sitk.RegionOfInterest(image, extent, start)
        if label is None:
            mask = _selected_mask(cropped, tuple(chosen_labels))
        else:
            mask = sitk.Cast(sitk.Equal(cropped, float(label)), sitk.sitkUInt8)
        cropped = sitk.Image()
        context = {} if label is None else {"label": label}
        label_step = _stage_runner(say, progress, 34, **context)
        meshed = pipeline.mesh_binary_mask(
            mask,
            path,
            settings=label_settings,
            cap_field_of_view=cap_field_of_view,
            allow_large_volume=allow_large_volume,
            step=label_step,
            log=say,
            warn=add_warning,
        )
        return LabelMesh(
            label=label,
            name=None if label is None else label_names.get(label),
            output_path=path,
            triangles=meshed.triangles,
            vertices=meshed.vertices,
            bounds_mm=meshed.bounds_mm,
            labelmap_components=meshed.labelmap_components,
            surface_components=meshed.surface_components,
            capped_field_of_view=meshed.capped_field_of_view,
            seconds=time.time() - before,
            quality=meshed.quality,
            surface=surface_provenance(label_settings),
        )

    meshes: list[LabelMesh] = []
    failed: list[dict[str, Any]] = []
    size = image.GetSize()
    spacing = image.GetSpacing()
    for position, label in enumerate(chosen_labels, start=1):
        name = label_names.get(label)
        path = str(directory / split_file_name(label, name, mesh_format))
        if progress is not None:
            progress(
                {
                    "event": "label_start",
                    "label": label,
                    "name": name,
                    "index": position,
                    "total": len(chosen_labels),
                }
            )
        say("label %d%s (%d of %d)" % (label, " %s" % name if name else "", position, len(chosen_labels)))
        try:
            label_settings = settings_for(label)
            item = stats[label]
            start, extent = _crop_region(
                size,
                item.bbox_index,
                item.bbox_size,
                crop_margin_voxels(label_settings, spacing),
            )
            mesh = mesh_region(label, start, extent, path, label_settings)
        except (ValueError, RuntimeError, OSError) as exc:
            failed.append({"label": label, "name": name, "error": str(exc)})
            add_warning("label %d%s failed: %s" % (label, " (%s)" % name if name else "", exc))
            if progress is not None:
                progress({"event": "label_failed", "label": label, "error": str(exc)})
            continue
        meshes.append(mesh)
        if progress is not None:
            progress(
                {
                    "event": "label_end",
                    "label": label,
                    "output": mesh.output_path,
                    "valid": bool(mesh.quality.get("valid")),
                    "seconds": mesh.seconds,
                }
            )
        if on_label is not None:
            on_label(mesh)

    combined: LabelMesh | None = None
    if combined_path is not None:
        chosen_stats = [stats[label] for label in chosen_labels]
        union_settings = combined_settings or settings_for(chosen_labels[0])
        start, extent = _union_region(chosen_stats)
        start, extent = _crop_region(
            size, start, extent, crop_margin_voxels(union_settings, spacing)
        )
        if progress is not None:
            progress({"event": "combined_start"})
        say("combined model of %d label(s)" % len(chosen_labels))
        try:
            combined = mesh_region(None, start, extent, combined_path, union_settings)
        except (ValueError, RuntimeError, OSError) as exc:
            failed.append({"label": None, "name": "combined", "error": str(exc)})
            add_warning("combined model failed: %s" % exc)
        else:
            if progress is not None:
                progress(
                    {
                        "event": "combined_end",
                        "output": combined.output_path,
                        "valid": bool(combined.quality.get("valid")),
                    }
                )

    provenance = {
        "input": {
            **volume_mod.provenance_for(volume),
            "input_kind": "labelmap",
            "foreground": _foreground_description(tuple(chosen_labels)),
        },
        "labels": {
            str(label): {
                "name": label_names.get(label),
                "voxels": stats[label].voxels,
                "volume_mm3": stats[label].volume_mm3,
                "touches_boundary": stats[label].touches_boundary,
            }
            for label in chosen_labels
        },
        "crop": "per-label bounding box plus smoothing and resampling margin",
        "allow_large_volume": bool(allow_large_volume),
        "coordinate_system": "SimpleITK physical space of the input labelmap",
    }
    return SplitResult(
        output_dir=str(directory),
        meshes=meshes,
        combined=combined,
        skipped=skipped,
        failed=failed,
        seconds=time.time() - started,
        warnings=warnings,
        provenance=provenance,
    )


# ------------------------------------------------------ N-way label fusion
@dataclass
class LabelFusionResult:
    """Result of registering and fusing two or more labelmaps."""

    output_path: str
    output_format: str
    compression: str
    pixel_type: str
    preserve_labels: bool
    grid_mm: float
    grid_size: tuple[int, int, int]
    grid_origin_mm: tuple[float, float, float]
    grid_direction: tuple[float, ...]
    registrations: list[registration.RegistrationResult]
    labels: dict[int, dict[str, Any]]
    seconds: float
    warnings: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def volume_fused_mm3(self) -> float:
        return sum(item["volume_mm3"] for item in self.labels.values())


def _label_pixel_type(max_label: int, binary: bool) -> tuple[int, str]:
    if binary or max_label <= 255:
        return sitk.sitkUInt8, "uint8"
    if max_label <= 65535:
        return sitk.sitkUInt16, "uint16"
    return sitk.sitkUInt32, "uint32"


def fuse_labels(
    fixed: VolumeCandidate,
    movings: Sequence[VolumeCandidate],
    output_path: str,
    *,
    grid_mm: float | None = None,
    labels: Sequence[int] | None = None,
    preserve_labels: bool = True,
    names: Mapping[int, str] | None = None,
    force: bool = False,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
    progress: ProgressSink | None = None,
) -> LabelFusionResult:
    """Register every moving labelmap to ``fixed`` and fuse them on one grid.

    Registration uses the union of the selected labels of each input, with the
    same quality gates as two-input fusion. With ``preserve_labels`` the output
    keeps label IDs (see :func:`fusion.fuse_label_fields`); otherwise it is a
    binary ``0``/``1`` union. ``grid_mm`` defaults to the finest input spacing
    when labels are preserved and to the binary fusion default otherwise.
    Label names from NIfTI label tables (and ``names``) are re-embedded in
    NIfTI outputs.
    """
    output = volume_output(output_path)
    if os.path.isdir(output_path):
        raise ValueError("volume output path is a directory: %s" % output_path)
    if not movings:
        raise ValueError("fusion needs at least one moving labelmap")
    inputs = [fixed, *movings]
    for index, first in enumerate(inputs):
        for second in inputs[index + 1 :]:
            if same_source(first, second):
                raise fusion.FusionError("two fusion inputs resolve to the same labelmap")
    if grid_mm is not None and (not np.isfinite(grid_mm) or grid_mm <= 0):
        raise ValueError("grid_mm must be finite and greater than zero")
    selection = normalize_labels(labels)
    started = time.time()

    def say(message: str) -> None:
        if log:
            log(message)

    step = _stage_runner(say, progress, 36)
    warnings: list[str] = []

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    add_warning(
        "subject identity is not verified; confirm that every input shows the "
        "same subject before using the fused labelmap"
    )
    add_warning(fusion.RIGID_REGISTRATION_WARNING)

    volumes: list[volume_mod.Volume] = []
    unions: list[sitk.Image] = []
    label_names: dict[int, str] = {}
    for index, candidate in enumerate(inputs):
        role = "fixed" if index == 0 else "moving %d" % index
        say("%-8s ID %d  %s  %s" % (role, candidate.id, candidate.format, candidate.source_name))
        volume = step(
            "load %s labelmap" % role,
            lambda candidate=candidate: _load_volume(
                candidate, allow_large_volume=allow_large_volume
            ),
        )
        for message in volume_mod.warnings_for(volume):
            add_warning(message)
        unions.append(_binary_mask(volume.image, selection))
        volumes.append(volume)
        for label, name in embedded_label_names(candidate).items():
            label_names.setdefault(label, name)
    if names:
        label_names.update({int(key): str(value) for key, value in names.items()})

    transforms: list[np.ndarray] = [np.eye(4)]
    registrations: list[registration.RegistrationResult] = []
    for index in range(1, len(inputs)):
        role = "moving %d" % index

        def register(index: int = index) -> registration.RegistrationResult:
            try:
                return registration.rigid_register(
                    unions[0], unions[index], log=lambda message: say("  " + message)
                )
            except registration.RegistrationError as exc:
                raise fusion.FusionError("%s: %s" % (inputs[index].source_name, exc)) from None

        registered = step("register %s to fixed" % role, register)
        try:
            fusion.check_registration(registered, force=force)
        except fusion.FusionError as exc:
            raise fusion.FusionError("%s: %s" % (inputs[index].source_name, exc)) from None
        for line in registered.summary().splitlines():
            say("  " + line.strip())
        registrations.append(registered)
        transforms.append(registered.transform)
    unions.clear()

    finest = min(min(volume.image.GetSpacing()) for volume in volumes)
    if grid_mm is None:
        grid_mm = float(finest) if preserve_labels else DEFAULT_FUSION_GRID_MM
    elif grid_mm > finest:
        add_warning(
            "the fused grid is %.2f mm but the finest input voxel is %.3f mm; "
            "structures thinner than the grid are lost" % (grid_mm, finest)
        )

    fused = step(
        "fuse label occupancy",
        lambda: fusion.fuse_label_fields(
            [volume.image for volume in volumes],
            transforms,
            float(grid_mm),
            labels=selection,
            binary=not preserve_labels,
            allow_large_volume=allow_large_volume,
            log=say,
            progress=progress,
        ),
    )
    provenance_inputs = [volume_mod.provenance_for(volume) for volume in volumes]
    volumes.clear()
    if not fused.voxels:
        raise fusion.FusionError("the fused labelmap contains no foreground voxels")

    max_label = max(fused.voxels)
    pixel_id, pixel_type = _label_pixel_type(max_label, not preserve_labels)
    image = sitk.Cast(sitk.GetImageFromArray(fused.labels), pixel_id)
    fused.labels = np.empty((0, 0, 0), dtype=np.uint32)
    image.SetSpacing((fused.grid_mm,) * 3)
    image.SetOrigin(fused.grid_origin_mm)
    image.SetDirection(tuple(np.eye(3).ravel()))
    owner = [image]
    image = sitk.Image()
    step(
        "validate and publish fused labelmap",
        lambda: volume_mod.write_verified_volume(
            owner[0], output_path, release_source=owner.clear
        ),
    )
    kept_names = {
        label: label_names[label]
        for label in fused.voxels
        if preserve_labels and label in label_names
    }
    if kept_names:
        step(
            "embed label names",
            lambda: labelnames.embed_label_names(output_path, kept_names),
        )

    voxel_mm3 = fused.grid_mm**3
    label_payload = {
        label: {
            "name": kept_names.get(label),
            "voxels": count,
            "volume_mm3": count * voxel_mm3,
        }
        for label, count in sorted(fused.voxels.items())
    }
    say(
        "fused %d label(s)  %.0f cm3"
        % (len(label_payload), sum(item["volume_mm3"] for item in label_payload.values()) / 1000.0)
    )
    result = LabelFusionResult(
        output_path=output_path,
        output_format=output.format,
        compression=output.compression,
        pixel_type=pixel_type,
        preserve_labels=preserve_labels,
        grid_mm=fused.grid_mm,
        grid_size=fused.grid_size,
        grid_origin_mm=fused.grid_origin_mm,
        grid_direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        registrations=registrations,
        labels=label_payload,
        seconds=time.time() - started,
        warnings=warnings,
    )
    result.provenance = {
        "output": {
            "kind": "labelmap" if preserve_labels else "binary labelmap",
            "format": output.format,
            "compression": output.compression,
            "pixel_type": pixel_type,
            "label_names_embedded": bool(kept_names) and output.format == "NIfTI",
        },
        "grid": {
            "size": list(fused.grid_size),
            "spacing_mm": [fused.grid_mm] * 3,
            "origin_mm": list(fused.grid_origin_mm),
            "direction": list(result.grid_direction),
        },
        "inputs": [
            {**record, "input_kind": "labelmap", "foreground": _foreground_description(selection)}
            for record in provenance_inputs
        ],
        "registrations": [
            fusion._registration_provenance(registered) for registered in registrations
        ],
        "fusion": {
            "method": "per-label antialiased occupancy, best label above 0.5",
            "labels_preserved": preserve_labels,
        },
        "forced": bool(force),
        "allow_large_volume": bool(allow_large_volume),
        "coordinate_system": (
            "axis-aligned isotropic lattice in the fixed labelmap's "
            "SimpleITK physical coordinate system"
        ),
    }
    return result
