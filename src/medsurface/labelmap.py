"""Convert or fuse externally produced discrete segmentation labelmaps."""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, replace
from typing import Any, Callable

import numpy as np
import SimpleITK as sitk

from . import merge as merge_mod
from . import pipeline, presets, surface
from . import volume as volume_mod
from .catalog import DicomSource, VolumeCandidate, same_source
from .defaults import DEFAULT_MERGE_GRID_MM

Logger = Callable[[str], None]


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
    """Use normal convert finishing while retaining every external mask shell."""
    return pipeline.surface_settings(
        presets.get("bone"),
        keep_largest_component=False,
    )


def resolve_surface_settings(
    *,
    resample_mm: float | None = None,
    smooth_iters: int | None = None,
    smooth_force: float | None = None,
    simplify_error_mm: float | None = None,
    post_smooth_iters: int | None = None,
) -> pipeline.SurfaceSettings:
    base = default_surface_settings()
    settings = replace(
        base,
        resample_mm=base.resample_mm if resample_mm is None else resample_mm,
        smooth_iters=base.smooth_iters if smooth_iters is None else smooth_iters,
        smooth_force=base.smooth_force if smooth_force is None else smooth_force,
        simplify_error_mm=(
            base.simplify_error_mm
            if simplify_error_mm is None
            else simplify_error_mm
        ),
        post_smooth_iters=(
            base.post_smooth_iters
            if post_smooth_iters is None
            else post_smooth_iters
        ),
    )
    pipeline.validate_surface_settings(settings)
    return settings


def _binary_mask(image: sitk.Image) -> sitk.Image:
    """Validate a discrete label volume and union every nonzero label."""
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

    if not has_foreground:
        raise ValueError("labelmap contains no nonzero foreground voxels")
    return sitk.Cast(sitk.NotEqual(image, 0), sitk.sitkUInt8)


def load(
    candidate: VolumeCandidate,
    *,
    allow_large_volume: bool = False,
) -> LoadedLabelmap:
    if isinstance(candidate.source, DicomSource):
        raise ValueError("labelmap input must be a NIfTI, NRRD, or MetaImage file, not DICOM")
    volume = volume_mod.load(candidate, allow_large_volume=allow_large_volume)
    mask = _binary_mask(volume.image)
    provenance = pipeline.volume_provenance(volume)
    provenance.update(
        {
            "input_kind": "labelmap",
            "foreground": "all nonzero voxels",
        }
    )
    return LoadedLabelmap(volume=volume, mask=mask, provenance=provenance)


def convert(
    candidate: VolumeCandidate,
    output_path: str,
    *,
    resample_mm: float | None = None,
    smooth_iters: int | None = None,
    smooth_force: float | None = None,
    simplify_error_mm: float | None = None,
    post_smooth_iters: int | None = None,
    cap_field_of_view: bool = True,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
) -> Result:
    settings = resolve_surface_settings(
        resample_mm=resample_mm,
        smooth_iters=smooth_iters,
        smooth_force=smooth_force,
        simplify_error_mm=simplify_error_mm,
        post_smooth_iters=post_smooth_iters,
    )
    surface.validate_output_path(output_path)
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
        "load labelmap",
        lambda: load(candidate, allow_large_volume=allow_large_volume),
    )
    for message in volume_mod.warnings_for(loaded.volume):
        add_warning(message)
    say(
        "labelmap %s  spacing %s mm  foreground nonzero"
        % (
            "x".join(str(value) for value in loaded.volume.size),
            " x ".join("%.3f" % value for value in loaded.volume.spacing),
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
        "surface": asdict(settings),
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


def merge(
    fixed: VolumeCandidate,
    moving: VolumeCandidate,
    output_path: str,
    *,
    grid_mm: float = DEFAULT_MERGE_GRID_MM,
    smooth_iters: int | None = None,
    smooth_force: float | None = None,
    simplify_error_mm: float | None = None,
    post_smooth_iters: int | None = None,
    force: bool = False,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
) -> merge_mod.MergeResult:
    settings = resolve_surface_settings(
        smooth_iters=smooth_iters,
        smooth_force=smooth_force,
        simplify_error_mm=simplify_error_mm,
        post_smooth_iters=post_smooth_iters,
    )
    surface.validate_output_path(output_path)
    if not np.isfinite(grid_mm) or grid_mm <= 0:
        raise ValueError("grid_mm must be finite and greater than zero")
    if same_source(fixed, moving):
        raise merge_mod.MergeError("fixed and moving inputs resolve to the same labelmap")
    started = time.time()

    def say(message: str) -> None:
        if log:
            log(message)

    def step(message: str, function):
        say("%s ..." % message)
        before = time.time()
        value = function()
        say("  %-36s %6.1fs" % (message, time.time() - before))
        return value

    warnings: list[str] = []

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    for message in merge_mod.check_compatible(fixed, moving):
        add_warning(message)
    add_warning(
        "labelmap contents are not verified; confirm that fixed and moving masks "
        "represent the same rigid structures before using the fused surface"
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

    fused = merge_mod.fuse_masks(
        fixed_loaded.mask,
        moving_loaded.mask,
        output_path,
        settings=settings,
        grid_mm=grid_mm,
        force=force,
        allow_large_volume=allow_large_volume,
        log=say,
        warn=warn,
    )
    warnings.extend(fused.warnings)
    registration = fused.registration
    provenance = {
        "fixed": fixed_loaded.provenance,
        "moving": moving_loaded.provenance,
        "surface": asdict(settings),
        "grid_mm": grid_mm,
        "surface_finishing": fused.surface_finishing,
        "transform_moving_to_fixed": registration.transform.tolist(),
        "rotation_deg": registration.rotation_deg,
        "registration": {
            "inlier_rms_mm": registration.inlier_rms_mm,
            "inlier_median_mm": registration.inlier_median_mm,
            "surface_overlap": registration.surface_overlap,
            "shared_fov_dice": registration.shared_fov_dice,
            "shared_fov_mm3": registration.shared_fov_mm3,
        },
        "coordinate_system": "SimpleITK physical space of the fixed labelmap",
        "forced": bool(force),
        "allow_large_volume": bool(allow_large_volume),
    }
    return merge_mod.MergeResult(
        output_path=output_path,
        triangles=fused.triangles,
        vertices=fused.vertices,
        bounds_mm=fused.bounds_mm,
        grid_mm=grid_mm,
        grid_size=fused.grid_size,
        registration=registration,
        volume_fixed_mm3=fused.volume_fixed_mm3,
        volume_moving_mm3=fused.volume_moving_mm3,
        volume_union_mm3=fused.volume_union_mm3,
        surface_components=fused.surface_components,
        seconds=time.time() - started,
        warnings=warnings,
        provenance=provenance,
        quality=fused.quality,
    )
