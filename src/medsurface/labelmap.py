"""Extract or fuse externally produced discrete segmentation labelmaps."""

from __future__ import annotations

import os
import time
from dataclasses import asdict, dataclass, replace
from typing import Any, Callable

import numpy as np
import SimpleITK as sitk

from . import fusion, pipeline, surface
from . import volume as volume_mod
from .catalog import DicomSource, VolumeCandidate, same_source
from .defaults import (
    DEFAULT_FUSION_GRID_MM,
    DEFAULT_LABELMAP_MASK_SMOOTH_MM,
    DEFAULT_LABELMAP_POST_SURFACE_SMOOTH_ITERS,
    DEFAULT_LABELMAP_SURFACE_SMOOTH_ITERS,
)
from .outputs import volume_output
from .presets import DestepSettings

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
        raise ValueError(
            "labelmap input must be a NIfTI, NRRD, or MetaImage file, not DICOM"
        )
    volume = volume_mod.load(candidate, allow_large_volume=allow_large_volume)
    mask = _binary_mask(volume.image)
    provenance = volume_mod.provenance_for(volume)
    provenance.update(
        {
            "input_kind": "labelmap",
            "foreground": "all nonzero voxels",
        }
    )
    return LoadedLabelmap(volume=volume, mask=mask, provenance=provenance)


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
) -> Result:
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

    for message in (
        pipeline.mask_smoothing_warning(settings.mask_smooth_mm),
        pipeline.destep_warning(settings.destep),
    ):
        if message:
            add_warning(message)

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
