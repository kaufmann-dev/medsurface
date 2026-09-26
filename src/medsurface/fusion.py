"""Rigid registration and binary foreground fusion for medical volumes."""

from __future__ import annotations

import math
import os
import time
from dataclasses import asdict, dataclass, field
from numbers import Real
from typing import Any, Callable

import numpy as np
import SimpleITK as sitk

from . import pipeline, registration, segment
from . import volume as volume_mod
from .catalog import VolumeCandidate, same_source
from .defaults import DEFAULT_FUSION_GRID_MM, MAX_VOXELS
from .outputs import volume_output
from .presets import Preset, validate_segmentation
from .registration import RegistrationResult

Logger = Callable[[str], None]

# Empirical gates selected from matching repeat studies and a deliberate
# unrelated-anatomy impostor. The weaker directional surface overlap is the
# strongest discriminator; the shared-field Dice and volume prevent apparently
# confident registration over too little common anatomy.
MIN_SURFACE_OVERLAP = 0.60
MIN_SHARED_FOV_DICE = 0.55
MIN_SHARED_FOV_MM3 = 20_000.0


class FusionError(RuntimeError):
    """Registration or fused-grid failure that makes publication unsafe."""


@dataclass
class FusionResult:
    output_path: str
    output_format: str
    compression: str
    grid_mm: float
    grid_size: tuple[int, int, int]
    grid_origin_mm: tuple[float, float, float]
    grid_direction: tuple[float, ...]
    registration: RegistrationResult
    foreground_fixed_voxels: int
    foreground_moving_voxels: int
    foreground_fused_voxels: int
    volume_fixed_mm3: float
    volume_moving_mm3: float
    volume_fused_mm3: float
    seconds: float
    warnings: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)


RIGID_REGISTRATION_WARNING = (
    "fusion uses rigid registration and is intended for matching non-deforming "
    "anatomy such as bone; anatomy that moved or deformed between acquisitions "
    "can produce a plausible but incorrect fusion"
)


def check_compatible(
    fixed: VolumeCandidate,
    moving: VolumeCandidate,
) -> list[str]:
    """Reject duplicate inputs and state the identity contract."""
    if same_source(fixed, moving):
        raise FusionError("fixed and moving inputs resolve to the same volume")
    return [
        "subject identity is not verified; confirm that fixed and moving volumes "
        "show the same subject before using the fused labelmap",
        RIGID_REGISTRATION_WARNING,
    ]


def _invalid_registration(detail: str) -> FusionError:
    return FusionError(
        "invalid registration result: %s; --force cannot override malformed "
        "registration diagnostics" % detail
    )


def _finite_registration_array(
    name: str,
    value: Any,
    shape: tuple[int, ...],
) -> np.ndarray:
    if not isinstance(value, np.ndarray) or value.shape != shape:
        raise _invalid_registration(
            "%s must be a numeric array with shape %s" % (name, shape)
        )
    if not np.issubdtype(value.dtype, np.number) or np.issubdtype(
        value.dtype, np.complexfloating
    ):
        raise _invalid_registration("%s must contain real numbers" % name)
    if not np.all(np.isfinite(value)):
        raise _invalid_registration("%s must contain only finite values" % name)
    return value


def _finite_registration_metric(name: str, value: Any) -> float:
    if not isinstance(value, Real) or isinstance(value, bool):
        raise _invalid_registration("%s must be a real number" % name)
    scalar = float(value)
    if not math.isfinite(scalar):
        raise _invalid_registration("%s must be finite" % name)
    return scalar


def _validate_registration_result(result: RegistrationResult) -> dict[str, float]:
    transform = _finite_registration_array("transform", result.transform, (4, 4))
    _finite_registration_array("fft_translation_mm", result.fft_translation_mm, (3,))
    _finite_registration_array("translation_mm", result.translation_mm, (3,))

    if not np.allclose(
        transform[3],
        np.array([0.0, 0.0, 0.0, 1.0]),
        rtol=0.0,
        atol=1e-6,
    ):
        raise _invalid_registration("transform is not homogeneous")

    rotation = transform[:3, :3]
    if np.any(np.abs(rotation) > 1.0 + 1e-5) or not np.allclose(
        rotation.T @ rotation,
        np.eye(3),
        rtol=0.0,
        atol=1e-5,
    ):
        raise _invalid_registration("transform is not rigid")
    if not math.isclose(float(np.linalg.det(rotation)), 1.0, rel_tol=0.0, abs_tol=1e-5):
        raise _invalid_registration("transform is not a proper rigid transform")

    metrics = {
        "rotation_deg": _finite_registration_metric(
            "rotation_deg", result.rotation_deg
        ),
        "inlier_rms_mm": _finite_registration_metric(
            "inlier_rms_mm", result.inlier_rms_mm
        ),
        "inlier_median_mm": _finite_registration_metric(
            "inlier_median_mm", result.inlier_median_mm
        ),
        "overlap_moving_in_fixed": _finite_registration_metric(
            "overlap_moving_in_fixed", result.overlap_moving_in_fixed
        ),
        "overlap_fixed_in_moving": _finite_registration_metric(
            "overlap_fixed_in_moving", result.overlap_fixed_in_moving
        ),
        "shared_fov_dice": _finite_registration_metric(
            "shared_fov_dice", result.shared_fov_dice
        ),
        "shared_fov_mm3": _finite_registration_metric(
            "shared_fov_mm3", result.shared_fov_mm3
        ),
    }

    if not 0.0 <= metrics["rotation_deg"] <= 180.0:
        raise _invalid_registration("rotation_deg must be between 0 and 180")
    for name in ("inlier_rms_mm", "inlier_median_mm", "shared_fov_mm3"):
        if metrics[name] < 0.0:
            raise _invalid_registration("%s must be non-negative" % name)
    for name in (
        "overlap_moving_in_fixed",
        "overlap_fixed_in_moving",
        "shared_fov_dice",
    ):
        if not 0.0 <= metrics[name] <= 1.0:
            raise _invalid_registration("%s must be between 0 and 1" % name)
    return metrics


def check_registration(result: RegistrationResult, force: bool = False) -> None:
    """Reject a transform which fails the established geometric gates."""
    metrics = _validate_registration_result(result)
    problems = []
    if metrics["shared_fov_mm3"] < MIN_SHARED_FOV_MM3:
        problems.append(
            "the scans share only %.0f cm3 of imaged space (need %.0f)"
            % (metrics["shared_fov_mm3"] / 1000.0, MIN_SHARED_FOV_MM3 / 1000.0)
        )
    surface_overlap = min(
        metrics["overlap_moving_in_fixed"],
        metrics["overlap_fixed_in_moving"],
    )
    if surface_overlap < MIN_SURFACE_OVERLAP:
        problems.append(
            "surfaces agree over only %.1f%% of the shared field of view "
            "(%.1f%% moving->fixed, %.1f%% fixed->moving; need %.0f%% both ways)"
            % (
                100 * surface_overlap,
                100 * metrics["overlap_moving_in_fixed"],
                100 * metrics["overlap_fixed_in_moving"],
                100 * MIN_SURFACE_OVERLAP,
            )
        )
    if metrics["shared_fov_dice"] < MIN_SHARED_FOV_DICE:
        problems.append(
            "foreground agreement in the shared field of view is %.3f (need %.2f)"
            % (metrics["shared_fov_dice"], MIN_SHARED_FOV_DICE)
        )
    if not problems or force:
        return
    raise FusionError(
        "registration failed its quality gates, so the two scans probably do not "
        "show the same anatomy:\n  - "
        + "\n  - ".join(problems)
        + "\n"
        + result.summary()
        + "\nPass --force to fuse them anyway."
    )


def _corners(image: sitk.Image) -> np.ndarray:
    size = np.asarray(image.GetSize()) - 1
    points = []
    for i in (0, size[0]):
        for j in (0, size[1]):
            for k in (0, size[2]):
                points.append(
                    image.TransformIndexToPhysicalPoint((int(i), int(j), int(k)))
                )
    return np.asarray(points, dtype=float)


def _common_grid(
    fixed_mask: sitk.Image,
    moving_mask: sitk.Image,
    transform: np.ndarray,
    grid_mm: float,
    *,
    allow_large_volume: bool = False,
) -> tuple[tuple[int, int, int], np.ndarray]:
    rotation, translation = transform[:3, :3], transform[:3, 3]
    points = np.vstack(
        [
            _corners(fixed_mask),
            (rotation @ _corners(moving_mask).T).T + translation,
        ]
    )
    scaled = points / grid_mm
    if not np.all(np.isfinite(scaled)):
        raise FusionError(
            "the fused grid coordinates overflow at %.4g mm; raise --grid-mm" % grid_mm
        )
    low_index = np.floor(scaled.min(0)) - 2
    high_index = np.ceil(scaled.max(0)) + 2
    planned = high_index - low_index + 1
    if not np.all(np.isfinite(planned)):
        raise FusionError(
            "the fused grid is too large at %.4g mm; raise --grid-mm" % grid_mm
        )
    if not allow_large_volume and np.any(planned > MAX_VOXELS):
        raise FusionError(
            "the fused grid is too large at %.4g mm; raise --grid-mm or pass "
            "--allow-large-volume" % grid_mm
        )
    size = (int(planned[0]), int(planned[1]), int(planned[2]))
    voxels = math.prod(size)
    if voxels > MAX_VOXELS and not allow_large_volume:
        raise FusionError(
            "the fused grid would hold %s voxels at %.4g mm, above the default "
            "limit of %s; raise --grid-mm or pass --allow-large-volume to attempt "
            "it (this may exhaust memory)" % (f"{voxels:,}", grid_mm, f"{MAX_VOXELS:,}")
        )
    origin = np.asarray(low_index * grid_mm, dtype=float)
    return size, origin


def _resample_field(
    image: sitk.Image,
    size: tuple[int, int, int],
    origin: np.ndarray,
    grid_mm: float,
    transform: sitk.Transform,
) -> sitk.Image:
    resample = sitk.ResampleImageFilter()
    resample.SetOutputSpacing((grid_mm, grid_mm, grid_mm))
    resample.SetSize([int(value) for value in size])
    resample.SetOutputOrigin([float(value) for value in origin])
    resample.SetOutputDirection(np.eye(3).ravel().tolist())
    resample.SetInterpolator(sitk.sitkLinear)
    resample.SetDefaultPixelValue(0.0)
    resample.SetTransform(transform)
    return resample.Execute(image)


def _registration_provenance(result: RegistrationResult) -> dict[str, Any]:
    return {
        "transform_moving_to_fixed": result.transform.tolist(),
        "fft_translation_mm": result.fft_translation_mm.tolist(),
        "rotation_deg": result.rotation_deg,
        "translation_mm": result.translation_mm.tolist(),
        "inlier_rms_mm": result.inlier_rms_mm,
        "inlier_median_mm": result.inlier_median_mm,
        "overlap_moving_in_fixed": result.overlap_moving_in_fixed,
        "overlap_fixed_in_moving": result.overlap_fixed_in_moving,
        "surface_overlap": result.surface_overlap,
        "shared_fov_dice": result.shared_fov_dice,
        "shared_fov_mm3": result.shared_fov_mm3,
    }


def _base_provenance(result: FusionResult) -> dict[str, Any]:
    return {
        "output": {
            "kind": "binary labelmap",
            "format": result.output_format,
            "compression": result.compression,
            "pixel_type": "uint8",
            "components": 1,
            "background": 0,
            "foreground": 1,
        },
        "grid": {
            "size": list(result.grid_size),
            "spacing_mm": [result.grid_mm] * 3,
            "origin_mm": list(result.grid_origin_mm),
            "direction": list(result.grid_direction),
        },
        "registration": _registration_provenance(result.registration),
    }


def fuse_masks(
    fixed_mask: sitk.Image,
    moving_mask: sitk.Image,
    output_path: str,
    *,
    grid_mm: float = DEFAULT_FUSION_GRID_MM,
    force: bool = False,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
    progress: Callable[[dict[str, Any]], None] | None = None,
) -> FusionResult:
    """Register, union, quantize, and publish two binary foreground masks."""
    output = volume_output(output_path)
    if os.path.isdir(output_path):
        raise ValueError("volume output path is a directory: %s" % output_path)
    if not math.isfinite(grid_mm) or grid_mm <= 0:
        raise ValueError("grid_mm must be finite and greater than zero")
    started = time.time()

    def say(message: str) -> None:
        if log:
            log(message)

    def step(message: str, function):
        say("%s ..." % message)
        if progress is not None:
            progress({"event": "stage_start", "stage": message})
        before = time.time()
        value = function()
        seconds = time.time() - before
        say("  %-36s %6.1fs" % (message, seconds))
        if progress is not None:
            progress({"event": "stage_end", "stage": message, "seconds": seconds})
        return value

    warnings: list[str] = []

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    say("registering ...")
    if progress is not None:
        progress({"event": "stage_start", "stage": "registering"})
    try:
        registered = registration.rigid_register(
            fixed_mask,
            moving_mask,
            log=lambda message: say("  " + message),
        )
    except registration.RegistrationError as exc:
        raise FusionError(str(exc)) from None
    check_registration(registered, force=force)
    if progress is not None:
        progress({"event": "stage_end", "stage": "registering"})
    for line in registered.summary().splitlines():
        say("  " + line.strip())

    finest = min(min(fixed_mask.GetSpacing()), min(moving_mask.GetSpacing()))
    if grid_mm > finest:
        add_warning(
            "the fused grid is %.2f mm but the finest input voxel is %.3f mm; "
            "structures thinner than the grid are lost. Lower --grid-mm to keep "
            "them, at cubic cost in memory." % (grid_mm, finest)
        )

    size, origin = step(
        "plan fused grid",
        lambda: _common_grid(
            fixed_mask,
            moving_mask,
            registered.transform,
            grid_mm,
            allow_large_volume=allow_large_volume,
        ),
    )
    say(
        "fused grid %s at %.2f mm isotropic (%.0f M voxels)"
        % (
            "x".join(str(value) for value in size),
            grid_mm,
            math.prod(size) / 1e6,
        )
    )
    identity = sitk.Transform(3, sitk.sitkIdentity)
    inverse = registration.inverse_transform(registered.transform)
    fixed_field = step(
        "resample fixed",
        lambda: _resample_field(
            segment.antialias_for_grid(fixed_mask, grid_mm),
            size,
            origin,
            grid_mm,
            identity,
        ),
    )
    moving_field = step(
        "resample moving",
        lambda: _resample_field(
            segment.antialias_for_grid(moving_mask, grid_mm),
            size,
            origin,
            grid_mm,
            inverse,
        ),
    )
    fused_field = step(
        "fuse occupancy fields",
        lambda: sitk.Maximum(fixed_field, moving_field),
    )

    def count_foreground() -> tuple[int, int, int]:
        fixed_count = int(
            np.count_nonzero(sitk.GetArrayViewFromImage(fixed_field) > 0.5)
        )
        moving_count = int(
            np.count_nonzero(sitk.GetArrayViewFromImage(moving_field) > 0.5)
        )
        fused_count = int(
            np.count_nonzero(sitk.GetArrayViewFromImage(fused_field) > 0.5)
        )
        return fixed_count, moving_count, fused_count

    fixed_count, moving_count, fused_count = step(
        "measure fused foreground",
        count_foreground,
    )
    if fused_count == 0:
        raise FusionError("the fused labelmap contains no foreground voxels")
    voxel_mm3 = grid_mm**3
    say(
        "foreground: fixed %.0f cm3 | moving %.0f cm3 | fused %.0f cm3"
        % (
            fixed_count * voxel_mm3 / 1000.0,
            moving_count * voxel_mm3 / 1000.0,
            fused_count * voxel_mm3 / 1000.0,
        )
    )

    binary = sitk.Cast(sitk.Greater(fused_field, 0.5), sitk.sitkUInt8)
    for key in binary.GetMetaDataKeys():
        binary.EraseMetaData(key)
    fixed_field = sitk.Image()
    moving_field = sitk.Image()
    fused_field = sitk.Image()
    fixed_mask = sitk.Image()
    moving_mask = sitk.Image()
    owner = [binary]
    binary = sitk.Image()
    step(
        "validate and publish fused labelmap",
        lambda: volume_mod.write_verified_volume(
            owner[0],
            output_path,
            release_source=owner.clear,
        ),
    )

    result = FusionResult(
        output_path=output_path,
        output_format=output.format,
        compression=output.compression,
        grid_mm=grid_mm,
        grid_size=size,
        grid_origin_mm=(float(origin[0]), float(origin[1]), float(origin[2])),
        grid_direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        registration=registered,
        foreground_fixed_voxels=fixed_count,
        foreground_moving_voxels=moving_count,
        foreground_fused_voxels=fused_count,
        volume_fixed_mm3=fixed_count * voxel_mm3,
        volume_moving_mm3=moving_count * voxel_mm3,
        volume_fused_mm3=fused_count * voxel_mm3,
        seconds=time.time() - started,
        warnings=warnings,
    )
    result.provenance = _base_provenance(result)
    return result


def _segmentation_settings(preset: Preset) -> dict[str, Any]:
    values = asdict(preset)
    keys = (
        "name",
        "description",
        "threshold",
        "threshold_unit",
        "threshold_max",
        "median_mm",
        "closing_mm",
        "opening_mm",
        "min_island_mm3",
        "keep_largest_island",
    )
    return {key: values[key] for key in keys}


def fuse(
    fixed: VolumeCandidate,
    moving: VolumeCandidate,
    preset: Preset,
    output_path: str,
    *,
    fixed_threshold: float | str | None = None,
    moving_threshold: float | str | None = None,
    grid_mm: float = DEFAULT_FUSION_GRID_MM,
    force: bool = False,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
) -> FusionResult:
    """Segment two intensity volumes independently and fuse their foreground."""
    volume_output(output_path)
    if os.path.isdir(output_path):
        raise ValueError("volume output path is a directory: %s" % output_path)
    if not math.isfinite(grid_mm) or grid_mm <= 0:
        raise ValueError("grid_mm must be finite and greater than zero")
    validate_segmentation(preset)
    for name, threshold in (
        ("fixed_threshold", fixed_threshold),
        ("moving_threshold", moving_threshold),
    ):
        if threshold not in (None, "auto") and not math.isfinite(float(threshold)):
            raise ValueError("%s must be finite or 'auto'" % name)
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

    for message in check_compatible(fixed, moving):
        add_warning(message)
    for message in pipeline.threshold_warnings(
        fixed, preset, fixed_threshold, "--fixed-threshold"
    ):
        add_warning(message)
    for message in pipeline.threshold_warnings(
        moving, preset, moving_threshold, "--moving-threshold"
    ):
        add_warning(message)

    say("fixed  ID %d  %s  %s" % (fixed.id, fixed.format, fixed.source_name))
    say("moving ID %d  %s  %s" % (moving.id, moving.format, moving.source_name))
    fixed_volume = step(
        "load fixed volume",
        lambda: volume_mod.load(fixed, allow_large_volume=allow_large_volume),
    )
    for message in volume_mod.warnings_for(fixed_volume):
        add_warning(message)
    moving_volume = step(
        "load moving volume",
        lambda: volume_mod.load(moving, allow_large_volume=allow_large_volume),
    )
    for message in volume_mod.warnings_for(moving_volume):
        add_warning(message)

    fixed_value, fixed_source = step(
        "resolve fixed threshold",
        lambda: pipeline.resolve_threshold(fixed_volume.image, preset, fixed_threshold),
    )
    moving_value, moving_source = step(
        "resolve moving threshold",
        lambda: pipeline.resolve_threshold(
            moving_volume.image, preset, moving_threshold
        ),
    )
    say(
        "threshold: fixed %.1f (%s), moving %.1f (%s)"
        % (fixed_value, fixed_source, moving_value, moving_source)
    )
    fixed_provenance = pipeline.source_provenance(
        fixed_volume, fixed_value, fixed_source
    )
    moving_provenance = pipeline.source_provenance(
        moving_volume, moving_value, moving_source
    )
    fixed_mask = step(
        "segment fixed",
        lambda: pipeline.build_mask(fixed_volume.image, preset, fixed_value, say),
    )
    fixed_volume.image = sitk.Image()
    moving_mask = step(
        "segment moving",
        lambda: pipeline.build_mask(moving_volume.image, preset, moving_value, say),
    )
    moving_volume.image = sitk.Image()

    result = fuse_masks(
        fixed_mask,
        moving_mask,
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
            "segmentation": _segmentation_settings(preset),
            "forced": bool(force),
            "allow_large_volume": bool(allow_large_volume),
            "coordinate_system": (
                "axis-aligned isotropic lattice in the fixed volume's "
                "SimpleITK physical coordinate system"
            ),
        }
    )
    return result


# ------------------------------------------------------- N-way label fusion
def common_grid_many(
    images: list[sitk.Image],
    transforms: list[np.ndarray],
    grid_mm: float,
    *,
    allow_large_volume: bool = False,
) -> tuple[tuple[int, int, int], np.ndarray]:
    """Axis-aligned grid covering every image after its moving->fixed transform."""
    points = np.vstack(
        [
            (transform[:3, :3] @ _corners(image).T).T + transform[:3, 3]
            for image, transform in zip(images, transforms)
        ]
    )
    scaled = points / grid_mm
    if not np.all(np.isfinite(scaled)):
        raise FusionError(
            "the fused grid coordinates overflow at %.4g mm; raise --grid-mm" % grid_mm
        )
    low_index = np.floor(scaled.min(0)) - 2
    high_index = np.ceil(scaled.max(0)) + 2
    planned = high_index - low_index + 1
    if not np.all(np.isfinite(planned)):
        raise FusionError("the fused grid is too large at %.4g mm; raise --grid-mm" % grid_mm)
    size = (int(planned[0]), int(planned[1]), int(planned[2]))
    voxels = math.prod(size)
    if voxels > MAX_VOXELS and not allow_large_volume:
        raise FusionError(
            "the fused grid would hold %s voxels at %.4g mm, above the default "
            "limit of %s; raise --grid-mm or pass --allow-large-volume to attempt "
            "it (this may exhaust memory)" % (f"{voxels:,}", grid_mm, f"{MAX_VOXELS:,}")
        )
    return size, np.asarray(low_index * grid_mm, dtype=float)


def _antialias_margin(image: sitk.Image, grid_mm: float) -> list[int]:
    """Input voxels needed around a crop for antialiasing plus linear support."""
    sigma_mm = segment._ANTIALIAS_SIGMA_FACTOR * grid_mm
    return [
        int(math.ceil(4.0 * sigma_mm / spacing)) + 2 for spacing in image.GetSpacing()
    ]


@dataclass
class LabelFieldResult:
    """Label-preserving occupancy fusion on one isotropic grid."""

    labels: np.ndarray  # (z, y, x)
    grid_size: tuple[int, int, int]
    grid_origin_mm: tuple[float, float, float]
    grid_mm: float
    voxels: dict[int, int]


def fuse_label_fields(
    images: list[sitk.Image],
    transforms: list[np.ndarray],
    grid_mm: float,
    *,
    labels: tuple[int, ...] | None = None,
    binary: bool = False,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    progress: Callable[[dict[str, Any]], None] | None = None,
) -> LabelFieldResult:
    """Stream every label of every input onto one grid and keep the best label.

    Each label is cropped to its bounding box, antialiased for the target grid,
    resampled with its input's transform, and compared against a running
    per-voxel best occupancy. A voxel keeps the label with the highest occupancy
    when that occupancy exceeds 0.5, which is exactly the union for one label and
    a fair split between touching labels. Memory is two grid-sized arrays
    regardless of how many labels or inputs are fused. With ``binary`` every
    selected label of every input is treated as foreground ``1``.
    """

    def say(message: str) -> None:
        if log:
            log(message)

    size, origin = common_grid_many(
        images, transforms, grid_mm, allow_large_volume=allow_large_volume
    )
    shape_zyx = (size[2], size[1], size[0])
    best_value = np.zeros(shape_zyx, dtype=np.float32)
    best_label = np.zeros(shape_zyx, dtype=np.uint32)
    say(
        "fused grid %s at %.2f mm isotropic (%.0f M voxels)"
        % ("x".join(str(value) for value in size), grid_mm, math.prod(size) / 1e6)
    )

    for input_index, (image, transform) in enumerate(zip(images, transforms)):
        discrete = sitk.Cast(image, sitk.sitkUInt32)
        shape = sitk.LabelShapeStatisticsImageFilter()
        shape.ComputePerimeterOff()
        shape.ComputeFeretDiameterOff()
        shape.ComputeOrientedBoundingBoxOff()
        shape.Execute(discrete)
        present = [int(label) for label in shape.GetLabels()]
        if labels is not None:
            present = [label for label in present if label in labels]
        if binary and present:
            groups: list[tuple[int, list[int]]] = [(1, present)]
        else:
            groups = [(label, [label]) for label in present]
        margin = _antialias_margin(image, grid_mm)
        inverse = registration.inverse_transform(transform)
        image_size = image.GetSize()
        for group_label, members in groups:
            boxes = [shape.GetBoundingBox(member) for member in members]
            low = [min(box[axis] for box in boxes) for axis in range(3)]
            high = [max(box[axis] + box[axis + 3] for box in boxes) for axis in range(3)]
            start = [max(0, low[axis] - margin[axis]) for axis in range(3)]
            stop = [min(image_size[axis], high[axis] + margin[axis]) for axis in range(3)]
            extent = [stop[axis] - start[axis] for axis in range(3)]
            if progress is not None:
                progress(
                    {
                        "event": "label_start",
                        "input": input_index,
                        "label": group_label,
                    }
                )
            cropped = sitk.RegionOfInterest(discrete, extent, start)
            if len(members) == 1:
                mask = sitk.Equal(cropped, members[0])
            else:
                values = sitk.GetArrayViewFromImage(cropped)
                selected = np.isin(values, members).astype(np.uint8)
                mask = sitk.GetImageFromArray(selected)
                mask.CopyInformation(cropped)
            occupancy = segment.antialias_for_grid(mask, grid_mm)
            corners = (transform[:3, :3] @ _corners(occupancy).T).T + transform[:3, 3]
            low_index = np.floor((corners.min(0) - origin) / grid_mm).astype(int) - 1
            high_index = np.ceil((corners.max(0) - origin) / grid_mm).astype(int) + 2
            low_index = np.maximum(low_index, 0)
            high_index = np.minimum(high_index, np.asarray(size))
            sub_size = high_index - low_index
            if np.any(sub_size <= 0):
                continue
            field = _resample_field(
                occupancy,
                (int(sub_size[0]), int(sub_size[1]), int(sub_size[2])),
                origin + low_index * grid_mm,
                grid_mm,
                inverse,
            )
            values = sitk.GetArrayViewFromImage(field)
            region = (
                slice(int(low_index[2]), int(high_index[2])),
                slice(int(low_index[1]), int(high_index[1])),
                slice(int(low_index[0]), int(high_index[0])),
            )
            current = best_value[region]
            better = values > current
            current[better] = values[better]
            best_label[region][better] = group_label
            if progress is not None:
                progress(
                    {"event": "label_end", "input": input_index, "label": group_label}
                )

    result = np.where(best_value > 0.5, best_label, 0).astype(np.uint32)
    del best_value, best_label
    found, counts = np.unique(result, return_counts=True)
    voxels = {int(label): int(count) for label, count in zip(found, counts) if label != 0}
    return LabelFieldResult(
        labels=result,
        grid_size=size,
        grid_origin_mm=(float(origin[0]), float(origin[1]), float(origin[2])),
        grid_mm=grid_mm,
        voxels=voxels,
    )
