"""Rigid registration and binary foreground fusion for medical volumes."""

from __future__ import annotations

import math
import time
from dataclasses import asdict, dataclass, field
from numbers import Real
from typing import Any

import numpy as np
import SimpleITK as sitk

from . import pipeline, registration, segment
from . import volume as volume_mod
from .catalog import VolumeCandidate, same_source
from .defaults import DEFAULT_FUSION_GRID_MM, MAX_VOXELS
from .outputs import volume_output
from .presets import Preset, validate_segmentation
from .registration import RegistrationResult, physical_corners, resample_axis_aligned
from .stages import Logger, ProgressSink, Reporter

# Empirical gates selected from matching repeat studies and a deliberate
# unrelated-anatomy impostor. The weaker directional surface overlap is the
# strongest discriminator; the shared-field Dice and volume prevent apparently
# confident registration over too little common anatomy.
MIN_SURFACE_OVERLAP = 0.60
MIN_SHARED_FOV_DICE = 0.55
MIN_SHARED_FOV_MM3 = 20_000.0


IDENTITY_DIRECTION = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)


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


def validate_grid_mm(grid_mm: float) -> None:
    if not math.isfinite(grid_mm) or grid_mm <= 0:
        raise ValueError("grid_mm must be finite and greater than zero")


def grid_warning(grid_mm: float, images: list[sitk.Image]) -> str | None:
    """Explain detail lost when the fused grid is coarser than an input voxel."""
    finest = min(min(image.GetSpacing()) for image in images)
    if grid_mm <= finest:
        return None
    return (
        "the fused grid is %.2f mm but the finest input voxel is %.3f mm; "
        "structures thinner than the grid are lost. Lower --grid-mm to keep "
        "them, at cubic cost in memory." % (grid_mm, finest)
    )


def registration_provenance(result: RegistrationResult) -> dict[str, Any]:
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
        "registration": registration_provenance(result.registration),
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
    progress: ProgressSink | None = None,
) -> FusionResult:
    """Register, union, quantize, and publish two binary foreground masks."""
    output = volume_output(output_path)
    validate_grid_mm(grid_mm)
    started = time.time()
    reporter = Reporter(log, warn, progress, 36)
    say, step = reporter.say, reporter.step

    def register() -> RegistrationResult:
        try:
            registered = registration.rigid_register(
                fixed_mask,
                moving_mask,
                log=lambda message: say("  " + message),
            )
        except registration.RegistrationError as exc:
            raise FusionError(str(exc)) from None
        check_registration(registered, force=force)
        return registered

    registered = step("register moving to fixed", register)
    for line in registered.summary().splitlines():
        say("  " + line.strip())

    coarse_grid = grid_warning(grid_mm, [fixed_mask, moving_mask])
    if coarse_grid is not None:
        reporter.warn(coarse_grid)

    fused = step(
        "fuse occupancy fields",
        lambda: fuse_label_fields(
            [fixed_mask, moving_mask],
            [np.eye(4), registered.transform],
            grid_mm,
            binary=True,
            allow_large_volume=allow_large_volume,
            log=say,
        ),
    )
    fixed_mask = sitk.Image()
    moving_mask = sitk.Image()
    fixed_count, moving_count = fused.input_voxels
    fused_count = fused.voxels.get(1, 0)
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

    owner = [label_field_image(fused, sitk.sitkUInt8)]
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
        grid_size=fused.grid_size,
        grid_origin_mm=fused.grid_origin_mm,
        grid_direction=IDENTITY_DIRECTION,
        registration=registered,
        foreground_fixed_voxels=fixed_count,
        foreground_moving_voxels=moving_count,
        foreground_fused_voxels=fused_count,
        volume_fixed_mm3=fixed_count * voxel_mm3,
        volume_moving_mm3=moving_count * voxel_mm3,
        volume_fused_mm3=fused_count * voxel_mm3,
        seconds=time.time() - started,
        warnings=reporter.warnings,
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
    progress: ProgressSink | None = None,
) -> FusionResult:
    """Segment two intensity volumes independently and fuse their foreground."""
    volume_output(output_path)
    validate_grid_mm(grid_mm)
    validate_segmentation(preset)
    for name, threshold in (
        ("fixed_threshold", fixed_threshold),
        ("moving_threshold", moving_threshold),
    ):
        if threshold not in (None, "auto") and not math.isfinite(float(threshold)):
            raise ValueError("%s must be finite or 'auto'" % name)
    started = time.time()
    reporter = Reporter(log, warn, progress, 36)
    say, step = reporter.say, reporter.step

    reporter.warn_all(check_compatible(fixed, moving))
    reporter.warn_all(
        pipeline.threshold_warnings(fixed, preset, fixed_threshold, "--fixed-threshold")
    )
    reporter.warn_all(
        pipeline.threshold_warnings(moving, preset, moving_threshold, "--moving-threshold")
    )

    say("fixed  ID %d  %s  %s" % (fixed.id, fixed.format, fixed.source_name))
    say("moving ID %d  %s  %s" % (moving.id, moving.format, moving.source_name))
    fixed_volume = step(
        "load fixed volume",
        lambda: volume_mod.load(fixed, allow_large_volume=allow_large_volume),
    )
    reporter.warn_all(volume_mod.warnings_for(fixed_volume))
    moving_volume = step(
        "load moving volume",
        lambda: volume_mod.load(moving, allow_large_volume=allow_large_volume),
    )
    reporter.warn_all(volume_mod.warnings_for(moving_volume))

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
        progress=progress,
    )
    result.seconds = time.time() - started
    result.warnings = reporter.warnings + result.warnings
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
def common_grid(
    images: list[sitk.Image],
    transforms: list[np.ndarray],
    grid_mm: float,
    *,
    allow_large_volume: bool = False,
) -> tuple[tuple[int, int, int], np.ndarray]:
    """Axis-aligned grid covering every image after its moving->fixed transform."""
    points = np.vstack(
        [physical_corners(image, transform) for image, transform in zip(images, transforms)]
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
    sigma_mm = segment.ANTIALIAS_SIGMA_FACTOR * grid_mm
    return [
        int(math.ceil(4.0 * sigma_mm / spacing)) + 2 for spacing in image.GetSpacing()
    ]


@dataclass
class LabelFieldResult:
    """Label-preserving occupancy fusion on one isotropic grid."""

    labels: np.ndarray  # (z, y, x); uint8 for a binary union, else uint32
    grid_size: tuple[int, int, int]
    grid_origin_mm: tuple[float, float, float]
    grid_mm: float
    voxels: dict[int, int]
    #: Grid voxels each input covers above half occupancy, in input order.
    input_voxels: list[int]


def label_field_image(fused: LabelFieldResult, pixel_id: int) -> sitk.Image:
    """Fused labels as an image on their grid; releases ``fused.labels``."""
    image = sitk.Cast(sitk.GetImageFromArray(fused.labels), pixel_id)
    fused.labels = np.empty((0, 0, 0), dtype=fused.labels.dtype)
    image.SetSpacing((fused.grid_mm,) * 3)
    image.SetOrigin(fused.grid_origin_mm)
    image.SetDirection(IDENTITY_DIRECTION)
    return image


def fuse_label_fields(
    images: list[sitk.Image],
    transforms: list[np.ndarray],
    grid_mm: float,
    *,
    labels: tuple[int, ...] | None = None,
    binary: bool = False,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    progress: ProgressSink | None = None,
) -> LabelFieldResult:
    """Stream every label of every input onto one grid.

    Each input's selected labels are cropped to their bounding box, antialiased
    for the target grid, and resampled with the input's transform. A voxel is
    foreground when the union occupancy of any input exceeds 0.5, exactly as in
    binary fusion, so touching labels never open cracks between them. With
    ``binary`` foreground is ``1``; otherwise it keeps the label with the highest
    individual occupancy. Memory is one grid-sized array for a binary union and
    three for labels, regardless of how many labels or inputs are fused.
    ``input_voxels`` counts, per input, the grid voxels its selected labels
    occupy above 0.5.
    """

    def say(message: str) -> None:
        if log:
            log(message)

    input_voxels: list[int] = []
    size, origin = common_grid(
        images, transforms, grid_mm, allow_large_volume=allow_large_volume
    )
    shape_zyx = (size[2], size[1], size[0])
    union_value = np.zeros(shape_zyx, dtype=np.float32)
    if not binary:
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
        if not present:
            input_voxels.append(0)
            continue
        margin = _antialias_margin(image, grid_mm)
        inverse = registration.inverse_transform(transform)
        image_size = image.GetSize()

        def resample(members: list[int]):
            """Occupancy of ``members`` on the sub-grid it covers, or ``None``."""
            boxes = [shape.GetBoundingBox(member) for member in members]
            low = [min(box[axis] for box in boxes) for axis in range(3)]
            high = [max(box[axis] + box[axis + 3] for box in boxes) for axis in range(3)]
            start = [max(0, low[axis] - margin[axis]) for axis in range(3)]
            stop = [min(image_size[axis], high[axis] + margin[axis]) for axis in range(3)]
            extent = [stop[axis] - start[axis] for axis in range(3)]
            cropped = sitk.RegionOfInterest(discrete, extent, start)
            if len(members) == 1:
                mask = sitk.Equal(cropped, members[0])
            else:
                selected = np.isin(sitk.GetArrayViewFromImage(cropped), members)
                mask = sitk.GetImageFromArray(selected.astype(np.uint8))
                mask.CopyInformation(cropped)
            occupancy = segment.antialias_for_grid(mask, grid_mm)
            corners = physical_corners(occupancy, transform)
            low_index = np.floor((corners.min(0) - origin) / grid_mm).astype(int) - 1
            high_index = np.ceil((corners.max(0) - origin) / grid_mm).astype(int) + 2
            low_index = np.maximum(low_index, 0)
            high_index = np.minimum(high_index, np.asarray(size))
            sub_size = high_index - low_index
            if np.any(sub_size <= 0):
                return None
            field = resample_axis_aligned(
                occupancy,
                (int(sub_size[0]), int(sub_size[1]), int(sub_size[2])),
                origin + low_index * grid_mm,
                grid_mm,
                inverse,
            )
            region = (
                slice(int(low_index[2]), int(high_index[2])),
                slice(int(low_index[1]), int(high_index[1])),
                slice(int(low_index[0]), int(high_index[0])),
            )
            return region, sitk.GetArrayFromImage(field)

        def report(event: str, label: int) -> None:
            if progress is not None:
                progress({"event": event, "input": input_index, "label": label})

        if binary:
            report("label_start", 1)
        sampled = resample(present)
        covered = 0
        if sampled is not None:
            region, values = sampled
            covered = int(np.count_nonzero(values > 0.5))
            np.maximum(union_value[region], values, out=union_value[region])
        input_voxels.append(covered)
        if binary:
            report("label_end", 1)
            continue
        for label in present:
            report("label_start", label)
            sampled = resample([label])
            if sampled is not None:
                region, values = sampled
                current = best_value[region]
                better = values > current
                current[better] = values[better]
                best_label[region][better] = label
            report("label_end", label)

    foreground = union_value > 0.5
    del union_value
    if binary:
        result = foreground.view(np.uint8)
        fused_voxels = int(np.count_nonzero(result))
        voxels = {1: fused_voxels} if fused_voxels else {}
    else:
        del best_value
        result = np.where(foreground, best_label, 0).astype(np.uint32)
        del best_label, foreground
        found, counts = np.unique(result, return_counts=True)
        voxels = {int(label): int(count) for label, count in zip(found, counts) if label != 0}
    return LabelFieldResult(
        labels=result,
        grid_size=size,
        grid_origin_mm=(float(origin[0]), float(origin[1]), float(origin[2])),
        grid_mm=grid_mm,
        voxels=voxels,
        input_voxels=input_voxels,
    )
