"""Fuse two scans of the same anatomy into one surface.

Two CT studies of one skull -- say a facial CT that stops mid-vault and a sinus CT
that stops at the maxilla -- together cover more than either alone. Fusing them
means registering one onto the other and unioning the bone.

The union happens on the **masks**, not the meshes. Boolean-unioning two shells
leaves a seam ridge wherever they disagree by a fraction of a millimetre; unioning
the solids and running marching cubes once gives a single continuous surface.

The common grid is **isotropic**. Resampling onto either scan's native grid would
force the other to inherit its slice stepping: merging a 0.3 mm sinus CT onto a
0.8 mm facial CT grid terraced the vault with the coarser scan's slice pitch.

Refusing bad input is the point of this module. Registration always returns *a*
transform -- FFT and ICP cannot fail, they can only converge somewhere useless --
so the answer must be checked, not trusted. Two unrelated scans produce a
confident-looking transform and a mesh made of two skulls stuck together at
random. The gates below exist to make that outcome an error instead.

Geometry cannot establish subject identity. Measured: a skull uniformly scaled
by 3%, which is well inside person-to-person variation, passes both geometric
gates (overlap 0.989, dice 0.675) because rigid registration parks it neatly on
top. The caller must therefore verify subject identity for every merge.
"""

from __future__ import annotations

import math
import time
from dataclasses import asdict, dataclass, field
from typing import Any, Callable

import numpy as np
import SimpleITK as sitk

from . import pipeline, registration, segment, surface
from . import volume as volume_mod
from .catalog import VolumeCandidate, same_source
from .defaults import DEFAULT_MERGE_GRID_MM, MAX_VOXELS
from .presets import Preset
from .presets import override as override_preset
from .presets import validate as validate_preset
from .registration import RegistrationResult

Logger = Callable[[str], None]

#: The surface stage is the preset's, unchanged. It is tempting to smooth a fused
#: grid more lightly, on the theory that an isotropic grid has no slice terracing
#: to remove. It has: the terracing is baked into each scan's *mask* by its own
#: slice pitch, long before anything is resampled. Resampling a 0.8 mm staircase
#: onto a 0.4 mm grid samples the staircase more finely; it does not flatten it.
#: A fused surface must be smoothed exactly as a single-scan one is, or it looks
#: visibly rougher than the scans it was built from.

# --- acceptance gates -------------------------------------------------------
# Empirically selected from four positive reconstructions/repeat studies of one
# skull and one deliberate impostor (a metal bar given the skull's identity):
#
#   pair                                   min overlap   dice    rms
#   2024 x 2023 (cross-study, cross-kernel)    0.922    0.855   0.123
#   2024 x 2023 J30s (soft kernel)             0.935    0.835   0.205
#   2024 x 2024 Hr40 3 mm                      0.969    0.852   0.187
#   2024 x 2024 sagittal reformat              0.991    0.932   0.090
#   bar phantom, unrelated anatomy             0.270    0.323   0.432
#
#: Symmetric surface overlap: the *weaker* of the two directions, each restricted
#: to the other scan's field of view. The single strongest discriminator.
MIN_SURFACE_OVERLAP = 0.60
#: Agreement of the two bone masks where both scans have data.
MIN_SHARED_FOV_DICE = 0.55
#: If the scans barely see the same space, there is nothing to register on --
#: and Dice over a few cubic centimetres proves nothing either way.
MIN_SHARED_FOV_MM3 = 20_000.0

# There is deliberately no gate on the ICP residual. It never discriminated: the
# impostor scored 0.43 mm and a 5%-oversized skull 0.41 mm, both well inside any
# plausible bound. ICP drives *some* residual down no matter what it is fitting.
# The residual is reported, because it is informative; it is not evidence.


class MergeError(RuntimeError):
    pass


@dataclass
class MergeResult:
    output_path: str
    triangles: int
    vertices: int
    bounds_mm: tuple[float, ...]
    grid_mm: float
    grid_size: tuple[int, int, int]
    registration: RegistrationResult
    volume_fixed_mm3: float
    volume_moving_mm3: float
    volume_union_mm3: float
    surface_components: int
    seconds: float
    warnings: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    quality: dict[str, Any] = field(default_factory=dict)


@dataclass
class MaskMergeResult:
    triangles: int
    vertices: int
    bounds_mm: tuple[float, ...]
    grid_size: tuple[int, int, int]
    registration: RegistrationResult
    volume_fixed_mm3: float
    volume_moving_mm3: float
    volume_union_mm3: float
    surface_components: int
    warnings: list[str]
    surface_finishing: dict[str, Any]
    quality: dict[str, Any]


RIGID_REGISTRATION_WARNING = (
    "merge uses rigid registration and is intended for matching non-deforming "
    "anatomy such as bone; anatomy that moved or deformed between acquisitions "
    "can produce a plausible but incorrect fusion"
)


def check_compatible(fixed: VolumeCandidate, moving: VolumeCandidate) -> list[str]:
    """Reject duplicate inputs and state the format-neutral identity contract."""
    if same_source(fixed, moving):
        raise MergeError("fixed and moving inputs resolve to the same volume")
    return [
        "subject identity is not verified; confirm that fixed and moving volumes "
        "show the same subject before using the fused surface",
        RIGID_REGISTRATION_WARNING,
    ]


def check_registration(result: RegistrationResult, force: bool = False) -> None:
    """Decide whether a transform describes the same anatomy, or garbage.

    Registration cannot report failure on its own: FFT always has a peak and ICP
    always converges to something. Every one of these numbers is near-perfect for
    a true match and hopeless for unrelated scans.
    """
    problems = []
    if result.shared_fov_mm3 < MIN_SHARED_FOV_MM3:
        problems.append(
            "the scans share only %.0f cm3 of imaged space (need %.0f)"
            % (result.shared_fov_mm3 / 1000.0, MIN_SHARED_FOV_MM3 / 1000.0)
        )
    if result.surface_overlap < MIN_SURFACE_OVERLAP:
        problems.append(
            "surfaces agree over only %.1f%% of the shared field of view "
            "(%.1f%% moving->fixed, %.1f%% fixed->moving; need %.0f%% both ways)"
            % (
                100 * result.surface_overlap,
                100 * result.overlap_moving_in_fixed,
                100 * result.overlap_fixed_in_moving,
                100 * MIN_SURFACE_OVERLAP,
            )
        )
    if result.shared_fov_dice < MIN_SHARED_FOV_DICE:
        problems.append(
            "foreground agreement in the shared field of view is %.3f (need %.2f)"
            % (result.shared_fov_dice, MIN_SHARED_FOV_DICE)
        )

    if not problems:
        return
    message = (
        "registration failed its quality gates, so the two scans probably do not "
        "show the same anatomy:\n  - "
        + "\n  - ".join(problems)
        + "\n"
        + result.summary()
    )
    if force:
        return
    raise MergeError(message + "\nPass --force to fuse them anyway.")


def _corners(image: sitk.Image) -> np.ndarray:
    size = np.asarray(image.GetSize()) - 1
    pts = []
    for i in (0, size[0]):
        for j in (0, size[1]):
            for k in (0, size[2]):
                pts.append(
                    image.TransformIndexToPhysicalPoint((int(i), int(j), int(k)))
                )
    return np.asarray(pts, dtype=float)


def _common_grid(
    fixed_mask,
    moving_mask,
    transform,
    grid_mm,
    *,
    allow_large_volume: bool = False,
):
    rot, trans = transform[:3, :3], transform[:3, 3]
    pts = np.vstack([_corners(fixed_mask), (rot @ _corners(moving_mask).T).T + trans])
    scaled = pts / grid_mm
    if not np.all(np.isfinite(scaled)):
        raise MergeError(
            "the fused grid coordinates overflow at %.4g mm; raise --grid-mm" % grid_mm
        )
    lo_index = np.floor(scaled.min(0)) - 2
    hi_index = np.ceil(scaled.max(0)) + 2
    planned = hi_index - lo_index + 1
    if not np.all(np.isfinite(planned)):
        raise MergeError(
            "the fused grid is too large at %.4g mm; raise --grid-mm" % grid_mm
        )
    if not allow_large_volume and np.any(planned > MAX_VOXELS):
        raise MergeError(
            "the fused grid is too large at %.4g mm; raise --grid-mm or pass "
            "--allow-large-volume" % grid_mm
        )
    size = tuple(int(value) for value in planned)
    voxels = math.prod(size)
    if voxels > MAX_VOXELS and not allow_large_volume:
        raise MergeError(
            "the fused grid would hold %s voxels at %.4g mm, above the default "
            "limit of %s; raise --grid-mm or pass --allow-large-volume to attempt "
            "it (this may exhaust memory)" % (f"{voxels:,}", grid_mm, f"{MAX_VOXELS:,}")
        )
    return size, lo_index * grid_mm


def _resample_field(image, size, origin, grid_mm, transform):
    r = sitk.ResampleImageFilter()
    r.SetOutputSpacing((grid_mm, grid_mm, grid_mm))
    r.SetSize([int(v) for v in size])
    r.SetOutputOrigin([float(v) for v in origin])
    r.SetOutputDirection(np.eye(3).ravel().tolist())
    r.SetInterpolator(sitk.sitkLinear)
    r.SetDefaultPixelValue(0.0)
    r.SetTransform(transform)
    return r.Execute(image)


def fuse_masks(
    fixed_mask: sitk.Image,
    moving_mask: sitk.Image,
    output_path: str,
    *,
    settings: pipeline.SurfaceSettings,
    grid_mm: float = DEFAULT_MERGE_GRID_MM,
    force: bool = False,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
) -> MaskMergeResult:
    """Register two binary masks, fuse their occupancy, and publish one mesh."""
    surface.validate_output_path(output_path)
    if not math.isfinite(grid_mm) or grid_mm <= 0:
        raise ValueError("grid_mm must be finite and greater than zero")
    pipeline.validate_surface_settings(settings)

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

    say("registering ...")
    try:
        reg = registration.rigid_register(
            fixed_mask,
            moving_mask,
            log=lambda message: say("  " + message),
        )
    except registration.RegistrationError as exc:
        raise MergeError(str(exc)) from None
    for line in reg.summary().splitlines():
        say("  " + line.strip() if line.startswith(" ") else "  " + line)
    check_registration(reg, force=force)

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
            reg.transform,
            grid_mm,
            allow_large_volume=allow_large_volume,
        ),
    )
    voxels = math.prod(size)
    say(
        "fused grid %s at %.2f mm isotropic (%.0f M voxels)"
        % ("x".join(str(value) for value in size), grid_mm, voxels / 1e6)
    )
    identity = sitk.Transform(3, sitk.sitkIdentity)
    inverse = registration.inverse_transform(reg.transform)

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
    fused = step(
        "fuse occupancy fields",
        lambda: sitk.Clamp(
            sitk.Maximum(fixed_field, moving_field),
            sitk.sitkFloat32,
            0.0,
            1.0,
        ),
    )

    if settings.mask_smooth_mm > 0:
        fused = step(
            "smooth fused mask occupancy field",
            lambda: segment.smooth_occupancy(fused, settings.mask_smooth_mm),
        )

    voxel_mm3 = grid_mm**3

    def measure_volumes() -> tuple[float, float, float]:
        return (
            float((sitk.GetArrayViewFromImage(fixed_field) > 0.5).sum()) * voxel_mm3,
            float((sitk.GetArrayViewFromImage(moving_field) > 0.5).sum()) * voxel_mm3,
            float((sitk.GetArrayViewFromImage(fused) > 0.5).sum()) * voxel_mm3,
        )

    fixed_volume_mm3, moving_volume_mm3, union_volume_mm3 = step(
        "measure fused volumes",
        measure_volumes,
    )
    say(
        "foreground: fixed %.0f cm3 | moving %.0f cm3 | fused %.0f cm3"
        % (
            fixed_volume_mm3 / 1000,
            moving_volume_mm3 / 1000,
            union_volume_mm3 / 1000,
        )
    )

    fused = segment.pad(fused, 1)
    affine = surface.index_to_physical(fused)
    poly = step(
        "marching cubes",
        lambda: surface.marching_cubes(fused, segment.ISO_OCCUPANCY),
    )
    if poly.topology.numValidFaces() == 0:
        raise MergeError("the fused volume produced no surface")
    say("  raw triangles %s" % f"{poly.topology.numValidFaces():,}")
    poly = step(
        "index -> fixed physical space",
        lambda: surface.transform(poly, affine),
    )

    finished = pipeline.finish_surface(
        poly,
        surface_smooth_iters=settings.surface_smooth_iters,
        simplify_error_mm=settings.simplify_error_mm,
        post_surface_smooth_iters=settings.post_surface_smooth_iters,
        keep_largest_component=settings.keep_largest_component,
        step=step,
        log=say,
    )
    for message in finished.warnings:
        add_warning(message)
    poly = finished.poly

    boundary, holes = step("check surface defects", lambda: surface.count_defects(poly))
    if boundary or holes:
        add_warning(
            "fused surface has %d boundary edge(s) and %d hole(s)" % (boundary, holes)
        )
    quality = step(
        "validate and publish mesh",
        lambda: surface.write_validated(poly, output_path),
    )
    return MaskMergeResult(
        triangles=int(poly.topology.numValidFaces()),
        vertices=int(poly.topology.numValidVerts()),
        bounds_mm=surface.bounds_mm(poly),
        grid_size=size,
        registration=reg,
        volume_fixed_mm3=fixed_volume_mm3,
        volume_moving_mm3=moving_volume_mm3,
        volume_union_mm3=union_volume_mm3,
        surface_components=finished.surface_components,
        warnings=warnings,
        surface_finishing=finished.provenance,
        quality=quality,
    )


def merge(
    fixed: VolumeCandidate,
    moving: VolumeCandidate,
    preset: Preset,
    output_path: str,
    fixed_threshold: float | str | None = None,
    moving_threshold: float | str | None = None,
    grid_mm: float = DEFAULT_MERGE_GRID_MM,
    mask_smooth_mm: float | None = None,
    surface_smooth_iters: int | None = None,
    simplify_error_mm: float | None = None,
    post_surface_smooth_iters: int | None = None,
    force: bool = False,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
) -> MergeResult:
    surface.validate_output_path(output_path)
    if not math.isfinite(grid_mm) or grid_mm <= 0:
        raise ValueError("grid_mm must be finite and greater than zero")
    preset = override_preset(
        preset,
        mask_smooth_mm=mask_smooth_mm,
        surface_smooth_iters=surface_smooth_iters,
        simplify_error_mm=simplify_error_mm,
        post_surface_smooth_iters=post_surface_smooth_iters,
    )
    validate_preset(preset)
    started = time.time()

    def say(msg: str) -> None:
        if log:
            log(msg)

    def step(label, fn):
        say("%s ..." % label)
        t = time.time()
        out = fn()
        say("  %-36s %6.1fs" % (label, time.time() - t))
        return out

    warnings: list[str] = []

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    def add_warnings(messages: list[str]) -> None:
        for message in messages:
            add_warning(message)

    add_warnings(check_compatible(fixed, moving))
    add_warnings(
        pipeline.threshold_warnings(
            fixed,
            preset,
            fixed_threshold,
            "--fixed-threshold",
        )
    )
    add_warnings(
        pipeline.threshold_warnings(
            moving,
            preset,
            moving_threshold,
            "--moving-threshold",
        )
    )
    smoothing_message = pipeline.mask_smoothing_warning(preset.mask_smooth_mm)
    if smoothing_message:
        add_warning(smoothing_message)

    say("fixed  ID %d  %s  %s" % (fixed.id, fixed.format, fixed.source_name))
    say("moving ID %d  %s  %s" % (moving.id, moving.format, moving.source_name))

    fixed_volume = step(
        "load fixed volume",
        lambda: volume_mod.load(fixed, allow_large_volume=allow_large_volume),
    )
    add_warnings(volume_mod.warnings_for(fixed_volume))
    moving_volume = step(
        "load moving volume",
        lambda: volume_mod.load(moving, allow_large_volume=allow_large_volume),
    )
    add_warnings(volume_mod.warnings_for(moving_volume))

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

    fixed_mask = step(
        "segment fixed",
        lambda: pipeline.build_mask(fixed_volume.image, preset, fixed_value),
    )
    moving_mask = step(
        "segment moving",
        lambda: pipeline.build_mask(moving_volume.image, preset, moving_value),
    )
    fused = fuse_masks(
        fixed_mask,
        moving_mask,
        output_path,
        settings=pipeline.surface_settings(preset),
        grid_mm=grid_mm,
        force=force,
        allow_large_volume=allow_large_volume,
        log=say,
        warn=warn,
    )
    warnings.extend(fused.warnings)

    provenance = {
        "fixed": pipeline.source_provenance(fixed_volume, fixed_value, fixed_source),
        "moving": pipeline.source_provenance(
            moving_volume, moving_value, moving_source
        ),
        "preset": asdict(preset),
        "grid_mm": grid_mm,
        "surface_finishing": fused.surface_finishing,
        "transform_moving_to_fixed": fused.registration.transform.tolist(),
        "rotation_deg": fused.registration.rotation_deg,
        "registration": {
            "inlier_rms_mm": fused.registration.inlier_rms_mm,
            "inlier_median_mm": fused.registration.inlier_median_mm,
            "surface_overlap": fused.registration.surface_overlap,
            "shared_fov_dice": fused.registration.shared_fov_dice,
            "shared_fov_mm3": fused.registration.shared_fov_mm3,
        },
        "coordinate_system": "SimpleITK physical space of the fixed volume",
        "forced": bool(force),
        "allow_large_volume": bool(allow_large_volume),
    }

    return MergeResult(
        output_path=output_path,
        triangles=fused.triangles,
        vertices=fused.vertices,
        bounds_mm=fused.bounds_mm,
        grid_mm=grid_mm,
        grid_size=fused.grid_size,
        registration=fused.registration,
        volume_fixed_mm3=fused.volume_fixed_mm3,
        volume_moving_mm3=fused.volume_moving_mm3,
        volume_union_mm3=fused.volume_union_mm3,
        surface_components=fused.surface_components,
        seconds=time.time() - started,
        warnings=warnings,
        provenance=provenance,
        quality=fused.quality,
    )
