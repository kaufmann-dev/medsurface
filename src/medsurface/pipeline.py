"""End-to-end medical image volume to surface extraction."""

from __future__ import annotations

import math
import time
from dataclasses import asdict, dataclass
from typing import Any, Callable, Sequence

import SimpleITK as sitk

from . import destep as destep_mod
from . import segment, surface
from . import volume as volume_mod
from .catalog import VolumeCandidate
from .defaults import DESTEP_CLAMP_WARNING_FRACTION, SURFACE_RELAX_FORCE
from .presets import DestepSettings, Preset, validate_destep
from .presets import validate as validate_preset
from .stages import ProgressSink, stage_runner

Logger = Callable[[str], None]
StepRunner = Callable[[str, Callable[[], Any]], Any]


@dataclass
class Result:
    output_path: str
    triangles: int
    vertices: int
    bounds_mm: tuple[float, ...]
    threshold_used: float
    threshold_source: str
    labelmap_components: int
    surface_components: int
    capped_field_of_view: bool
    seconds: float
    warnings: list[str]
    provenance: dict[str, Any]
    quality: dict[str, Any]


@dataclass
class SurfaceFinish:
    poly: Any
    surface_components: int
    warnings: list[str]
    provenance: dict[str, Any]


@dataclass(frozen=True)
class SurfaceSettings:
    """Mask-to-mesh controls shared by intensity and labelmap extraction."""

    resample_mm: float
    mask_smooth_mm: float
    surface_smooth_iters: int
    simplify_error_mm: float
    post_surface_smooth_iters: int
    keep_largest_component: bool
    destep: DestepSettings | None = None


@dataclass
class MaskSurfaceResult:
    triangles: int
    vertices: int
    bounds_mm: tuple[float, ...]
    labelmap_components: int
    surface_components: int
    capped_field_of_view: bool
    surface_finishing: dict[str, Any]
    quality: dict[str, Any]


def surface_settings(
    preset: Preset, *, keep_largest_component: bool | None = None
) -> SurfaceSettings:
    """Extract only the mask-to-mesh portion of an extraction preset."""
    return SurfaceSettings(
        resample_mm=preset.resample_mm,
        mask_smooth_mm=preset.mask_smooth_mm,
        surface_smooth_iters=preset.surface_smooth_iters,
        simplify_error_mm=preset.simplify_error_mm,
        post_surface_smooth_iters=preset.post_surface_smooth_iters,
        keep_largest_component=(
            preset.keep_largest_component
            if keep_largest_component is None
            else keep_largest_component
        ),
        destep=preset.destep,
    )


def validate_surface_settings(settings: SurfaceSettings) -> None:
    nonnegative = {
        "resample_mm": settings.resample_mm,
        "mask_smooth_mm": settings.mask_smooth_mm,
        "simplify_error_mm": settings.simplify_error_mm,
    }
    for name, value in nonnegative.items():
        if not math.isfinite(value) or value < 0:
            raise ValueError("%s must be finite and non-negative" % name)
    if (
        not isinstance(settings.surface_smooth_iters, int)
        or settings.surface_smooth_iters < 0
    ):
        raise ValueError("surface_smooth_iters must be a non-negative integer")
    if (
        not isinstance(settings.post_surface_smooth_iters, int)
        or settings.post_surface_smooth_iters < 0
    ):
        raise ValueError("post_surface_smooth_iters must be a non-negative integer")
    if settings.destep is not None:
        validate_destep(settings.destep)


def mask_smoothing_warning(mask_smooth_mm: float) -> str | None:
    """Explain the topology-changing tradeoff of occupancy-mask smoothing."""
    if mask_smooth_mm <= 0:
        return None
    return (
        "mask smoothing uses a Gaussian sigma of %.2f mm before meshing; "
        "it can round boundaries, merge narrow gaps, or erase structures near "
        "this scale. Use --mask-smooth-mm 0 to disable it" % mask_smooth_mm
    )


def destep_warning(destep: DestepSettings | None) -> str | None:
    """Explain that stair-step fairing deliberately moves anatomy."""
    if destep is None:
        return None
    return (
        "destep fairing (%s region, %d iterations) moves the surface up to %.2f mm; "
        "it can flatten shallow anatomy such as sutures in the faired region"
        % (destep.region, destep.iterations, destep.max_displacement_mm)
    )


def _destep_clamp_warning(stats: Any) -> str | None:
    """Flag fairing that mostly ran into its displacement clamp."""
    if not stats.safeguard.accepted:
        return None
    if stats.clamped_vertex_fraction < DESTEP_CLAMP_WARNING_FRACTION:
        return None
    return (
        "destep fairing reached the %.2f mm displacement limit on %.1f%% of "
        "vertices, which can leave clamp seams; lower --destep-iters or raise "
        "--destep-max-mm"
        % (stats.max_displacement_mm, 100.0 * stats.clamped_vertex_fraction)
    )


def _smoothing_warnings(
    stage: str,
    smoothing: surface.SmoothingSafeguard,
) -> list[str]:
    if not smoothing.accepted:
        return [
            "%s was discarded because %s; kept the valid input surface"
            % (stage, smoothing.rejection_reason)
        ]
    if not (
        smoothing.initial_self_intersecting_faces or smoothing.initial_disoriented_faces
    ):
        return []

    avoided = []
    if smoothing.initial_self_intersecting_faces:
        avoided.append(
            "%s self-intersecting face(s)"
            % f"{smoothing.initial_self_intersecting_faces:,}"
        )
    if smoothing.initial_disoriented_faces:
        avoided.append(
            "%s disoriented face(s)" % f"{smoothing.initial_disoriented_faces:,}"
        )
    return [
        "%s kept %s vertices at their input positions to prevent %s; all %d "
        "requested iterations were retained elsewhere"
        % (
            stage,
            f"{smoothing.protected_vertices:,}",
            " and ".join(avoided),
            smoothing.requested_iterations,
        )
    ]


def finish_surface(
    poly,
    *,
    surface_smooth_iters: int,
    simplify_error_mm: float,
    post_surface_smooth_iters: int,
    keep_largest_component: bool,
    step: StepRunner,
    log: Logger,
    destep: DestepSettings | None = None,
) -> SurfaceFinish:
    """Shared, intersection-safe finishing for intensity and labelmap extraction."""
    warnings = []

    poly, pre_smoothing = step(
        "relax surface",
        lambda: surface.smooth_safely(
            poly,
            surface_smooth_iters,
            SURFACE_RELAX_FORCE,
        ),
    )
    warnings.extend(_smoothing_warnings("surface relaxation", pre_smoothing))

    surface_components = surface.component_count(poly)
    if keep_largest_component:
        poly, surface_components = step(
            "largest component",
            lambda: surface.largest_component(poly),
        )
        log("  surface shells: %d (kept 1)" % surface_components)
    else:
        log("  surface shells: %d (kept all)" % surface_components)

    if simplify_error_mm > 0:
        poly, decimation = step(
            "simplify",
            lambda: surface.decimate_safely(poly, simplify_error_mm),
        )
    else:
        poly, decimation = surface.decimate_safely(poly, simplify_error_mm)
    log(
        "  simplify error %.3f mm; MeshLib estimate %.4f mm; triangles %s; "
        "repair attempts %d; protected source faces %s"
        % (
            decimation.simplify_error_mm,
            decimation.error_introduced_mm,
            f"{decimation.actual_faces:,}",
            decimation.repair_attempts,
            f"{decimation.protected_input_faces:,}",
        )
    )
    if decimation.attempted and not decimation.accepted:
        warnings.append(
            "simplification at %.3f mm was discarded because %s; kept the valid "
            "%s-triangle surface"
            % (
                decimation.simplify_error_mm,
                decimation.rejection_reason,
                f"{decimation.actual_faces:,}",
            )
        )
    elif decimation.repair_attempts:
        avoided_defects = []
        if decimation.initial_self_intersecting_faces:
            avoided_defects.append(
                "%s self-intersecting face(s)"
                % f"{decimation.initial_self_intersecting_faces:,}"
            )
        if decimation.initial_disoriented_faces:
            avoided_defects.append(
                "%s disoriented face(s)" % f"{decimation.initial_disoriented_faces:,}"
            )
        warnings.append(
            "simplification at %.3f mm protected %s source faces (%.2f%%) from collapse "
            "to prevent %s; MeshLib estimates %.4f mm "
            "introduced error and the valid result contains %s triangles"
            % (
                decimation.simplify_error_mm,
                f"{decimation.protected_input_faces:,}",
                100.0 * decimation.protected_input_face_fraction,
                " and ".join(avoided_defects),
                decimation.error_introduced_mm,
                f"{decimation.actual_faces:,}",
            )
        )

    if post_surface_smooth_iters > 0:
        poly, post_smoothing = step(
            "finish surface",
            lambda: surface.smooth_safely(
                poly,
                post_surface_smooth_iters,
                SURFACE_RELAX_FORCE,
            ),
        )
        log(
            "  post-smoothing displacement rms %.4f mm, max %.4f mm"
            % (
                post_smoothing.rms_displacement_mm,
                post_smoothing.max_displacement_mm,
            )
        )
    else:
        poly, post_smoothing = surface.smooth_safely(
            poly,
            post_surface_smooth_iters,
            SURFACE_RELAX_FORCE,
        )
    warnings.extend(
        _smoothing_warnings("post-simplification relaxation", post_smoothing)
    )

    destep_stats = None
    if destep is not None:
        settings = destep
        poly, destep_stats = step(
            "destep fairing",
            lambda: destep_mod.destep_safely(poly, settings),
        )
        log(
            "  destep region %s: faired %.1f%%, frozen %.1f%%; displacement "
            "rms %.4f mm, p95 %.4f mm, max %.4f mm; volume %+.2f%%"
            % (
                destep_stats.region,
                100.0 * destep_stats.faired_vertex_fraction,
                100.0 * destep_stats.frozen_vertex_fraction,
                destep_stats.safeguard.rms_displacement_mm,
                destep_stats.p95_displacement_mm,
                destep_stats.safeguard.max_displacement_mm,
                destep_stats.volume_change_percent,
            )
        )
        warnings.extend(_smoothing_warnings("destep fairing", destep_stats.safeguard))
        clamp_warning = _destep_clamp_warning(destep_stats)
        if clamp_warning:
            warnings.append(clamp_warning)

    return SurfaceFinish(
        poly=poly,
        surface_components=surface_components,
        warnings=warnings,
        provenance={
            "pre_smoothing": asdict(pre_smoothing),
            "decimation": asdict(decimation),
            "post_smoothing": asdict(post_smoothing),
            "destep": asdict(destep_stats) if destep_stats is not None else None,
        },
    )


def mesh_binary_mask(
    binary: sitk.Image,
    output_path: str,
    *,
    settings: SurfaceSettings,
    cap_field_of_view: bool,
    allow_large_volume: bool,
    step: StepRunner,
    log: Logger,
    warn: Logger,
    lattice_origin: Sequence[float] | None = None,
) -> MaskSurfaceResult:
    """Extract, finish, validate, and publish one already-binary mask.

    ``lattice_origin`` anchors resampling to an uncropped source's lattice.
    """
    validate_surface_settings(settings)
    surface.validate_output_path(output_path)

    def count_components() -> int:
        label_stats = sitk.LabelShapeStatisticsImageFilter()
        label_stats.Execute(sitk.ConnectedComponent(binary))
        return len(label_stats.GetLabels())

    labelmap_components = step("analyse components", count_components)
    touches = volume_mod.touches_boundary(binary)
    if touches and not cap_field_of_view:
        warn(
            "segmented foreground reaches the input boundary and --no-cap was given, "
            "so the surface will be left open there"
        )
    if touches and cap_field_of_view:
        warn(
            "segmented foreground reaches the input boundary; the opening has been "
            "capped flat. Missing anatomy cannot be recovered."
        )

    if settings.resample_mm > 0:
        native = min(binary.GetSpacing())
        if settings.resample_mm > native:
            warn(
                "--resample-mm %.2f is coarser than the native %.3f mm voxel, so "
                "structures thinner than the target voxel can be erased. Use "
                "--simplify-error-mm to reduce triangles without changing the grid."
                % (settings.resample_mm, native)
            )
        grid = step(
            "resample isotropic",
            lambda: segment.resample_isotropic(
                binary,
                settings.resample_mm,
                log,
                pad_border=settings.mask_smooth_mm <= 0 and cap_field_of_view,
                allow_large_volume=allow_large_volume,
                lattice_origin=lattice_origin,
            ),
        )
        isovalue = segment.ISO_OCCUPANCY
    else:
        grid = binary
        isovalue = 0.5

    if settings.mask_smooth_mm > 0:
        grid = step(
            "smooth mask occupancy field",
            lambda: segment.smooth_occupancy(grid, settings.mask_smooth_mm),
        )
        isovalue = segment.ISO_OCCUPANCY

    if cap_field_of_view and not (
        settings.resample_mm > 0 and settings.mask_smooth_mm <= 0
    ):
        grid = segment.pad(grid, 1)

    affine = surface.index_to_physical(grid)
    poly = step("marching cubes", lambda: surface.marching_cubes(grid, isovalue))
    raw_faces = int(poly.topology.numValidFaces())
    log("  raw triangles: %s" % f"{raw_faces:,}")
    if raw_faces == 0:
        raise ValueError("marching cubes produced no triangles")
    poly = step("index -> physical space", lambda: surface.transform(poly, affine))

    finished = finish_surface(
        poly,
        surface_smooth_iters=settings.surface_smooth_iters,
        simplify_error_mm=settings.simplify_error_mm,
        post_surface_smooth_iters=settings.post_surface_smooth_iters,
        keep_largest_component=settings.keep_largest_component,
        step=step,
        log=log,
        destep=settings.destep,
    )
    for message in finished.warnings:
        warn(message)
    poly = finished.poly
    quality = step(
        "validate and publish mesh",
        lambda: surface.write_validated(poly, output_path),
    )

    return MaskSurfaceResult(
        triangles=int(poly.topology.numValidFaces()),
        vertices=int(poly.topology.numValidVerts()),
        bounds_mm=surface.bounds_mm(poly),
        labelmap_components=labelmap_components,
        surface_components=finished.surface_components,
        capped_field_of_view=bool(touches and cap_field_of_view),
        surface_finishing=finished.provenance,
        quality=quality,
    )


def build_mask(
    image: sitk.Image,
    preset: Preset,
    threshold: float,
    log: Callable[[str], None] | None = None,
) -> sitk.Image:
    """Threshold and clean a volume into a binary bone mask.

    Shared by ``extract`` and ``fuse`` so both independently produce the same
    cleaned foreground mask from one intensity volume.

    """
    binary = segment.binarize(image, threshold, preset.threshold_max)
    binary = segment.islands(
        binary, preset.keep_largest_island, preset.min_island_mm3, log
    )
    binary = segment.median(binary, preset.median_mm, log)
    if preset.opening_mm > 0:
        binary = segment.opening(binary, preset.opening_mm, log)
    binary = segment.closing(binary, preset.closing_mm, log)
    binary = segment.islands(
        binary, preset.keep_largest_island, preset.min_island_mm3, log
    )
    return binary


def resolve_threshold(
    image: sitk.Image, preset: Preset, override: float | str | None
) -> tuple[float, str]:
    if override == "auto":
        return segment.auto_threshold(image), "otsu"
    if override is not None:
        return float(override), "explicit"

    if preset.threshold == "auto":
        return segment.auto_threshold(image), "otsu"
    return float(preset.threshold), "preset:%s" % preset.name


def threshold_warnings(
    candidate: VolumeCandidate,
    preset: Preset,
    override: float | str | None,
    option_name: str = "--threshold",
) -> list[str]:
    if override is not None or preset.threshold_unit != "HU":
        return []
    if volume_mod.has_calibrated_hu(candidate):
        return []
    return [
        "preset %r applies a %g HU threshold, but HU calibration cannot be verified "
        "for %s; the value will be applied to stored intensities. Use %s "
        "with an intentional stored value or --preset auto for Otsu."
        % (preset.name, float(preset.threshold), candidate.source_name, option_name)
    ]


def source_provenance(
    vol: volume_mod.Volume,
    threshold: float,
    threshold_source: str,
) -> dict[str, Any]:
    record = volume_mod.provenance_for(vol)
    record.update(
        {
            "hu_calibration": (
                "verified"
                if volume_mod.has_calibrated_hu(vol.candidate)
                else "unverified"
            ),
            "threshold": threshold,
            "threshold_source": threshold_source,
        }
    )
    return record


def extract(
    candidate: VolumeCandidate,
    preset: Preset,
    output_path: str,
    threshold: float | str | None = None,
    cap_field_of_view: bool = True,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
    progress: ProgressSink | None = None,
) -> Result:
    validate_preset(preset)
    surface.validate_output_path(output_path)
    t0 = time.time()

    def say(msg: str) -> None:
        if log:
            log(msg)

    step = stage_runner(say, progress, 34)

    warnings: list[str] = []

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    def add_warnings(messages: list[str]) -> None:
        for message in messages:
            add_warning(message)

    add_warnings(threshold_warnings(candidate, preset, threshold))
    for message in (
        mask_smoothing_warning(preset.mask_smooth_mm),
        destep_warning(preset.destep),
    ):
        if message:
            add_warning(message)
    vol = step(
        "load volume",
        lambda: volume_mod.load(candidate, allow_large_volume=allow_large_volume),
    )
    add_warnings(volume_mod.warnings_for(vol))
    lo, hi = step("measure intensity range", vol.intensity_range)
    say(
        "volume %s  spacing %s mm  intensity %.0f..%.0f"
        % (
            "x".join(str(v) for v in vol.size),
            " x ".join("%.3f" % s for s in vol.spacing),
            lo,
            hi,
        )
    )

    value, source = step(
        "resolve threshold",
        lambda: resolve_threshold(vol.image, preset, threshold),
    )
    unit = "HU" if volume_mod.has_calibrated_hu(candidate) else "intensity"
    say("threshold %.1f %s (%s)" % (value, unit, source))
    if value > hi:
        raise ValueError(
            "threshold %.1f is above the brightest voxel (%.1f); nothing would be "
            "segmented" % (value, hi)
        )

    binary = step("segment", lambda: build_mask(vol.image, preset, value, say))
    meshed = mesh_binary_mask(
        binary,
        output_path,
        settings=surface_settings(preset),
        cap_field_of_view=cap_field_of_view,
        allow_large_volume=allow_large_volume,
        step=step,
        log=say,
        warn=add_warning,
    )

    provenance = {
        "input": source_provenance(vol, value, source),
        "preset": asdict(preset),
        "surface_finishing": meshed.surface_finishing,
        "capped_field_of_view": meshed.capped_field_of_view,
        "allow_large_volume": bool(allow_large_volume),
        "coordinate_system": "SimpleITK physical space of the input volume",
    }

    return Result(
        output_path=output_path,
        triangles=meshed.triangles,
        vertices=meshed.vertices,
        bounds_mm=meshed.bounds_mm,
        threshold_used=value,
        threshold_source=source,
        labelmap_components=meshed.labelmap_components,
        surface_components=meshed.surface_components,
        capped_field_of_view=meshed.capped_field_of_view,
        seconds=time.time() - t0,
        warnings=warnings,
        provenance=provenance,
        quality=meshed.quality,
    )
