"""End-to-end medical image volume to surface conversion."""

from __future__ import annotations

import math
import time
from dataclasses import asdict, dataclass
from typing import Any, Callable

import SimpleITK as sitk

from . import segment, surface
from . import volume as volume_mod
from .catalog import DicomSource, VolumeCandidate
from .defaults import SURFACE_RELAX_FORCE, SURFACE_RELAX_ITERS
from .presets import Preset
from .presets import validate as validate_preset

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
    """Mask-to-mesh controls shared by conversion and fusion workflows."""

    resample_mm: float
    smooth_mm: float
    simplify_error_mm: float
    keep_largest_component: bool


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


def surface_settings(preset: Preset, *, keep_largest_component: bool | None = None) -> SurfaceSettings:
    """Extract only the mask-to-mesh portion of a conversion preset."""
    return SurfaceSettings(
        resample_mm=preset.resample_mm,
        smooth_mm=preset.smooth_mm,
        simplify_error_mm=preset.simplify_error_mm,
        keep_largest_component=(
            preset.keep_largest_component
            if keep_largest_component is None
            else keep_largest_component
        ),
    )


def validate_surface_settings(settings: SurfaceSettings) -> None:
    nonnegative = {
        "resample_mm": settings.resample_mm,
        "smooth_mm": settings.smooth_mm,
        "simplify_error_mm": settings.simplify_error_mm,
    }
    for name, value in nonnegative.items():
        if not math.isfinite(value) or value < 0:
            raise ValueError("%s must be finite and non-negative" % name)


def smoothing_warning(smooth_mm: float) -> str | None:
    """Explain the intentional geometry tradeoff of physical smoothing."""
    if smooth_mm <= 0:
        return None
    return (
        "surface smoothing uses a Gaussian sigma of %.2f mm before meshing; "
        "it can round boundaries, merge narrow gaps, or erase structures near "
        "this scale. Use --smooth-mm 0 to disable it" % smooth_mm
    )


def finish_surface(
    poly,
    *,
    relax_surface: bool,
    simplify_error_mm: float,
    keep_largest_component: bool,
    step: StepRunner,
    log: Logger,
) -> SurfaceFinish:
    """Shared, intersection-safe finishing for conversion and fusion."""
    warnings = []

    poly, smoothing = step(
        "relax surface",
        lambda: surface.smooth_safely(
            poly,
            SURFACE_RELAX_ITERS if relax_surface else 0,
            SURFACE_RELAX_FORCE,
        ),
    )
    if smoothing.initial_self_intersecting_faces:
        warnings.append(
            "surface relaxation kept %s vertices at their pre-relaxation positions "
            "to prevent %d self-intersecting face(s); all %d requested iterations "
            "were retained elsewhere"
            % (
                f"{smoothing.protected_vertices:,}",
                smoothing.initial_self_intersecting_faces,
                smoothing.requested_iterations,
            )
        )

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
        warnings.append(
            "simplification at %.3f mm protected %s source faces (%.2f%%) from collapse "
            "to prevent %s self-intersecting candidate faces; MeshLib estimates %.4f mm "
            "introduced error and the collision-free result contains %s triangles"
            % (
                decimation.simplify_error_mm,
                f"{decimation.protected_input_faces:,}",
                100.0 * decimation.protected_input_face_fraction,
                f"{decimation.initial_self_intersecting_faces:,}",
                decimation.error_introduced_mm,
                f"{decimation.actual_faces:,}",
            )
        )

    return SurfaceFinish(
        poly=poly,
        surface_components=surface_components,
        warnings=warnings,
        provenance={
            "smoothing": asdict(smoothing),
            "decimation": asdict(decimation),
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
) -> MaskSurfaceResult:
    """Extract, finish, validate, and publish one already-binary mask."""
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
                pad_border=settings.smooth_mm <= 0 and cap_field_of_view,
                allow_large_volume=allow_large_volume,
            ),
        )
        isovalue = segment.ISO_OCCUPANCY
    else:
        grid = binary
        isovalue = 0.5

    if settings.smooth_mm > 0:
        grid = step(
            "smooth occupancy field",
            lambda: segment.smooth_occupancy(grid, settings.smooth_mm),
        )
        isovalue = segment.ISO_OCCUPANCY

    if cap_field_of_view and not (
        settings.resample_mm > 0 and settings.smooth_mm <= 0
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
        relax_surface=settings.smooth_mm > 0,
        simplify_error_mm=settings.simplify_error_mm,
        keep_largest_component=settings.keep_largest_component,
        step=step,
        log=log,
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


def build_mask(image: sitk.Image, preset: Preset, threshold: float,
               log: Callable[[str], None] | None = None) -> sitk.Image:
    """Threshold and clean a volume into a binary bone mask.

    Shared by ``convert`` and ``merge`` so the two can never drift apart: a fused
    surface must be built from exactly the mask a single-scan conversion would
    have produced.

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


def volume_provenance(vol: volume_mod.Volume) -> dict[str, Any]:
    """Format-neutral provenance shared by intensity and labelmap inputs."""
    candidate = vol.candidate
    if isinstance(candidate.source, DicomSource):
        path = str(candidate.source.catalog_path)
    else:
        path = str(candidate.source.path)
    record: dict[str, Any] = {
        "id": candidate.id,
        "format": candidate.format,
        "path": path,
        "source": candidate.source_name,
        "size": list(vol.size),
        "pixel_type": vol.image.GetPixelIDTypeAsString(),
        "spacing_mm": list(vol.spacing),
        "origin_mm": list(vol.image.GetOrigin()),
        "direction": list(vol.image.GetDirection()),
        "modality": candidate.modality,
        "description": candidate.description,
        "plane": candidate.plane,
    }
    if candidate.dicom is not None:
        series = candidate.dicom
        record["dicom"] = {
            "series_uid": series.uid,
            "series_orientation_part": [series.part, series.n_parts],
            "series_number": series.series_number,
            "description": series.description,
            "convolution_kernel": list(series.kernel_values),
            "slices": series.n_slices,
            "image_type": list(series.image_type),
            "rescale_type": series.rescale_type,
            "multi_energy_ct_acquisition": series.multi_energy_ct_acquisition,
            "hu_calibration_consistent": series.hu_calibration_consistent,
        }
    return record


def source_provenance(
    vol: volume_mod.Volume,
    threshold: float,
    threshold_source: str,
) -> dict[str, Any]:
    record = volume_provenance(vol)
    record.update(
        {
            "hu_calibration": (
                "verified" if volume_mod.has_calibrated_hu(vol.candidate) else "unverified"
            ),
            "threshold": threshold,
            "threshold_source": threshold_source,
        }
    )
    return record


def convert(
    candidate: VolumeCandidate,
    preset: Preset,
    output_path: str,
    threshold: float | str | None = None,
    cap_field_of_view: bool = True,
    allow_large_volume: bool = False,
    log: Logger | None = None,
    warn: Logger | None = None,
) -> Result:
    validate_preset(preset)
    surface.validate_output_path(output_path)
    t0 = time.time()

    def say(msg: str) -> None:
        if log:
            log(msg)

    def step(msg: str, fn):
        say("%s ..." % msg)
        t = time.time()
        out = fn()
        say("  %-34s %6.1fs" % (msg, time.time() - t))
        return out

    warnings: list[str] = []

    def add_warning(message: str) -> None:
        warnings.append(message)
        if warn:
            warn(message)

    def add_warnings(messages: list[str]) -> None:
        for message in messages:
            add_warning(message)

    add_warnings(threshold_warnings(candidate, preset, threshold))
    smoothing_message = smoothing_warning(preset.smooth_mm)
    if smoothing_message:
        add_warning(smoothing_message)
    vol = step(
        "load volume",
        lambda: volume_mod.load(candidate, allow_large_volume=allow_large_volume),
    )
    add_warnings(volume_mod.warnings_for(vol))
    lo, hi = step("measure intensity range", vol.intensity_range)
    say("volume %s  spacing %s mm  intensity %.0f..%.0f"
        % ("x".join(str(v) for v in vol.size),
           " x ".join("%.3f" % s for s in vol.spacing), lo, hi))

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
