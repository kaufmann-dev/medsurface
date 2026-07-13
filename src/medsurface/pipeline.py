"""End-to-end medical image volume to surface conversion."""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass
from typing import Any, Callable

import SimpleITK as sitk

from . import segment, surface, volume as volume_mod
from .catalog import DicomSource, VolumeCandidate
from .presets import Preset

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


def finish_surface(
    poly,
    *,
    smooth_iters: int,
    smooth_force: float,
    simplify_error_mm: float,
    post_smooth_iters: int,
    keep_largest_component: bool,
    step: StepRunner,
    log: Logger,
) -> SurfaceFinish:
    """Shared, intersection-safe finishing for conversion and fusion."""
    warnings = []

    poly, initial = step(
        "smooth",
        lambda: surface.smooth_safely(poly, smooth_iters, smooth_force),
    )
    if initial.initial_self_intersecting_faces:
        warnings.append(
            "initial smoothing kept %s vertices at their pre-smooth positions "
            "to prevent %d self-intersecting face(s); all %d requested iterations "
            "were retained elsewhere"
            % (
                f"{initial.protected_vertices:,}",
                initial.initial_self_intersecting_faces,
                initial.requested_iterations,
            )
        )

    surface_components = 1
    if keep_largest_component:
        poly, surface_components = step(
            "largest component",
            lambda: surface.largest_component(poly),
        )
        log("  surface shells: %d (kept 1)" % surface_components)

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

    poly, final = step(
        "post-smooth",
        lambda: surface.smooth_safely(poly, post_smooth_iters, smooth_force),
    )
    if final.initial_self_intersecting_faces:
        warnings.append(
            "post-smoothing kept %s vertices at their pre-smooth positions "
            "to prevent %d self-intersecting face(s); all %d requested iterations "
            "were retained elsewhere"
            % (
                f"{final.protected_vertices:,}",
                final.initial_self_intersecting_faces,
                final.requested_iterations,
            )
        )

    return SurfaceFinish(
        poly=poly,
        surface_components=surface_components,
        warnings=warnings,
        provenance={
            "initial_smoothing": asdict(initial),
            "decimation": asdict(decimation),
            "post_smoothing": asdict(final),
        },
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


def source_provenance(
    vol: volume_mod.Volume,
    threshold: float,
    threshold_source: str,
) -> dict[str, Any]:
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
        "hu_calibration": "verified" if volume_mod.has_calibrated_hu(candidate) else "unverified",
        "threshold": threshold,
        "threshold_source": threshold_source,
    }
    if candidate.dicom is not None:
        series = candidate.dicom
        record["dicom"] = {
            "series_uid": series.uid,
            "series_orientation_part": [series.part, series.n_parts],
            "series_number": series.series_number,
            "description": series.description,
            "convolution_kernel": series.kernel,
            "slices": series.n_slices,
        }
    return record


def convert(
    candidate: VolumeCandidate,
    preset: Preset,
    output_path: str,
    threshold: float | str | None = None,
    cap_field_of_view: bool = True,
    log: Logger | None = None,
) -> Result:
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

    vol = step("load volume", lambda: volume_mod.load(candidate))
    warnings = volume_mod.warnings_for(vol)
    warnings.extend(threshold_warnings(candidate, preset, threshold))
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

    def count_components() -> int:
        label_stats = sitk.LabelShapeStatisticsImageFilter()
        label_stats.Execute(sitk.ConnectedComponent(binary))
        return len(label_stats.GetLabels())

    labelmap_components = step("analyse components", count_components)

    touches = volume_mod.touches_boundary(binary)
    if touches and not cap_field_of_view:
        warnings.append(
            "anatomy reaches the edge of the scanned volume and --no-cap was given, "
            "so the surface will be left open there"
        )
    if touches and cap_field_of_view:
        warnings.append(
            "anatomy is truncated by the scanner's field of view; the opening has "
            "been capped flat. The missing anatomy cannot be recovered."
        )

    if cap_field_of_view:
        binary = segment.pad(binary, 1)

    if preset.resample_mm > 0:
        native = min(vol.spacing)
        if preset.resample_mm > native:
            warnings.append(
                "--resample-mm %.2f is coarser than the native %.3f mm voxel, so structures "
                "thinner than the target voxel are erased. On a head CT, resampling to 0.6 mm "
                "reopened 257 pores that the morphological closing had sealed, and terraced "
                "the vault. Use --simplify-error-mm to shed triangles without touching the grid."
                % (preset.resample_mm, native)
            )
        grid = step("resample isotropic",
                    lambda: segment.resample_isotropic(binary, preset.resample_mm, say))
        isovalue = segment.ISO_OCCUPANCY
    else:
        grid = binary
        isovalue = 0.5

    affine = surface.index_to_physical(grid)
    poly = step("marching cubes", lambda: surface.marching_cubes(grid, isovalue))
    raw_faces = int(poly.topology.numValidFaces())
    say("  raw triangles: %s" % f"{raw_faces:,}")
    if raw_faces == 0:
        raise ValueError("marching cubes produced no triangles")

    finished = finish_surface(
        poly,
        smooth_iters=preset.smooth_iters,
        smooth_force=preset.smooth_force,
        simplify_error_mm=preset.simplify_error_mm,
        post_smooth_iters=preset.post_smooth_iters,
        keep_largest_component=preset.keep_largest_component,
        step=step,
        log=say,
    )
    poly = finished.poly
    surface_components = finished.surface_components
    warnings.extend(finished.warnings)
    poly = step("index -> physical space", lambda: surface.transform(poly, affine))

    quality = step("validate and publish mesh", lambda: surface.write_validated(poly, output_path))

    provenance = {
        "input": source_provenance(vol, value, source),
        "preset": asdict(preset),
        "surface_finishing": finished.provenance,
        "capped_field_of_view": bool(touches and cap_field_of_view),
        "coordinate_system": "SimpleITK physical space of the input volume",
    }

    return Result(
        output_path=output_path,
        triangles=int(poly.topology.numValidFaces()),
        vertices=int(poly.topology.numValidVerts()),
        bounds_mm=surface.bounds_mm(poly),
        threshold_used=value,
        threshold_source=source,
        labelmap_components=labelmap_components,
        surface_components=surface_components,
        capped_field_of_view=bool(touches and cap_field_of_view),
        seconds=time.time() - t0,
        warnings=warnings,
        provenance=provenance,
        quality=quality,
    )
