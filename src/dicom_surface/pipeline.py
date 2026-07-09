"""End-to-end DICOM -> surface conversion."""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass
from typing import Any, Callable

import SimpleITK as sitk

from . import segment, surface, volume as volume_mod
from .presets import Preset, accepts_modality
from .series import Series

Logger = Callable[[str], None]


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


class ModalityMismatch(ValueError):
    pass


def build_mask(image: sitk.Image, preset: Preset, threshold: float,
               log: Callable[[str], None] | None = None) -> sitk.Image:
    """Threshold and clean a volume into a binary bone mask.

    Shared by ``convert`` and ``merge`` so the two can never drift apart: a fused
    surface must be built from exactly the mask a single-scan conversion would
    have produced.
    """
    binary = segment.binarize(image, threshold, preset.threshold_max)
    binary = segment.islands(binary, preset.keep_largest_island, preset.min_island_mm3, log)
    binary = segment.median(binary, preset.median_mm, log)
    if preset.opening_mm > 0:
        binary = segment.opening(binary, preset.opening_mm, log)
    binary = segment.closing(binary, preset.closing_mm, log)
    binary = segment.islands(binary, preset.keep_largest_island, preset.min_island_mm3, log)
    return binary


def resolve_threshold(
    image: sitk.Image, series: Series, preset: Preset, override: float | None
) -> tuple[float, str]:
    if override is not None:
        return float(override), "explicit"

    if preset.threshold == "auto":
        return segment.auto_threshold(image, series.modality), "otsu"

    if not accepts_modality(preset, series.modality):
        raise ModalityMismatch(
            "preset %r uses a Hounsfield-unit threshold (%g HU) which is only "
            "meaningful for CT, but this series is %s. Use --preset auto for an "
            "Otsu threshold, or pass --threshold with an explicit intensity."
            % (preset.name, preset.threshold, series.modality)
        )
    return float(preset.threshold), "preset:%s" % preset.name


def convert(
    series: Series,
    preset: Preset,
    output_path: str,
    threshold: float | None = None,
    cap_field_of_view: bool = True,
    log: Logger | None = None,
) -> Result:
    t0 = time.time()

    def say(msg: str) -> None:
        if log:
            log(msg)

    def step(msg: str, fn):
        t = time.time()
        out = fn()
        say("  %-34s %6.1fs" % (msg, time.time() - t))
        return out

    vol = volume_mod.load(series)
    warnings = volume_mod.warnings_for(vol)
    lo, hi = vol.intensity_range()
    say("volume %s  spacing %s mm  intensity %.0f..%.0f"
        % ("x".join(str(v) for v in vol.size),
           " x ".join("%.3f" % s for s in vol.spacing), lo, hi))

    value, source = resolve_threshold(vol.image, series, preset, threshold)
    unit = "HU" if volume_mod.has_calibrated_hu(series) else "intensity"
    say("threshold %.1f %s (%s)" % (value, unit, source))
    if value > hi:
        raise ValueError(
            "threshold %.1f is above the brightest voxel (%.1f); nothing would be "
            "segmented" % (value, hi)
        )

    binary = step("segment", lambda: build_mask(vol.image, preset, value, say))

    label_stats = sitk.LabelShapeStatisticsImageFilter()
    label_stats.Execute(sitk.ConnectedComponent(binary))
    labelmap_components = len(label_stats.GetLabels())

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
                "the vault. Use --target-faces to shed triangles without touching geometry."
                % (preset.resample_mm, native)
            )
        grid = step("resample isotropic",
                    lambda: segment.resample_isotropic(binary, preset.resample_mm, say))
        isovalue = segment.ISO_OCCUPANCY
    else:
        grid = binary
        isovalue = 0.5

    affine = surface.index_to_physical(grid)
    vtk_img = surface.to_vtk_image(grid)

    poly = step("marching cubes", lambda: surface.marching_cubes(vtk_img, isovalue))
    say("  raw triangles: %s" % f"{poly.GetNumberOfPolys():,}")
    if poly.GetNumberOfPolys() == 0:
        raise ValueError("marching cubes produced no triangles")

    poly = step("smooth", lambda: surface.smooth(poly, preset.smooth_iters, preset.passband))

    surface_components = 1
    if preset.keep_largest_component:
        poly, surface_components = step("largest component",
                                        lambda: surface.largest_component(poly))
        say("  surface shells: %d (kept 1)" % surface_components)

    if preset.target_faces > 0 and poly.GetNumberOfPolys() > preset.target_faces:
        before = surface.count_defects(poly)
        poly = step("decimate", lambda: surface.decimate(poly, preset.target_faces))
        after = surface.count_defects(poly)
        if after > before:
            warnings.append(
                "decimation to %s triangles introduced %d boundary and %d non-manifold "
                "edge(s); the mesh is no longer watertight. Raise --target-faces, or run "
                "'dicom-surface repair'."
                % (f"{preset.target_faces:,}", after[0] - before[0], after[1] - before[1])
            )

    poly = step("post-smooth",
                lambda: surface.smooth(poly, preset.post_smooth_iters, preset.passband))
    poly = step("index -> patient space (LPS)", lambda: surface.transform(poly, affine))
    poly = step("normals", lambda: surface.compute_normals(poly))

    surface.write(poly, output_path)

    provenance = {
        "series_uid": series.uid,
        "series_ident": series.ident,
        "series_orientation_part": [series.part, series.n_parts],
        "series_description": series.description,
        "series_number": series.series_number,
        "modality": series.modality,
        "convolution_kernel": series.kernel,
        "slices": series.n_slices,
        "spacing_mm": list(vol.spacing),
        "preset": asdict(preset),
        "threshold": value,
        "threshold_source": source,
        "capped_field_of_view": bool(touches and cap_field_of_view),
        "coordinate_system": "LPS (DICOM patient space)",
    }

    return Result(
        output_path=output_path,
        triangles=int(poly.GetNumberOfPolys()),
        vertices=int(poly.GetNumberOfPoints()),
        bounds_mm=surface.bounds_mm(poly),
        threshold_used=value,
        threshold_source=source,
        labelmap_components=labelmap_components,
        surface_components=surface_components,
        capped_field_of_view=bool(touches and cap_field_of_view),
        seconds=time.time() - t0,
        warnings=warnings,
        provenance=provenance,
    )
