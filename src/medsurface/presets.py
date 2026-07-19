"""Named segmentation and surface-finishing parameter bundles."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Union

Threshold = Union[float, str]  # a number (intensity) or the string "auto"


@dataclass(frozen=True)
class Preset:
    name: str
    description: str
    threshold: Threshold
    #: ``HU`` for calibrated CT presets; ``auto`` when computed from the image.
    threshold_unit: str
    #: Upper intensity bound; None = +inf. Useful to exclude metal.
    threshold_max: float | None = None
    #: Kernel extents in mm (total width, not radius). 0 disables the step.
    median_mm: float = 0.0
    closing_mm: float = 0.0
    opening_mm: float = 0.0
    #: Drop connected foreground blobs smaller than this physical volume.
    min_island_mm3: float = 0.0
    keep_largest_island: bool = True
    #: Isotropic voxel size for the surface grid, mm. 0 = keep the native grid.
    #: Cheap and fast, but it destroys structures thinner than the target voxel:
    #: on a head CT, 0.6 mm reopened 257 pores that the closing had sealed.
    #: Prefer simplify_error_mm to control triangle count without coarsening the grid.
    resample_mm: float = 0.0
    #: Gaussian sigma applied to the segmented occupancy mask in physical mm.
    mask_smooth_mm: float = 0.0
    #: Topology-preserving MeshLib relaxation iterations before simplification.
    surface_smooth_iters: int = 20
    #: MeshLib estimated surface-deviation/QEM limit in model millimetres.
    #: This is not a certified Hausdorff bound. 0 disables simplification.
    simplify_error_mm: float = 0.0
    #: Final topology-preserving relaxation iterations after simplification.
    post_surface_smooth_iters: int = 0
    keep_largest_component: bool = True


PRESETS: dict[str, Preset] = {
    "bone": Preset(
        name="bone",
        description="Cortical bone from CT. Balanced: denoises without erasing teeth or sutures.",
        threshold=300.0,
        threshold_unit="HU",
        median_mm=1.0,
        closing_mm=2.4,
        min_island_mm3=50.0,
        mask_smooth_mm=0.0,
        surface_smooth_iters=20,
        simplify_error_mm=0.25,
        post_surface_smooth_iters=40,
    ),
    "teeth": Preset(
        name="teeth",
        description="Enamel and dense dentin only. Minimal morphology to keep cusps sharp.",
        threshold=1200.0,
        threshold_unit="HU",
        median_mm=0.6,
        closing_mm=0.6,
        min_island_mm3=5.0,
        keep_largest_island=False,
        keep_largest_component=False,
        mask_smooth_mm=0.0,
        surface_smooth_iters=10,
        simplify_error_mm=0.12,
        post_surface_smooth_iters=0,
    ),
    "skin": Preset(
        name="skin",
        description="Outer skin/air boundary from CT. Closes generously to bridge hair and noise.",
        threshold=-300.0,
        threshold_unit="HU",
        median_mm=1.4,
        closing_mm=3.2,
        min_island_mm3=500.0,
        mask_smooth_mm=0.0,
        surface_smooth_iters=25,
        simplify_error_mm=0.35,
        post_surface_smooth_iters=10,
    ),
    "auto": Preset(
        name="auto",
        description="Otsu threshold. Use for MR, CBCT, ultrasound or any uncalibrated intensity.",
        threshold="auto",
        threshold_unit="auto",
        median_mm=1.0,
        closing_mm=2.0,
        min_island_mm3=50.0,
        mask_smooth_mm=0.0,
        surface_smooth_iters=20,
        simplify_error_mm=0.25,
        post_surface_smooth_iters=40,
    ),
}


def get(name: str) -> Preset:
    try:
        return PRESETS[name]
    except KeyError:
        raise KeyError(
            "unknown preset %r; available: %s" % (name, ", ".join(sorted(PRESETS)))
        ) from None


def validate(preset: Preset) -> None:
    """Reject invalid processing values before image or mesh allocation begins."""
    validate_segmentation(preset)
    nonnegative = {
        "resample_mm": preset.resample_mm,
        "mask_smooth_mm": preset.mask_smooth_mm,
        "simplify_error_mm": preset.simplify_error_mm,
    }
    for name, value in nonnegative.items():
        if not math.isfinite(value) or value < 0:
            raise ValueError("%s must be finite and non-negative" % name)
    if (
        not isinstance(preset.surface_smooth_iters, int)
        or preset.surface_smooth_iters < 0
    ):
        raise ValueError("surface_smooth_iters must be a non-negative integer")
    if (
        not isinstance(preset.post_surface_smooth_iters, int)
        or preset.post_surface_smooth_iters < 0
    ):
        raise ValueError("post_surface_smooth_iters must be a non-negative integer")


def validate_segmentation(preset: Preset) -> None:
    """Validate only fields used to turn intensities into a binary mask."""
    for name, value in {
        "median_mm": preset.median_mm,
        "closing_mm": preset.closing_mm,
        "opening_mm": preset.opening_mm,
        "min_island_mm3": preset.min_island_mm3,
    }.items():
        if not math.isfinite(value) or value < 0:
            raise ValueError("%s must be finite and non-negative" % name)
    for name, threshold_value in (
        ("threshold", preset.threshold),
        ("threshold_max", preset.threshold_max),
    ):
        if (
            threshold_value is not None
            and threshold_value != "auto"
            and not math.isfinite(float(threshold_value))
        ):
            raise ValueError("%s must be finite" % name)


def override(preset: Preset, **kwargs) -> Preset:
    """Apply non-None overrides on top of a preset."""
    changes = {k: v for k, v in kwargs.items() if v is not None}
    return replace(preset, **changes) if changes else preset
