"""Named parameter bundles.

A preset is a *starting point*, not a guarantee. Hounsfield-unit thresholds are
only meaningful for CT (where 0 HU is water and -1000 HU is air); on MR, CBCT or
any modality without a calibrated rescale the same number is nonsense, so those
presets declare which modalities they accept and the CLI refuses to apply them
elsewhere.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Union

Threshold = Union[float, str]  # a number (intensity) or the string "auto"


@dataclass(frozen=True)
class Preset:
    name: str
    description: str
    #: Modalities this preset's threshold is valid for. Empty = any.
    modalities: tuple[str, ...]
    threshold: Threshold
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
    #: Surface stage.
    smooth_iters: int = 20
    #: MeshLib relaxation strength per iteration.
    smooth_force: float = 0.1
    #: MeshLib estimated surface-deviation/QEM limit in model millimetres.
    #: This is not a certified Hausdorff bound. 0 disables simplification.
    simplify_error_mm: float = 0.0
    post_smooth_iters: int = 0
    keep_largest_component: bool = True


_HU = ("CT",)

PRESETS: dict[str, Preset] = {
    "bone": Preset(
        name="bone",
        description="Cortical bone from CT. Balanced: denoises without erasing teeth or sutures.",
        modalities=_HU,
        threshold=300.0,
        median_mm=1.0,
        closing_mm=2.4,
        min_island_mm3=50.0,
        smooth_iters=20,
        simplify_error_mm=0.25,
        # 25, not 12. Measured on a head CT: 12 leaves visible slice terracing,
        # 25 removes it for 0.03 mm of extra mean displacement -- against a 0.8 mm
        # slice pitch whose stair-step amplitude is ~0.4 mm. Past 25 the returns
        # collapse (40 iterations buy 0.9 degrees for another 0.024 mm of RMS).
        post_smooth_iters=25,
    ),
    "bone-detail": Preset(
        name="bone-detail",
        description="Bone at full resolution, minimally smoothed. For measurement. Huge files.",
        modalities=_HU,
        threshold=300.0,
        median_mm=0.6,
        closing_mm=1.2,
        min_island_mm3=20.0,
        # Deliberately light. Smoothing displaces the surface, and this preset
        # exists for people who would rather see the scanner's stair-steps than
        # have a filter move their geometry: mean displacement here is ~0.02 mm
        # against ~0.085 mm for the `bone` preset.
        smooth_iters=8,
        simplify_error_mm=0.0,
    ),
    "teeth": Preset(
        name="teeth",
        description="Enamel and dense dentin only. Minimal morphology to keep cusps sharp.",
        modalities=_HU,
        threshold=1200.0,
        median_mm=0.6,
        closing_mm=0.6,
        min_island_mm3=5.0,
        keep_largest_island=False,
        keep_largest_component=False,
        smooth_iters=10,
        simplify_error_mm=0.12,
    ),
    "skin": Preset(
        name="skin",
        description="Outer skin/air boundary from CT. Closes generously to bridge hair and noise.",
        modalities=_HU,
        threshold=-300.0,
        median_mm=1.4,
        closing_mm=3.2,
        min_island_mm3=500.0,
        smooth_iters=25,
        simplify_error_mm=0.35,
        post_smooth_iters=10,
    ),
    "auto": Preset(
        name="auto",
        description="Otsu threshold. Use for MR, CBCT, ultrasound or any uncalibrated intensity.",
        modalities=(),
        threshold="auto",
        median_mm=1.0,
        closing_mm=2.0,
        min_island_mm3=50.0,
        smooth_iters=20,
        simplify_error_mm=0.25,
        post_smooth_iters=25,
    ),
}


def get(name: str) -> Preset:
    try:
        return PRESETS[name]
    except KeyError:
        raise KeyError(
            "unknown preset %r; available: %s" % (name, ", ".join(sorted(PRESETS)))
        ) from None
def override(preset: Preset, **kwargs) -> Preset:
    """Apply non-None overrides on top of a preset."""
    changes = {k: v for k, v in kwargs.items() if v is not None}
    return replace(preset, **changes) if changes else preset


def accepts_modality(preset: Preset, modality: str) -> bool:
    return not preset.modalities or modality in preset.modalities
