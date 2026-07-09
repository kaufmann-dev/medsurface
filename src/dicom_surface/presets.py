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
    #: Prefer target_faces to control triangle count.
    resample_mm: float = 0.0
    #: Surface stage.
    smooth_iters: int = 20
    passband: float = 0.1
    #: Decimation to a triangle budget. 0 = off. Topology-preserving.
    target_faces: int = 0
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
        target_faces=600_000,
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
        target_faces=0,
    ),
    "bone-print": Preset(
        name="bone-print",
        description="Smooth, watertight, low-poly bone for 3D printing. Sacrifices fine detail.",
        modalities=_HU,
        threshold=350.0,
        median_mm=1.4,
        closing_mm=3.2,
        min_island_mm3=200.0,
        smooth_iters=30,
        target_faces=250_000,
        post_smooth_iters=15,
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
        target_faces=300_000,
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
        target_faces=400_000,
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
        target_faces=600_000,
        post_smooth_iters=25,
    ),
}


@dataclass(frozen=True)
class PrintProfile:
    """What a printer needs, as distinct from what the anatomy is.

    A :class:`Preset` says which tissue to extract and how finely. A profile says
    how to make that tissue survive a printer: seal the pores that would print as
    fragile holes, thicken the walls that are thinner than the machine can lay
    down, and drop the fragments too small to handle.

    Every field is a no-op at zero, so ``ANATOMICAL`` -- all zeros -- is a true
    identity. There is no "printing disabled" branch to drift out of sync with
    the enabled one.
    """

    name: str
    description: str
    #: Kernel extent (diameter) the closing is raised to, if the preset's is smaller.
    closing_mm: float = 0.0
    #: Outward growth of every surface, a radius. Inflates outer dimensions too.
    thicken_mm: float = 0.0
    #: Floor on the island filter: fragments below this are unprintable anyway.
    min_island_mm3: float = 0.0
    #: The printer's minimum feature size. Diagnostic only -- never enforced.
    min_feature_mm: float = 0.0


#: The identity profile: geometry stays faithful to the scan.
ANATOMICAL = PrintProfile(
    name="anatomical",
    description="No printability changes. Geometry faithful to the scan.",
)

PRINT_PROFILES: dict[str, PrintProfile] = {
    "anatomical": ANATOMICAL,
    "resin": PrintProfile(
        name="resin",
        description="Seals pores and thickens walls a little. Keeps fine detail.",
        closing_mm=3.2,
        thicken_mm=0.4,
        min_island_mm3=100.0,
        min_feature_mm=0.6,
    ),
    "fdm": PrintProfile(
        name="fdm",
        description="Seals hard and thickens for a nozzle. Sacrifices fine detail.",
        closing_mm=4.8,
        thicken_mm=1.2,
        min_island_mm3=200.0,
        min_feature_mm=1.2,
    ),
}


def get(name: str) -> Preset:
    try:
        return PRESETS[name]
    except KeyError:
        raise KeyError(
            "unknown preset %r; available: %s" % (name, ", ".join(sorted(PRESETS)))
        ) from None


def get_print_profile(name: str) -> PrintProfile:
    try:
        return PRINT_PROFILES[name]
    except KeyError:
        raise KeyError(
            "unknown print profile %r; available: %s"
            % (name, ", ".join(sorted(PRINT_PROFILES)))
        ) from None


def override(preset: Preset, **kwargs) -> Preset:
    """Apply non-None overrides on top of a preset."""
    changes = {k: v for k, v in kwargs.items() if v is not None}
    return replace(preset, **changes) if changes else preset


def accepts_modality(preset: Preset, modality: str) -> bool:
    return not preset.modalities or modality in preset.modalities
