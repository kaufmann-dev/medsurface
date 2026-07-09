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
        post_smooth_iters=12,
    ),
    "bone-detail": Preset(
        name="bone-detail",
        description="Bone at full mesh resolution. No decimation. Very large files.",
        modalities=_HU,
        threshold=300.0,
        median_mm=0.6,
        closing_mm=1.2,
        min_island_mm3=20.0,
        smooth_iters=15,
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
        post_smooth_iters=12,
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
