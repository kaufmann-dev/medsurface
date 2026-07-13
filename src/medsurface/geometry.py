"""Voxel geometry helpers.

The single most important function here is :func:`kernel_radius_voxels`. Medical
volumes are almost never isotropic -- a head CT is typically ~0.3 mm in-plane and
0.8-1.0 mm between slices -- so a filter size expressed in millimetres maps to a
*different* number of voxels on every axis.

Getting the rounding wrong here is a silent data-quality bug: rounding a 1.0 mm
request against 0.8 mm slices yields a 3-voxel kernel spanning 2.4 mm, i.e. 2.4x
the requested size, which erases real anatomy while appearing to honour the
parameter. We therefore define the contract as:

    the realised kernel extent never exceeds the requested extent

and pick the largest radius satisfying ``(2r + 1) * spacing <= mm``.
"""

from __future__ import annotations

import math
from typing import Sequence

#: Guards against binary-float artefacts in the exact-fit case. 2.4 / 0.8 is
#: 2.9999999999999996, not 3.0, which would floor a legitimate radius of 1 down
#: to 0 and silently drop the filter on that axis.
_EPS = 1e-9


def kernel_radius_voxels(mm: float, spacing: Sequence[float]) -> list[int]:
    """Half-width, in voxels per axis, of a kernel of ``mm`` total extent.

    Whenever the radius is non-zero the realised extent ``(2r + 1) * spacing`` is
    guaranteed ``<= mm``. A radius of 0 means the kernel is a single voxel on that
    axis -- an identity filter -- which is the correct answer when the requested
    extent is smaller than one voxel.
    """
    if mm <= 0:
        return [0 for _ in spacing]
    radii = []
    for s in spacing:
        if s <= 0:
            raise ValueError("non-positive spacing: %r" % (spacing,))
        r = math.floor((mm / s - 1.0) / 2.0 + _EPS)
        radii.append(max(0, int(r)))
    return radii


def kernel_extent_mm(radii: Sequence[int], spacing: Sequence[float]) -> list[float]:
    """Realised kernel extent in mm, given per-axis radii."""
    return [(2 * r + 1) * s for r, s in zip(radii, spacing)]


def voxel_volume_mm3(spacing: Sequence[float]) -> float:
    v = 1.0
    for s in spacing:
        v *= s
    return v


def mm3_to_voxels(mm3: float, spacing: Sequence[float]) -> int:
    """Convert a physical volume to a voxel count for the given spacing."""
    if mm3 <= 0:
        return 0
    return max(1, int(round(mm3 / voxel_volume_mm3(spacing))))


def is_anisotropic(spacing: Sequence[float], tol: float = 1.05) -> bool:
    lo, hi = min(spacing), max(spacing)
    return hi / lo > tol
