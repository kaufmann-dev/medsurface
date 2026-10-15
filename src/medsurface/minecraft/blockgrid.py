"""Optional cubic-block resampling and anatomical orientation of labelmaps.

Without these options one labelmap voxel becomes one Minecraft block and the
index axes are neither flipped nor resampled. Medical scans often have
anisotropic voxels (for example 0.7 x 0.7 x 5 mm), which makes such a world look
squashed, and their index axes can mirror the anatomy once mapped onto
Minecraft's axes. The functions here build a new label grid of cubic blocks with
matching SimpleITK geometry, which the rest of the pipeline treats as the input.

Physical coordinates are ITK's LPS millimetres. Arrays are indexed ``[i, j, k]``
here (SimpleITK index order); SimpleITK arrays are ``[k, j, i]``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
import SimpleITK as sitk

#: Largest supersampling factor per axis when a block covers many voxels.
MAX_SUPERSAMPLING = 6
#: Automatic block sizes, in millimetres, tried from finest to coarsest.
AUTO_BLOCK_LADDER_MM = (
    0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0,
)
#: Automatic block size keeps the vertical axis at or below this many blocks.
AUTO_MAX_VERTICAL_BLOCKS = 320.0
#: Automatic block size keeps the whole grid at or below this many blocks.
AUTO_MAX_BLOCKS = 64_000_000.0
#: NIfTI-1 stores dimensions as signed 16-bit integers.
MAX_AXIS_BLOCKS = 32767

#: Target index axes for anatomical output in LPS: +i = Right, +j = Posterior,
#: +k = Superior. The world places i on Minecraft X, j on Z, and k on Y, so the
#: model stands upright and is not mirrored.
ANATOMICAL_FRAME = np.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])


@dataclass(frozen=True)
class TransformOptions:
    #: Block size in mm, ``"auto"``, or ``None`` to keep the voxel grid.
    block_mm: float | str | None = None
    orientation: str = "voxel"
    fill_threshold: float = 0.5
    crop: bool = False

    @property
    def is_identity(self) -> bool:
        return self.block_mm is None and self.orientation == "voxel" and not self.crop


@dataclass(frozen=True)
class TransformRecord:
    """How a world or resampled volume was derived from its source grid."""

    orientation: str
    block_mm: float | None
    fill_threshold: float | None
    cropped: bool
    source_size: tuple[int, int, int]
    source_spacing: tuple[float, float, float]
    size: tuple[int, int, int]
    spacing: tuple[float, float, float]
    origin: tuple[float, float, float]
    direction: tuple[float, ...]

    def as_dict(self) -> dict[str, Any]:
        return {
            "orientation": self.orientation,
            "block_mm": self.block_mm,
            "fill_threshold": self.fill_threshold,
            "cropped": self.cropped,
            "source_size": list(self.source_size),
            "source_spacing": list(self.source_spacing),
            "size": list(self.size),
            "spacing": list(self.spacing),
            "origin": list(self.origin),
            "direction": list(self.direction),
        }


def parse_block_mm(value: str) -> float | str:
    """Parse ``--block-mm``: a positive number of millimetres or ``auto``."""
    if value.strip().lower() == "auto":
        return "auto"
    try:
        millimetres = float(value)
    except ValueError:
        raise ValueError("--block-mm expects a size in mm or 'auto', got %r" % value) from None
    if not math.isfinite(millimetres) or millimetres <= 0:
        raise ValueError("--block-mm must be a positive number of millimetres")
    return millimetres


def validate_options(options: TransformOptions) -> None:
    if options.orientation not in ("voxel", "anatomical"):
        raise ValueError("orientation must be voxel or anatomical")
    if not (math.isfinite(options.fill_threshold) and 0 < options.fill_threshold <= 1):
        raise ValueError("--fill-threshold must be greater than 0 and at most 1")
    if options.block_mm is None and options.crop:
        raise ValueError("--crop requires --block-mm")


# ---------------------------------------------------------------- affines
def source_affine(image: sitk.Image) -> np.ndarray:
    """3x4 affine from continuous index ``(i, j, k)`` to LPS millimetres."""
    direction = np.array(image.GetDirection(), dtype=float).reshape(3, 3)
    affine = np.zeros((3, 4))
    affine[:, :3] = direction * np.array(image.GetSpacing(), dtype=float)
    affine[:, 3] = image.GetOrigin()
    if not np.all(np.isfinite(affine)) or abs(np.linalg.det(affine[:, :3])) < 1e-12:
        raise ValueError("the image orientation matrix is singular or not finite")
    return affine


def _invert(affine: np.ndarray) -> np.ndarray:
    inverse = np.linalg.inv(affine[:, :3])
    return np.hstack([inverse, -(inverse @ affine[:, 3])[:, None]])


def _compose(outer: np.ndarray, inner: np.ndarray) -> np.ndarray:
    result = outer[:, :3] @ inner
    result[:, 3] += outer[:, 3]
    return result


def _voxel_frame(source: np.ndarray) -> np.ndarray:
    """Unit directions of the source index axes (as columns)."""
    linear = source[:, :3]
    return linear / np.linalg.norm(linear, axis=0)


def _image(voxels_ijk: np.ndarray, affine: np.ndarray, dtype: np.dtype) -> sitk.Image:
    """SimpleITK image from ``[i, j, k]`` labels and a 3x4 index-to-LPS affine."""
    spacing = np.linalg.norm(affine[:, :3], axis=0)
    image = sitk.GetImageFromArray(np.ascontiguousarray(voxels_ijk.transpose(2, 1, 0).astype(dtype)))
    image.SetSpacing(tuple(float(value) for value in spacing))
    image.SetOrigin(tuple(float(value) for value in affine[:, 3]))
    image.SetDirection(tuple(float(value) for value in (affine[:, :3] / spacing).reshape(-1)))
    return image


# ------------------------------------------------------------ target grid
def _source_bounds(labels_kji: np.ndarray, crop: bool) -> tuple[np.ndarray, np.ndarray]:
    """Continuous source index box: the volume, or the labelled box plus one voxel."""
    size = np.array(labels_kji.shape[::-1], dtype=float)
    if not crop:
        return np.full(3, -0.5), size - 0.5
    foreground = labels_kji != 0
    if not foreground.any():
        raise ValueError("--crop found no nonzero labels to crop to")
    low = np.zeros(3)
    high = np.zeros(3)
    for axis, others in ((0, (0, 1)), (1, (0, 2)), (2, (1, 2))):
        present = np.flatnonzero(foreground.any(axis=others))
        low[axis] = present[0]
        high[axis] = present[-1]
    return low - 1.5, high + 1.5


def _frame_extent(
    source: np.ndarray, frame: np.ndarray, low: np.ndarray, high: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    corners = np.array(
        [[(high if corner >> axis & 1 else low)[axis] for axis in range(3)] for corner in range(8)]
    )
    points = corners @ source[:, :3].T + source[:, 3]
    coordinates = points @ frame
    return coordinates.min(axis=0), coordinates.max(axis=0)


def _auto_block_size(
    spacing: tuple[float, ...], extent: tuple[np.ndarray, np.ndarray], orientation: str
) -> float:
    finest = min(spacing)
    lengths = extent[1] - extent[0]
    # Anatomical output is always vertical along k; voxel output may use any
    # axis as vertical, so the tallest axis is a safe bound.
    vertical = lengths[2] if orientation == "anatomical" else float(lengths.max())
    for candidate in AUTO_BLOCK_LADDER_MM:
        if candidate + 1e-9 < finest:
            continue
        blocks = float(np.prod(np.ceil(lengths / candidate)))
        if vertical / candidate <= AUTO_MAX_VERTICAL_BLOCKS and blocks <= AUTO_MAX_BLOCKS:
            return candidate
    return AUTO_BLOCK_LADDER_MM[-1]


def _target_grid(
    frame: np.ndarray, extent: tuple[np.ndarray, np.ndarray], block_mm: float
) -> tuple[tuple[int, int, int], np.ndarray]:
    size = []
    first_center = np.zeros(3)
    for axis in range(3):
        length = extent[1][axis] - extent[0][axis]
        count = max(math.ceil(length / block_mm - 1e-6), 1)
        if count > MAX_AXIS_BLOCKS:
            raise ValueError(
                "a %g mm block grid needs %d blocks on one axis; at most %d are supported"
                % (block_mm, count, MAX_AXIS_BLOCKS)
            )
        size.append(count)
        first_center[axis] = extent[0][axis] + (length - count * block_mm) / 2 + block_mm / 2
    affine = np.zeros((3, 4))
    affine[:, :3] = frame * block_mm
    affine[:, 3] = frame @ first_center
    return (size[0], size[1], size[2]), affine


# ------------------------------------------------------------- resampling
def _resample_labels(
    labels_kji: np.ndarray,
    source: np.ndarray,
    target: np.ndarray,
    size: tuple[int, int, int],
    threshold: float,
) -> np.ndarray:
    """Coverage-based label resampling; returns ``[i, j, k]`` labels.

    Each block is supersampled; every sample weights its eight neighbouring
    voxels trilinearly. A block is filled when labels cover at least
    ``threshold`` of it, and takes the label covering the most of it (ties go to
    the smaller label).
    """
    to_source = _compose(_invert(source), target)
    steps = np.linalg.norm(to_source[:, :3], axis=0)
    factors = [min(max(math.ceil(step), 1), MAX_SUPERSAMPLING) for step in steps]
    offsets = [
        ((a + 0.5) / factors[0] - 0.5, (b + 0.5) / factors[1] - 0.5, (c + 0.5) / factors[2] - 0.5)
        for c in range(factors[2])
        for b in range(factors[1])
        for a in range(factors[0])
    ]
    samples = float(len(offsets))

    present = np.unique(labels_kji)
    present = present[present != 0]
    output = np.zeros(int(np.prod(size)), dtype=np.int64)
    if not len(present):
        return output.reshape(size, order="F")
    columns = len(present) + 1
    ranks = (np.searchsorted(present, labels_kji) + 1).astype(np.int32)
    ranks[labels_kji == 0] = 0
    ranks_flat = ranks.reshape(-1)
    source_size = np.array(labels_kji.shape[::-1], dtype=np.int64)
    corner_offsets = [((corner >> 0) & 1, (corner >> 1) & 1, (corner >> 2) & 1) for corner in range(8)]

    total = len(output)
    batch = int(min(max((1 << 24) // columns, 4096), 1 << 20))
    for start in range(0, total, batch):
        flat = np.arange(start, min(start + batch, total), dtype=np.int64)
        index = np.stack(
            [flat % size[0], (flat // size[0]) % size[1], flat // (size[0] * size[1])], axis=1
        ).astype(float)
        count = len(flat)
        weights = np.zeros(count * columns)
        row_base = np.arange(count, dtype=np.int64) * columns
        for offset in offsets:
            position = (index + offset) @ to_source[:, :3].T + to_source[:, 3]
            base = np.floor(position)
            fraction = position - base
            base_int = base.astype(np.int64)
            for corner in corner_offsets:
                corner_index = base_int + corner
                weight = np.prod(
                    np.where(np.array(corner, dtype=bool), fraction, 1.0 - fraction), axis=1
                )
                valid = np.all((corner_index >= 0) & (corner_index < source_size), axis=1)
                valid &= weight > 0
                lookup = np.where(
                    valid,
                    corner_index[:, 0]
                    + source_size[0] * (corner_index[:, 1] + source_size[1] * corner_index[:, 2]),
                    0,
                )
                rank = np.where(valid, ranks_flat[lookup], 0)
                weights += np.bincount(
                    row_base + rank, weights=np.where(valid, weight, 0.0), minlength=count * columns
                )
        per_label = weights.reshape(count, columns)[:, 1:]
        foreground = per_label.sum(axis=1)
        filled = foreground / samples + 1e-9 >= threshold
        best = present[np.argmax(per_label, axis=1)]
        output[start : start + count] = np.where(filled, best, 0)
    return output.reshape(size, order="F")


def _permute_to_anatomy(
    labels_kji: np.ndarray, source: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Without resampling, anatomical orientation only permutes and flips axes."""
    directions = _voxel_frame(source)
    used = [False, False, False]
    mapping: list[tuple[int, float]] = []
    for target_axis in range(3):
        best_axis = -1
        best_dot = 0.0
        for source_axis in range(3):
            if used[source_axis]:
                continue
            dot = float(directions[:, source_axis] @ ANATOMICAL_FRAME[:, target_axis])
            if best_axis < 0 or abs(dot) > abs(best_dot):
                best_axis, best_dot = source_axis, dot
        used[best_axis] = True
        mapping.append((best_axis, -1.0 if best_dot < 0 else 1.0))

    source_size = labels_kji.shape[::-1]
    index_map = np.zeros((3, 4))
    voxels = labels_kji.transpose(2, 1, 0)
    for target_axis, (source_axis, sign) in enumerate(mapping):
        index_map[source_axis, target_axis] = sign
        if sign < 0:
            index_map[source_axis, 3] = source_size[source_axis] - 1
            voxels = np.flip(voxels, axis=source_axis)
    permuted = np.transpose(voxels, [source_axis for source_axis, _ in mapping])
    return permuted, _compose(source, index_map)


def transform(image: sitk.Image, options: TransformOptions) -> tuple[sitk.Image, TransformRecord]:
    """Apply optional cropping, anatomical orientation, and cubic resampling.

    Label values and the pixel type never change.
    """
    validate_options(options)
    if options.block_mm is None and options.orientation == "voxel":
        raise ValueError("nothing to transform: pass --block-mm and/or --orient anatomical")
    source = source_affine(image)
    labels = sitk.GetArrayViewFromImage(image)
    dtype = labels.dtype

    block_mm: float | None = None
    if options.block_mm is None:
        voxels, affine = _permute_to_anatomy(labels, source)
    else:
        frame = ANATOMICAL_FRAME if options.orientation == "anatomical" else _voxel_frame(source)
        low, high = _source_bounds(labels, options.crop)
        extent = _frame_extent(source, frame, low, high)
        block_mm = (
            _auto_block_size(image.GetSpacing(), extent, options.orientation)
            if options.block_mm == "auto"
            else float(options.block_mm)
        )
        size, affine = _target_grid(frame, extent, block_mm)
        voxels = _resample_labels(labels, source, affine, size, options.fill_threshold)

    derived = _image(voxels, affine, dtype)
    record = TransformRecord(
        orientation=options.orientation,
        block_mm=block_mm,
        fill_threshold=options.fill_threshold if block_mm is not None else None,
        cropped=options.crop,
        source_size=tuple(image.GetSize()),  # type: ignore[arg-type]
        source_spacing=tuple(image.GetSpacing()),  # type: ignore[arg-type]
        size=tuple(derived.GetSize()),  # type: ignore[arg-type]
        spacing=tuple(derived.GetSpacing()),  # type: ignore[arg-type]
        origin=tuple(derived.GetOrigin()),  # type: ignore[arg-type]
        direction=tuple(derived.GetDirection()),
    )
    return derived, record
