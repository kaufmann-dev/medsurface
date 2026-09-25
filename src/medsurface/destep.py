"""Masked Taubin fairing that removes broad stair-step ripples.

CT slice terraces on smooth anatomy such as a skull dome can survive as shallow
ripples several millimetres long. Voxel-scale mask smoothing and one-ring mesh
relaxation are too local to remove them, and stronger global smoothing erases
thin anatomy. This stage instead runs many uniform-Laplacian Taubin iterations
under a smooth per-vertex weight: weight 0 vertices never move, so detail kept
outside the faired region is bit-identical.

Every vertex displacement is clamped, faces folded by the fairing return to
their input positions, and the existing smoothing safeguard then removes any
remaining self-intersections or inconsistent winding. Connectivity never
changes, so a watertight input stays watertight.
"""

from __future__ import annotations

from dataclasses import dataclass

import meshlib.mrmeshpy as mrmeshpy
import numpy as np
from scipy import sparse

from . import surface
from .defaults import (
    DESTEP_AUTO_DETAIL_RADIUS_MM,
    DESTEP_AUTO_FULL_RADIUS_MM,
    DESTEP_AUTO_GUARD_RINGS,
    DESTEP_AUTO_SPREAD_PASSES,
    DESTEP_LAMBDA,
    DESTEP_MU,
    DESTEP_UNFOLD_DOT,
)
from .presets import DestepSettings, validate_destep

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
_MAX_UNFOLD_ROUNDS = 5


@dataclass(frozen=True)
class DestepStats:
    region: str
    iterations: int
    max_displacement_mm: float
    faired_vertex_fraction: float
    frozen_vertex_fraction: float
    clamped_vertices: int
    clamped_vertex_fraction: float
    unfolded_vertices: int
    p95_displacement_mm: float
    volume_change_percent: float
    safeguard: surface.SmoothingSafeguard


def adjacency(vertex_count: int, faces: np.ndarray) -> sparse.csr_matrix:
    """Symmetric vertex adjacency of the unique triangle edges."""
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    return sparse.csr_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(vertex_count, vertex_count)
    )


def neighbour_average(adjacent: sparse.csr_matrix) -> sparse.csr_matrix:
    """Operator that replaces each vertex value with its neighbours' mean."""
    degree = np.asarray(adjacent.sum(axis=1)).ravel()
    degree[degree == 0] = 1
    return sparse.csr_matrix(sparse.diags(1.0 / degree) @ adjacent)


def masked_taubin(
    vertices: np.ndarray,
    average: sparse.csr_matrix,
    weights: np.ndarray,
    iterations: int,
) -> np.ndarray:
    """Alternate lambda/mu uniform-Laplacian steps scaled per vertex."""
    result = vertices.copy()
    lambda_step = (weights * DESTEP_LAMBDA)[:, None]
    mu_step = (weights * DESTEP_MU)[:, None]
    for iteration in range(iterations):
        step = lambda_step if iteration % 2 == 0 else mu_step
        result += step * (average @ result - result)
    return result


def _cosine_ramp(t: np.ndarray) -> np.ndarray:
    t = np.clip(t, 0.0, 1.0)
    return 0.5 - 0.5 * np.cos(np.pi * t)


def band_weights(vertices: np.ndarray, settings: DestepSettings) -> np.ndarray:
    """Ramp from 0 at the frozen coordinate to 1 at the fully faired one."""
    assert settings.full_mm is not None and settings.frozen_mm is not None
    coordinate = vertices[:, _AXIS_INDEX[settings.axis]]
    weights = _cosine_ramp(
        (coordinate - settings.frozen_mm) / (settings.full_mm - settings.frozen_mm)
    )
    if not len(weights) or weights.max() <= 0:
        low, high = (
            (float(coordinate.min()), float(coordinate.max()))
            if len(coordinate)
            else (0.0, 0.0)
        )
        raise ValueError(
            "destep band from frozen %s=%.3f mm to full %s=%.3f mm does not reach "
            "the surface, whose %s range is %.3f to %.3f mm"
            % (
                settings.axis,
                settings.frozen_mm,
                settings.axis,
                settings.full_mm,
                settings.axis,
                low,
                high,
            )
        )
    return weights


def _dilate(selected: np.ndarray, adjacent: sparse.csr_matrix, rings: int) -> np.ndarray:
    grown = selected.copy()
    for _ in range(rings):
        grown |= (adjacent @ grown.astype(np.float64)) > 0
    return grown


def auto_weights(
    vertices: np.ndarray,
    faces: np.ndarray,
    adjacent: sparse.csr_matrix,
    average: sparse.csr_matrix,
    iterations: int,
) -> np.ndarray:
    """Fair broad surfaces and freeze tightly curved detail.

    Curvature is estimated on an unmasked faired copy, where terraces no longer
    create artificial edges. Each vertex takes the largest normal change per
    millimetre along its edges, which also catches rims and saddles that mean
    curvature misses, and the estimate is averaged over its neighbourhood.
    """
    faired = masked_taubin(vertices, average, np.ones(len(vertices)), iterations)
    normals = surface.vertex_normals_from_arrays(faired, faces)
    edges = adjacent.tocoo()
    lengths = np.linalg.norm(faired[edges.row] - faired[edges.col], axis=1)
    edge_curvature = np.linalg.norm(
        normals[edges.row] - normals[edges.col], axis=1
    ) / np.maximum(lengths, 1e-12)
    curvature = np.zeros(len(vertices))
    np.maximum.at(curvature, edges.row, edge_curvature)
    for _ in range(DESTEP_AUTO_SPREAD_PASSES):
        curvature = 0.5 * (curvature + average @ curvature)
    radius = 1.0 / np.maximum(curvature, 1e-12)

    weights = _cosine_ramp(
        (radius - DESTEP_AUTO_DETAIL_RADIUS_MM)
        / (DESTEP_AUTO_FULL_RADIUS_MM - DESTEP_AUTO_DETAIL_RADIUS_MM)
    )
    detail = radius <= DESTEP_AUTO_DETAIL_RADIUS_MM
    weights[_dilate(detail, adjacent, DESTEP_AUTO_GUARD_RINGS)] = 0.0
    for _ in range(DESTEP_AUTO_GUARD_RINGS):
        weights = 0.5 * (weights + average @ weights)
    weights[detail] = 0.0
    return weights


def _unit_face_normals(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    triangles = vertices[faces]
    normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    lengths = np.linalg.norm(normals, axis=1)
    normals[lengths > 0] /= lengths[lengths > 0, None]
    return normals


def revert_folds(
    vertices: np.ndarray,
    original: np.ndarray,
    faces: np.ndarray,
    adjacent: sparse.csr_matrix,
) -> tuple[np.ndarray, int]:
    """Restore vertices of faces whose normal flipped, plus one ring."""
    result = vertices.copy()
    original_normals = _unit_face_normals(original, faces)
    comparable = np.linalg.norm(original_normals, axis=1) > 0
    reverted = np.zeros(len(result), dtype=bool)
    for _ in range(_MAX_UNFOLD_ROUNDS):
        normals = _unit_face_normals(result, faces)
        alignment = np.einsum("ij,ij->i", normals, original_normals)
        folded = comparable & (alignment < DESTEP_UNFOLD_DOT)
        if not folded.any():
            break
        selected = np.zeros(len(result), dtype=bool)
        selected[np.unique(faces[folded])] = True
        selected = _dilate(selected, adjacent, 1)
        result[selected] = original[selected]
        reverted |= selected
    return result, int(np.count_nonzero(reverted))


def destep_safely(
    mesh: mrmeshpy.Mesh,
    settings: DestepSettings,
) -> tuple[mrmeshpy.Mesh, DestepStats]:
    """Fair the selected region and keep only a clean, clamped result."""
    validate_destep(settings)
    vertices, faces = surface.to_arrays(mesh)
    adjacent = adjacency(len(vertices), faces)
    average = neighbour_average(adjacent)

    if settings.region == "band":
        weights = band_weights(vertices, settings)
    elif settings.region == "auto":
        weights = auto_weights(
            vertices,
            faces,
            adjacent,
            average,
            settings.iterations,
        )
    else:
        weights = np.ones(len(vertices))

    faired = masked_taubin(vertices, average, weights, settings.iterations)
    delta = faired - vertices
    distance = np.linalg.norm(delta, axis=1)
    clamped = distance > settings.max_displacement_mm
    delta[clamped] *= (settings.max_displacement_mm / distance[clamped])[:, None]
    candidate_vertices, unfolded = revert_folds(
        vertices + delta, vertices, faces, adjacent
    )

    result, safeguard = surface.protect_smoothed_surface(
        mesh,
        surface.from_arrays(candidate_vertices, faces),
        settings.iterations,
    )
    result_vertices = surface.to_arrays(result)[0]
    moved = np.linalg.norm(result_vertices - vertices, axis=1)
    volume_before = abs(float(mesh.volume()))
    volume_after = abs(float(result.volume()))
    return result, DestepStats(
        region=settings.region,
        iterations=settings.iterations,
        max_displacement_mm=float(settings.max_displacement_mm),
        faired_vertex_fraction=float(np.mean(weights > 0.99)) if len(weights) else 0.0,
        frozen_vertex_fraction=float(np.mean(weights < 0.01)) if len(weights) else 0.0,
        clamped_vertices=int(np.count_nonzero(clamped)),
        clamped_vertex_fraction=float(np.mean(clamped)) if len(clamped) else 0.0,
        unfolded_vertices=unfolded,
        p95_displacement_mm=float(np.percentile(moved, 95)) if len(moved) else 0.0,
        volume_change_percent=(
            100.0 * (volume_after - volume_before) / volume_before
            if volume_before > 0
            else 0.0
        ),
        safeguard=safeguard,
    )
