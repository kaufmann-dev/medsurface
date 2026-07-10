"""Rigid registration of two scans of the same anatomy.

This implementation does not inspect ``FrameOfReferenceUID`` or assume that two
series already align. DICOM patient coordinates are patient-oriented (LPS), but
origin, pose, and head tilt can differ between acquisitions, so ``merge`` always
registers the moving scan.

The search is global-then-local:

1. **Exhaustive translation** by FFT cross-correlation of the two bone masks on a
   coarse lattice. Every integer translation is scored at once, so the translation
   needs no initial guess. Rotation is *not* searched here.
2. **Point-to-plane ICP** on the surfaces to recover rotation and refine. Rotation
   therefore still depends on ICP's basin of attraction: measured on a synthetic
   phantom, it recovers 3 to 40 degrees to under 1 degree and fails beyond about
   60. Point-to-point ICP slides far too slowly on smooth bone; on a real skull
   pair it was still descending after 60 iterations at 0.93 mm RMS, where
   point-to-plane converged to 0.176 mm in 40.

Correspondences are *trimmed*. The scans overlap only partially -- one holds the
cranial vault, the other the mandible -- so a plain least-squares fit would drag
the non-overlapping points into correspondence and skew the result.

Nothing here decides whether a registration is *good*; it reports diagnostics and
lets the caller set the bar. See :class:`RegistrationResult`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import SimpleITK as sitk
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree
from vtk.util import numpy_support  # noqa: N813

from . import segment, surface

#: Voxel size of the lattice used for the exhaustive translation search. Coarse
#: on purpose: it only has to land within ICP's capture radius.
LATTICE_MM = 2.0

#: A source point further than this from the target surface is assumed to have no
#: counterpart (the scans cover different anatomy) rather than to be misaligned.
CORRESPONDENCE_TOL_MM = 4.0

#: Fraction of source points kept in the first, robust fitting pass. The rest are
#: assumed to have no counterpart.
#:
#: A second pass then re-fits with the trim widened to the overlap actually
#: measured. Trimming is necessary when the scans overlap partially, but harmful
#: when they overlap fully: the discarded correspondences are the ones furthest
#: from the rotation axis, which is where the rotational signal lives. Held at
#: 0.45, a 5 degree misalignment of two identical volumes converges 4.49 degrees
#: off; widening the trim brings the same case to 0.44.
TRIM_KEEP = 0.45

#: Never trim beyond this, so a stray outlier cannot dominate the final fit.
TRIM_KEEP_MAX = 0.95


class RegistrationError(RuntimeError):
    pass


@dataclass
class RegistrationResult:
    """Transform mapping scan B into scan A's frame, with quality diagnostics."""

    transform: np.ndarray  # 4x4, homogeneous
    fft_translation_mm: np.ndarray
    rotation_deg: float
    translation_mm: np.ndarray
    #: Point-to-plane residual over the trimmed inliers.
    inlier_rms_mm: float
    inlier_median_mm: float
    #: Fraction of the moving surface, *restricted to where the fixed scan has
    #: data*, that found a counterpart -- and the same measured the other way.
    #:
    #: Both directions are needed. A one-directional measure is scale-free: a
    #: small dense object buried in a large one finds a nearby surface almost
    #: everywhere. A bar phantom scored 96.4% against a skull. Measured the other
    #: way -- how much of the skull's surface the bar accounts for -- it collapses.
    overlap_moving_in_fixed: float
    overlap_fixed_in_moving: float
    #: Dice of the two bone masks, computed only where both scans actually have
    #: data. Insensitive to how much anatomy each scan uniquely covers.
    shared_fov_dice: float
    shared_fov_mm3: float

    @property
    def surface_overlap(self) -> float:
        """The weaker of the two directions. This is what should be gated on."""
        return min(self.overlap_moving_in_fixed, self.overlap_fixed_in_moving)

    def summary(self) -> str:
        return (
            "rotation %.2f deg, translation [%.1f %.1f %.1f] mm\n"
            "  point-to-plane rms %.3f mm (median %.3f mm)\n"
            "  surface overlap    %.1f%% moving->fixed, %.1f%% fixed->moving\n"
            "  dice in shared fov %.3f  (%.0f cm3 seen by both)"
            % (
                self.rotation_deg, *self.translation_mm,
                self.inlier_rms_mm, self.inlier_median_mm,
                100.0 * self.overlap_moving_in_fixed,
                100.0 * self.overlap_fixed_in_moving,
                self.shared_fov_dice, self.shared_fov_mm3 / 1000.0,
            )
        )


# --------------------------------------------------------------------- helpers
def _world_lattice(mask: sitk.Image, mm: float) -> tuple[np.ndarray, np.ndarray]:
    """Resample onto an axis-aligned world grid so lattices can be correlated."""
    size = np.asarray(mask.GetSize())
    corners = []
    for i in (0, size[0] - 1):
        for j in (0, size[1] - 1):
            for k in (0, size[2] - 1):
                corners.append(mask.TransformIndexToPhysicalPoint((int(i), int(j), int(k))))
    corners = np.asarray(corners, dtype=float)

    lo = np.floor(corners.min(0) / mm) * mm - mm
    hi = np.ceil(corners.max(0) / mm) * mm + mm
    out_size = np.maximum(1, np.round((hi - lo) / mm).astype(int) + 1)

    r = sitk.ResampleImageFilter()
    r.SetOutputSpacing((mm, mm, mm))
    r.SetSize([int(v) for v in out_size])
    r.SetOutputOrigin([float(v) for v in lo])
    r.SetOutputDirection(np.eye(3).ravel().tolist())
    r.SetInterpolator(sitk.sitkNearestNeighbor)
    r.SetDefaultPixelValue(0)
    arr = sitk.GetArrayFromImage(r.Execute(mask)).astype(np.float32)  # (z, y, x)
    return arr, lo


def _fft_translation(a: np.ndarray, oa: np.ndarray,
                     b: np.ndarray, ob: np.ndarray, mm: float) -> np.ndarray:
    """World translation t such that B + t best overlaps A, over all integer lags."""
    if not a.any() or not b.any():
        raise RegistrationError("one of the masks is empty; nothing to register")
    corr = fftconvolve(a, b[::-1, ::-1, ::-1], mode="full")
    lag = np.array(np.unravel_index(int(np.argmax(corr)), corr.shape)) - (np.array(b.shape) - 1)
    return (np.asarray(oa) - np.asarray(ob)) + mm * lag[::-1]  # (z,y,x) -> (x,y,z)


def _matched_fraction(src, tree, transform, max_dist) -> float:
    moved = (transform[:3, :3] @ src.T).T + transform[:3, 3]
    dist, _ = tree.query(moved, workers=-1)
    return float((dist < max_dist).mean())


def _icp_point_to_plane(src, tgt, tgt_normals, tgt_tree, transform,
                        iterations=80, keep=TRIM_KEEP, max_dist=CORRESPONDENCE_TOL_MM):
    rot, trans = transform[:3, :3].copy(), transform[:3, 3].copy()
    for _ in range(iterations):
        p = (rot @ src.T).T + trans
        d, idx = tgt_tree.query(p, workers=-1)
        ok = (d < max_dist) & (d <= max(np.quantile(d, keep), 1e-9))
        if ok.sum() < 200:
            break

        pi, qi, ni = p[ok], tgt[idx[ok]], tgt_normals[idx[ok]]
        # linearised rigid step: minimise sum(((p + w x p + dt) - q) . n)^2
        a = np.hstack([np.cross(pi, ni), ni])
        b = -np.einsum("ij,ij->i", pi - qi, ni)
        sol, *_ = np.linalg.lstsq(a, b, rcond=None)
        w, dt = sol[:3], sol[3:]

        theta = float(np.linalg.norm(w))
        if theta < 1e-12:
            drot = np.eye(3)
        else:
            k = w / theta
            kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
            drot = np.eye(3) + math.sin(theta) * kx + (1 - math.cos(theta)) * (kx @ kx)

        rot = drot @ rot
        trans = drot @ trans + dt

    out = np.eye(4)
    out[:3, :3], out[:3, 3] = rot, trans
    return out


def surface_points(mask: sitk.Image, smooth_iters: int = 10):
    """Vertices and vertex normals of a mask's surface, in patient coordinates."""
    padded = segment.pad(mask, 1)
    poly = surface.marching_cubes(surface.to_vtk_image(padded), 0.5)
    if poly.GetNumberOfPolys() == 0:
        raise RegistrationError("mask has no surface; the threshold left no foreground")
    poly = surface.smooth(poly, smooth_iters, 0.1)
    poly = surface.transform(poly, surface.index_to_physical(padded))
    poly = surface.compute_normals(poly)

    pts = numpy_support.vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float64)
    normals = numpy_support.vtk_to_numpy(poly.GetPointData().GetNormals()).astype(np.float64)
    return pts, normals


def _inside_fov(image: sitk.Image, points: np.ndarray) -> np.ndarray:
    """Which physical points fall inside the image's acquired volume."""
    direction = np.asarray(image.GetDirection()).reshape(3, 3)
    origin = np.asarray(image.GetOrigin())
    spacing = np.asarray(image.GetSpacing())
    index = (np.linalg.inv(direction) @ (points - origin).T).T / spacing
    upper = np.asarray(image.GetSize()) - 1
    return np.all((index >= 0) & (index <= upper), axis=1)


def _directional_overlap(points, tree, fov_image, tol) -> float:
    """Fraction of `points` inside `fov_image` that have a neighbour within `tol`.

    Points outside the other scan's field of view are excluded rather than
    counted as misses: a scan legitimately images anatomy the other never saw.
    """
    inside = _inside_fov(fov_image, points)
    if inside.sum() < 100:
        return 0.0
    dist, _ = tree.query(points[inside], workers=-1)
    return float((dist < tol).mean())


def _coverage(mask: sitk.Image) -> sitk.Image:
    ones = sitk.Image(mask.GetSize(), sitk.sitkUInt8)
    ones.CopyInformation(mask)
    return ones + 1


def _resample_binary(image, size, origin, spacing, transform) -> np.ndarray:
    """Resample a 0/1 image and re-threshold at half occupancy.

    Linear, not nearest-neighbour. When the output lattice is commensurate with
    the input (1.5 mm samples of a 1.0 mm grid), nearest-neighbour lands on exact
    midpoints, and a 1e-15 perturbation of the transform flips the tie. Dice then
    drops from 1.000 to 0.916 for a mesh registered against *itself*.
    """
    r = sitk.ResampleImageFilter()
    r.SetOutputSpacing([float(spacing)] * 3)
    r.SetSize([int(v) for v in size])
    r.SetOutputOrigin([float(v) for v in origin])
    r.SetOutputDirection(np.eye(3).ravel().tolist())
    r.SetInterpolator(sitk.sitkLinear)
    r.SetDefaultPixelValue(0.0)
    r.SetTransform(transform)
    out = r.Execute(sitk.Cast(image, sitk.sitkFloat32))
    return sitk.GetArrayViewFromImage(out) > 0.5


def inverse_transform(t: np.ndarray) -> sitk.AffineTransform:
    """SimpleITK resampling maps OUTPUT points to INPUT points, so it wants T^-1."""
    rot_t = t[:3, :3].T
    aff = sitk.AffineTransform(3)
    aff.SetMatrix([float(v) for v in rot_t.ravel()])
    aff.SetTranslation([float(v) for v in (-rot_t @ t[:3, 3])])
    return aff


def _shared_fov_dice(mask_a, mask_b, t, mm=1.5):
    """Dice restricted to the region both scans actually imaged.

    Plain Dice over the union would punish a perfectly good registration just
    because one scan covers anatomy the other never saw.
    """
    corners = []
    for mask, xform in ((mask_a, None), (mask_b, t)):
        size = np.asarray(mask.GetSize()) - 1
        for i in (0, size[0]):
            for j in (0, size[1]):
                for k in (0, size[2]):
                    p = np.asarray(mask.TransformIndexToPhysicalPoint((int(i), int(j), int(k))))
                    corners.append(p if xform is None else xform[:3, :3] @ p + xform[:3, 3])
    corners = np.asarray(corners)
    lo = np.floor(corners.min(0) / mm) * mm - mm
    hi = np.ceil(corners.max(0) / mm) * mm + mm
    size = np.maximum(1, np.round((hi - lo) / mm).astype(int) + 1)

    ident = sitk.Transform(3, sitk.sitkIdentity)
    inv = inverse_transform(t)

    aa = _resample_binary(mask_a, size, lo, mm, ident)
    bb = _resample_binary(mask_b, size, lo, mm, inv)
    ca = _resample_binary(_coverage(mask_a), size, lo, mm, ident)
    cb = _resample_binary(_coverage(mask_b), size, lo, mm, inv)
    both = ca & cb

    shared_mm3 = float(both.sum()) * mm ** 3
    aa, bb = aa & both, bb & both
    denom = float(aa.sum() + bb.sum())
    dice = 0.0 if denom == 0 else 2.0 * float((aa & bb).sum()) / denom
    return dice, shared_mm3


# ----------------------------------------------------------------------- api
def rigid_register(mask_a: sitk.Image, mask_b: sitk.Image,
                   samples: int = 60000, seed: int = 0,
                   log=None) -> RegistrationResult:
    """Rigidly map ``mask_b`` into ``mask_a``'s frame."""

    def say(msg):
        if log:
            log(msg)

    say("exhaustive translation search (FFT) ...")
    a_lat, a_origin = _world_lattice(mask_a, LATTICE_MM)
    b_lat, b_origin = _world_lattice(mask_b, LATTICE_MM)
    t0_translation = _fft_translation(a_lat, a_origin, b_lat, b_origin, LATTICE_MM)
    say("  translation %s mm" % np.round(t0_translation, 1))

    transform = np.eye(4)
    transform[:3, 3] = t0_translation

    say("ICP ...")
    tgt, tgt_normals = surface_points(mask_a)
    src, _ = surface_points(mask_b)

    rng = np.random.default_rng(seed)
    if len(src) > samples:
        src = src[rng.choice(len(src), samples, replace=False)]

    tree = cKDTree(tgt)
    # Pass 1: trim hard, so non-overlapping anatomy cannot drag the fit.
    transform = _icp_point_to_plane(src, tgt, tgt_normals, tree, transform, keep=TRIM_KEEP)
    # Pass 2: widen the trim to the overlap we actually found, so a fully
    # overlapping pair keeps the correspondences that carry the rotation.
    matched = _matched_fraction(src, tree, transform, CORRESPONDENCE_TOL_MM)
    keep = min(TRIM_KEEP_MAX, max(TRIM_KEEP, 0.9 * matched))
    say("  matched %.0f%% of moving surface -> refitting with trim %.2f" % (100 * matched, keep))
    transform = _icp_point_to_plane(src, tgt, tgt_normals, tree, transform, keep=keep)

    rot, trans = transform[:3, :3], transform[:3, 3]
    moved = (rot @ src.T).T + trans
    dist, idx = tree.query(moved, workers=-1)

    # Symmetric overlap, each direction restricted to the other scan's field of
    # view. See RegistrationResult.overlap_moving_in_fixed for why.
    forward = _directional_overlap(moved, tree, mask_a, CORRESPONDENCE_TOL_MM)

    src_tree = cKDTree(src)
    tgt_sample = tgt if len(tgt) <= samples else tgt[rng.choice(len(tgt), samples, replace=False)]
    inv = np.linalg.inv(transform)
    tgt_in_b = (inv[:3, :3] @ tgt_sample.T).T + inv[:3, 3]
    backward = _directional_overlap(tgt_in_b, src_tree, mask_b, CORRESPONDENCE_TOL_MM)

    inliers = dist <= max(np.quantile(dist, TRIM_KEEP), 1e-9)
    resid = np.abs(np.einsum("ij,ij->i", moved[inliers] - tgt[idx[inliers]],
                             tgt_normals[idx[inliers]]))
    rms = float(np.sqrt((resid ** 2).mean())) if resid.size else float("inf")
    median = float(np.median(resid)) if resid.size else float("inf")

    say("computing shared-field-of-view agreement ...")
    dice, shared_mm3 = _shared_fov_dice(mask_a, mask_b, transform)

    angle = math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(rot) - 1) / 2))))
    return RegistrationResult(
        transform=transform,
        fft_translation_mm=t0_translation,
        rotation_deg=angle,
        translation_mm=trans,
        inlier_rms_mm=rms,
        inlier_median_mm=median,
        overlap_moving_in_fixed=forward,
        overlap_fixed_in_moving=backward,
        shared_fov_dice=dice,
        shared_fov_mm3=shared_mm3,
    )
