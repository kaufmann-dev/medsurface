"""Threshold, clean and prepare a binary label volume for surface extraction."""

from __future__ import annotations

import math

import numpy as np
import SimpleITK as sitk

from .geometry import (
    dilation_extent_mm,
    dilation_radius_voxels,
    kernel_extent_mm,
    kernel_radius_voxels,
    mm3_to_voxels,
)

FOREGROUND = 1


def otsu_threshold(values: np.ndarray, nbins: int = 512) -> float:
    """Otsu's threshold over a 1-D sample of intensities.

    Computed on a masked sample rather than the whole volume: in a CT the air
    background is so dominant that whole-volume Otsu simply separates air from
    everything else, which is never the boundary of interest.
    """
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("no finite voxels to threshold")
    lo, hi = float(values.min()), float(values.max())
    if hi <= lo:
        return lo

    hist, edges = np.histogram(values, bins=nbins, range=(lo, hi))
    hist = hist.astype(np.float64)
    centers = (edges[:-1] + edges[1:]) / 2.0

    total = hist.sum()
    if total == 0:
        return lo
    w0 = np.cumsum(hist)
    w1 = total - w0
    valid = (w0 > 0) & (w1 > 0)
    if not np.any(valid):
        return lo

    sum_total = float((hist * centers).sum())
    sum0 = np.cumsum(hist * centers)
    mean0 = np.divide(sum0, w0, out=np.zeros_like(sum0), where=w0 > 0)
    mean1 = np.divide(sum_total - sum0, w1, out=np.zeros_like(sum0), where=w1 > 0)

    between = w0 * w1 * (mean0 - mean1) ** 2
    between[~valid] = -1.0
    return float(centers[int(np.argmax(between))])


def auto_threshold(image: sitk.Image, modality: str) -> float:
    """Otsu over the voxels plausibly containing signal."""
    a = sitk.GetArrayViewFromImage(image).ravel()
    if modality == "CT":
        # Exclude air and the reconstruction-circle padding.
        sample = a[a > -400]
    else:
        sample = a[a > np.percentile(a, 50)] if a.size else a
    if sample.size == 0:
        sample = a
    # Subsample: Otsu on a histogram does not need 80M voxels.
    if sample.size > 4_000_000:
        step = sample.size // 4_000_000 + 1
        sample = sample[::step]
    return otsu_threshold(sample.astype(np.float64))


def binarize(image: sitk.Image, lower: float, upper: float | None) -> sitk.Image:
    """Threshold to a 0/1 label volume.

    When no upper bound is given we use the image's own maximum rather than a
    sentinel like float32-max: ITK casts the bound to the image's pixel type, and
    on an int16 CT that sentinel overflows to a negative number, making the upper
    bound smaller than the lower one.
    """
    if upper is None:
        mm_filter = sitk.MinimumMaximumImageFilter()
        mm_filter.Execute(image)
        hi = float(mm_filter.GetMaximum())
    else:
        hi = float(upper)

    lo = float(lower)
    if lo > hi:
        raise ValueError(
            "threshold %.1f is above the brightest voxel (%.1f); nothing would be "
            "segmented" % (lo, hi)
        )
    return sitk.BinaryThreshold(
        image, lowerThreshold=lo, upperThreshold=hi,
        insideValue=FOREGROUND, outsideValue=0,
    )


def _radii(mm: float, spacing) -> tuple[list[int], list[float], bool]:
    r = kernel_radius_voxels(mm, spacing)
    extent = kernel_extent_mm(r, spacing)
    return r, extent, any(x > 0 for x in r)


def median(image: sitk.Image, mm: float, log=None) -> sitk.Image:
    """Despeckle. Removes isolated noise voxels without rounding real edges."""
    if mm <= 0:
        return image
    r, extent, active = _radii(mm, image.GetSpacing())
    if not active:
        if log:
            log("median %.2f mm: kernel collapses to 1 voxel on every axis, skipping" % mm)
        return image
    if log:
        log("median  %-11s -> %s voxels = %s mm"
            % ("%.2f mm" % mm, r, _fmt(extent)))
    return sitk.BinaryMedian(image, radius=r, foregroundValue=FOREGROUND, backgroundValue=0)


def closing(image: sitk.Image, mm: float, log=None) -> sitk.Image:
    """Seal pores smaller than the kernel (trabecular bone, thin cortical gaps)."""
    if mm <= 0:
        return image
    r, extent, active = _radii(mm, image.GetSpacing())
    if not active:
        if log:
            log("closing %.2f mm: kernel collapses to 1 voxel on every axis, skipping" % mm)
        return image
    if log:
        log("closing %-11s -> %s voxels = %s mm" % ("%.2f mm" % mm, r, _fmt(extent)))
    return sitk.BinaryMorphologicalClosing(
        image, kernelRadius=r, kernelType=sitk.sitkBall,
        foregroundValue=FOREGROUND, safeBorder=True,
    )


def opening(image: sitk.Image, mm: float, log=None) -> sitk.Image:
    """Break thin bridges and shave protrusions thinner than the kernel."""
    if mm <= 0:
        return image
    r, extent, active = _radii(mm, image.GetSpacing())
    if not active:
        if log:
            log("opening %.2f mm: kernel collapses to 1 voxel on every axis, skipping" % mm)
        return image
    if log:
        log("opening %-11s -> %s voxels = %s mm" % ("%.2f mm" % mm, r, _fmt(extent)))
    return sitk.BinaryMorphologicalOpening(
        image, kernelRadius=r, kernelType=sitk.sitkBall,
        foregroundValue=FOREGROUND,
    )


def dilate(image: sitk.Image, mm: float, log=None) -> sitk.Image:
    """Grow the foreground outward by at least ``mm``.

    ``mm`` is a radius, not a kernel extent -- the distance the surface moves --
    which is why this uses :func:`geometry.dilation_radius_voxels` and not the
    flooring converter every other filter here uses. Anisotropic voxels make the
    realised growth exceed the request on the coarse axis.

    This is the primitive. For printing, use :func:`thicken`, which applies it
    only where the material is actually too thin.
    """
    if mm <= 0:
        return image
    spacing = image.GetSpacing()
    radii = dilation_radius_voxels(mm, spacing)
    if log:
        log("dilate %-11s -> %s voxels = grows %s mm"
            % ("%.2f mm" % mm, radii, _fmt(dilation_extent_mm(radii, spacing))))
    return sitk.BinaryDilate(
        image, kernelRadius=radii, kernelType=sitk.sitkBall,
        foregroundValue=FOREGROUND,
    )


def thin_mask(image: sitk.Image, min_feature_mm: float) -> sitk.Image:
    """The material a ball of ``min_feature_mm`` diameter cannot reach.

    That ball fits exactly where an opening by radius ``min_feature_mm / 2``
    survives, so the thin material is what the mask keeps and the opening throws
    away. No distance transform, no local-thickness estimator: an opening *is*
    the definition.
    """
    if min_feature_mm <= 0:
        empty = sitk.Image(image.GetSize(), sitk.sitkUInt8)
        empty.CopyInformation(image)
        return empty

    radii = dilation_radius_voxels(min_feature_mm / 2.0, image.GetSpacing())
    opened = sitk.BinaryMorphologicalOpening(
        image, kernelRadius=radii, kernelType=sitk.sitkBall,
        foregroundValue=FOREGROUND,
    )
    return sitk.And(sitk.Cast(image, sitk.sitkUInt8), sitk.Not(opened))


def count_foreground(image: sitk.Image) -> int:
    """Number of non-zero voxels.

    Exists so that no caller ever writes ``GetArrayViewFromImage(some_filter(...))``.
    That returns a numpy *view* into the image's buffer, and the temporary image is
    freed the moment the call returns, leaving the view pointing at released memory.
    It reads fine on a small test volume and segfaults on an 80M-voxel scan. Binding
    the image to a parameter keeps it alive for the duration of the count.
    """
    return int(np.count_nonzero(sitk.GetArrayViewFromImage(image)))


def thin_fraction(image: sitk.Image, min_feature_mm: float) -> float:
    """Fraction of the material that is thinner than ``min_feature_mm``.

    Reported, never enforced. A printer's minimum feature size is a property of
    the printer, not of the anatomy, and the right response to a thin orbital
    floor is a decision, not an automatic edit.
    """
    if min_feature_mm <= 0:
        return 0.0
    total = count_foreground(image)
    if total == 0:
        return 0.0
    return count_foreground(thin_mask(image, min_feature_mm)) / total


def _thick_core(image: sitk.Image, min_feature_mm: float) -> sitk.Image:
    """Everything a ball of ``min_feature_mm`` diameter can reach: the opening."""
    radii = dilation_radius_voxels(min_feature_mm / 2.0, image.GetSpacing())
    return sitk.BinaryMorphologicalOpening(
        image, kernelRadius=radii, kernelType=sitk.sitkBall,
        foregroundValue=FOREGROUND,
    )


def thicken(image: sitk.Image, thicken_mm: float, min_feature_mm: float = 0.0,
            log=None) -> sitk.Image:
    """Grow only the material that is too thin to print, and leave the rest alone.

    Dilating everything is simpler, and it is what a naive print-prep step does,
    but it is dimensionally wrong. A cranial vault is 5 mm of solid bone and needs
    nothing; inflating it moves the model's outer surface for no benefit. Only the
    paper-thin structures -- orbital floor, ethmoid, nasal septum -- need material.

    Selecting the thin set is subtler than it looks. ``mask \\ opening(mask, r)``
    is the textbook answer and it is wrong here: an opening is the union of the
    balls it contains, so it cannot reach into a sharp convex corner, and *every
    surface of a voxelised object is locally sharp*. That set is a speckle over
    the whole surface, and dilating it inflates the entire model -- a voxelised
    sphere grows its bounding box by the full thickening radius.

    So thin material is material that no sufficiently thick region can reach::

        core = opening(mask, feature / 2)          # everything thick enough
        thin = mask \\ dilate(core, thicken_mm)     # beyond the core's reach
        out  = mask | dilate(thin, thicken_mm)

    A solid block is entirely within its own core's reach and does not move. A
    one-voxel sheet has no core at all and is grown everywhere.

    With no minimum feature size, material thinner than twice the thickening
    radius is treated as thin -- what asking to "thicken by 1 mm" implies.
    """
    if thicken_mm <= 0:
        return image

    feature = min_feature_mm if min_feature_mm > 0 else 2.0 * thicken_mm
    binary = sitk.Cast(image, sitk.sitkUInt8)

    reachable = dilate(_thick_core(binary, feature), thicken_mm)
    thin = sitk.And(binary, sitk.Not(reachable))

    thin_voxels = count_foreground(thin)
    if thin_voxels == 0:
        if log:
            log("thicken %.2f mm: nothing thinner than %.2f mm, nothing to do"
                % (thicken_mm, feature))
        return image

    total = count_foreground(binary)
    spacing = image.GetSpacing()
    radii = dilation_radius_voxels(thicken_mm, spacing)
    if log:
        log("thicken %-11s -> %s voxels = grows %s mm, applied to the %.1f%% of "
            "material out of reach of bone thicker than %.2f mm"
            % ("%.2f mm" % thicken_mm, radii, _fmt(dilation_extent_mm(radii, spacing)),
               100.0 * thin_voxels / max(total, 1), feature))

    return sitk.Or(binary, dilate(thin, thicken_mm))


def islands(
    image: sitk.Image,
    keep_largest: bool,
    min_mm3: float,
    log=None,
) -> sitk.Image:
    """Drop small connected components; optionally keep only the largest.

    Note this operates on the *labelmap*. It removes free-floating specks, but it
    does not remove enclosed internal cavities -- those are background, not
    foreground, and never appear as separate foreground components. Removing them
    is a mesh-level operation (see ``surface.largest_component``).
    """
    if not keep_largest and min_mm3 <= 0:
        return image

    min_voxels = mm3_to_voxels(min_mm3, image.GetSpacing()) if min_mm3 > 0 else 0
    cc = sitk.ConnectedComponent(image)
    relabelled = sitk.RelabelComponent(cc, minimumObjectSize=min_voxels, sortByObjectSize=True)

    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(relabelled)
    n = len(stats.GetLabels())
    if n == 0:
        raise ValueError(
            "thresholding left no foreground; the threshold is probably above the "
            "brightest voxel in the volume"
        )

    if keep_largest:
        if log:
            log("islands: %d component(s) >= %.0f mm3, keeping largest" % (n, min_mm3))
        return sitk.BinaryThreshold(relabelled, 1, 1, FOREGROUND, 0)

    if log:
        log("islands: keeping %d component(s) >= %.0f mm3" % (n, min_mm3))
    return sitk.BinaryThreshold(relabelled, 1, n, FOREGROUND, 0)


#: Isolevel of the fractional-occupancy field. A voxel of 1.0 is fully inside,
#: 0.0 fully outside, so the surface lies at 0.5 -- the same convention as
#: marching cubes on the raw 0/1 labelmap, which keeps volumes consistent
#: between the resampled and native paths.
ISO_OCCUPANCY = 0.5

#: Gaussian sigma, as a fraction of the target voxel size, applied before
#: downsampling. This is the anti-aliasing pre-filter: point-sampling a binary
#: mask onto a coarser grid makes thin structures blink in and out depending on
#: where the samples happen to land. A symmetric kernel leaves the 0.5 level of
#: a straight edge exactly where it was, so it costs almost no volume
#: (measured: -0.26% on a 20 mm sphere at 1 mm sampling).
_ANTIALIAS_SIGMA_FACTOR = 0.3


def antialias_for_grid(image: sitk.Image, grid_mm: float) -> sitk.Image:
    """Blur only the axes that are about to be downsampled.

    Point-sampling a mask onto a coarser grid makes thin structures blink in and
    out depending on where the samples land. Axes being *up*sampled need no
    pre-filter, and blurring them would only cost detail. SimpleITK rejects a
    vector sigma containing zeros, so the axes are filtered one at a time.
    """
    out = sitk.Cast(image, sitk.sitkFloat32)
    for axis, spacing in enumerate(image.GetSpacing()):
        if grid_mm <= spacing:
            continue
        blur = sitk.RecursiveGaussianImageFilter()
        blur.SetDirection(axis)
        blur.SetSigma(_ANTIALIAS_SIGMA_FACTOR * grid_mm)
        out = blur.Execute(out)
    return out


def resample_isotropic(binary: sitk.Image, mm: float, log=None) -> sitk.Image:
    """Resample the mask onto an isotropic grid as a fractional-occupancy field.

    This is the topology-safe way to control triangle count. Decimating the mesh
    afterwards tears thin structures -- a skull's orbital walls are one voxel
    thick, and edge collapses weld their opposite faces together, silently
    opening a closed surface. Resampling cannot do that, because marching cubes
    always returns a manifold surface.

    A signed distance field would be the textbook choice here, but ITK's Maurer
    transform quantises distance to voxel centres, placing its zero level half a
    voxel inside the true boundary. Measured on a 20 mm sphere that shrinks the
    volume by 6%. Smoothed occupancy keeps the error under 0.5%.
    """
    field = antialias_for_grid(binary, mm)

    size = [
        max(1, int(math.ceil(n * s / mm)))
        for n, s in zip(field.GetSize(), field.GetSpacing())
    ]
    r = sitk.ResampleImageFilter()
    r.SetOutputSpacing((mm, mm, mm))
    r.SetSize(size)
    r.SetOutputOrigin(field.GetOrigin())
    r.SetOutputDirection(field.GetDirection())
    r.SetInterpolator(sitk.sitkLinear)
    r.SetDefaultPixelValue(0.0)  # outside
    out = r.Execute(field)

    # The coarser grid may clip the object at its outer face; a border of
    # background keeps the isosurface closed.
    out = sitk.ConstantPad(out, [1, 1, 1], [1, 1, 1], 0.0)

    if log:
        log("resample %.2f mm isotropic -> %s voxels"
            % (mm, "x".join(str(v) for v in out.GetSize())))
    return out


def pad(image: sitk.Image, width: int = 1) -> sitk.Image:
    """Surround the volume with background.

    Marching cubes does not close the surface where foreground meets the edge of
    the image. Anatomy truncated by the scanner's field of view therefore yields
    an *open* mesh. One voxel of background caps it flat. SimpleITK shifts the
    image origin accordingly, so physical coordinates stay correct.
    """
    if width <= 0:
        return image
    w = [int(width)] * image.GetDimension()
    return sitk.ConstantPad(image, w, w, 0.0)


def _fmt(values) -> str:
    return " x ".join("%.2f" % v for v in values)
