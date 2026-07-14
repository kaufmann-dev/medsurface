"""Threshold, clean and prepare a binary label volume for surface extraction."""

from __future__ import annotations

import math

import numpy as np
import SimpleITK as sitk

from .defaults import MAX_VOXELS
from .geometry import (
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


def auto_threshold(image: sitk.Image) -> float:
    """Format-neutral Otsu over voxels plausibly containing signal."""
    a = sitk.GetArrayViewFromImage(image).ravel()
    finite = a[np.isfinite(a)]
    sample = finite[finite > np.percentile(finite, 50)] if finite.size else finite
    if sample.size == 0:
        sample = finite
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


def smooth_occupancy(image: sitk.Image, sigma_mm: float) -> sitk.Image:
    """Regularize a binary or fractional occupancy boundary in physical space.

    The 0.5 isosurface of a straight boundary remains centered while voxel-scale
    terraces are attenuated. Curved boundaries, narrow gaps, and structures near
    the configured sigma can move, merge, or disappear; callers surface that
    tradeoff to users instead of treating this as lossless anti-aliasing.
    """
    if not math.isfinite(sigma_mm) or sigma_mm < 0:
        raise ValueError("occupancy smoothing sigma must be finite and non-negative")
    if sigma_mm == 0:
        return image
    real = sitk.Cast(image, sitk.sitkFloat32)
    if min(image.GetSize()) < 4:
        return sitk.DiscreteGaussian(
            real,
            variance=[float(sigma_mm) ** 2] * image.GetDimension(),
            useImageSpacing=True,
        )
    return sitk.SmoothingRecursiveGaussian(real, float(sigma_mm))


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


def resample_isotropic(
    binary: sitk.Image,
    mm: float,
    log=None,
    pad_border: bool = True,
    *,
    allow_large_volume: bool = False,
) -> sitk.Image:
    """Resample the mask onto an isotropic grid as a fractional-occupancy field.

    Resampling can reduce triangle count but can erase thin structures and change
    components, cavities, tunnels, or genus. Padding usually enables a closed
    isosurface, but output validity is measured rather than assumed.

    A signed distance field would be the textbook choice here, but ITK's Maurer
    transform quantises distance to voxel centres, placing its zero level half a
    voxel inside the true boundary. Measured on a 20 mm sphere that shrinks the
    volume by 6%. Smoothed occupancy keeps the error under 0.5%.
    """
    if not math.isfinite(mm) or mm <= 0:
        raise ValueError("resample spacing must be finite and greater than zero")

    size = [
        max(1, int(math.ceil(n * s / mm)))
        for n, s in zip(binary.GetSize(), binary.GetSpacing())
    ]
    voxels = math.prod(size)
    if voxels > MAX_VOXELS and not allow_large_volume:
        raise ValueError(
            "resampled grid would hold %s voxels at %.4g mm, above the default "
            "limit of %s; raise --resample-mm or pass --allow-large-volume "
            "to attempt it (this may exhaust memory)"
            % (f"{voxels:,}", mm, f"{MAX_VOXELS:,}")
        )
    field = antialias_for_grid(binary, mm)
    r = sitk.ResampleImageFilter()
    r.SetOutputSpacing((mm, mm, mm))
    r.SetSize(size)
    r.SetOutputOrigin(field.GetOrigin())
    r.SetOutputDirection(field.GetDirection())
    r.SetInterpolator(sitk.sitkLinear)
    r.SetDefaultPixelValue(0.0)  # outside
    out = r.Execute(field)

    # The coarser grid may clip the object at its outer face; a border of
    # background keeps the isosurface closed. Callers that go on to mesh the
    # result want this. Callers that go on to run more morphology do not: the
    # border would make the mask stop touching the volume edge, and the
    # field-of-view truncation warning would silently never fire.
    if pad_border:
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
