"""Fuse two scans of the same anatomy into one surface.

Two CT studies of one skull -- say a facial CT that stops mid-vault and a sinus CT
that stops at the maxilla -- together cover more than either alone. Fusing them
means registering one onto the other and unioning the bone.

The union happens on the **masks**, not the meshes. Boolean-unioning two shells
leaves a seam ridge wherever they disagree by a fraction of a millimetre; unioning
the solids and running marching cubes once gives a single continuous surface.

The common grid is **isotropic**. Resampling onto either scan's native grid would
force the other to inherit its slice stepping: merging a 0.3 mm sinus CT onto a
0.8 mm facial CT grid terraced the vault with the coarser scan's slice pitch.

Refusing bad input is the point of this module. Registration always returns *a*
transform -- FFT and ICP cannot fail, they can only converge somewhere useless --
so the answer must be checked, not trusted. Two unrelated scans produce a
confident-looking transform and a mesh made of two skulls stuck together at
random. The gates below exist to make that outcome an error instead.
"""

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import pydicom
import SimpleITK as sitk

from . import pipeline, registration, segment, surface, volume as volume_mod
from .presets import Preset
from .registration import RegistrationResult
from .series import Series

Logger = Callable[[str], None]

#: Voxel size of the fused grid. Fine enough to keep sub-millimetre bone, coarse
#: enough that a whole head fits in memory: a head at 0.4 mm is ~170M voxels.
DEFAULT_GRID_MM = 0.4

#: Surface defaults for the fused grid. Lighter than `convert`'s on purpose: an
#: isotropic grid has no slice terracing to smooth away, so heavy smoothing would
#: only cost detail.
DEFAULT_SMOOTH_ITERS = 6
DEFAULT_PASSBAND = 0.18
DEFAULT_TARGET_FACES = 2_000_000
DEFAULT_POST_SMOOTH_ITERS = 6

# --- acceptance gates -------------------------------------------------------
# Calibrated against four true pairs and one deliberate impostor (a metal bar
# phantom given the skull's own patient identity):
#
#   pair                                   min overlap   dice    rms
#   2024 x 2023 (cross-study, cross-kernel)    0.922    0.855   0.123
#   2024 x 2023 J30s (soft kernel)             0.935    0.835   0.205
#   2024 x 2024 Hr40 3 mm                      0.969    0.852   0.187
#   2024 x 2024 sagittal reformat              0.991    0.932   0.090
#   bar phantom, unrelated anatomy             0.270    0.323   0.432
#
#: Symmetric surface overlap: the *weaker* of the two directions, each restricted
#: to the other scan's field of view. The single strongest discriminator.
MIN_SURFACE_OVERLAP = 0.60
#: Agreement of the two bone masks where both scans have data.
MIN_SHARED_FOV_DICE = 0.55
#: A sanity bound, not a discriminator: the impostor's residual (0.43 mm) sits
#: comfortably inside it. ICP drives *some* residual down no matter what it is
#: fitting, which is exactly why overlap and Dice carry the decision.
MAX_INLIER_RMS_MM = 1.0
#: If the scans barely see the same space, there is nothing to register on.
MIN_SHARED_FOV_MM3 = 20_000.0


class MergeError(RuntimeError):
    pass


@dataclass
class MergeResult:
    output_path: str
    triangles: int
    vertices: int
    bounds_mm: tuple[float, ...]
    grid_mm: float
    grid_size: tuple[int, int, int]
    registration: RegistrationResult
    volume_a_mm3: float
    volume_b_mm3: float
    volume_union_mm3: float
    surface_components: int
    seconds: float
    warnings: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)


def _patient_fingerprint(path: str) -> str | None:
    """A hash of the patient identifiers. Compared, never stored or shown."""
    fields = _patient_fields(path)
    if fields is None or not any(fields.values()):
        return None
    joined = "|".join(fields[k] for k in sorted(fields))
    return hashlib.sha256(joined.encode("utf-8", "replace")).hexdigest()[:16]


def _normalise_name(value: str) -> str:
    """`Doe^Jane` and `DOE  JANE ` are the same person.

    DICOM writes names as caret-separated components, and institutions disagree
    about case and padding.
    """
    return " ".join(str(value).replace("^", " ").split()).casefold()


def _patient_fields(path: str) -> dict[str, str] | None:
    try:
        ds = pydicom.dcmread(path, stop_before_pixels=True)
    except Exception:  # noqa: BLE001
        return None
    return {
        "id": str(getattr(ds, "PatientID", "")).strip(),
        "issuer": str(getattr(ds, "IssuerOfPatientID", "")).strip(),
        "name": _normalise_name(getattr(ds, "PatientName", "")),
        "birth_date": str(getattr(ds, "PatientBirthDate", "")).strip(),
        "sex": str(getattr(ds, "PatientSex", "")).strip().upper(),
    }


@dataclass
class PatientMatch:
    verdict: str  # "same" | "different" | "unknown"
    matched: list[str]
    conflicts: list[str]


def compare_patients(path_a: str, path_b: str) -> PatientMatch:
    """Decide whether two instances describe the same person.

    ``PatientID`` equality is *not* the test. Medical record numbers are scoped to
    the issuing institution -- which is why DICOM carries ``IssuerOfPatientID`` --
    so the same person legitimately has different IDs at different hospitals. Two
    real studies of one skull differed in both ``PatientID`` (7 vs 15 characters)
    and ``PatientName`` formatting, while agreeing exactly on birth date and sex.

    Conflicting demographics are treated as proof of different people; a matching
    ID, or a matching name plus birth date, as proof of the same one. Anything
    else is ``unknown`` and left to the caller.

    Only field *names* ever appear in the result. Values are never returned,
    logged, or raised.
    """
    fa, fb = _patient_fields(path_a), _patient_fields(path_b)
    if fa is None or fb is None or not any(fa.values()) or not any(fb.values()):
        return PatientMatch("unknown", [], [])

    def both(key):
        return bool(fa[key]) and bool(fb[key])

    def agree(key):
        return both(key) and fa[key] == fb[key]

    def conflict(key):
        return both(key) and fa[key] != fb[key]

    matched = [k for k in ("id", "name", "birth_date", "sex") if agree(k)]

    # Demographics are the strong evidence: they are not institution-scoped.
    hard_conflicts = [k for k in ("birth_date", "sex", "name") if conflict(k)]
    if hard_conflicts:
        return PatientMatch("different", matched, hard_conflicts)

    # An ID match only proves identity if the issuers do not disagree.
    if agree("id") and not conflict("issuer"):
        return PatientMatch("same", matched, [])

    if agree("name") and (agree("birth_date") or not both("birth_date")):
        return PatientMatch("same", matched, [])

    # IDs differ and nothing corroborates that this is one person.
    if conflict("id"):
        return PatientMatch("different", matched, ["id"])

    return PatientMatch("unknown", matched, [])


def check_compatible(a: Series, b: Series, force: bool = False) -> list[str]:
    """Refuse pairs that cannot sensibly be fused. Returns non-fatal warnings."""
    warnings: list[str] = []

    if a.uid == b.uid and a.part == b.part:
        raise MergeError("both inputs resolve to the same series (%s)" % a.ident)

    if a.modality != b.modality:
        msg = ("modalities differ (%s vs %s); intensities are not comparable"
               % (a.modality, b.modality))
        if not force:
            raise MergeError(msg + ". Pass --force if you really mean it.")
        warnings.append(msg)

    if not a.files or not b.files:
        return warnings

    match = compare_patients(a.files[0], b.files[0])
    if match.verdict == "different":
        fields = ", ".join(match.conflicts)
        msg = ("the two series appear to be different patients (%s %s; values not "
               "shown)" % (fields, "differ" if len(match.conflicts) > 1 else "differs"))
        if not force:
            raise MergeError(msg + ". Pass --force if these really are the same body.")
        warnings.append(msg)
    elif match.verdict == "unknown":
        warnings.append("cannot verify that both scans are of the same person "
                        "(no comparable patient identifiers)")

    return warnings


def check_registration(result: RegistrationResult, force: bool = False) -> None:
    """Decide whether a transform describes the same anatomy, or garbage.

    Registration cannot report failure on its own: FFT always has a peak and ICP
    always converges to something. Every one of these numbers is near-perfect for
    a true match and hopeless for unrelated scans.
    """
    problems = []
    if result.shared_fov_mm3 < MIN_SHARED_FOV_MM3:
        problems.append(
            "the scans share only %.0f cm3 of imaged space (need %.0f)"
            % (result.shared_fov_mm3 / 1000.0, MIN_SHARED_FOV_MM3 / 1000.0))
    if result.surface_overlap < MIN_SURFACE_OVERLAP:
        problems.append(
            "surfaces agree over only %.1f%% of the shared field of view "
            "(%.1f%% moving->fixed, %.1f%% fixed->moving; need %.0f%% both ways)"
            % (100 * result.surface_overlap,
               100 * result.overlap_moving_in_fixed,
               100 * result.overlap_fixed_in_moving,
               100 * MIN_SURFACE_OVERLAP))
    if result.inlier_rms_mm > MAX_INLIER_RMS_MM:
        problems.append(
            "point-to-plane residual %.2f mm exceeds %.2f mm"
            % (result.inlier_rms_mm, MAX_INLIER_RMS_MM))
    if result.shared_fov_dice < MIN_SHARED_FOV_DICE:
        problems.append(
            "bone agreement in the shared field of view is %.3f (need %.2f)"
            % (result.shared_fov_dice, MIN_SHARED_FOV_DICE))

    if not problems:
        return
    message = (
        "registration failed its quality gates, so the two scans probably do not "
        "show the same anatomy:\n  - " + "\n  - ".join(problems)
        + "\n" + result.summary()
    )
    if force:
        return
    raise MergeError(message + "\nPass --force to fuse them anyway.")


def _antialias(image: sitk.Image, grid_mm: float) -> sitk.Image:
    """Blur only the axes being downsampled, or thin bone aliases away.

    A vector sigma with zeros in it is rejected by SimpleITK, so the axes are
    filtered one at a time.
    """
    out = sitk.Cast(image, sitk.sitkFloat32)
    for axis, spacing in enumerate(image.GetSpacing()):
        if grid_mm <= spacing:
            continue
        blur = sitk.RecursiveGaussianImageFilter()
        blur.SetDirection(axis)
        blur.SetSigma(segment._ANTIALIAS_SIGMA_FACTOR * grid_mm)
        out = blur.Execute(out)
    return out


def _corners(image: sitk.Image) -> np.ndarray:
    size = np.asarray(image.GetSize()) - 1
    pts = []
    for i in (0, size[0]):
        for j in (0, size[1]):
            for k in (0, size[2]):
                pts.append(image.TransformIndexToPhysicalPoint((int(i), int(j), int(k))))
    return np.asarray(pts, dtype=float)


def _common_grid(mask_a, mask_b, transform, grid_mm):
    rot, trans = transform[:3, :3], transform[:3, 3]
    pts = np.vstack([_corners(mask_a), (rot @ _corners(mask_b).T).T + trans])
    lo = np.floor(pts.min(0) / grid_mm) * grid_mm - 2 * grid_mm
    hi = np.ceil(pts.max(0) / grid_mm) * grid_mm + 2 * grid_mm
    size = np.ceil((hi - lo) / grid_mm).astype(int) + 1
    return tuple(int(v) for v in size), lo


def _resample_field(image, size, origin, grid_mm, transform):
    r = sitk.ResampleImageFilter()
    r.SetOutputSpacing((grid_mm, grid_mm, grid_mm))
    r.SetSize([int(v) for v in size])
    r.SetOutputOrigin([float(v) for v in origin])
    r.SetOutputDirection(np.eye(3).ravel().tolist())
    r.SetInterpolator(sitk.sitkLinear)
    r.SetDefaultPixelValue(0.0)
    r.SetTransform(transform)
    return r.Execute(image)


def merge(
    series_a: Series,
    series_b: Series,
    preset: Preset,
    output_path: str,
    threshold: float | None = None,
    grid_mm: float = DEFAULT_GRID_MM,
    smooth_iters: int = DEFAULT_SMOOTH_ITERS,
    passband: float = DEFAULT_PASSBAND,
    target_faces: int = DEFAULT_TARGET_FACES,
    post_smooth_iters: int = DEFAULT_POST_SMOOTH_ITERS,
    force: bool = False,
    log: Logger | None = None,
) -> MergeResult:
    started = time.time()

    def say(msg: str) -> None:
        if log:
            log(msg)

    def step(label, fn):
        t = time.time()
        out = fn()
        say("  %-36s %6.1fs" % (label, time.time() - t))
        return out

    warnings = check_compatible(series_a, series_b, force=force)

    say("fixed  %s  %s  (%d slices)" % (series_a.ident, series_a.label(), series_a.n_slices))
    say("moving %s  %s  (%d slices)" % (series_b.ident, series_b.label(), series_b.n_slices))

    vol_a = volume_mod.load(series_a)
    vol_b = volume_mod.load(series_b)
    warnings.extend(volume_mod.warnings_for(vol_a))
    warnings.extend(volume_mod.warnings_for(vol_b))

    value_a, source_a = pipeline.resolve_threshold(vol_a.image, series_a, preset, threshold)
    value_b, source_b = pipeline.resolve_threshold(vol_b.image, series_b, preset, threshold)
    say("threshold: fixed %.1f (%s), moving %.1f (%s)" % (value_a, source_a, value_b, source_b))

    mask_a = step("segment fixed", lambda: pipeline.build_mask(vol_a.image, preset, value_a))
    mask_b = step("segment moving", lambda: pipeline.build_mask(vol_b.image, preset, value_b))

    say("registering ...")
    reg = registration.rigid_register(mask_a, mask_b, log=lambda m: say("  " + m))
    for line in reg.summary().splitlines():
        say("  " + line.strip() if line.startswith(" ") else "  " + line)
    check_registration(reg, force=force)

    finest = min(min(vol_a.spacing), min(vol_b.spacing))
    if grid_mm > finest:
        warnings.append(
            "the fused grid is %.2f mm but the finest input voxel is %.3f mm; "
            "structures thinner than the grid are lost. Lower --grid-mm to keep "
            "them, at cubic cost in memory." % (grid_mm, finest)
        )

    size, origin = _common_grid(mask_a, mask_b, reg.transform, grid_mm)
    voxels = int(np.prod(size))
    say("fused grid %s at %.2f mm isotropic (%.0f M voxels)"
        % ("x".join(str(v) for v in size), grid_mm, voxels / 1e6))
    if voxels > 800e6:
        raise MergeError(
            "the fused grid would hold %.0f M voxels at %.2f mm; raise --grid-mm"
            % (voxels / 1e6, grid_mm))

    identity = sitk.Transform(3, sitk.sitkIdentity)
    inverse = registration.inverse_transform(reg.transform)

    field_a = step("resample fixed",
                   lambda: _resample_field(_antialias(mask_a, grid_mm), size, origin,
                                           grid_mm, identity))
    field_b = step("resample moving",
                   lambda: _resample_field(_antialias(mask_b, grid_mm), size, origin,
                                           grid_mm, inverse))

    # Union of occupancy, not of labels: keeps the sub-voxel boundary each scan
    # carries, so the fused surface is not quantised to the grid.
    fused = sitk.Clamp(sitk.Maximum(field_a, field_b), sitk.sitkFloat32, 0.0, 1.0)

    voxel_mm3 = grid_mm ** 3
    vol_a_mm3 = float((sitk.GetArrayViewFromImage(field_a) > 0.5).sum()) * voxel_mm3
    vol_b_mm3 = float((sitk.GetArrayViewFromImage(field_b) > 0.5).sum()) * voxel_mm3
    vol_u_mm3 = float((sitk.GetArrayViewFromImage(fused) > 0.5).sum()) * voxel_mm3
    say("bone: fixed %.0f cm3 | moving %.0f cm3 | fused %.0f cm3"
        % (vol_a_mm3 / 1000, vol_b_mm3 / 1000, vol_u_mm3 / 1000))

    fused = segment.pad(fused, 1)

    affine = surface.index_to_physical(fused)
    poly = step("marching cubes",
                lambda: surface.marching_cubes(surface.to_vtk_image(fused),
                                               segment.ISO_OCCUPANCY))
    if poly.GetNumberOfPolys() == 0:
        raise MergeError("the fused volume produced no surface")
    say("  raw triangles %s" % f"{poly.GetNumberOfPolys():,}")

    poly = step("smooth", lambda: surface.smooth(poly, smooth_iters, passband))
    poly, shells = step("largest component", lambda: surface.largest_component(poly))
    say("  surface shells %d (kept 1)" % shells)
    poly = step("decimate", lambda: surface.decimate(poly, target_faces))
    poly = step("post-smooth", lambda: surface.smooth(poly, post_smooth_iters, passband))
    poly = surface.transform(poly, affine)
    poly = surface.compute_normals(poly)

    boundary, nonmanifold = surface.count_defects(poly)
    if boundary or nonmanifold:
        warnings.append("fused surface has %d boundary and %d non-manifold edge(s)"
                        % (boundary, nonmanifold))

    surface.write(poly, output_path)

    provenance = {
        "fixed": {"uid": series_a.uid, "ident": series_a.ident,
                  "description": series_a.description, "slices": series_a.n_slices,
                  "spacing_mm": list(vol_a.spacing), "threshold": value_a},
        "moving": {"uid": series_b.uid, "ident": series_b.ident,
                   "description": series_b.description, "slices": series_b.n_slices,
                   "spacing_mm": list(vol_b.spacing), "threshold": value_b},
        "grid_mm": grid_mm,
        "transform_moving_to_fixed": reg.transform.tolist(),
        "rotation_deg": reg.rotation_deg,
        "registration": {
            "inlier_rms_mm": reg.inlier_rms_mm,
            "inlier_median_mm": reg.inlier_median_mm,
            "surface_overlap": reg.surface_overlap,
            "shared_fov_dice": reg.shared_fov_dice,
            "shared_fov_mm3": reg.shared_fov_mm3,
        },
        "coordinate_system": "LPS (DICOM patient space of the fixed series)",
        "forced": bool(force),
    }

    return MergeResult(
        output_path=output_path,
        triangles=int(poly.GetNumberOfPolys()),
        vertices=int(poly.GetNumberOfPoints()),
        bounds_mm=surface.bounds_mm(poly),
        grid_mm=grid_mm,
        grid_size=size,
        registration=reg,
        volume_a_mm3=vol_a_mm3,
        volume_b_mm3=vol_b_mm3,
        volume_union_mm3=vol_u_mm3,
        surface_components=shells,
        seconds=time.time() - started,
        warnings=warnings,
        provenance=provenance,
    )
