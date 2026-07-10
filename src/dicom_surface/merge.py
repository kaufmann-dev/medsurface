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

Geometry alone cannot do it. Measured: a skull uniformly scaled by 3%, which is
well inside person-to-person variation, passes both geometric gates (overlap
0.989, dice 0.675) because rigid registration parks it neatly on top. So the
demographics are checked too -- not out of bureaucracy, but because the geometry
demonstrably cannot tell two similar bodies apart. Scans with no identifiers at
all (de-identified data) warn and proceed.
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field
from typing import Any, Callable

import numpy as np
import pydicom
import SimpleITK as sitk

from . import pipeline, registration, segment, surface, volume as volume_mod
from .defaults import DEFAULT_MERGE_GRID_MM
from .presets import ANATOMICAL, Preset, PrintProfile
from .registration import RegistrationResult
from .series import Series

Logger = Callable[[str], None]

#: The surface stage is the preset's, unchanged. It is tempting to smooth a fused
#: grid more lightly, on the theory that an isotropic grid has no slice terracing
#: to remove. It has: the terracing is baked into each scan's *mask* by its own
#: slice pitch, long before anything is resampled. Resampling a 0.8 mm staircase
#: onto a 0.4 mm grid samples the staircase more finely; it does not flatten it.
#: A fused surface must be smoothed exactly as a single-scan one is, or it looks
#: visibly rougher than the scans it was built from.

# --- acceptance gates -------------------------------------------------------
# Empirically selected from four positive reconstructions/repeat studies of one
# skull and one deliberate impostor (a metal bar given the skull's identity):
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
#: If the scans barely see the same space, there is nothing to register on --
#: and Dice over a few cubic centimetres proves nothing either way.
MIN_SHARED_FOV_MM3 = 20_000.0

# There is deliberately no gate on the ICP residual. It never discriminated: the
# impostor scored 0.43 mm and a 5%-oversized skull 0.41 mm, both well inside any
# plausible bound. ICP drives *some* residual down no matter what it is fitting.
# The residual is reported, because it is informative; it is not evidence.


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

    Conflicting demographics are treated as evidence of different people; a
    matching ID, or a matching name plus birth date, as sufficient corroboration
    that the studies come from the same person. Neither is proof -- identifiers
    are re-issued, pseudonymised and mistyped -- so `--force` exists. Anything
    else is ``unknown`` and left to the caller.

    This is the only place in the package that reads patient-identifying fields,
    and it reads them solely to compare them. Only field *names* ever appear in
    the result. Values are never returned, logged, raised, or written to
    provenance.
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
        raise MergeError("both inputs resolve to the same series (%s)" % a.uid)

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
        if match.conflicts == ["id"]:
            msg += (". If these are de-identified studies of one person that were "
                    "given different pseudonyms, pass --force")
        if not force:
            raise MergeError(msg + ". Pass --force if these really are the same body.")
        warnings.append(msg)
    elif match.verdict == "unknown":
        # De-identified data lands here, warns, and proceeds.
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
    grid_mm: float = DEFAULT_MERGE_GRID_MM,
    smooth_iters: int | None = None,
    passband: float | None = None,
    target_faces: int | None = None,
    post_smooth_iters: int | None = None,
    print_profile: PrintProfile = ANATOMICAL,
    force: bool = False,
    log: Logger | None = None,
) -> MergeResult:
    started = time.time()

    # Fall back to the preset, so a fused surface is finished exactly as a
    # single-scan one is.
    smooth_iters = preset.smooth_iters if smooth_iters is None else smooth_iters
    passband = preset.passband if passband is None else passband
    target_faces = preset.target_faces if target_faces is None else target_faces
    post_smooth_iters = (preset.post_smooth_iters if post_smooth_iters is None
                         else post_smooth_iters)

    def say(msg: str) -> None:
        if log:
            log(msg)

    def step(label, fn):
        say("%s ..." % label)
        t = time.time()
        out = fn()
        say("  %-36s %6.1fs" % (label, time.time() - t))
        return out

    warnings = check_compatible(series_a, series_b, force=force)

    say("fixed  DICOM #%s  %s  (%d slices)"
        % (series_a.series_number if series_a.series_number is not None else "-",
           series_a.label(), series_a.n_slices))
    say("moving DICOM #%s  %s  (%d slices)"
        % (series_b.series_number if series_b.series_number is not None else "-",
           series_b.label(), series_b.n_slices))

    vol_a = step("load fixed DICOM volume", lambda: volume_mod.load(series_a))
    vol_b = step("load moving DICOM volume", lambda: volume_mod.load(series_b))
    warnings.extend(volume_mod.warnings_for(vol_a))
    warnings.extend(volume_mod.warnings_for(vol_b))

    value_a, source_a = step(
        "resolve fixed threshold",
        lambda: pipeline.resolve_threshold(vol_a.image, series_a, preset, threshold),
    )
    value_b, source_b = step(
        "resolve moving threshold",
        lambda: pipeline.resolve_threshold(vol_b.image, series_b, preset, threshold),
    )
    say("threshold: fixed %.1f (%s), moving %.1f (%s)" % (value_a, source_a, value_b, source_b))

    # Register on anatomy, always. Thickening both scans would inflate Dice and
    # surface overlap -- the very numbers used by the empirical gates --
    # so a print profile would quietly make `merge` easier to fool.
    mask_a = step("segment fixed", lambda: pipeline.build_mask(vol_a.image, preset, value_a))
    mask_b = step("segment moving", lambda: pipeline.build_mask(vol_b.image, preset, value_b))

    say("registering ...")
    reg = registration.rigid_register(mask_a, mask_b, log=lambda m: say("  " + m))
    for line in reg.summary().splitlines():
        say("  " + line.strip() if line.startswith(" ") else "  " + line)
    check_registration(reg, force=force)

    if print_profile != ANATOMICAL:
        say("print profile %s: %s" % (print_profile.name, print_profile.description))
        # Profile each scan before fractional-occupancy fusion. This preserves
        # each scan's post-profile boundary, but the complete transform is not
        # distributive: closing, island filtering, grid conversion, and thin-set
        # selection can differ from profiling the fused mask. Controlled overlap
        # cases show that this order can over-thicken shared structures; changing
        # it requires a separate geometry task because fusion-first closing can
        # also seal an anatomical inter-scan gap.
        mask_a = step("re-segment fixed for printing",
                      lambda: pipeline.build_mask(vol_a.image, preset, value_a,
                                                  print_profile, say))
        mask_b = step("re-segment moving for printing",
                      lambda: pipeline.build_mask(vol_b.image, preset, value_b,
                                                  print_profile, say))
        # Per scan: the two rarely share a voxel grid, so each gets its own
        # printability grid and each can fall short of the request differently.
        for label, mask, spacing in (("fixed", mask_a, vol_a.spacing),
                                     ("moving", mask_b, vol_b.spacing)):
            for message in pipeline.grid_warnings(mask, print_profile, spacing):
                warnings.append("%s scan: %s" % (label, message))

    finest = min(min(vol_a.spacing), min(vol_b.spacing))
    if grid_mm > finest:
        warnings.append(
            "the fused grid is %.2f mm but the finest input voxel is %.3f mm; "
            "structures thinner than the grid are lost. Lower --grid-mm to keep "
            "them, at cubic cost in memory." % (grid_mm, finest)
        )

    size, origin = step(
        "plan fused grid",
        lambda: _common_grid(mask_a, mask_b, reg.transform, grid_mm),
    )
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
                   lambda: _resample_field(segment.antialias_for_grid(mask_a, grid_mm), size, origin,
                                           grid_mm, identity))
    field_b = step("resample moving",
                   lambda: _resample_field(segment.antialias_for_grid(mask_b, grid_mm), size, origin,
                                           grid_mm, inverse))

    # Union of occupancy, not of labels: keeps the sub-voxel boundary each scan
    # carries, so the fused surface is not quantised to the grid.
    fused = step(
        "fuse occupancy fields",
        lambda: sitk.Clamp(sitk.Maximum(field_a, field_b), sitk.sitkFloat32, 0.0, 1.0),
    )

    voxel_mm3 = grid_mm ** 3

    def measure_volumes() -> tuple[float, float, float]:
        return (
            float((sitk.GetArrayViewFromImage(field_a) > 0.5).sum()) * voxel_mm3,
            float((sitk.GetArrayViewFromImage(field_b) > 0.5).sum()) * voxel_mm3,
            float((sitk.GetArrayViewFromImage(fused) > 0.5).sum()) * voxel_mm3,
        )

    vol_a_mm3, vol_b_mm3, vol_u_mm3 = step("measure fused volumes", measure_volumes)
    say("bone: fixed %.0f cm3 | moving %.0f cm3 | fused %.0f cm3"
        % (vol_a_mm3 / 1000, vol_b_mm3 / 1000, vol_u_mm3 / 1000))

    if print_profile.min_feature_mm > 0:
        solid = sitk.BinaryThreshold(fused, segment.ISO_OCCUPANCY, 1e9, 1, 0)
        warnings.extend(step("thin-material check",
                             lambda: pipeline.thin_material_warning(solid, print_profile)))

    fused = segment.pad(fused, 1)

    affine = surface.index_to_physical(fused)
    poly = step("marching cubes",
                lambda: surface.marching_cubes(surface.to_vtk_image(fused),
                                               segment.ISO_OCCUPANCY))
    if poly.GetNumberOfPolys() == 0:
        raise MergeError("the fused volume produced no surface")
    say("  raw triangles %s" % f"{poly.GetNumberOfPolys():,}")

    finished = pipeline.finish_surface(
        poly,
        smooth_iters=smooth_iters,
        passband=passband,
        target_faces=target_faces,
        post_smooth_iters=post_smooth_iters,
        keep_largest_component=True,
        step=step,
        log=say,
    )
    poly = finished.poly
    shells = finished.surface_components
    warnings.extend(finished.warnings)
    poly = step("index -> patient space (LPS)", lambda: surface.transform(poly, affine))
    poly = step("normals", lambda: surface.compute_normals(poly))

    boundary, nonmanifold = step("check surface defects", lambda: surface.count_defects(poly))
    if boundary or nonmanifold:
        warnings.append("fused surface has %d boundary and %d non-manifold edge(s)"
                        % (boundary, nonmanifold))

    step("write mesh", lambda: surface.write(poly, output_path))

    provenance = {
        "fixed": {"uid": series_a.uid, "series_number": series_a.series_number,
                  "series_orientation_part": [series_a.part, series_a.n_parts],
                  "description": series_a.description, "slices": series_a.n_slices,
                  "spacing_mm": list(vol_a.spacing), "threshold": value_a},
        "moving": {"uid": series_b.uid, "series_number": series_b.series_number,
                   "series_orientation_part": [series_b.part, series_b.n_parts],
                   "description": series_b.description, "slices": series_b.n_slices,
                   "spacing_mm": list(vol_b.spacing), "threshold": value_b},
        "grid_mm": grid_mm,
        "print_profile": asdict(print_profile),
        "surface_finishing": finished.provenance,
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
