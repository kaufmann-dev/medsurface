"""Command line interface."""

from __future__ import annotations

import argparse
import json
import sys

from dataclasses import replace

from . import presets as presets_mod
from .merge import DEFAULT_GRID_MM
from .presets import PRINT_PROFILES
from . import pipeline, series as series_mod, validate as validate_mod
from .presets import PRESETS


def _log(msg: str) -> None:
    print(msg, flush=True)


def _quiet(msg: str) -> None:  # noqa: ARG001
    pass


def _warn(msg: str) -> None:
    print("warning: %s" % msg, file=sys.stderr, flush=True)


def _resolve_print_profile(args: argparse.Namespace):
    """Compose the print profile, letting explicit flags beat it.

    ``build_mask`` takes ``max(preset.closing_mm, profile.closing_mm)``, so an
    explicit ``--closing-mm`` has to be written to *both* sides or the profile
    would silently win whenever it asks for more.
    """
    profile = presets_mod.get_print_profile(args.print_profile)

    if getattr(args, "thicken_mm", None) is not None:
        profile = replace(profile, thicken_mm=args.thicken_mm)
    if args.closing_mm is not None:
        profile = replace(profile, closing_mm=args.closing_mm)
    if args.min_island_mm3 is not None:
        profile = replace(profile, min_island_mm3=args.min_island_mm3)
    return profile


# --------------------------------------------------------------------- list
def cmd_list(args: argparse.Namespace) -> int:
    found = series_mod.discover(args.dicom_dir)
    if not found:
        print("no DICOM instances found under %s" % args.dicom_dir, file=sys.stderr)
        return 1

    if args.json:
        payload = [
            {
                "ident": s.ident,
                "uid": s.uid,
                "part": s.part,
                "n_parts": s.n_parts,
                "series_number": s.series_number,
                "modality": s.modality,
                "description": s.description,
                "slices": s.n_slices,
                "rows": s.rows,
                "columns": s.columns,
                "pixel_spacing": list(s.pixel_spacing) if s.pixel_spacing else None,
                "slice_spacing": s.slice_spacing,
                "plane": s.plane,
                "kernel": s.kernel,
                "sharp_kernel": s.sharp_kernel,
                "spacing_uniform": s.spacing_uniform,
                "spacing_spread_mm": s.spacing_spread_mm,
                "localizer": s.is_localizer,
                "usable": s.usable,
                "unusable_reason": s.unusable_reason,
            }
            for s in found
        ]
        print(json.dumps(payload, indent=2))
        return 0

    best = series_mod.rank([s for s in found if s.usable])
    recommended = best[0] if best else None

    header = ("#", "MOD", "DESCRIPTION", "SLICES", "VOXEL mm", "PLANE", "NOTES")
    print("%-6s %-4s %-32s %6s %-22s %-9s %s" % header)
    print("-" * 114)
    split_seen = False
    for s in found:
        voxel = "-"
        if s.pixel_spacing and s.slice_spacing:
            voxel = "%.3f x %.3f x %.3f" % (s.pixel_spacing[0], s.pixel_spacing[1], s.slice_spacing)
        notes = []
        if recommended is not None and s is recommended:
            notes.append("<- default")
        reason = s.unusable_reason
        if reason:
            notes.append(reason)
        if s.n_parts > 1:
            split_seen = True
            notes.append("orientation %d of %d in this UID" % (s.part, s.n_parts))
        if s.sharp_kernel:
            notes.append("sharp kernel %s" % s.kernel)
        print("%-6s %-4s %-32s %6d %-22s %-9s %s" % (
            s.ident,
            s.modality,
            (s.description or "(none)")[:32],
            s.n_slices,
            voxel,
            s.plane if s.usable else "-",
            ", ".join(notes),
        ))
    print()
    if split_seen:
        print("Some SeriesInstanceUIDs hold more than one orientation and were split;")
        print("select those with their dotted ident, e.g. --series 1021.1")
        print()
    print("Convert the default with:  dicom-surface convert %s -o out.stl" % args.dicom_dir)
    return 0


# ------------------------------------------------------------------ convert
def cmd_convert(args: argparse.Namespace) -> int:
    log = _quiet if args.quiet else _log

    found = series_mod.discover(args.dicom_dir)
    if not found:
        print("no DICOM instances found under %s" % args.dicom_dir, file=sys.stderr)
        return 1

    try:
        chosen = series_mod.select(found, args.series)
    except ValueError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    log("series %s  %s  (%d slices)" % (chosen.ident, chosen.label(), chosen.n_slices))

    preset = presets_mod.get(args.preset)
    preset = presets_mod.override(
        preset,
        median_mm=args.median_mm,
        closing_mm=args.closing_mm,
        opening_mm=args.opening_mm,
        min_island_mm3=args.min_island_mm3,
        resample_mm=args.resample_mm,
        smooth_iters=args.smooth_iters,
        passband=args.passband,
        target_faces=args.target_faces,
        post_smooth_iters=args.post_smooth_iters,
        keep_largest_island=False if args.all_islands else None,
        keep_largest_component=False if args.all_components else None,
    )

    threshold = None
    if args.threshold is not None:
        if args.threshold == "auto":
            preset = presets_mod.override(preset, threshold="auto")
        else:
            try:
                threshold = float(args.threshold)
            except ValueError:
                print("error: --threshold must be a number or 'auto'", file=sys.stderr)
                return 2

    try:
        result = pipeline.convert(
            series=chosen,
            preset=preset,
            output_path=args.output,
            threshold=threshold,
            cap_field_of_view=not args.no_cap,
            print_profile=_resolve_print_profile(args),
            log=log,
        )
    except pipeline.ModalityMismatch as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2
    except ValueError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 1

    for w in result.warnings:
        _warn(w)

    log("")
    log("wrote %s" % result.output_path)
    log("  triangles %s   vertices %s   %.1fs"
        % (f"{result.triangles:,}", f"{result.vertices:,}", result.seconds))

    report = None
    if not args.no_validate:
        report = validate_mod.validate(result.output_path,
                                       self_intersections=args.self_intersections)
        log("")
        log("quality:")
        log(validate_mod.summarise(report))
        if not report["watertight"]:
            _warn("output is not watertight; try 'dicom-surface repair' "
                  "or relax --closing-mm")

    if args.json:
        payload = {
            "result": {
                "output": result.output_path,
                "triangles": result.triangles,
                "vertices": result.vertices,
                "bounds_mm": list(result.bounds_mm),
                "seconds": result.seconds,
                "capped_field_of_view": result.capped_field_of_view,
                "labelmap_components": result.labelmap_components,
                "surface_components": result.surface_components,
                "warnings": result.warnings,
            },
            "provenance": result.provenance,
            "quality": report,
        }
        with open(args.json, "w") as fh:
            json.dump(payload, fh, indent=2)
        log("wrote %s" % args.json)

    return 0


# -------------------------------------------------------------------- merge
def cmd_merge(args: argparse.Namespace) -> int:
    from . import merge as merge_mod

    log = _quiet if args.quiet else _log

    dir_b = args.dicom_dir_b or args.dicom_dir_a
    found_a = series_mod.discover(args.dicom_dir_a)
    found_b = found_a if dir_b == args.dicom_dir_a else series_mod.discover(dir_b)
    if not found_a or not found_b:
        print("no DICOM instances found", file=sys.stderr)
        return 1

    try:
        chosen_a = series_mod.select(found_a, args.series_a)
        chosen_b = series_mod.select(found_b, args.series_b)
    except ValueError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    preset = presets_mod.get(args.preset)
    preset = presets_mod.override(
        preset,
        median_mm=args.median_mm,
        closing_mm=args.closing_mm,
        min_island_mm3=args.min_island_mm3,
    )

    threshold = None
    if args.threshold is not None:
        if args.threshold == "auto":
            preset = presets_mod.override(preset, threshold="auto")
        else:
            try:
                threshold = float(args.threshold)
            except ValueError:
                print("error: --threshold must be a number or 'auto'", file=sys.stderr)
                return 2

    try:
        result = merge_mod.merge(
            series_a=chosen_a,
            series_b=chosen_b,
            preset=preset,
            output_path=args.output,
            threshold=threshold,
            grid_mm=args.grid_mm,
            smooth_iters=args.smooth_iters,
            passband=args.passband,
            target_faces=args.target_faces,
            post_smooth_iters=args.post_smooth_iters,
            print_profile=_resolve_print_profile(args),
            force=args.force,
            log=log,
        )
    except merge_mod.MergeError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 3
    except (pipeline.ModalityMismatch, ValueError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    for w in result.warnings:
        _warn(w)

    log("")
    log("wrote %s" % result.output_path)
    log("  triangles %s   vertices %s   %.1fs"
        % (f"{result.triangles:,}", f"{result.vertices:,}", result.seconds))

    report = None
    if not args.no_validate:
        report = validate_mod.validate(result.output_path,
                                       self_intersections=args.self_intersections)
        log("")
        log("quality:")
        log(validate_mod.summarise(report))
        if not report["watertight"]:
            _warn("fused mesh is not watertight; try 'dicom-surface repair'")

    if args.json:
        payload = {
            "result": {
                "output": result.output_path,
                "triangles": result.triangles,
                "vertices": result.vertices,
                "bounds_mm": list(result.bounds_mm),
                "grid_mm": result.grid_mm,
                "grid_size": list(result.grid_size),
                "volume_fixed_mm3": result.volume_a_mm3,
                "volume_moving_mm3": result.volume_b_mm3,
                "volume_fused_mm3": result.volume_union_mm3,
                "seconds": result.seconds,
                "warnings": result.warnings,
            },
            "provenance": result.provenance,
            "quality": report,
        }
        with open(args.json, "w") as fh:
            json.dump(payload, fh, indent=2)
        log("wrote %s" % args.json)

    return 0


# ----------------------------------------------------------------- validate
def cmd_validate(args: argparse.Namespace) -> int:
    report = validate_mod.validate(args.mesh, self_intersections=args.self_intersections)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(args.mesh)
        print(validate_mod.summarise(report))
    return 0 if report["watertight"] else 1


# -------------------------------------------------------------------- repair
def cmd_repair(args: argparse.Namespace) -> int:
    from . import repair as repair_mod

    try:
        stats = repair_mod.repair(args.mesh, args.output, log=_log if not args.quiet else None)
    except repair_mod.RepairUnavailable as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    print("wrote %s" % args.output)
    report = validate_mod.validate(args.output, self_intersections=args.self_intersections)
    print(validate_mod.summarise(report))
    if args.json:
        print(json.dumps({"repair": stats, "quality": report}, indent=2))
    return 0 if report["watertight"] else 1


# ---------------------------------------------------------------- presets
def cmd_presets(_args: argparse.Namespace) -> int:
    print("PRESETS  --  what tissue to extract, and how finely  (--preset)")
    print()
    for name in sorted(PRESETS):
        p = PRESETS[name]
        modality = ", ".join(p.modalities) if p.modalities else "any"
        thr = p.threshold if isinstance(p.threshold, str) else "%g" % p.threshold
        print("  %-12s  [%s]" % (name, modality))
        print("    %s" % p.description)
        triangles = f"{p.target_faces:,}" if p.target_faces else "all"
        print("    threshold=%s  median=%.1fmm  closing=%.1fmm  smooth=%d  triangles=%s"
              % (thr, p.median_mm, p.closing_mm, p.smooth_iters, triangles))
        print()

    print("PRINT PROFILES  --  what a printer needs  (--print-profile)")
    print()
    print("  Composes with any preset: raises its closing and island filter, and")
    print("  thickens walls. Orthogonal to `bone-print`, which is a triangle budget.")
    print()
    for name in sorted(PRINT_PROFILES):
        p = PRINT_PROFILES[name]
        print("  %-12s" % name)
        print("    %s" % p.description)
        if p == presets_mod.ANATOMICAL:
            print("    (the default: no geometric changes at all)")
        else:
            print("    closing>=%.1fmm  thicken=%.1fmm radius  islands>=%.0fmm3  "
                  "min feature=%.1fmm"
                  % (p.closing_mm, p.thicken_mm, p.min_island_mm3, p.min_feature_mm))
        print()
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="dicom-surface",
        description="Turn a DICOM series into a watertight 3D surface mesh.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    pl = sub.add_parser("list", help="show every series in a DICOM directory")
    pl.add_argument("dicom_dir")
    pl.add_argument("--json", action="store_true", help="machine-readable output")
    pl.set_defaults(func=cmd_list)

    pp = sub.add_parser("presets", help="describe the built-in presets")
    pp.set_defaults(func=cmd_presets)

    pc = sub.add_parser("convert", help="extract a surface mesh from a DICOM series")
    pc.add_argument("dicom_dir")
    pc.add_argument("-o", "--output", required=True, help="output .stl/.ply/.obj/.vtp")
    pc.add_argument("--series", help="series number, UID, or description substring")
    pc.add_argument("--preset", default="bone", choices=sorted(PRESETS))
    pc.add_argument("--threshold", help="intensity (HU for CT) or 'auto' for Otsu")
    pc.add_argument("--median-mm", type=float, help="despeckle kernel extent, mm")
    pc.add_argument("--closing-mm", type=float, help="pore-sealing kernel extent, mm")
    pc.add_argument("--opening-mm", type=float, help="bridge-breaking kernel extent, mm")
    pc.add_argument("--min-island-mm3", type=float, help="drop blobs smaller than this")
    pc.add_argument("--print-profile", default="anatomical", choices=sorted(PRINT_PROFILES),
                    help="prepare the mesh for a printer (default: %(default)s, which "
                         "makes no printability changes)")
    pc.add_argument("--thicken-mm", type=float,
                    help="grow every surface outward by this radius, thickening walls. "
                         "Also enlarges the model's outer dimensions.")
    pc.add_argument("--all-islands", action="store_true",
                    help="keep every labelmap island, not just the largest")
    pc.add_argument("--all-components", action="store_true",
                    help="keep every surface shell, including internal cavities")
    pc.add_argument("--resample-mm", type=float,
                    help="isotropic voxel size for the surface grid, mm (0 = native). "
                         "Fast and low-memory, but ERASES structures thinner than the "
                         "target voxel; prefer --target-faces")
    pc.add_argument("--smooth-iters", type=int, help="windowed-sinc iterations")
    pc.add_argument("--passband", type=float, help="windowed-sinc passband (lower = smoother)")
    pc.add_argument("--target-faces", type=int,
                    help="decimate to this many triangles, preserving topology (0 = off)")
    pc.add_argument("--post-smooth-iters", type=int, help="smoothing after decimation")
    pc.add_argument("--no-cap", action="store_true",
                    help="do not close the surface where anatomy leaves the field of view")
    pc.add_argument("--no-validate", action="store_true")
    pc.add_argument("--self-intersections", action="store_true",
                    help="also count self-intersecting faces (needs the 'quality' extra)")
    pc.add_argument("--json", help="write results and provenance to this JSON file")
    pc.add_argument("-q", "--quiet", action="store_true")
    pc.set_defaults(func=cmd_convert)

    pm = sub.add_parser(
        "merge",
        help="fuse two scans of the same anatomy into one surface",
        description="Rigidly register two DICOM series of the same anatomy and fuse "
                    "them into a single watertight surface. Refuses pairs that do "
                    "not pass registration quality gates.",
    )
    pm.add_argument("dicom_dir_a", help="fixed scan (defines the output coordinate frame)")
    pm.add_argument("dicom_dir_b", nargs="?",
                    help="moving scan; omit to fuse two series from the first directory")
    pm.add_argument("-o", "--output", required=True)
    pm.add_argument("--series-a", help="series ident in the fixed scan")
    pm.add_argument("--series-b", help="series ident in the moving scan")
    pm.add_argument("--preset", default="bone", choices=sorted(PRESETS))
    pm.add_argument("--threshold", help="intensity (HU for CT) or 'auto'; applies to both")
    pm.add_argument("--median-mm", type=float)
    pm.add_argument("--closing-mm", type=float)
    pm.add_argument("--min-island-mm3", type=float)
    pm.add_argument("--print-profile", default="anatomical", choices=sorted(PRINT_PROFILES),
                    help="prepare the fused mesh for a printer (default: %(default)s). "
                         "Registration always runs on unmodified anatomy.")
    pm.add_argument("--thicken-mm", type=float,
                    help="grow every surface outward by this radius, thickening walls")
    pm.add_argument("--grid-mm", type=float, default=DEFAULT_GRID_MM,
                    help="isotropic voxel size of the fused grid (default %(default)s). "
                         "Finer keeps thinner bone, at cubic memory cost.")
    pm.add_argument("--smooth-iters", type=int, help="default: from the preset")
    pm.add_argument("--passband", type=float, help="default: from the preset")
    pm.add_argument("--target-faces", type=int, help="default: from the preset")
    pm.add_argument("--post-smooth-iters", type=int, help="default: from the preset")
    pm.add_argument("--force", action="store_true",
                    help="fuse even if the scans look like different patients or the "
                         "registration fails its quality gates")
    pm.add_argument("--no-validate", action="store_true")
    pm.add_argument("--self-intersections", action="store_true")
    pm.add_argument("--json", help="write results and provenance to this JSON file")
    pm.add_argument("-q", "--quiet", action="store_true")
    pm.set_defaults(func=cmd_merge)

    pv = sub.add_parser("validate", help="report mesh quality")
    pv.add_argument("mesh")
    pv.add_argument("--json", action="store_true")
    pv.add_argument("--self-intersections", action="store_true")
    pv.set_defaults(func=cmd_validate)

    pr = sub.add_parser("repair", help="make a non-watertight mesh watertight (needs 'repair' extra)")
    pr.add_argument("mesh")
    pr.add_argument("-o", "--output", required=True)
    pr.add_argument("--json", action="store_true")
    pr.add_argument("--self-intersections", action="store_true")
    pr.add_argument("-q", "--quiet", action="store_true")
    pr.set_defaults(func=cmd_repair)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
