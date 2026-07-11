# Technical reference

This document describes implementation details and limitations that are useful
for auditing results or contributing to `dicom-surface`. Start with the project
[README](../README.md) for installation and first use, or the [user
guide](user-guide.md) for presets, input limitations, and safety.

## Command-line architecture

The console script points directly to a Typer application. Typer owns command
dispatch, Rich-formatted help, required input path validation, typed preset
choices, and usage errors. The root callback prints help and
returns success when no command is supplied. Shell completion options are not
installed by the application.

The CLI module imports only Typer, Rich, defaults, and preset registries while
constructing commands. MeshLib, registration, validation, and processing
modules are imported inside the command that needs them. Root help, command
help, and malformed invocations therefore finish without loading the processing
pipeline. Expected selection, modality, file, and safety failures are concise;
unexpected programming exceptions remain visible, with traceback locals
hidden.

Human output uses shared Rich stdout and stderr consoles with terminal color
detection. Series, preset, and quality results use responsive
tables. Dynamic paths, UIDs, descriptions, and error text are treated as plain
text rather than Rich markup. Long operations announce a stage before they
begin. Interactive terminals show an indeterminate spinner, current stage, and
per-stage elapsed time. The interactive renderer runs in a small helper process
so native MeshLib calls that hold Python's interpreter lock cannot freeze its
animation. Redirected output receives persistent ANSI-free stage lines. No
percentage is shown because the processing libraries do not expose a reliable
completed-work total. Normal progress can be suppressed with `--quiet`, while
warnings and failures remain visible.

JSON is a separate plain-output contract. `list --json`, `validate --json`, and
`repair --json` write only JSON to stdout. `convert --json FILE` and
`merge --json FILE` write JSON only to the requested file. Human output never
shares the JSON destination, and JSON contains no ANSI control sequences.
Warnings raised by pydicom during discovery are captured, deduplicated with
their occurrence counts preserved, and rendered concisely on stderr without
Python source locations.

## Series discovery and selection

Discovery sorts stacks by DICOM `SeriesNumber`, complete SeriesInstanceUID, and
orientation part, then assigns unique 1-based row IDs. Those IDs are stable for
unchanged directory contents but intentionally local to one discovery result.
They are not written into processing provenance. Provenance uses the complete
UID, DICOM SeriesNumber, and orientation-part metadata instead.

`--series`, `--series-a`, and `--series-b` resolve a displayed row ID, a complete
SeriesInstanceUID, or a case-insensitive description substring. DICOM
`SeriesNumber` is not a selector because it need not be unique. When a UID or
description matches several orientation stacks, selection fails with the row
IDs that disambiguate it.

## Processing pipeline

`convert` performs these stages:

1. Discover DICOM instances, group them by series UID and orientation, and
   assign deterministic row IDs after sorting.
2. Order slices by `ImagePositionPatient` projected onto the slice normal.
3. Load the ordered stack with SimpleITK and resolve the intensity threshold.
4. Apply island filtering, median filtering, optional opening, and closing.
5. Pad field-of-view boundaries, extract the 0.5 isosurface with MeshLib marching
   cubes, smooth, select surface components, and optionally simplify.
6. Transform vertices into DICOM patient LPS coordinates; validate the mesh in
   memory, write and validate a temporary file, then atomically publish it.

The output is normally a closed surface because the mask is padded with
background before extraction. When anatomy touches the scan boundary, the cap is
flat and a warning explains that missing anatomy was not recovered.

## DICOM compatibility

pydicom reads headers for discovery. SimpleITK's GDCM-backed reader loads pixel
data. The loader is for classic image series, not arbitrary DICOM objects.

| case                                        | behavior                                                                                               |
| ------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| Classic single-frame stacks                 | Supported and covered by generated uncompressed CT tests                                               |
| Enhanced multi-frame objects                | Unsupported; one file is counted as one instance and per-frame geometry is not parsed                  |
| Compressed transfer syntaxes                | Delegated to codecs in the installed SimpleITK/GDCM build; not tested here                             |
| `RescaleSlope` / `RescaleIntercept`         | Applied by SimpleITK; CT presets rely on calibrated HU                                                 |
| Pixel padding and `MONOCHROME1`             | No explicit project handling; untested                                                                 |
| Consistent oblique stacks                   | Supported; slice normal and direction cosines are preserved                                            |
| Gantry tilt or nonparallel slices           | No correction or complete geometry validation; unsupported without independent checks                  |
| Irregular spacing                           | Warned and regularized when moderate; rejected when spread exceeds `max(0.1 mm, 0.5 × median spacing)` |
| Duplicate positions                         | Not rejected; unsupported                                                                              |
| Missing geometry tags or varying dimensions | May fail or load incorrectly; unsupported                                                              |
| Localizers                                  | Rejected when `ImageType` contains `LOCALIZER`; short stacks are also rejected                         |
| Reformats                                   | A consistent classic stack can load; axial data is preferred after voxel-volume ranking                |
| Vendor mosaics                              | Unsupported                                                                                            |

`FrameOfReferenceUID` is not read. `merge` always registers the moving scan and
does not assume cross-series coordinates already align.

When `ConvolutionKernel` is present, a case-insensitive warning heuristic uses
the literal pattern
`(?:^|[^0-9])(?:[BHUY]r?|BONE|EDGE|LUNG)\s*_?([6-9]\d)`. It recognizes names
such as `Hr68` and `B70f`; vendor naming outside that pattern can be missed or
misclassified. A match only produces a sharp-kernel warning. It does not change
series ranking, thresholds, or processing.

## Morphology and physical units

Median, opening, and closing parameters are kernel extents in millimetres. Each
axis is floored to the largest integer radius whose realized `(2r+1) × spacing`
does not exceed the request. A zero radius is an identity operation on that axis.

Fractional occupancy is linearly interpolated and thresholded at 0.5. Resampling
can erase thin structures or change components, cavities, tunnels, and genus. It
does not preserve anatomical topology.

## Surface extraction and finishing

The isosurface implementation is MeshLib `marchingCubes` at 0.5. SimpleITK arrays
are transposed from z/y/x to x/y/z before extraction, then shifted half a voxel
to preserve the established sample-coordinate convention. Mesh validity is
measured rather than assumed.

Smoothing uses MeshLib `relaxKeepVolume`. The preset controls initial iterations,
relaxation force, and post-simplification iterations. The `--smooth-force` option
replaces the former backend-specific passband. Smoothing moves surfaces, and the
project does not provide a general deviation bound. Every requested
iteration runs before self-intersection detection. If smoothing makes
non-adjacent faces collide, vertices in those collision patches return to their
pre-smooth positions. The protected set grows by topological rings until the
mesh is collision-free. Smoothing is therefore retained globally instead of
reducing the iteration count for the whole surface.

Simplification uses MeshLib's quadric edge-collapse implementation with
`DecimateStrategy.MinimizeError`. `--simplify-error-mm` sets the estimated
surface-deviation/QEM limit in model millimetres; this is not a certified
Hausdorff bound, and `0` disables simplification. A candidate is accepted only
when it preserves the component/hole/Euler signature, does not increase boundary
or non-manifold edges, and has no self-intersections. If a candidate
contains collisions, the colliding triangles are projected back onto the
pre-decimation source mesh. Four-ring source neighborhoods around those
locations are excluded from collapse and decimation restarts. Up to eight local
protection passes are allowed. MeshLib supplies every collision check and marks
both faces from each colliding pair.

This keeps simplification active outside small unsafe patches. If topology or
manifold checks fail, a collision patch cannot be mapped, or all protection
passes are exhausted, the valid pre-decimation mesh is retained and the command
reports a warning. The resulting face count and MeshLib's introduced-error
estimate are reported. Conversion and merging use this same finishing path.

## Merge registration and gates

Registration uses thresholded and morphologically cleaned anatomical masks:

1. Resample each mask onto a separate axis-aligned 2.0 mm world lattice.
2. Find the best integer 3-D translation with FFT cross-correlation.
3. Extract smoothed surface vertices and normals.
4. Run two point-to-plane ICP passes. The first keeps the closest 45%; the
   second widens to `clamp(0.9 × matched_fraction, 0.45, 0.95)`.

Moving vertices are sampled to at most 60,000 with seed 0. ICP uses a 4 mm
correspondence limit and up to 80 iterations per pass.

| gate                      | definition                                                                                                           |  threshold |
| ------------------------- | -------------------------------------------------------------------------------------------------------------------- | ---------: |
| Symmetric surface overlap | Weaker of moving→fixed and fixed→moving fractions within the other scan's field of view and 4 mm of a surface vertex |       0.60 |
| Shared-field Dice         | Binary-mask Dice on a 1.5 mm world grid, restricted to common acquired coverage                                      |       0.55 |
| Shared field of view      | Common coverage volume on that 1.5 mm grid                                                                           | 20,000 mm³ |

Point-to-plane RMS and median are reported but not gated. The thresholds were
selected on a small exploratory development set, not clinically calibrated.

## Patient comparison

`merge` reads `PatientID`, `IssuerOfPatientID`, `PatientName`,
`PatientBirthDate`, and `PatientSex` from the first instance of each series.

- A present name, birth-date, or sex conflict refuses the merge.
- Matching IDs corroborate identity unless both issuers are present and differ.
- A matching normalized name corroborates identity when birth dates match or
  are not both present.
- Differing IDs without earlier corroboration refuse the merge.
- Missing or unverifiable identity warns and proceeds.

`--force` overrides modality mismatch, a different-patient verdict, and failed
registration gates. It does not override selecting the same series twice,
loading errors, or an excessive fused grid.

Patient-field values are not logged or stored in provenance. Series UIDs,
descriptions, paths, and derived anatomy can still be identifying.

## Validation and repair

Validation reports MeshLib-imported triangle/vertex counts, components,
watertightness, winding, holes, boundary edges, disoriented faces, genus when
defined, volume when closed, bounding box, and self-intersecting faces. The same
complete check runs after `convert`, `merge`, and `repair`, and for `validate`.

A report is valid only when MeshLib imports the mesh as watertight, consistently
wound, and enclosing a volume, with zero holes, boundary edges, disoriented
faces, and self-intersecting faces. Multiple closed components are allowed.
MeshLib can normalize raw face configurations while loading, so validation does
not expose separate non-manifold or degenerate-face counters. A failed
self-intersection measurement is invalid. `convert` and `merge` have no
validation bypass: they validate both the
in-memory mesh and the serialized temporary file, then atomically replace the
requested destination only with a valid output.

STL, PLY, and OBJ are loaded directly by MeshLib. VTP is not supported. All
metrics describe the MeshLib-imported representation, and self-intersections
count the unique faces in MeshLib collision pairs.

`repair` uses MeshLib to unite vertices within 1e-6, fix multiple edges,
decimate degeneracies, and fill holes. It writes a new file; validation never
mutates the input it measures. A repaired file remains on disk if it does not
pass the complete validation contract.

## Dependencies and development

The repository uses uv 0.11.28 for dependency resolution, environments, command
execution, and CI. Typer 0.21 defines the CLI and Rich 14 renders human terminal
output. `uv.lock` covers Python 3.10–3.13; local development defaults to Python
3.12.

```sh
uv sync --locked
uv run pytest -q
```

Hatchling is the PEP 517 build backend, and uv orchestrates the workflow.
MeshLib is pinned to `3.1.3.297` for consistent collision and simplification
behavior; no OpenGL system dependency is required.

## Verification

Synthetic regression tests cover series grouping, geometry, segmentation,
surface extraction, registration, validation, and repair. Run them with
`uv run pytest -q`.
