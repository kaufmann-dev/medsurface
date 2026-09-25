# Technical reference

This document describes implementation details and limitations useful for
auditing results or contributing to `medsurface`. Start with the project
[README](../README.md) for installation and first use, or the [user
guide](user-guide.md) for operational guidance and safety.

## Command-line architecture

The console script points directly to a Typer application. Typer owns command
dispatch, Rich-formatted help, required input-path validation, typed preset
choices, and usage errors. The root callback prints help successfully when no
command is supplied; `--version` prints the installed package version without
loading processing modules. Shell completion options are not installed by the
application.

The public operations are deliberately separate:

- `convert` preserves one image while changing only its storage format.
- `extract` segments one intensity volume and produces a surface mesh.
- `fuse` segments, registers, and unions two intensity volumes into a binary
  labelmap.
- `labelmap extract` produces a surface from an existing labelmap.
- `labelmap fuse` registers and unions two existing labelmaps into a binary
  labelmap.

Fusion never dispatches to mesh output. A fused labelmap must be passed
explicitly to `labelmap extract`, which keeps volume publication independently
inspectable and gives extraction its normal labelmap defaults.

The CLI module imports Typer, Rich, defaults, output classification, and preset
registries while constructing commands. MeshLib, registration, validation, and
processing modules are imported inside the command that needs them. Root help,
command help, and malformed invocations therefore finish without loading heavy
processing modules. Expected selection, volume, file, and registration failures
are concise; unexpected programming exceptions remain visible with traceback
locals hidden.

Human output uses shared Rich stdout and stderr consoles. Dynamic paths, UIDs,
descriptions, and error text are plain text rather than Rich markup. Suggested
commands are formatted from argument vectors: POSIX uses shell quoting and
Windows uses `list2cmdline`. Volume discovery and successful fusion therefore
emit copyable, shell-safe extraction hints even when paths contain spaces or
metacharacters, unless `--quiet` suppresses normal output.

Long operations announce a stage before they begin. Interactive terminals show
an indeterminate spinner and elapsed time in a helper process, so native MeshLib
work that holds Python's interpreter lock does not freeze animation. Live
warnings and errors clear and redraw the display. The helper ignores `SIGINT`;
the main process owns cancellation, reaps the helper, cleans temporary resources,
prints one `Cancelled.` diagnostic, and exits with status 130. Redirected output
receives persistent ANSI-free stage lines. `--quiet` suppresses normal progress
but not warnings or failures.

## JSON and output transactions

`list --json`, `validate --json`, and `repair --json` write plain JSON to stdout.
Conversion, extraction, and fusion accept `--json FILE`; the human stream never
shares that destination and JSON contains no ANSI control sequences.

Report families are operation-specific:

- Conversion records output format/compression, dimensions, pixel type,
  components, geometry, duration, warnings, metadata policy, and input
  provenance.
- Extraction records surface counts, bounds, mask/surface components, finishing
  provenance, validation quality, duration, and warnings.
- Fusion records output format/compression, scalar `uint8` type, grid geometry,
  fixed/moving/fused foreground counts and volumes, complete registration data,
  segmentation settings, duration, warnings, and both inputs' provenance. It
  has no triangle, surface-finishing, or mesh-quality fields.

Reports are written, flushed, and atomically replaced from a same-directory
temporary file. For conversion and fusion, the CLI also stages a same-filesystem
rollback sibling for any existing primary destination, using a hard link when
available and a copy otherwise. Any volume-output or report failure restores the
old primary file, or removes a newly created one; success removes the rollback
sibling.

Before pixel loading, every processing command rejects primary/report aliases
to each other and to every discovered input through lexical identity, symlinks,
or hard links. Protected paths include all selected and unselected discovered
volume files, every DICOM instance, detached headers, and their referenced
payloads. Unsupported output suffixes are also rejected before discovery.

Warnings raised by pydicom during discovery are captured and deduplicated with
occurrence counts while omitting Python source locations.

## Volume discovery and selection

A direct supported file contributes one catalog candidate. A directory is
walked recursively; classic DICOM instances are grouped into physical stacks and
supported NIfTI, NRRD, and MetaImage headers each contribute a candidate. The
combined list is sorted deterministically and assigned unique 1-based IDs. IDs
are stable for unchanged contents but local to one discovery result.

The common record carries format, source, modality, description, size, spacing,
direction, origin, pixel type, component count, plane, and usability. DICOM adds
UID, orientation part, `SeriesNumber`, convolution-kernel values, HU evidence,
and spacing diagnostics. Multi-valued kernels are arrays in JSON/provenance and
comma-separated in human output. Direction matrices provide plane information
for file inputs.

Only displayed integer IDs are accepted by `--volume`, `--fixed-volume`, and
`--moving-volume`. UIDs, `SeriesNumber`, descriptions, and paths are not
selectors.

Selection without an ID follows four rules: one usable candidate is automatic;
multiple usable DICOM candidates with one shared modality retain DICOM ranking;
DICOM candidates spanning modalities require an explicit ID; every other
multi-volume catalog also requires an ID. Missing DICOM modality is its own
`unknown` value. `fuse` builds or reuses a catalog for each positional input and
applies the rules independently.

The labelmap group is narrower. `labelmap extract` accepts one direct NIfTI,
NRRD, or MetaImage file. `labelmap fuse` accepts one such file for each role.
Directories, DICOM, catalog IDs, presets, and label selectors are rejected
because label selection belongs to the program that produced the segmentation.

## Strict conversion and verified volume publication

`volume.convert` performs these stages:

1. Classify the requested atomic output before loading pixels.
2. Load exactly one already-selected candidate under the source voxel limit,
   requesting source metadata unless `--strip-metadata` is active.
3. Capture output/result geometry and format-neutral source provenance.
4. Publish through `write_verified_volume`.

No thresholding, casting, normalization, resampling, reorientation, or
segmentation occurs. The preserved core contract is loaded voxel bytes, scalar
pixel ID/type, component count, image dimension, voxel dimensions, spacing,
origin, and direction.

Metadata is operation-specific. Strict conversion preserves it by default;
extraction and fusion load only pixels and geometry because their derived
outputs do not carry source metadata. For DICOM preservation,
`ImageSeriesReader` orders the already-discovered series files, updates the
per-slice metadata dictionary, loads private tags, and promotes valid metadata
from the representative first slice onto the resulting 3-D image. Individual
values that cannot be encoded as UTF-8 are omitted with a warning rather than
lossily repaired. File readers retain metadata exposed by their ImageIO only
when preservation is requested.

The verified writer compares preserved source metadata values to readback and
warns when the destination cannot represent every entry. Metadata loss does not
relax or fail the core image contract. With `--strip-metadata`, source metadata
is never requested from DICOM and is removed from file inputs; required
destination-format headers are still synthesized by ImageIO.

One centralized output classifier defines the storage contract:

| suffix    | format    | compression | atomic storage |
| --------- | --------- | ----------- | -------------- |
| `.nii`    | NIfTI     | none        | single file    |
| `.nii.gz` | NIfTI     | gzip        | single file    |
| `.nrrd`   | NRRD      | gzip        | single file    |
| `.mha`    | MetaImage | zlib        | single file    |

Suffix matching is case-insensitive. The temporary file always uses the
normalized lowercase compound suffix because ITK's writer dispatch is
case-sensitive for NIfTI and can interpret uppercase `.MHA` as detached
MetaImage. Explicit ImageIO selection allows uppercase results to be read back
as later inputs. `.nhdr` and `.mhd` remain input-only.

`write_verified_volume` implements the publication contract:

1. Capture a SHA-256 digest by iterating array planes without retaining another
   complete volume, plus primitive type/dimension/geometry facts and metadata.
2. Create a same-directory temporary sibling with the normalized final suffix.
3. Write with suffix-selected compression at level 9 for compressed formats.
4. Invoke the ownership-release callback and drop the writer's source reference
   before reading the temporary file, avoiding two complete output buffers.
5. Compare dimension, voxel dimensions, pixel ID, components, and digest exactly.
   Compare NRRD and MetaImage geometry with zero relative tolerance and absolute
   tolerance `1e-5`. NIfTI geometry additionally accepts at most one float32 ULP
   per field because NIfTI-1 stores its affine and spacing at that precision.
6. Inspect serialized headers to verify the requested compression contract.
7. Release readback and atomically replace the destination.

All failure and cancellation paths remove the temporary file. An existing
destination is untouched until the final replacement. A format that coerces a
core property—for example, a NIfTI writer orthogonalizing a sheared direction—is
rejected after readback rather than silently accepted. Sub-micrometre NIfTI
rounding within one representable float32 step is not treated as a geometry
change.

## Extraction pipelines

### Intensity extraction

`pipeline.extract` resolves the preset/override threshold, measures the loaded
range, and builds a binary mask. `pipeline.build_mask` applies thresholding,
island filtering, median filtering, optional opening, closing, and a final
island pass. The mask is optionally resampled and smoothed, padded for
field-of-view capping, converted to a 0.5 isosurface, transformed into SimpleITK
physical coordinates, finished, validated, and atomically written as STL, PLY,
or OBJ.

When foreground touches the source boundary, default padding creates a flat cap
and emits a warning. `--no-cap` leaves the surface open there. Missing anatomy
is never reconstructed.

### Labelmap extraction

`labelmap.extract` starts at the binary-mask boundary:

1. Discover and load one direct self-describing image.
2. Require finite, non-negative, integer-valued voxels and nonempty foreground.
3. Cast every nonzero value to one shared `uint8` foreground mask.
4. Apply the labelmap surface settings, extract, finish, validate, and publish.

It skips intensity thresholds, Otsu, morphology, mask-island removal, and label
interpretation. Integer-valued floating-point labelmaps are accepted;
fractional probability maps are rejected. Binary and multilabel inputs follow
the same path.

Labelmap surface defaults are native grid, `0.8 mm` occupancy smoothing, 20
pre-simplification relaxation iterations, `0.25 mm` simplification error, zero
post-simplification iterations, and all surface components retained. These
settings are resolved only when extraction begins; fusion does not embed them.

Every input axis must contain at least four samples because the shared recursive
Gaussian and processing contract have no alternate small-image path.

## File-volume compatibility

SimpleITK reads headers during discovery and pixels only after selection.
Discovery resolves detached MetaImage/NRRD payload references and marks a header
unusable if any payload is missing or unreadable. Candidates must be scalar,
real-valued 3-D images with at least four voxels per axis, finite origin and
direction, finite positive spacing, and a nonsingular direction matrix.

The selected header dimensions are multiplied with Python integers before pixel
loading. Missing dimensions are unusable. Source volumes above 500 million
voxels are refused by default; the same limit applies to planned extraction and
fusion grids. `--allow-large-volume` bypasses only these count guards.

| format    | accepted inputs   | atomic outputs    | behavior                                            |
| --------- | ----------------- | ----------------- | --------------------------------------------------- |
| NIfTI     | `.nii`, `.nii.gz` | `.nii`, `.nii.gz` | Single-file input and output                        |
| NRRD      | `.nrrd`, `.nhdr`  | `.nrrd`           | Detached input requires every referenced payload    |
| MetaImage | `.mha`, `.mhd`    | `.mha`            | Detached input requires every referenced payload    |
| HDF5      | `.h5`, `.hdf5`    | none              | Not accepted                                        |
| Raw       | `.raw`, `.bin`    | none              | Not accepted as an independent self-described image |

NRRD and MetaImage can retain arbitrary metadata keys their readers expose.
NIfTI metadata is limited to fields supported by its ImageIO. Plane is always
derived from geometry.

## DICOM compatibility

pydicom reads headers for discovery. SimpleITK's GDCM-backed reader loads pixel
data. The loader targets classic image series, not arbitrary DICOM objects.

| case                                        | behavior                                                                                               |
| ------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| Classic single-frame stacks                 | Supported and covered by generated uncompressed CT tests                                               |
| Enhanced multi-frame objects                | Unsupported; one file is counted as one instance and per-frame geometry is not parsed                  |
| Compressed transfer syntaxes                | Delegated to codecs in the installed SimpleITK/GDCM build                                              |
| `RescaleSlope` / `RescaleIntercept`         | Applied by SimpleITK; both must be present and consistent before HU can be verified                    |
| Pixel padding and `MONOCHROME1`             | No explicit project handling; untested                                                                 |
| Consistent oblique stacks                   | Supported; slice normal and direction cosines are preserved                                            |
| Gantry tilt or nonparallel slices           | No correction or complete geometry validation; unsupported without independent checks                  |
| Irregular spacing                           | Warned and regularized when moderate; rejected when spread exceeds `max(0.1 mm, 0.5 × median spacing)` |
| Duplicate positions                         | Not rejected; unsupported                                                                              |
| Missing geometry tags or varying dimensions | Missing dimensions are rejected; varying dimensions need independent checks                            |
| Localizers                                  | Rejected when `ImageType` contains `LOCALIZER`; short stacks are also rejected                         |
| Reformats                                   | A consistent classic stack can load; axial data is preferred after voxel-volume ranking                |
| Vendor mosaics                              | Unsupported                                                                                            |

`FrameOfReferenceUID` is not used. Fusion registers moving foreground rather
than assuming that cross-series coordinates already agree.

HU calibration is conservative. Every slice must carry the same finite rescale
slope/intercept. An explicit `RescaleType` must be `HU`. Without it, the series
is accepted as HU only when `ImageType` begins with `ORIGINAL`, it is not a
localizer, and `MultienergyCTAcquisition` is not `YES`. Derived or multienergy
CT can still verify itself with `RescaleType=HU`; otherwise HU presets warn.

Each `ConvolutionKernel` value is independently checked by a case-insensitive
heuristic for sharp/edge-enhancing names. A match warns but does not change
selection, thresholds, or processing.

## Morphology and physical units

Median, opening, and closing parameters are kernel extents in millimetres. Each
axis is floored to the largest integer radius whose realized `(2r+1) × spacing`
does not exceed the request. A zero radius is an identity on that axis.
Intensity extraction and intensity fusion expose the same segmentation
overrides. Fusion applies them independently to fixed and moving masks.

Labelmap commands do not run morphology or mask-island filtering. Their input is
treated as the completed segmentation, and every disconnected nonzero region is
preserved until grid resampling or extraction controls change it.

Fractional occupancy is linearly interpolated and classified at `> 0.5` for
fusion. Resampling can erase thin structures or alter components, cavities,
tunnels, and genus. Planned non-finite/non-positive spacing and grids above the
voxel ceiling are rejected before allocation.

## Surface extraction and finishing

MeshLib `marchingCubes` extracts the 0.5 isosurface. SimpleITK arrays are
transposed from z/y/x to x/y/z and shifted half a voxel to preserve sample
coordinates. The index-to-physical affine transforms vertices before finishing.
When its linear component has a negative determinant, triangle order is reversed
so a left-handed image direction does not produce inward winding.

`--mask-smooth-mm` applies one recursive Gaussian to the segmented occupancy
mask before marching cubes. Intensity presets default to zero; external
labelmaps default to `0.8 mm`. Smoothing can move boundaries, join narrow gaps,
or remove small structures and always produces a warning when enabled.

`--mesh-smooth-iters` controls fixed-force MeshLib `relaxKeepVolume` before
simplification. `--post-mesh-smooth-iters` controls the same guarded relaxation
afterward. If a pass creates collisions or inconsistent winding, vertices in
unsafe neighborhoods return to their input positions. The protected area grows
by topological rings until the mesh is clean; the whole pass is discarded if
local protection cannot succeed. Zero disables only that pass.

`--simplify-error-mm` supplies MeshLib's estimated QEM/deviation limit in
physical model millimetres. It is not a certified Hausdorff bound. Candidates
must preserve component/hole/Euler signatures, avoid increasing boundary or
non-manifold edges, and contain no self-intersections or disoriented faces.
Unsafe source neighborhoods are protected and decimation retries up to eight
times. If no candidate passes, the valid pre-decimation surface is retained.

`--destep auto|all|band` adds optional final masked Taubin fairing after
post-simplification relaxation. It runs alternating `λ 0.5` and `μ −0.53`
uniform-Laplacian steps for `--destep-iters` iterations, each scaled by a
per-vertex weight, so weight-zero vertices stay bit-identical. `band` weights
follow a cosine ramp between `--destep-frozen-mm` and `--destep-full-mm` on
one axis and must reach the surface. `all` uses weight one. `auto` fairs an
unmasked copy, averages its vertex normals over five neighbour passes, and
takes each vertex's largest normal change per millimetre along its edges as
its curvature. Radii below `4 mm` are detail. Connected detail smaller than
`20 mm²` is an isolated speck and is faired, because freezing it would pin the
surrounding ripples. Other detail is frozen and returns to full fairing over
`1–4 mm` of geodesic edge-path distance. Detail clusters of at least `100 mm²`
also freeze everything within `10 mm` and fully fair beyond `20 mm`, so smooth
patches between facial features stay unchanged. The constants were tuned on a
fused CT skull labelmap, where `auto` matched a manually banded result on the
outer vault while leaving the orbital rims, face, and dentition untouched.
Displacements are clamped to `--destep-max-mm`.
Faces whose normal rotates past a dot product of `0.2` return, with one ring,
to their input positions. The relaxation safeguard then protects any remaining
self-intersections or disoriented faces, or keeps the unfaired surface.
Provenance `surface_finishing.destep` records region, iterations, faired and
frozen fractions, clamped and unfolded vertices, displacement, volume change,
and the safeguard result; it is `null` when fairing is off. The method follows
a masked-fairing experiment on a CT skull dome, where 600 iterations cleared
8–16 mm ripples with at most 1 mm displacement.

`--components all|largest` applies only to extraction. Omitted intensity
extraction follows the preset: bone, skin, and auto retain the largest shell;
teeth retains all. Labelmap extraction defaults to all. Fusion has no mask
smoothing, mesh smoothing, simplification, capping, or component-selection
parameters and does not invoke a surface-output path.

## Fusion registration and gates

`fusion.fuse` supplies independently thresholded and cleaned intensity masks.
`labelmap.fuse` supplies independently validated all-nonzero masks.
`fusion.fuse_masks` owns their common registration, grid, quantization, and
publication path.

Registration proceeds as follows:

1. Resample each mask onto its own axis-aligned 2.0 mm world lattice.
2. Find the best integer 3-D translation by FFT cross-correlation.
3. Extract temporary smoothed surface vertices/normals for registration metrics.
4. Run two point-to-plane ICP passes. The first keeps the closest 45%; the second
   keeps `clamp(0.9 × matched_fraction, 0.45, 0.95)`.

Moving vertices are sampled to at most 60,000 with seed 0. ICP uses a 4 mm
correspondence limit and up to 80 iterations per pass. Temporary surface work is
internal to registration only; fused-output generation never extracts a mesh.

| gate                      | definition                                                                                                           |  threshold |
| ------------------------- | -------------------------------------------------------------------------------------------------------------------- | ---------: |
| Symmetric surface overlap | Weaker of moving→fixed and fixed→moving fractions within the other scan's field of view and 4 mm of a surface vertex |       0.60 |
| Shared-field Dice         | Binary-mask Dice on a 1.5 mm world grid, restricted to common acquired coverage                                      |       0.55 |
| Shared field of view      | Common coverage volume on that 1.5 mm grid                                                                           | 20,000 mm³ |

Point-to-plane RMS and median are reported but not gated. Thresholds came from a
small exploratory development set, not clinical calibration. `--force` bypasses
only gate failure.

After registration, grid planning transforms moving bounds into the fixed
physical coordinate system, combines them with fixed bounds, adds a two-voxel
margin, and creates an identity-direction isotropic lattice at `--grid-mm`
spacing. Both masks are antialiased for that grid and linearly resampled. The
voxelwise maximum combines occupancies; values strictly greater than `0.5` are
cast to `uint8` foreground `1`, with background `0`.

Fixed, moving, and fused foreground counts are measured under the same
quantization. Volumes equal count × `grid_mm³`. Empty fusion is an error. Derived
metadata is cleared before verified volume publication. Positive source label
IDs and their collisions are intentionally collapsed to one class.

The grid must have finite positive spacing and no coordinate/integer overflow.
Plans above 500 million voxels are rejected unless `--allow-large-volume` is
given. A grid coarser than the finest input voxel emits a feature-loss warning.

## Fusion identity contract

Neither fusion command can establish subject identity. The intensity path does
not compare modalities. Every operation warns the user to confirm both inputs.
The labelmap path additionally warns that masks must represent matching rigid
structures. Different positive IDs are equivalent.

The same file, hard-linked file, or DICOM UID/orientation part cannot fill both
roles. `--force` does not override duplicate selection, catalog ambiguity,
loading failures, or voxel limits. `--allow-large-volume` overrides only source
and processing-grid count ceilings.

DICOM UIDs, descriptions, paths, JSON provenance, and derived anatomy can remain
identifying even though catalog discovery does not read patient identity fields.

## Validation and repair

Mesh validation reports imported triangle/vertex counts, components,
watertightness, winding, holes, boundary edges, disoriented faces, genus when
defined, closed volume, bounding box, and self-intersecting faces. It runs after
both extraction commands, repair, and explicit validation.

A mesh is valid only when MeshLib imports it as watertight, consistently wound,
and enclosing a volume, with zero holes, boundary edges, disoriented faces, and
self-intersecting faces. Multiple closed components are allowed. Metrics describe
the MeshLib-imported representation. STL, PLY, and OBJ are supported; VTP is not.

Extraction and repair validate the in-memory mesh and the serialized temporary
file, then replace the destination only with a valid result. There is no
validation bypass. `validate` never mutates its input.

Volume conversion and fusion do not use mesh validation. Their voxel/geometry,
compression, and atomic replacement checks are the verified-volume contract
described above.

`repair` unites vertices within `1e-6`, fixes multiple edges, decimates
degeneracies, and fills holes. An invalid result is discarded and an existing
destination is preserved.

## Dependencies and development

The repository uses uv 0.11.28 for dependency resolution, environments, command
execution, and CI. Typer 0.21 defines the CLI and Rich 14 renders terminal
output. `uv.lock` covers Python 3.11–3.14; local development defaults to Python
3.12.

```sh
uv sync --locked
uv run ruff check .
uv run mypy
uv run pytest -q
uv build
```

Hatchling is the PEP 517 backend. MeshLib is pinned to `3.1.3.297` for consistent
collision and simplification behavior; no OpenGL system dependency is required.
CI runs linting, source type checks, package builds, an installed-command smoke
test, and the full suite on Python 3.11 through 3.14.

## Verification

Synthetic regression tests cover catalog selection, every supported conversion
input/output family, metadata policy, verified atomic publication, binary
fusion, registration gates, complete fusion-to-extraction workflows, surface
extraction, validation, and repair. Run them with `uv run pytest -q`.
