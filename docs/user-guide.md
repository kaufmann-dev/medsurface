# User guide

[Project README](../README.md) · [Choosing a volume](#choosing-a-volume) ·
[Strict volume conversion](#strict-volume-conversion) ·
[Extracting an intensity surface](#extracting-an-intensity-surface) ·
[Using an external labelmap](#using-an-external-labelmap) ·
[Fusion](#fusion) ·
[Input limitations](#input-requirements-and-limitations) ·
[Safety and privacy](#safety-and-privacy) ·
[Technical reference](technical-reference.md)

`medsurface` separates storage conversion, segmentation-based fusion, and
surface extraction. This distinction keeps derived labelmaps inspectable and
editable before a mesh is created. See the [technical
reference](technical-reference.md) for algorithms, metrics, and dependency
details.

## Safety and privacy

`medsurface` is not validated for diagnosis, treatment planning, or other
clinical decisions. Independently review source images, binary labelmaps, and
generated surfaces.

- Segmentation thresholds, morphology, resampling, registration, smoothing,
  field-of-view capping, and scan resolution can change or omit anatomy.
- A structurally valid mesh and a storage-verified volume can still be
  anatomically wrong.
- Anatomy cut off by the scan is capped flat during surface extraction by
  default; missing anatomy cannot be recovered.
- Fusion does not inspect patient identity fields or establish subject
  identity. Confirm that both inputs show the same subject.
- Fusion uses rigid registration and is intended for matching non-deforming
  anatomy such as bone. Motion or deformation can produce a plausible but
  incorrect union.
- Before `labelmap fuse`, confirm that both masks represent the same selected
  structures. Label values and names are not compared.
- `--force` bypasses failed registration-quality gates only. It does not bypass
  duplicate-input, selection, loading, malformed-image, or voxel-limit errors.
- `--allow-large-volume` bypasses the default voxel-count guard. It does not make
  an operation memory-safe; allocation failure or operating-system termination
  remains possible.
- Conversion preserves source metadata best effort by default. Use
  `--strip-metadata` when source metadata must be removed, and independently
  inspect the result.
- Fused labelmaps do not copy source image metadata, but paths, UIDs,
  descriptions, registration data, and derived anatomy in logs or optional JSON
  provenance can remain identifying.
- Primary outputs and JSON reports cannot alias each other or any discovered
  medical input through a lexical path, symlink, or hard link. Protected inputs
  include every DICOM instance, detached header, and detached payload.

Treat source data, outputs, logs, and provenance according to the same privacy
rules as other personal health data.

## Choosing a volume

Run `medsurface list INPUT` immediately before conversion, intensity extraction,
or fusion when an input may contain more than one volume. `INPUT` may be a
supported file or a directory tree. Directory discovery combines classic DICOM
series with `.nii`, `.nii.gz`, `.nrrd`, `.nhdr`, `.mha`, and `.mhd` files in one
ID space.

The table shows modality, description, slice count, voxel spacing, anatomical
plane, and usability when available. DICOM also shows `SeriesNumber` as
`DICOM #`. Plane comes from stored geometry rather than a text label.

Every row has one positive integer `ID`, which is the only selector:

```sh
medsurface convert scans/ --volume 2 -o selected.nrrd
medsurface extract scans/ --volume 2 -o selected.stl
medsurface fuse scans-a/ scans-b/ \
  --fixed-volume 1 --moving-volume 3 -o fused.nii.gz
```

Selection without an ID is safe only in these cases:

- A catalog with exactly one usable candidate selects it.
- A DICOM-only catalog whose usable entries share one modality uses the
  established ranking: usable stacks, smaller voxels, axial plane, then slice
  count.
- A DICOM-only catalog spanning modalities requires an ID.
- Any other catalog with multiple usable volumes requires an ID.

`list` prints a shell-safe `medsurface extract` hint for the automatic default,
or a selector template when no default exists. IDs are deterministic for
unchanged contents but local to one discovery result; run `list` again after
files change. `list --json` reports the integer `id` and nests DICOM-specific
fields under `dicom`.

To fuse two volumes found under one directory, repeat the path and choose two
different IDs:

```sh
medsurface fuse scans/ scans/ \
  --fixed-volume 1 --moving-volume 2 -o fused.mha
```

UIDs, `SeriesNumber`, descriptions, and paths are metadata, not selectors.

## Strict volume conversion

`convert` changes storage format for exactly one selected image:

```sh
medsurface convert INPUT -o OUTPUT \
  [--volume ID] [--strip-metadata] \
  [--allow-large-volume] [--json REPORT.json] [-q]
```

It preserves the image loaded by SimpleITK:

- voxel values, without thresholding or casting;
- scalar pixel type and component count;
- image dimension and voxel dimensions;
- spacing, origin, and direction.

It does not resample, normalize, reorient, segment, or extract a surface.
Conversion output is always one atomic file:

| output format | suffixes          | compression |
| ------------- | ----------------- | ----------- |
| NIfTI         | `.nii`, `.nii.gz` | none, gzip  |
| NRRD          | `.nrrd`           | gzip        |
| MetaImage     | `.mha`            | zlib        |

Detached `.nhdr` and `.mhd` files are accepted as inputs but never outputs,
because a header and payload cannot be committed with one atomic replacement.
Suffix matching is case-insensitive.

By default, source metadata is written when the destination format can represent
it. Metadata models differ: an entry supported by NRRD or MetaImage may have no
NIfTI equivalent. Such loss produces a warning but does not weaken the image
contract. `--strip-metadata` removes source metadata before writing while the
writer still creates required format headers. DICOM conversion copies available
metadata from the representative slice when preservation is requested.

Publication is strict:

1. The image is written to a temporary sibling using suffix-selected lossless
   compression.
2. The source pixel buffer is released before the temporary file is loaded.
3. Readback must match the source voxel digest, pixel type, components,
   dimensions, and physical geometry within `1e-5`.
4. The destination is atomically replaced only after verification.

If a target format cannot represent the core image contract, conversion fails,
removes the temporary file, and leaves an existing destination unchanged. A
same-format conversion is supported when input and output are different files;
in-place conversion and aliases to any discovered input are rejected.

Example conversions:

```sh
medsurface convert dicom-study/ -o study.nii.gz
medsurface convert segmentation.nhdr -o segmentation.mha
medsurface convert scan.nrrd -o anonymized.nii --strip-metadata
```

With `--json FILE`, the report records the output path, format, compression,
dimensions, pixel type, components, spacing, origin, direction, duration,
warnings, metadata policy, and selected-volume provenance.

## Extracting an intensity surface

`extract` segments one intensity volume and creates an STL, PLY, or OBJ mesh:

```sh
medsurface extract scans/ --preset bone -o skull.stl
medsurface extract scans/ --preset auto -o uncalibrated.stl
medsurface extract scans/ --threshold 250 -o bone-250hu.stl
```

### Choosing a preset

A preset supplies segmentation and surface-finishing defaults. Run
`medsurface presets` to inspect the installed values.

| preset  | use it for                                              | threshold | median | closing | island floor | mask smoothing | pre-mesh smoothing | simplify error | post-mesh smoothing |
| ------- | ------------------------------------------------------- | --------: | -----: | ------: | -----------: | -------------: | -----------------: | -------------: | ------------------: |
| `bone`  | General CT bone models                                  |    300 HU | 1.0 mm |  2.4 mm |       50 mm³ |            off |      20 iterations |        0.25 mm |       40 iterations |
| `teeth` | Enamel and dense dentin; keeps separate teeth           |  1,200 HU | 0.6 mm |  0.6 mm |        5 mm³ |            off |      10 iterations |        0.12 mm |                 off |
| `skin`  | Outer skin surface from CT                              |   −300 HU | 1.4 mm |  3.2 mm |      500 mm³ |            off |      25 iterations |        0.35 mm |       10 iterations |
| `auto`  | MR, CBCT, ultrasound, or other uncalibrated intensities |      Otsu | 1.0 mm |  2.0 mm |       50 mm³ |            off |      20 iterations |        0.25 mm |       40 iterations |

Numeric thresholds are inclusive lower bounds. `auto` calculates Otsu's
threshold. HU presets are treated as calibrated HU only when a DICOM CT series
has complete, consistent rescale evidence. Use an intentional `--threshold` or
`--preset auto` for uncalibrated values.

`--median-mm`, `--opening-mm`, `--closing-mm`, `--min-island-mm3`, and
`--all-islands` control segmentation cleanup. Surface-specific controls are:

- `--resample-mm` for the surface grid;
- `--mask-smooth-mm` for physical Gaussian occupancy smoothing;
- `--mesh-smooth-iters` and `--post-mesh-smooth-iters` for guarded MeshLib
  relaxation before and after simplification;
- `--simplify-error-mm` for MeshLib's estimated QEM deviation limit;
- `--components all|largest` for output shells;
- `--destep auto|all|band` for optional final stair-step fairing;
- `--no-cap` to leave field-of-view openings uncapped.

Set a numeric control to `0` to disable its stage. Mask smoothing can move
boundaries, close gaps, join regions, or erase thin structures. Simplification's
error is an estimate for that step, not a certified Hausdorff bound.

The `teeth` preset keeps all surviving islands and surface components; `bone`,
`skin`, and `auto` keep the largest by default. An explicit component option
overrides the preset.

## Using an external labelmap

Use `labelmap extract` when another program already segmented the image:

```sh
medsurface labelmap extract segmentation.nii.gz -o surface.stl
```

The command requires one direct NIfTI, NRRD, or MetaImage file. DICOM,
directories, volume selectors, and presets are absent. Voxel values must be
finite, non-negative integers with at least one nonzero voxel. Integer-valued
floating-point files are accepted; fractional probability maps are not. By
default every nonzero value becomes one foreground class, so a multilabel file
creates one mesh with multiple structures.

Labelmap extraction skips intensity thresholds, Otsu, morphology, and mask-island
removal. Its independent surface defaults are native grid, `0.8 mm` mask
smoothing, 20 pre-simplification relaxation iterations, a `0.25 mm`
simplification limit, no post-simplification relaxation, and all surviving
surface components. Use the same extraction-only surface options to override
those defaults.

For example, select several structures in one TotalSegmentator multilabel file,
then extract them:

```sh
TotalSegmentator -i scan.nii.gz -o selected.nii.gz \
  --ml --roi_subset skull vertebrae_C1
medsurface labelmap extract selected.nii.gz -o selected.stl
```

TotalSegmentator is not installed or run by medsurface.

### Selecting labels

`--labels` restricts the foreground to some label IDs. Lists and inclusive
ranges are accepted:

```sh
medsurface labelmap extract segmentations.nii.gz --labels 25-50 -o spine.stl
```

Selecting every label present gives exactly the same mesh as omitting
`--labels`.

### One mesh per label

`--split DIR` writes one validated mesh per label into `DIR`, named
`NNN_name.stl` from the label ID and its name:

```sh
medsurface labelmap extract segmentations.nii.gz --labels 25-50,91 \
  --split parts/ -o parts/combined.stl --json parts.json
```

- Names come from `--label-names names.json` (an object such as
  `{"91": "skull"}`) or, when absent, from a label table embedded in a NIfTI
  input the way TotalSegmentator writes it. Unnamed labels become
  `NNN_label-NNN.stl`.
- With `--split`, `-o` is optional. When given, it receives one combined mesh of
  every selected label, extracted from their union, and its extension sets the
  format of the per-label files. Without `-o` the files are `.stl`.
- Each label is cropped to its bounding box plus a margin covering mask
  smoothing and resampling, so each mesh matches what a separate
  `--labels ID` extraction would produce while running much faster. A label that
  reaches the image boundary is capped exactly as it would be uncropped.
- Selected labels without voxels are skipped with a warning. A label that fails
  is reported and does not stop the others; the command then exits `1`.
- The JSON report lists every mesh with its quality report, the combined mesh,
  skipped and failed labels, and per-label voxel counts, volumes, and boundary
  contact.

## Reducing stair-step ripples

Slice terraces from CT or MR can remain visible as broad, shallow ripples on
smooth anatomy even after mask smoothing and mesh relaxation. Those stages are
too local for ripples several millimetres long, and stronger global smoothing
erases thin anatomy. `--destep` adds one final, optional fairing stage to
`extract` and `labelmap extract`. It is off unless requested:

```sh
medsurface labelmap extract skull.nrrd -o skull.stl --destep auto
```

Choose the faired region:

| region | faired vertices                                                                   |
| ------ | --------------------------------------------------------------------------------- |
| `auto` | Broad surfaces; detail such as teeth and rims, and the area around it, stays put  |
| `band` | A coordinate band you choose along one model axis                                 |
| `all`  | Every vertex, including detail                                                    |

`band` needs `--destep-full-mm` and `--destep-frozen-mm`, in model millimetres
along `--destep-axis` (default `z`). Vertices at or beyond the frozen
coordinate never move, those at or beyond the full coordinate are fully faired,
and a smooth ramp joins them. Their order selects the faired side. To choose
values, read the surface's `bbox_min` and `bbox_max` from
`medsurface validate model.stl --json`, then place the ramp with margin between
the rippled area and the detail you must keep:

```sh
medsurface extract scans/ -o skull.stl --destep band \
  --destep-axis z --destep-full-mm -555 --destep-frozen-mm -600
```

`--destep-iters` (default 600) sets how far the fairing reaches and
`--destep-max-mm` (default 1 mm) limits every vertex's displacement. Raise the
iterations until the bands clear, and inspect the result with smooth shading
and low-angle light; soft lighting and decimated previews hide the bands. The
reach is measured in triangles, not millimetres, so meshes with smaller
triangles need more iterations. Too many iterations slowly inflate the surface
until vertices meet the displacement limit; the command warns when more than
10% of vertices reach it.

Fairing deliberately changes the faired anatomy. It flattens shallow features
such as sutures along with the ripples, and `all` also rounds teeth and edges up
to the displacement limit. `auto` is a curvature heuristic. It freezes tightly
curved detail, keeps the surroundings of large detail such as the face frozen
for about 10 mm, and fairs broad surfaces such as the cranial vault. Tiny
isolated bumps are faired with their surroundings. Check the faired and frozen
percentages it logs, and use `band` when it protects too much or too little. Fairing that would fold or intersect the surface is reverted locally;
if no clean result is possible, the unfaired surface is kept with a warning.

## Fusion

Fusion publishes a labelmap volume and never writes a mesh. Valid destinations
are `.nii`, `.nii.gz`, `.nrrd`, and `.mha`. The result is binary `uint8` unless
`labelmap fuse --preserve-labels` keeps label IDs.

### Fuse intensity volumes

`fuse` independently thresholds and cleans each selected intensity input before
registration:

```sh
medsurface fuse fixed.nii.gz moving.mha -o fused.nrrd \
  --fixed-threshold auto --moving-threshold 420
```

It accepts the preset's segmentation controls: threshold, median, opening,
closing, island floor, and island policy. The preset's surface grid, mask
smoothing, component selection, mesh smoothing, simplification, and capping
settings are irrelevant and have no fusion flags.

### Fuse labelmaps

`labelmap fuse` validates both masks, maps every nonzero value to foreground,
then registers and unions them:

```sh
medsurface labelmap fuse fixed-selected.nii.gz moving-selected.nii.gz \
  -o fused.mha
```

Different positive IDs are equivalent. Source IDs, semantic names, and metadata
are not preserved in the derived labelmap.

### Fuse several labelmaps and keep their labels

`labelmap fuse` accepts more than one moving labelmap. Each one is registered
to the fixed labelmap with the same quality gates. `--preserve-labels` keeps
label IDs instead of writing a binary union:

```sh
medsurface labelmap fuse head.nii.gz neck.nii.gz chest.nii.gz \
  --preserve-labels -o body.nii.gz --json body.json
medsurface labelmap extract body.nii.gz --split parts/ -o parts/body.stl
```

- Registration uses the union of each input's labels (or of `--labels`).
- Every label of every input is antialiased and resampled onto the shared grid.
  A voxel keeps the label with the highest occupancy when that occupancy is
  above `0.5`. For a single label this is the ordinary union; where two labels
  touch, the boundary is split by their occupancies.
- Without `--grid-mm`, the grid uses the finest input spacing. Memory stays at
  two grid-sized arrays however many labels and inputs are fused.
- The output is `uint8`, or `uint16`/`uint32` when label IDs need it. NIfTI
  outputs carry the input label table (and `--label-names`) so `--split`
  names its files.
- Without `--preserve-labels`, every selected label is foreground `1` and the
  output is a binary union on a `0.4 mm` grid unless `--grid-mm` is given.

### Shared registration and grid

Both commands:

1. Rigidly register moving foreground to fixed foreground.
2. Enforce shared-field volume, symmetric overlap, and Dice gates unless
   `--force` is supplied.
3. Antialias and linearly resample both masks onto an axis-aligned isotropic
   `--grid-mm` lattice in the fixed input's physical coordinate system.
4. Take the voxelwise maximum of the fractional occupancies and classify values
   strictly greater than `0.5` as foreground `1`; all others are `0`.
5. Verify and atomically publish the result through the same volume writer used
   by strict conversion.

The default grid spacing is `0.4 mm`. A coarser grid can erase structures; a
finer grid has cubic memory cost. Fusion does not smooth, pad, run marching
cubes, select surface components, simplify, cap, or validate a mesh.

After success, the CLI prints a shell-safe extraction hint unless `--quiet`
suppresses normal output. The complete workflow remains explicit:

```sh
medsurface fuse fixed-scan/ moving-scan/ -o fused.nii.gz
# Inspect or edit fused.nii.gz.
medsurface labelmap extract fused.nii.gz -o fused.stl
```

The labelmap path is identical:

```sh
medsurface labelmap fuse fixed-mask.nrrd moving-mask.mha -o fused.nrrd
# Inspect or edit fused.nrrd.
medsurface labelmap extract fused.nrrd -o fused.stl
```

The later extraction uses normal labelmap defaults. No intensity-preset surface
settings are embedded in or transferred through the fused volume.

With `--json FILE`, `fuse` reports output format/compression, scalar type,
components, grid geometry, fixed/moving/fused foreground voxel counts and
volumes, complete registration metrics and transform, warnings, segmentation
settings, and fixed/moving provenance. `labelmap fuse` reports the same output
and grid fields, per-label voxel counts and volumes, one registration per
moving input, and per-input provenance. Neither contains mesh-quality or
surface-finishing fields.

## Input requirements and limitations

Supported file inputs are NIfTI (`.nii`, `.nii.gz`), NRRD (`.nrrd`, `.nhdr`),
and MetaImage (`.mha`, `.mhd`). Detached headers must retain every referenced
payload. HDF5, headerless raw data, and `.bin` are not accepted as independent
volumes because their geometry and voxel interpretation are not self-describing.

Inputs must be real-valued, scalar, three-dimensional images with at least four
voxels on each axis, finite origin and direction values, positive finite
spacing, and a nonsingular direction matrix. A supported file with a readable
but invalid header appears in `list` with its reason.

Source volumes and planned processing/fusion grids above 500 million voxels are
refused before pixel reading or allocation. This is a coarse emergency ceiling,
not a peak-memory estimate. Prefer a coarser upstream volume or fusion grid;
expert users may pass `--allow-large-volume` at the risk of allocation failure.

Classic single-frame DICOM stacks have additional limitations:

- Enhanced multi-frame objects and vendor mosaics are unsupported.
- Compressed pixels require a codec in the installed SimpleITK/GDCM build.
- Consistent oblique stacks are supported; gantry tilt, nonparallel slices,
  duplicate positions, missing geometry, and varying dimensions are not fully
  validated.
- Localizers and short stacks are rejected. Moderate spacing irregularity is
  warned about and regularized; severe inconsistency is rejected.
- Pixel padding and `MONOCHROME1` have no project-specific handling.
- `FrameOfReferenceUID` is not used. Fusion always registers moving foreground
  rather than trusting cross-series coordinates.
- Sampling interval and partial-volume effects limit recoverable thin anatomy.

See the [DICOM compatibility table](technical-reference.md#dicom-compatibility)
for exact behavior.

## Long-running commands

`list`, `convert`, `extract`, `fuse`, both labelmap commands, `validate`, and
`repair` announce the current stage. Interactive terminals show an indeterminate
spinner and elapsed time; redirected output receives persistent ANSI-free lines.
No percentage is invented when processing libraries do not expose one.

Use `--quiet` with processing commands to suppress normal progress. Warnings and
errors remain visible.

Programs driving `labelmap extract` or `labelmap fuse` can pass
`--progress json` to receive JSON lines on stderr instead of the human display:
`stage_start`/`stage_end` (with `seconds`), `label_start`/`label_end`/
`label_failed` (with the label ID, position, and total), `combined_start`/
`combined_end`, `status`, `warning`, and `error` events. The human summary on
stdout is unchanged and can still be suppressed with `--quiet`. `list --json`, `validate --json`, and `repair --json`
write machine-readable JSON to stdout. Conversion, extraction, and fusion use
`--json FILE` for a separate report.

Press `Ctrl+C` once to cancel. Temporary files are cleaned before control returns
to the terminal. Native image or mesh work may delay cancellation until control
returns to Python.

Numeric options reject non-finite or out-of-range values. Unsupported output
extensions are rejected before discovery or pixel loading. Fusion commands do
not expose surface-only options, so there is no output-dependent processing
branch.

## Understanding validation

Mesh extraction and repair validate both the in-memory surface and its serialized
temporary file. A mesh is accepted only when MeshLib imports it as watertight,
consistently wound, and enclosing a volume, with no holes, boundary edges,
disoriented faces, or self-intersecting faces. Multiple closed components are
allowed. Validation describes the MeshLib-imported representation.

Strict conversion and fusion use a different image contract: the temporary
volume is read back and compared for voxel digest, scalar pixel type, component
count, dimensions, and spacing/origin/direction within `1e-5`. Compression is
also checked against the suffix. These checks establish faithful storage, not
anatomical correctness.

Conversion, extraction, and fusion accept `--json FILE`; validation and repair
offer plain `--json` on stdout. Reports are written to temporary siblings and
atomically replaced. Output and report paths must differ and cannot overwrite
any discovered input. Conversion and fusion make primary publication
transactional with a requested report: report failure restores an existing
primary output or removes a newly created one.

Surface safeguards protect both relaxation passes and simplification against
self-intersections and inconsistent winding. Unsafe local vertices remain at
their valid input positions; a stage is discarded when local protection cannot
produce a clean result. Occupancy smoothing occurs before extraction and can
still change anatomy or topology even when the resulting mesh is valid.

Structural and storage checks do not establish manufacturability, dimensional
accuracy, anatomical truth, or fitness for a clinical purpose.
