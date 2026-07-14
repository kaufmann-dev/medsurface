# User guide

[Project README](../README.md) · [Choosing a volume](#choosing-a-volume) ·
[Choosing a preset](#choosing-a-preset) ·
[Using an external labelmap](#using-an-external-labelmap) ·
[Input limitations](#input-requirements-and-limitations) ·
[Long-running commands](#long-running-commands) ·
[Safety and privacy](#safety-and-privacy) · [Technical reference](technical-reference.md)

This guide explains the user-visible choices and limitations of
`medsurface`. See the [technical reference](technical-reference.md) for the
processing algorithms, metrics, thresholds, and dependency details.

## Safety and privacy

`medsurface` is not validated for diagnosis, treatment planning, or other
clinical decisions. Independently review both the source images and generated
surface.

- Segmentation thresholds, smoothing, field-of-view capping, and scan resolution
  can change or omit anatomy.
- A valid and watertight mesh can still be anatomically wrong.
- Anatomy cut off by the scan is capped flat by default; missing anatomy cannot
  be recovered from the input.
- `merge` does not read identity fields or establish subject identity. Confirm
  that both inputs show the same subject before merging them.
- Both merge commands use rigid registration and are intended for matching
  non-deforming anatomy such as bone. Movement or deformation between
  acquisitions can produce a plausible but incorrect fusion.
- Before `labelmap merge`, also confirm that both masks contain the same selected
  structures. Label values and structure names are not inspected.
- `--force` can bypass failed registration-quality checks and produce a
  plausible-looking but incorrect fusion. It does not bypass selection, loading,
  duplicate-input, or voxel-count errors.
- `--allow-large-volume` bypasses the default voxel-count guard for source,
  resampled, and fused grids. It does not make the operation memory-safe; the
  process or operating system may still terminate when memory is exhausted.
- Series UIDs, descriptions, file paths, and derived anatomy can remain
  identifying even when demographic fields have been removed.
- Treat source data, logs, provenance, and output meshes according to the same
  privacy rules as other personal health data.
- Mesh and JSON report paths cannot alias any discovered medical image input,
  including DICOM instances and detached NRRD or MetaImage payloads.

## Choosing a volume

Run `medsurface list INPUT` immediately before conversion when an input can
contain more than one volume. `INPUT` may be a supported image file or a
directory tree. Directory discovery combines all classic DICOM series and all
`.nii`, `.nii.gz`, `.nrrd`, `.nhdr`, `.mha`, and `.mhd` files into one list. For
example, a directory containing a DICOM stack and `mask.nrrd` produces two rows
in one ID space; neither source type hides the other.

The table retains useful imaging metadata across formats: modality, description,
slice count, voxel spacing, and anatomical plane. DICOM also shows
`SeriesNumber` as `DICOM #`. File-format metadata is used when it is present,
and plane is derived from the stored direction matrix. A missing value is shown
as `-`.

Every row has a unique positive integer `ID`. It is the only selector:

```sh
medsurface convert scans/ --volume 2 -o model.stl
medsurface merge scans-a/ scans-b/ \
  --fixed-volume 1 --moving-volume 3 -o merged.stl
```

Selection stays automatic in the unambiguous cases:

- A direct supported file or any catalog with exactly one usable volume is
  selected automatically.
- A directory whose usable entries are all DICOM with the same modality retains
  the established DICOM ranking: usable stacks first, then smaller voxels,
  axial plane, and finally slice count.
- A DICOM-only catalog spanning multiple modalities has no implicit winner.
  `list` explains the ambiguity and shows every modality; pass the intended ID.
- A catalog with multiple usable entries that includes a file volume has no
  implicit winner. Run `list` and pass its ID.

To merge two volumes found under the same directory, repeat the path and choose
two IDs:

```sh
medsurface merge scans/ scans/ \
  --fixed-volume 1 --moving-volume 2 -o merged.stl
```

`DICOM #` is metadata only because `SeriesNumber` need not be unique. A DICOM
UID can also contain more than one orientation; discovery gives each orientation
its own row. UIDs, descriptions, and filenames are never accepted as selectors.

IDs are deterministic for unchanged contents but local to one discovery result.
Run `list` again after adding, removing, or replacing files instead of reusing an
older ID. In `list --json`, the selector is integer `id`; DICOM-specific values
remain nested under `dicom`.

## Choosing a preset

A preset supplies the segmentation and mesh-finishing defaults. List the
installed values at any time with `medsurface presets`.

| preset  | use it for                                              | threshold | median | closing | island floor | smoothing iterations (initial + final) | simplify error |
| ------- | ------------------------------------------------------- | --------: | -----: | ------: | -----------: | -------------------------------------: | -------------: |
| `bone`  | General CT bone models                                  |    300 HU | 1.0 mm |  2.4 mm |       50 mm³ |                                20 + 40 |        0.25 mm |
| `teeth` | Enamel and dense dentin; keeps separate teeth           |  1,200 HU | 0.6 mm |  0.6 mm |        5 mm³ |                                 10 + 0 |        0.12 mm |
| `skin`  | Outer skin surface from CT                              |   −300 HU | 1.4 mm |  3.2 mm |      500 mm³ |                                25 + 10 |        0.35 mm |
| `auto`  | MR, CBCT, ultrasound, or other uncalibrated intensities |      Otsu | 1.0 mm |  2.0 mm |       50 mm³ |                                20 + 40 |        0.25 mm |

Numeric thresholds are inclusive lower bounds. `auto` calculates a
format-neutral Otsu threshold from the volume instead of assuming calibrated
Hounsfield units. The `bone`, `teeth`, and `skin` values are treated as HU only
when a DICOM CT series has a complete, consistent rescale transform and either
an explicit `HU` rescale type or original non-multienergy CT image metadata.
Other inputs, derived CT without explicit units, inconsistent series, and
ambiguous multienergy CT receive a warning. Use an intentional numeric
`--threshold` or `--preset auto` when values are not calibrated HU.

MeshLib smoothing uses the preset's iteration count and relaxation force. Use
`--smooth-iters` and `--smooth-force` to override them; the built-in presets use
force `0.1`, selected to match the preceding surface finish on the two reference
CT studies.

`teeth` keeps every mask island and surface component that survives its size
floor. The other presets keep only the largest component. `--simplify-error-mm`
sets MeshLib's estimated surface-deviation/QEM limit in model millimetres; it is
not a certified Hausdorff bound. `0` disables simplification. The resulting
triangle count is an outcome, not a target. Simplification protects small
source-surface neighborhoods when collapsing them would create
self-intersections. If no candidate can preserve topology and mesh validity,
the valid higher-resolution surface is retained with a warning.

```sh
medsurface convert scans/ --preset teeth -o teeth.stl
medsurface convert scans/ --preset auto -o uncalibrated.stl
medsurface convert scans/ --threshold 250 -o bone-250hu.stl
```

For `merge`, thresholds belong to their input roles and may differ:

```sh
medsurface merge fixed.nii.gz moving.mha \
  --fixed-threshold auto --moving-threshold 420 -o merged.stl
```

An explicit CLI value overrides the corresponding preset value. Run
`medsurface convert --help` or `medsurface merge --help` for the complete set of
overrides. Both commands accept the shared median, opening, closing, island and
component selection, smoothing, and simplification controls.

## Using an external labelmap

Use `labelmap convert` when TotalSegmentator or another program has already
segmented the image:

```sh
medsurface labelmap convert segmentation.nii.gz -o surface.stl
```

This starts after segmentation. It loads the labelmap, treats zero as background
and every nonzero value as foreground, then extracts, finishes, validates, and
writes one surface mesh. It does not threshold image intensities, run Otsu,
apply median/opening/closing filters, remove mask islands, or interpret label
numbers. Separate nonzero regions are all retained, so one multilabel file can
produce one STL containing multiple structures.

The input must be one direct NIfTI, NRRD, or MetaImage file. Its voxels must be
finite, non-negative integers; integer-valued floating-point images are
accepted, but probability maps and fractional labels are not. Directories,
DICOM, selectors, presets, and structure-name flags are deliberately absent.
The surface uses the normal `bone` finishing defaults and retains every surface
component. Surface controls such as `--smooth-iters`, `--smooth-force`, and
`--simplify-error-mm` remain available.

TotalSegmentator is not installed or run by medsurface. Its default output is a
directory containing one binary `.nii.gz` file per structure; pass any one of
those files directly. To put several selected structures in one mesh, ask
TotalSegmentator for one multilabel file and select the structures there:

```sh
TotalSegmentator -i scan.nii.gz -o selected.nii.gz \
  --ml --roi_subset skull vertebrae_C1
medsurface labelmap convert selected.nii.gz -o selected.stl
```

All nonzero labels in `selected.nii.gz` become one mesh. Selecting every
structure is the same workflow without `--roi_subset`, although a whole-body
surface can be large and contain many disconnected components.

Use `labelmap merge` only for matching segmentations from two acquisitions:

```sh
medsurface labelmap merge fixed-selected.nii.gz moving-selected.nii.gz \
  -o fused.stl
```

It rigidly registers the moving foreground to the fixed foreground, checks the
registration, unions the masks on the shared fusion grid, and creates one mesh
in the fixed labelmap's physical coordinate system. Confirm that both inputs
belong to the same subject and contain the same selected, non-deforming
structures. `--force` bypasses only failed registration-quality gates. Do not
use this command merely to combine separate structure files from one scan;
create one multilabel file upstream instead.

## Input requirements and limitations

Supported file volumes are NIfTI (`.nii`, `.nii.gz`), NRRD (`.nrrd`, `.nhdr`),
and MetaImage (`.mha`, `.mhd`). `.mha` stores header and voxels together;
`.mhd` references a separate payload, so both files must remain together. NRRD
has the equivalent single-file `.nrrd` and paired-header `.nhdr` forms. HDF5,
headerless raw data, and `.bin` are not accepted because their geometry and
voxel interpretation are not self-describing.

File inputs must be real-valued, scalar, three-dimensional images with at least
two voxels on each axis, finite origin/direction values, positive finite spacing,
and a nonsingular direction matrix. Supported-format files that violate these
constraints can still appear in `list` with a reason when their header is
readable. A detached `.mhd` or `.nhdr` whose referenced payload is absent or
unreadable is listed as unusable instead of being selected and failing later.

Source volumes and planned conversion or merge grids above 500 million voxels
are refused before their corresponding pixel read or allocation. This is a
coarse emergency ceiling rather than a memory guarantee: pixel types and
processing stages use different amounts of memory per voxel. Prefer resampling
to a coarser spacing. Expert users can pass `--allow-large-volume` to any
`convert` or `merge` command to bypass every voxel-count ceiling at the risk of
swapping, an allocation failure, or an operating-system OOM termination.

Classic single-frame DICOM stacks additionally have these limitations:

- Enhanced multi-frame DICOM objects and vendor mosaic formats are unsupported.
- Compressed pixel data works only when the installed SimpleITK/GDCM build has a
  codec for its transfer syntax; compressed inputs are not covered by this
  project's tests.
- Consistent oblique stacks are supported. Gantry tilt, nonparallel slices,
  duplicate slice positions, missing geometry, and varying image dimensions are
  not fully validated and should be treated as unsupported without independent
  checks.
- Labelled localizers and short stacks are rejected. Moderately irregular slice
  spacing is warned about and regularized; severe inconsistency is rejected.
- Pixel padding and `MONOCHROME1` have no project-specific handling and are
  untested.
- `FrameOfReferenceUID` is not used. `merge` always registers its second scan
  instead of assuming that coordinates agree across series.
- Through-plane resolution, orientation, and partial-volume effects limit what
  thin anatomy can be inferred. A feature smaller than the sampling interval is
  not necessarily absent, but its geometry cannot be recovered reliably.

The [DICOM compatibility table](technical-reference.md#dicom-compatibility)
gives the exact behavior for each known case.

## Long-running commands

`list`, both `convert` commands, both `merge` commands, `validate`, and `repair`
show the current stage while they work. On an interactive terminal, an
indeterminate spinner and elapsed time make activity visible without inventing
a percentage that the processing libraries cannot measure. When output is
redirected, the same stage changes are written as persistent plain-text lines
without animation or ANSI control sequences.

Use `--quiet` with either `convert`, either `merge`, or `repair` to suppress
normal progress; warnings and failures remain visible. Warnings are printed
when they become known, so identity, rigid-registration, structure-content, and
unverified-HU warnings appear before loading and processing rather than after an
output has been written. Machine-readable `list --json`, `validate --json`, and
`repair --json` suppress progress so stdout contains only JSON.

Numeric processing options reject non-finite and out-of-range values as usage
errors. Unsupported mesh output extensions are rejected before discovery or
image processing. Grid and resampling allocations are also bounded before the
image toolkit is asked to allocate them; increase the requested voxel spacing
if the planned volume is too large.

## Understanding validation

Conversion, merging, and repair validate the in-memory surface, write a
temporary file in the destination directory, validate that serialized file,
and atomically publish it only when both checks pass. A failed operation leaves
an existing destination unchanged. `validate` runs the same checks without
changing its input.

Both `convert` commands and both `merge` commands accept `--json FILE` and
publish the report atomically. The report and mesh must be different files, and
neither may overwrite a discovered image header, detached payload, or DICOM
instance.

A report is valid only when MeshLib imports the mesh as watertight, consistently
wound, and enclosing a volume, with no holes, boundary edges, disoriented faces,
or self-intersecting faces. Multiple closed components are allowed. MeshLib may
normalize unsupported raw face configurations while loading; the report describes
the imported mesh rather than exposing separate raw non-manifold or degenerate
face counters. A failed self-intersection measurement is invalid.

Before writing, conversion and merging also guard both smoothing stages and
simplification against self-intersections. Smoothing still runs every requested
iteration; only vertices in collision neighborhoods retain their pre-smooth
positions. These safeguards are reported as warnings and recorded in JSON
provenance when they are used.

These checks establish mesh structure, not anatomical correctness,
manufacturability, dimensional accuracy, or fitness for a clinical purpose.
