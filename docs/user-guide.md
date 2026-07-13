# User guide

[Project README](../README.md) · [Choosing a volume](#choosing-a-volume) ·
[Choosing a preset](#choosing-a-preset) ·
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
- `--force` can bypass failed registration-quality checks and produce a
  plausible-looking but incorrect fusion. It does not bypass selection, loading,
  duplicate-input, or size errors.
- Series UIDs, descriptions, file paths, and derived anatomy can remain
  identifying even when demographic fields have been removed.
- Treat source data, logs, provenance, and output meshes according to the same
  privacy rules as other personal health data.

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
- A directory whose usable entries are all DICOM retains the established DICOM
  ranking: usable stacks first, then smaller voxels, axial plane, and finally
  slice count.
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
Hounsfield units. The `bone`, `teeth`, and `skin` preset values are verified as
HU only for DICOM CT. On every other input they are applied to stored values and
produce a warning; use an intentional numeric `--threshold` or `--preset auto`
when those values are not calibrated HU.

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
`medsurface convert --help` for the complete set of overrides.

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
readable.

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

`list`, `convert`, `merge`, `validate`, and `repair` show the current stage while
they work. On an interactive terminal, an indeterminate spinner and elapsed
time make activity visible without inventing a percentage that the processing
libraries cannot measure. When output is redirected, the same stage changes are
written as persistent plain-text lines without animation or ANSI control
sequences.

Use `--quiet` with `convert`, `merge`, or `repair` to suppress normal progress;
warnings and failures remain visible. Machine-readable `list --json`,
`validate --json`, and `repair --json` suppress progress so stdout contains only
JSON.

## Understanding validation

Conversion and merging validate the in-memory surface, write a temporary file
in the destination directory, validate that serialized file, and atomically
publish it only when both checks pass. `repair` validates the file it writes;
`validate` runs the same checks without changing its input.

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
