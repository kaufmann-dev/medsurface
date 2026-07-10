# dicom-surface

[Install](#install) · [Quick start](#quick-start) · [Commands](#commands) ·
[Printing](#prepare-a-model-for-printing) · [Merging](#merge-two-scans) ·
[Safety](#safety-and-privacy) · [Technical documentation](docs/README.md)

Turn a DICOM image series into an STL, PLY, OBJ, or VTP surface mesh from the
command line.

## Install

Python 3.10 or newer and [uv](https://docs.astral.sh/uv/) are required.

```sh
uv tool install dicom-surface
```

On Linux, PyMeshLab needs system OpenGL libraries. Many desktop distributions
already provide them. On Debian or Ubuntu, install them with:

```sh
sudo apt install libgl1 libopengl0
```

## Quick start

First inspect the available image series:

```sh
dicom-surface list ~/scans/head-ct
```

Then convert the recommended series:

```sh
dicom-surface convert ~/scans/head-ct -o skull.stl
```

The command prints mesh-quality results after writing the file. A normal result
should be watertight, consistently wound, and have zero boundary and
non-manifold edges.

To select a different series, use the number, dotted identifier, UID, or part of
its description shown by `list`:

```sh
dicom-surface convert ~/scans/head-ct --series 6 -o skull.stl
```

## Commands

| command                                        | purpose                                             |
| ---------------------------------------------- | --------------------------------------------------- |
| `dicom-surface list DICOM_DIR`                 | Show every image series and the recommended default |
| `dicom-surface presets`                        | Show tissue presets and print profiles              |
| `dicom-surface convert DICOM_DIR -o MODEL.stl` | Convert one series to a surface mesh                |
| `dicom-surface merge DIR_A DIR_B -o MODEL.stl` | Register and combine two scans of the same person   |
| `dicom-surface validate MODEL.stl`             | Report mesh quality without changing the file       |
| `dicom-surface repair MODEL.stl -o FIXED.stl`  | Repair an open or non-manifold mesh                 |

Run `dicom-surface COMMAND --help` for every option.

Output format follows the extension: `.stl`, `.ply`, `.obj`, or `.vtp`.
Coordinates are millimetres in DICOM patient LPS space. STL does not store that
coordinate-system label.

## Choose what to extract

| preset        | use it for                                                | threshold | median | closing | island floor | smoothing iterations (initial + final) | triangle target |
| ------------- | --------------------------------------------------------- | --------: | -----: | ------: | -----------: | -------------------------------------: | --------------: |
| `bone`        | General CT bone models                                    |    300 HU | 1.0 mm |  2.4 mm |       50 mm³ |                                20 + 25 |         600,000 |
| `bone-detail` | Maximum detail and measurement; produces very large files |    300 HU | 0.6 mm |  1.2 mm |       20 mm³ |                                  8 + 0 |             off |
| `teeth`       | Enamel and dense dentin; keeps separate teeth             |  1,200 HU | 0.6 mm |  0.6 mm |        5 mm³ |                                 10 + 0 |         300,000 |
| `skin`        | Outer skin surface from CT                                |   −300 HU | 1.4 mm |  3.2 mm |      500 mm³ |                                25 + 10 |         400,000 |
| `auto`        | MR, CBCT, ultrasound, or other uncalibrated intensities   |      Otsu | 1.0 mm |  2.0 mm |       50 mm³ |                                20 + 25 |         600,000 |

Numeric thresholds are inclusive lower bounds; `auto` calculates an Otsu
threshold from the scan. The triangle target is a requested face-count budget,
not a guaranteed exact count. `teeth` keeps every surviving mask island and
surface component; the other presets keep only the largest.

Examples:

```sh
dicom-surface convert scans/ --preset teeth -o teeth.stl
dicom-surface convert scans/ --threshold 250 -o bone-250hu.stl
```

CT Hounsfield-unit presets are refused on non-CT data unless you provide an
explicit threshold. Use `--preset auto` when intensities are not calibrated HU.

## Prepare a model for printing

Print profiles seal small pores, remove tiny fragments, and add material to
selected thin mask regions:

```sh
dicom-surface convert scans/ --print-profile resin -o resin-skull.stl
dicom-surface convert scans/ --print-profile fdm -o fdm-skull.stl
```

| profile      | intended use                      | closing floor | island floor | feature target |
| ------------ | --------------------------------- | ------------: | -----------: | -------------: |
| `anatomical` | No print-profile changes; default |     unchanged |    unchanged |           none |
| `resin`      | Fine-detail resin printing        |        3.2 mm |      100 mm³ |         0.6 mm |
| `fdm`        | FDM/nozzle printing               |        4.8 mm |      200 mm³ |         1.2 mm |

Closing and island floors raise the preset values only when the profile value
is larger.

These are mask-processing targets, not guarantees about final STL thickness or
the manufactured part. Inspect the result in your slicer. Short thin structures
next to thick anatomy are a documented limitation; see the
[technical reference](docs/technical-reference.md#print-profile-behavior).

## Merge two scans

Use `merge` only for scans of the same person and anatomy:

```sh
dicom-surface merge facial-ct/ sinus-ct/ \
  --series-a 6 --series-b 2 -o skull.stl
```

The first scan defines the output coordinate frame. The second is rigidly
registered before both masks are combined. Identity and registration-quality
checks can refuse unsafe input.

`--force` overrides patient, modality, and registration checks. It can create a
plausible-looking but incorrect model, so use it only after independently
confirming the scans belong together.

## Validate and repair

Validation never repairs the file it measures:

```sh
dicom-surface validate model.stl
dicom-surface validate model.stl --self-intersections
```

Self-intersection checking is slower and opt-in. Watertightness alone does not
rule out self-intersections.

Repair creates a new file and validates it afterward:

```sh
dicom-surface repair broken.stl -o repaired.stl
```

## Input limitations

The tested input path is a classic single-frame DICOM image stack. Enhanced
multi-frame DICOM and vendor mosaic formats are not supported. Compressed pixel
data depends on the codecs included with SimpleITK/GDCM.

Use `list` before converting. It rejects short stacks, labelled localizers, and
severely inconsistent spacing, but it cannot prove that every DICOM geometry is
correct. See the [compatibility table](docs/technical-reference.md#dicom-compatibility).

## Safety and privacy

- Segmentation thresholds, smoothing, capping, print preparation, and scan
  resolution can all change or omit anatomy.
- A watertight mesh can still be anatomically wrong or self-intersecting.
- Field-of-view truncation is capped flat by default and cannot recover missing
  anatomy.
- `merge` compares selected patient fields but does not establish identity.
- Series UIDs, descriptions, paths, and derived anatomy can remain identifying.
- Treat every output mesh as personal health data.

## Development

```sh
uv sync --locked
uv run pytest -q
```

The project pins uv 0.11.28, uses Python 3.12 for local development, and tests
Python 3.10–3.13 in CI.
