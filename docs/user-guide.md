# User guide

[Project README](../README.md) · [Choosing a series](#choosing-a-series) ·
[Choosing a preset](#choosing-a-preset) ·
[Print profiles](#print-profiles) · [Input limitations](#input-requirements-and-limitations) ·
[Long-running commands](#long-running-commands) ·
[Safety and privacy](#safety-and-privacy) · [Technical reference](technical-reference.md)

This guide explains the user-visible choices and limitations of
`dicom-surface`. See the [technical reference](technical-reference.md) for the
processing algorithms, metrics, thresholds, and dependency details.

## Safety and privacy

`dicom-surface` is not validated for diagnosis, treatment planning, or other
clinical decisions. Independently review both the source images and generated
surface.

- Segmentation thresholds, smoothing, field-of-view capping, print preparation,
  and scan resolution can change or omit anatomy.
- A valid and watertight mesh can still be anatomically wrong.
- Anatomy cut off by the scan is capped flat by default; missing anatomy cannot
  be recovered from the input.
- `merge` compares selected patient fields but does not establish identity.
- `--force` can bypass patient, modality, and registration checks and produce a
  plausible-looking but incorrect fusion.
- Series UIDs, descriptions, file paths, and derived anatomy can remain
  identifying even when demographic fields have been removed.
- Treat source data, logs, provenance, and output meshes according to the same
  privacy rules as other personal health data.

## Choosing a series

Run `dicom-surface list DICOM_DIR` immediately before conversion. Discovery is
recursive and can list series from mixed patient or study directory trees. The
`default` status marks the series automatic selection would use.

Every displayed stack has a unique, 1-based `ID`. Pass that ID, a complete
SeriesInstanceUID, or a case-insensitive description substring to `--series`:

```sh
dicom-surface convert scans/ --series 1 -o model.stl
dicom-surface convert scans/ --series 1.2.840.113619.2.55.3.604688435.123 -o model.stl
dicom-surface convert scans/ --series "thin axial" -o model.stl
```

For `merge`, the same syntax applies independently to `--series-a` and
`--series-b`. A row ID refers only to the corresponding directory's latest
discovery result.

`DICOM #` is the source file's `SeriesNumber`. Several unrelated UIDs can carry
the same number, so it is shown as metadata but is not a selector. A UID can
also contain more than one orientation; discovery gives each orientation its
own row ID and records its orientation-part number. A complete UID or
description that matches several rows is rejected as ambiguous and reports the
row IDs to choose from.

Row IDs are deterministic for unchanged contents but are local to a discovery
result. Run `list` again after adding, removing, or replacing files instead of
reusing an older ID. In `list --json`, the selector is the integer `id`; `uid`,
`series_number`, `part`, and `n_parts` remain separate metadata.

## Choosing a preset

A preset supplies the segmentation and mesh-finishing defaults. List the
installed values at any time with `dicom-surface presets`.

| preset        | use it for                                                | threshold | median | closing | island floor | smoothing iterations (initial + final) | triangle target |
| ------------- | --------------------------------------------------------- | --------: | -----: | ------: | -----------: | -------------------------------------: | --------------: |
| `bone`        | General CT bone models                                    |    300 HU | 1.0 mm |  2.4 mm |       50 mm³ |                                20 + 25 |         600,000 |
| `bone-detail` | Maximum detail and measurement; produces very large files |    300 HU | 0.6 mm |  1.2 mm |       20 mm³ |                                  8 + 0 |             off |
| `teeth`       | Enamel and dense dentin; keeps separate teeth             |  1,200 HU | 0.6 mm |  0.6 mm |        5 mm³ |                                 10 + 0 |         300,000 |
| `skin`        | Outer skin surface from CT                                |   −300 HU | 1.4 mm |  3.2 mm |      500 mm³ |                                25 + 10 |         400,000 |
| `auto`        | MR, CBCT, ultrasound, or other uncalibrated intensities   |      Otsu | 1.0 mm |  2.0 mm |       50 mm³ |                                20 + 25 |         600,000 |

Numeric thresholds are inclusive lower bounds. `auto` calculates an Otsu
threshold from the scan instead of assuming calibrated Hounsfield units. CT
presets with HU thresholds are refused on non-CT data unless you supply an
explicit `--threshold`.

`teeth` keeps every mask island and surface component that survives its size
floor. The other presets keep only the largest component. Triangle targets are
requested face-count budgets; constrained meshes may finish above them, and
meshes already below the target are not enlarged.

```sh
dicom-surface convert scans/ --preset teeth -o teeth.stl
dicom-surface convert scans/ --preset auto -o uncalibrated.stl
dicom-surface convert scans/ --threshold 250 -o bone-250hu.stl
```

An explicit CLI value overrides the corresponding preset value. Run
`dicom-surface convert --help` for the complete set of overrides.

## Print profiles

A print profile raises selected morphology settings after the tissue preset is
resolved. It does not replace the preset or change the intensity threshold.

| profile      | intended use                      | closing floor | island floor | feature target |
| ------------ | --------------------------------- | ------------: | -----------: | -------------: |
| `anatomical` | No print-profile changes; default |     unchanged |    unchanged |           none |
| `resin`      | Fine-detail resin printing        |        3.2 mm |      100 mm³ |         0.6 mm |
| `fdm`        | FDM/nozzle printing               |        4.8 mm |      200 mm³ |         1.2 mm |

Closing and island floors raise the preset values only when the profile value
is larger. The feature target selects thin mask regions and adds material around
them; there is no separate user-facing thickening distance.

```sh
dicom-surface convert scans/ --print-profile resin -o resin-skull.stl
dicom-surface convert scans/ --print-profile fdm -o fdm-skull.stl
```

Profile values are mask-processing targets, not guarantees about final mesh or
manufactured wall thickness. Scan sampling, interpolation, surface extraction,
smoothing, decimation, and the printer can all change the realized result.
Inspect the final mesh in a slicer. Short thin features close to thick anatomy
may not be selected for thickening; the [technical
reference](technical-reference.md#print-profile-behavior) describes this known
limitation.

## Input requirements and limitations

The tested input path is a classic single-frame image stack with consistent
geometry.

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

Conversion, merging, and repair validate the file they write. `validate` runs
the same checks without changing its input.

A report is valid only when the mesh is watertight, consistently wound, an
enclosed volume, and has no boundary edges, non-manifold edge uses, degenerate
faces, or self-intersecting faces. Multiple closed components are allowed. A
failed or unavailable self-intersection measurement makes the result incomplete
and therefore invalid.

These checks establish mesh structure, not anatomical correctness,
manufacturability, dimensional accuracy, or fitness for a clinical purpose.
