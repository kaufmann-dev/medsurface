# dicom-surface

[Install](#install) · [Quick start](#quick-start) · [Commands](#commands) ·
[User guide][user-guide] · [Technical reference][technical-reference]

Turn a DICOM image series into an STL, PLY, or OBJ surface mesh from the
command line.

> **Safety:** `dicom-surface` is not validated for diagnosis, treatment
> planning, or other clinical decisions. Segmentation and mesh processing can
> change or omit anatomy. Review the source images and output independently.
> See [Safety and privacy][safety].

## Install

Python 3.10 or newer and [uv](https://docs.astral.sh/uv/) are required.

```sh
uv tool install .
```

Mesh finishing and validation use the pinned MeshLib Python package; no system
OpenGL libraries are required.

## Quick start

First inspect the available image series:

```sh
dicom-surface list ~/scans/head-ct
```

Then convert the recommended series:

```sh
dicom-surface convert ~/scans/head-ct -o skull.stl
```

The default `bone` preset extracts CT voxels at or above 300 HU. Use the unique
row ID in the `ID` column, a complete SeriesInstanceUID, or part of the
description printed by `list` to select a different series:

```sh
dicom-surface convert ~/scans/head-ct --series 1 -o skull.stl
```

Read [Choosing a preset][presets] before converting other tissues or non-CT
data. For a printable model, also read [Print profiles][profiles].

## Commands

| command                                        | purpose                                             |
| ---------------------------------------------- | --------------------------------------------------- |
| `dicom-surface list DICOM_DIR`                 | Show every image series and the recommended default |
| `dicom-surface presets`                        | Show tissue presets and print profiles              |
| `dicom-surface convert DICOM_DIR -o MODEL.stl` | Convert one series to a surface mesh                |
| `dicom-surface merge DIR_A DIR_B -o MODEL.stl` | Register and combine two scans of the same person   |
| `dicom-surface validate MODEL.stl`             | Report mesh quality without changing the file       |
| `dicom-surface repair MODEL.stl -o FIXED.stl`  | Repair an open or non-manifold mesh                 |

Run `dicom-surface COMMAND --help` for every option. Output format follows the
extension: `.stl`, `.ply`, or `.obj`.

## Common workflows

Choose a different tissue preset or an explicit threshold:

```sh
dicom-surface convert scans/ --preset teeth -o teeth.stl
dicom-surface convert scans/ --threshold 250 -o bone-250hu.stl
```

Prepare selected thin regions for resin or FDM printing:

```sh
dicom-surface convert scans/ --print-profile resin -o resin-skull.stl
```

Merge two scans after confirming that they show the same person and anatomy:

```sh
dicom-surface merge facial-ct/ sinus-ct/ \
  --series-a 1 --series-b 1 -o skull.stl
```

The first scan defines the output coordinate frame. Identity and registration
checks can refuse unsafe input. `--force` overrides those checks and can create
a plausible-looking but incorrect model.

## Validate and repair

`convert` and `merge` validate the written mesh by default. Inspect an existing
mesh without changing it:

```sh
dicom-surface validate model.stl
```

Repair writes and then validates a separate output file:

```sh
dicom-surface repair broken.stl -o repaired.stl
```

A valid mesh can still be anatomically wrong. Read the [input limitations][input-limitations]
and the [validation explanation][validation] before relying on an output.

[user-guide]: https://github.com/kaufmann-dev/dicom-surface/blob/main/docs/user-guide.md
[technical-reference]: https://github.com/kaufmann-dev/dicom-surface/blob/main/docs/technical-reference.md
[safety]: https://github.com/kaufmann-dev/dicom-surface/blob/main/docs/user-guide.md#safety-and-privacy
[presets]: https://github.com/kaufmann-dev/dicom-surface/blob/main/docs/user-guide.md#choosing-a-preset
[profiles]: https://github.com/kaufmann-dev/dicom-surface/blob/main/docs/user-guide.md#print-profiles
[input-limitations]: https://github.com/kaufmann-dev/dicom-surface/blob/main/docs/user-guide.md#input-requirements-and-limitations
[validation]: https://github.com/kaufmann-dev/dicom-surface/blob/main/docs/user-guide.md#understanding-validation
