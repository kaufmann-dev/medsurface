# medsurface

[Install](#install) · [Quick start](#quick-start) · [Commands](#commands) ·
[User guide][user-guide] · [Technical reference][technical-reference]

Turn a medical image volume into an STL, PLY, or OBJ surface mesh from the
command line. DICOM series, NIfTI, NRRD, and MetaImage inputs are supported.

> **Safety:** `medsurface` is not validated for diagnosis, treatment
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

First inspect the available volumes:

```sh
medsurface list ~/scans/head-ct
```

Then convert the recommended volume:

```sh
medsurface convert ~/scans/head-ct -o skull.stl
```

The default `bone` preset extracts values at or above 300. It is treated as HU
only when DICOM CT metadata provides sufficient calibration evidence; other
inputs receive a warning. Use the unique integer ID printed by `list` to select
a different volume. A DICOM directory spanning multiple modalities, such as CT
and MR, deliberately has no automatic default because those volumes are not
comparable acquisitions:

```sh
medsurface convert ~/scans/head-ct --volume 1 -o skull.stl
```

A direct `.nii`, `.nii.gz`, `.nrrd`, `.nhdr`, `.mha`, or `.mhd` file needs no
selector:

```sh
medsurface convert segmentation-input.nrrd -o surface.stl
```

Read [Choosing a preset][presets] before converting other tissues or
uncalibrated data.

## Commands

| command                                      | purpose                                               |
| -------------------------------------------- | ----------------------------------------------------- |
| `medsurface list INPUT`                      | Show every supported volume and any automatic default |
| `medsurface presets`                         | Show the available tissue presets                     |
| `medsurface convert INPUT -o MODEL.stl`      | Convert one volume to a surface mesh                  |
| `medsurface merge FIXED MOVING -o MODEL.stl` | Register and combine two volumes of the same subject  |
| `medsurface validate MODEL.stl`              | Report mesh quality without changing the file         |
| `medsurface repair MODEL.stl -o FIXED.stl`   | Repair an open or non-manifold mesh                   |

Run `medsurface COMMAND --help` for every option and `medsurface --version` for
the installed version. Output format follows the extension: `.stl`, `.ply`, or
`.obj`; unsupported output extensions are rejected before image processing.

## Common workflows

Choose a different tissue preset or an explicit threshold:

```sh
medsurface convert scans/ --preset teeth -o teeth.stl
medsurface convert scans/ --threshold 250 -o bone-250hu.stl
medsurface convert scan.mha --preset auto -o automatic.stl
```

Merge two scans after confirming that they show the same person and anatomy:

```sh
medsurface merge facial-ct/ sinus-ct/ \
  --fixed-volume 1 --moving-volume 1 -o skull.stl
```

The first input defines the output coordinate frame. `merge` cannot verify
subject identity, so confirm it yourself. Registration-quality checks can
refuse geometrically unsafe input; `--force` overrides only those checks and can
create a plausible-looking but incorrect model.

## Validate and repair

`convert` and `merge` validate the written mesh by default. Inspect an existing
mesh without changing it:

```sh
medsurface validate model.stl
```

Repair writes and then validates a separate output file:

```sh
medsurface repair broken.stl -o repaired.stl
```

A valid mesh can still be anatomically wrong. Read the [input limitations][input-limitations]
and the [validation explanation][validation] before relying on an output.

[user-guide]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/user-guide.md
[technical-reference]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/technical-reference.md
[safety]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/user-guide.md#safety-and-privacy
[presets]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/user-guide.md#choosing-a-preset
[input-limitations]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/user-guide.md#input-requirements-and-limitations
[validation]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/user-guide.md#understanding-validation
