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

A direct `.nii`, `.nii.gz`, `.nrrd`, `.nhdr`, `.mha`, or `.mhd` intensity
volume needs no selector:

```sh
medsurface convert scan.nrrd -o surface.stl
```

Read [Choosing a preset][presets] before converting other tissues or
uncalibrated data.

If another tool already created a segmentation labelmap, skip medsurface's
thresholding and cleanup stages:

```sh
medsurface labelmap convert segmentation.nii.gz -o surface.stl
```

Labelmap commands smooth voxel-scale terracing with a `0.8 mm` physical
Gaussian by default. Use `--smooth-mm 0` to preserve the supplied boundary as
literally as possible; smoothing can round or erase structures near its scale.
Every labelmap axis must contain at least four voxels, including when smoothing
is disabled.

## Commands

| command                                                         | purpose                                                        |
| --------------------------------------------------------------- | -------------------------------------------------------------- |
| `medsurface list INPUT`                                         | Show every supported volume and any automatic default          |
| `medsurface presets`                                            | Show the available tissue presets                              |
| `medsurface convert INPUT -o MODEL.stl`                         | Segment one intensity volume and create a surface mesh         |
| `medsurface merge FIXED MOVING -o MODEL.stl`                    | Segment, register, and combine two volumes of the same subject |
| `medsurface labelmap convert MASK -o MODEL.stl`                 | Create one surface from all nonzero labels in one mask         |
| `medsurface labelmap merge FIXED_MASK MOVING_MASK -o MODEL.stl` | Register two matching labelmaps and create their fused surface |
| `medsurface validate MODEL.stl`                                 | Report mesh quality without changing the file                  |
| `medsurface repair MODEL.stl -o FIXED.stl`                      | Repair an open or non-manifold mesh                            |

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
subject identity, so confirm it yourself. Its rigid registration is intended
for matching non-deforming anatomy such as bone. Registration-quality checks
can refuse geometrically unsafe input; `--force` overrides only those checks
and can create a plausible-looking but incorrect model.

Use a segmentation made by another tool, including multiple structures in one
mesh. For example, TotalSegmentator normally writes one binary NIfTI per
structure; `--ml` instead writes one multilabel NIfTI. Select the desired
structures there, then medsurface unions every nonzero label:

```sh
TotalSegmentator -i scan.nii.gz -o selected.nii.gz \
  --ml --roi_subset skull vertebrae_C1
medsurface labelmap convert selected.nii.gz -o selected.stl
```

TotalSegmentator remains a separate, optional program. A single per-structure
file from its default output can also be passed directly. To fuse matching
segmentations from two acquisitions, use:

```sh
medsurface labelmap merge fixed-selected.nii.gz moving-selected.nii.gz \
  -o fused.stl
```

Both labelmaps must contain the same selected structures. This command performs
rigid registration, so it is not a shortcut for combining separate structure
files from one acquisition; create one multilabel input upstream for that case.
Physical smoothing is applied after fusion, so it also attenuates voxel terraces
and small boundary disagreements introduced by the union.

## Validate and repair

Every `convert` and `merge` command validates the written mesh. Inspect an
existing mesh without changing it:

```sh
medsurface validate model.stl
```

Repair validates both its in-memory result and serialized temporary file, then
atomically publishes a separate output:

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
