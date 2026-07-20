# medsurface

[Install](#install) · [Quick start](#quick-start) · [Commands](#commands) ·
[User guide][user-guide] · [Technical reference][technical-reference]

Convert medical volumes without changing their image data, extract STL/PLY/OBJ
surface meshes, or register two segmentations into an editable binary labelmap.
DICOM series, NIfTI, NRRD, and MetaImage inputs are supported.

> **Safety:** `medsurface` is not validated for diagnosis, treatment planning,
> or other clinical decisions. Segmentation, registration, and mesh processing
> can change or omit anatomy. Review source images and every output independently.
> See [Safety and privacy][safety].

## Install

Python 3.11 or newer and [uv](https://docs.astral.sh/uv/) are required.

```sh
uv tool install .
```

## Quick start

First inspect the available volumes:

```sh
medsurface list ~/scans/head-ct
```

Then extract the recommended volume's surface:

```sh
medsurface extract ~/scans/head-ct -o skull.stl
```

The default `bone` preset extracts stored values at or above 300. It is treated
as HU only when DICOM CT metadata supplies sufficient calibration evidence;
other inputs receive a warning. Use the unique integer ID printed by `list` when
the catalog has no safe automatic choice:

```sh
medsurface extract ~/scans/head-ct --volume 1 -o skull.stl
```

A direct `.nii`, `.nii.gz`, `.nrrd`, `.nhdr`, `.mha`, or `.mhd` volume needs no
selector. To change only its storage format, use strict conversion:

```sh
medsurface convert scan.mhd -o scan.nrrd
```

`convert` preserves loaded voxel values, scalar pixel type, component count,
dimensions, spacing, origin, and direction. It does not threshold, cast,
resample, or reorient the image. Source metadata is retained when the output
format can represent it; add `--strip-metadata` to remove it while keeping
required format headers.

If another tool already created a segmentation labelmap, start at surface
extraction instead of intensity segmentation:

```sh
medsurface labelmap extract segmentation.nii.gz -o surface.stl
```

## Commands

| command                                                       | purpose                                                            |
| ------------------------------------------------------------- | ------------------------------------------------------------------ |
| `medsurface list INPUT`                                       | Show every supported volume and any automatic default              |
| `medsurface presets`                                          | Show the available tissue presets                                  |
| `medsurface convert INPUT -o VOLUME`                          | Preserve one image while changing its single-file storage format   |
| `medsurface extract INPUT -o MODEL.stl`                       | Segment one intensity volume and extract a surface mesh            |
| `medsurface fuse FIXED MOVING -o LABELMAP`                    | Segment, register, and union two volumes into a binary labelmap    |
| `medsurface labelmap extract MASK -o MODEL.stl`               | Extract one surface from all nonzero labels in one mask            |
| `medsurface labelmap fuse FIXED_MASK MOVING_MASK -o LABELMAP` | Register and union two matching labelmaps into one binary labelmap |
| `medsurface validate MODEL.stl`                               | Report mesh quality without changing the file                      |
| `medsurface repair MODEL.stl -o FIXED.stl`                    | Repair an open or non-manifold mesh                                |

Run `medsurface COMMAND --help` for every option and `medsurface --version` for
the installed version. Mesh destinations are `.stl`, `.ply`, or `.obj`.
Conversion and fusion destinations are atomic single-file `.nii`, `.nii.gz`,
`.nrrd`, or `.mha` volumes. Detached `.nhdr` and `.mhd` remain input-only.

## Common workflows

### Convert a volume without processing it

Convert one selected catalog volume and record the operation:

```sh
medsurface convert scans/ --volume 2 -o selected.nii.gz \
  --json selected-conversion.json
```

`.nii` is uncompressed, `.nii.gz` uses gzip, `.nrrd` uses gzip, and `.mha` uses
zlib. The volume is written to a temporary sibling, read back, and checked for
exact voxel data, type, components, dimensions, and physical geometry before it
atomically replaces the destination. An unrepresentable core image contract is
an error; metadata preservation remains best effort.

### Extract an intensity surface

Choose a preset or explicit threshold when producing a mesh:

```sh
medsurface extract scans/ --preset teeth -o teeth.stl
medsurface extract scans/ --threshold 250 -o bone-250hu.stl
medsurface extract scan.mha --preset auto -o automatic.stl
```

Surface controls such as `--mask-smooth-mm`, `--mesh-smooth-iters`,
`--simplify-error-mm`, `--post-mesh-smooth-iters`, and `--components` belong to
`extract` commands only.

### Fuse intensity volumes, then extract

Confirm that both scans show the same subject and matching non-deforming anatomy,
then create an inspectable binary union before meshing it:

```sh
medsurface fuse facial-ct/ sinus-ct/ \
  --fixed-volume 1 --moving-volume 1 -o skull-union.nrrd
medsurface labelmap extract skull-union.nrrd -o skull.stl
```

`fuse` independently segments each input, rigidly registers moving foreground to
fixed foreground, and publishes scalar `uint8` values `0` and `1` on an
axis-aligned isotropic grid. It never writes a mesh or applies surface controls.
After success it prints a shell-safe `medsurface labelmap extract` command hint,
unless normal output is suppressed with `--quiet`.
Registration-quality checks can refuse unsafe geometry; `--force` overrides only
those gates and can produce a plausible but incorrect labelmap.

### Fuse external labelmaps, then extract

TotalSegmentator and similar tools remain separate optional programs. Every
nonzero label is treated as one foreground class:

```sh
TotalSegmentator -i scan.nii.gz -o selected.nii.gz \
  --ml --roi_subset skull vertebrae_C1
medsurface labelmap extract selected.nii.gz -o selected.stl
```

To register matching masks from two acquisitions, publish the binary union,
optionally edit it, and explicitly extract its surface:

```sh
medsurface labelmap fuse fixed-selected.nii.gz moving-selected.nii.gz \
  -o fused.mha --json fused.json
# Edit fused.mha if needed.
medsurface labelmap extract fused.mha -o fused.stl
```

Different positive source IDs are equivalent and are not preserved. Fusion does
not transfer source metadata or intensity-preset surface settings to the result;
the later `labelmap extract` uses its normal labelmap defaults.

## Validate and repair

Every mesh-producing extraction validates its serialized mesh before atomic
publication. Inspect an existing mesh without changing it:

```sh
medsurface validate model.stl
```

Repair validates both its in-memory result and serialized temporary file, then
atomically publishes a separate output:

```sh
medsurface repair broken.stl -o repaired.stl
```

Primary outputs and optional JSON reports cannot alias each other or any
discovered input, including DICOM instances and detached payloads. Conversion
and fusion roll their primary destination back if requested report publication
fails. A valid mesh or verified volume can still be anatomically wrong; read the
[input limitations][input-limitations] and [validation
explanation][validation].

[user-guide]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/user-guide.md
[technical-reference]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/technical-reference.md
[safety]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/user-guide.md#safety-and-privacy
[input-limitations]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/user-guide.md#input-requirements-and-limitations
[validation]: https://github.com/kaufmann-dev/medsurface/blob/main/docs/user-guide.md#understanding-validation
