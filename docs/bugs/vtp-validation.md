# VTP output failed default validation

Fixed: 2026-07-10 07:39:38 UTC (+0000)

Baseline commit: `f6f90b067f8b628592cfc69ac7a5bc8886a5f5cf`

## Symptom

`dicom-surface convert ... -o model.vtp` wrote a valid VTP file and then failed
its default validation step with:

```text
NotImplementedError: file_type 'vtp' not supported
```

`dicom-surface validate model.vtp` failed for the same reason.

## Cause

VTK wrote VTP, but validation sent every format to Trimesh. Trimesh 4.12.2 does
not load VTP. The collision checker must therefore receive VTP triangles from
memory rather than loading the file itself.

## Fix

Validation loads VTP through `vtkXMLPolyDataReader`, triangulates it with
`vtkTriangleFilter`, and converts the resulting arrays into the same MeshLib
representation used for every quality metric and self-intersection check.

Regression tests cover ordinary VTP validation and the VTP self-intersection
path. STL, PLY, and OBJ continue to use their existing Trimesh loaders.
