# VTP output failed default validation

Status: fixed.

## Symptom

`dicom-surface convert ... -o model.vtp` wrote a valid VTP file and then failed
its default validation step with:

```text
NotImplementedError: file_type 'vtp' not supported
```

`dicom-surface validate model.vtp` failed for the same reason.

## Cause

VTK wrote VTP, but validation sent every format to Trimesh. Trimesh 4.12.2 does
not load VTP. PyMeshLab also cannot load VTP directly, so the optional
self-intersection path had the same format gap.

## Fix

Validation now loads VTP through `vtkXMLPolyDataReader`, triangulates it with
`vtkTriangleFilter`, and passes the resulting arrays into the existing Trimesh
metric pipeline. For self-intersection checks, the in-memory vertices and faces
are passed to PyMeshLab instead of asking PyMeshLab to load the VTP file.

Regression tests cover ordinary VTP validation and the VTP self-intersection
path. STL, PLY, and OBJ continue to use their existing Trimesh loaders.
