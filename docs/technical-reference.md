# Technical reference

This document describes implementation details and limitations that are useful
for auditing results or contributing to `dicom-surface`. The project
[README](../README.md) is the user guide.

## Processing pipeline

`convert` performs these stages:

1. Discover DICOM instances and group them by series UID and orientation.
2. Order slices by `ImagePositionPatient` projected onto the slice normal.
3. Load the ordered stack with SimpleITK and resolve the intensity threshold.
4. Apply island filtering, median filtering, optional opening, and closing.
5. Optionally apply a print profile on a selected isotropic mask grid.
6. Pad field-of-view boundaries, extract the 0.5 isosurface with VTK Flying
   Edges, smooth, select surface components, and optionally decimate.
7. Transform vertices into DICOM patient LPS coordinates, write the requested
   format, and validate the written file.

The output is normally a closed surface because the mask is padded with
background before extraction. When anatomy touches the scan boundary, the cap is
flat and a warning explains that missing anatomy was not recovered.

## DICOM compatibility

pydicom reads headers for discovery. SimpleITK's GDCM-backed reader loads pixel
data. The loader is for classic image series, not arbitrary DICOM objects.

| case | behavior |
|---|---|
| Classic single-frame stacks | Supported and covered by generated uncompressed CT tests |
| Enhanced multi-frame objects | Unsupported; one file is counted as one instance and per-frame geometry is not parsed |
| Compressed transfer syntaxes | Delegated to codecs in the installed SimpleITK/GDCM build; not tested here |
| `RescaleSlope` / `RescaleIntercept` | Applied by SimpleITK; CT presets rely on calibrated HU |
| Pixel padding and `MONOCHROME1` | No explicit project handling; untested |
| Consistent oblique stacks | Supported; slice normal and direction cosines are preserved |
| Gantry tilt or nonparallel slices | No correction or complete geometry validation; unsupported without independent checks |
| Irregular spacing | Warned and regularized when moderate; rejected when spread exceeds `max(0.1 mm, 0.5 × median spacing)` |
| Duplicate positions | Not rejected; unsupported |
| Missing geometry tags or varying dimensions | May fail or load incorrectly; unsupported |
| Localizers | Rejected when `ImageType` contains `LOCALIZER`; short stacks are also rejected |
| Reformats | A consistent classic stack can load; axial data is preferred after voxel-volume ranking |
| Vendor mosaics | Unsupported |

`FrameOfReferenceUID` is not read. `merge` always registers the moving scan and
does not assume cross-series coordinates already align.

When `ConvolutionKernel` is present, a case-insensitive warning heuristic uses
the literal pattern
`(?:^|[^0-9])(?:[BHUY]r?|BONE|EDGE|LUNG)\s*_?([6-9]\d)`. It recognizes names
such as `Hr68` and `B70f`; vendor naming outside that pattern can be missed or
misclassified. A match only produces a sharp-kernel warning. It does not change
series ranking, thresholds, or processing.

## Morphology and physical units

Median, opening, and closing parameters are kernel extents in millimetres. Each
axis is floored to the largest integer radius whose realized `(2r+1) × spacing`
does not exceed the request. A zero radius is an identity operation on that axis.

Selective print-profile dilation is different: its radius is rounded up so a
positive request cannot silently become zero. The mask is first resampled to an
isotropic printability grid when needed to prevent a thick source slice from
turning the dilation ball into a large ellipsoid.

Fractional occupancy is linearly interpolated and thresholded at 0.5. Resampling
can erase thin structures or change components, cavities, tunnels, and genus. It
does not preserve anatomical topology.

## Print-profile behavior

Profiles add a closing target, an island-volume floor, and one minimum-feature
target:

| profile | closing | island floor | feature target |
|---|---:|---:|---:|
| `anatomical` | 0 | 0 | 0 |
| `resin` | 3.2 mm | 100 mm³ | 0.6 mm |
| `fdm` | 4.8 mm | 200 mm³ | 1.2 mm |

There is no separate thickening distance. For target `T`, production computes:

```text
r         = T / 2
core      = opening(mask, r)
reachable = dilate(core, r)
thin      = mask AND NOT reachable
output    = mask OR dilate(thin, r)
```

This avoids inflating large solid regions, but it has a known blind spot: a
short fragile feature inside the thick core's dilation reach is excluded from
`thin`. Controlled FDM-grid fin and bridge cases remained 0.3 mm thick and did
not trigger the global 5% thin-material warning. No final-mesh thickness
measurement exists, so profile values are targets rather than exact STL or
manufacturing guarantees.

## Surface extraction and finishing

The isosurface implementation is `vtkFlyingEdges3D` at 0.5. Padding normally
closes the volume boundary, but Flying Edges can emit degenerate triangles and
mesh validity is measured rather than assumed.

Smoothing uses `vtkWindowedSincPolyDataFilter`. The preset controls initial
iterations, passband, and post-decimation iterations. Smoothing moves surfaces,
and the project does not provide a general deviation bound.

Decimation uses PyMeshLab's
`meshing_decimation_quadric_edge_collapse` with boundary, normal, and topology
preservation enabled. `targetfacenum` is a requested budget, not a guarantee
that every constrained mesh can reach it. Meshes already under budget are left
unchanged. Boundary and non-manifold edges are compared before and after
decimation.

## Merge registration and gates

Registration uses unmodified anatomical masks even when a print profile is
requested:

1. Resample each mask onto a separate axis-aligned 2.0 mm world lattice.
2. Find the best integer 3-D translation with FFT cross-correlation.
3. Extract smoothed surface vertices and normals.
4. Run two point-to-plane ICP passes. The first keeps the closest 45%; the
   second widens to `clamp(0.9 × matched_fraction, 0.45, 0.95)`.

Moving vertices are sampled to at most 60,000 with seed 0. ICP uses a 4 mm
correspondence limit and up to 80 iterations per pass.

| gate | definition | threshold |
|---|---|---:|
| Symmetric surface overlap | Weaker of moving→fixed and fixed→moving fractions within the other scan's field of view and 4 mm of a surface vertex | 0.60 |
| Shared-field Dice | Binary-mask Dice on a 1.5 mm world grid, restricted to common acquired coverage | 0.55 |
| Shared field of view | Common coverage volume on that 1.5 mm grid | 20,000 mm³ |

Point-to-plane RMS and median are reported but not gated. The thresholds were
selected on a small exploratory development set, not clinically calibrated.

After registration passes, the current print merge profiles each scan before
resampling fractional occupancies and taking their maximum. Only dilation
distributes over union; the whole profile does not. Controlled offset-plate
cases showed over-thickening, while fusion-first closing could incorrectly seal
an inter-scan gap. Fusion-first morphology could also operate on a fused grid of
up to 800 million voxels instead of separate printability grids capped at 300
million voxels each, increasing peak memory. The production order remains
unchanged.

## Patient comparison

`merge` reads `PatientID`, `IssuerOfPatientID`, `PatientName`,
`PatientBirthDate`, and `PatientSex` from the first instance of each series.

- A present name, birth-date, or sex conflict refuses the merge.
- Matching IDs corroborate identity unless both issuers are present and differ.
- A matching normalized name corroborates identity when birth dates match or
  are not both present.
- Differing IDs without earlier corroboration refuse the merge.
- Missing or unverifiable identity warns and proceeds.

`--force` overrides modality mismatch, a different-patient verdict, and failed
registration gates. It does not override selecting the same series twice,
loading errors, or an excessive fused grid.

Patient-field values are not logged or stored in provenance. Series UIDs,
descriptions, paths, and derived anatomy can still be identifying.

## Validation and repair

Validation reports triangle/vertex counts, components, watertightness, winding,
boundary edges, non-manifold edge uses, degenerate faces, genus when defined,
volume when closed, bounding box, and self-intersecting faces. The same complete
check runs after `convert`, `merge`, and `repair`, and for `validate`.

A report is valid only when the mesh is watertight, consistently wound, a valid
enclosed volume, and has zero boundary edges, non-manifold edge uses, degenerate
faces, and self-intersecting faces. Multiple closed components are allowed. A
failed or unavailable self-intersection measurement makes validation incomplete
and therefore invalid. Commands return exit status 1 for an invalid report;
`convert` and `merge` retain `--no-validate` as an explicit full bypass.

STL facets do not encode shared topology. Trimesh processing welds positions at
digits derived from its merge tolerance before topology is calculated. Duplicate
faces are not removed and inconsistent winding is reported rather than repaired.

STL, PLY, and OBJ are loaded by Trimesh. VTP is loaded and triangulated with VTK,
then passed through the same Trimesh metric pipeline. VTP self-intersection checks
pass the in-memory triangle mesh to PyMeshLab because PyMeshLab cannot load VTP
directly.

`repair` uses MeshLib to unite vertices within 1e-6, fix multiple edges,
decimate degeneracies, and fill holes. It writes a new file; validation never
mutates the input it measures. A repaired file remains on disk if it does not
pass the complete validation contract.

## Dependencies and development

The repository uses uv 0.11.28 for dependency resolution, environments, command
execution, and CI. `uv.lock` covers Python 3.10–3.13; local development defaults
to Python 3.12.

```sh
uv sync --locked
uv run pytest -q
```

Hatchling is the PEP 517 build backend, and uv orchestrates the workflow.

On Linux, the PyMeshLab wheel expects `libGL.so.1`, and its meshing plugin used
for decimation expects `libOpenGL.so.0`. Distribution package names differ; on
Debian and Ubuntu they are provided by `libgl1` and `libopengl0`. CI installs
both before `uv sync --locked`.

## Verification

Synthetic regression tests cover series grouping, geometry, segmentation,
surface extraction, registration, validation, VTP, repair, and print-profile
primitives. Run them with `uv run pytest -q`.
