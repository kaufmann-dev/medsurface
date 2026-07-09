# dicom-surface

Turn a DICOM series into a **watertight** 3D surface mesh.

```console
$ dicom-surface list ~/scans/head-ct
#    MOD  DESCRIPTION                      SLICES VOXEL mm               PLANE     NOTES
--------------------------------------------------------------------------------------------
1    CT   Topogramm 0,60 sag Tr20 MPR           1 -                      -         localizer
2    CT   GS nativ 3,00 ax Hr40 A3 MPR         92 0.345 x 0.345 x 2.000  axial
6    CT   GS nativ 1,00 ax Hr68 A1 MPR        231 0.315 x 0.315 x 0.800  axial     <- default, sharp kernel Hr68
8    CT   GS nativ 1,00 sag Hr68 A1 MPR       242 0.315 x 0.315 x 0.800  sagittal  sharp kernel Hr68
501  CT   Patientenprotokoll                    1 -                      -         not an image series

$ dicom-surface convert ~/scans/head-ct -o skull.stl --preset bone
...
quality:
  triangles           1,916,624
  components          1
  watertight          yes
  boundary edges      0
  non-manifold edges  0
  degenerate faces    0
  volume              384295 mm3
```

Every conversion is validated before it is handed back. If the mesh is not
watertight, the tool says so rather than letting you find out in a slicer.

## Install

```sh
pip install dicom-surface                 # core
pip install 'dicom-surface[repair]'       # + repair of meshes from elsewhere
pip install 'dicom-surface[quality]'      # + self-intersection counting
```

Python 3.10+. No 3D Slicer, no MeshLab, no system packages.

## Usage

```sh
dicom-surface list    DICOM_DIR                    # what's in this folder?
dicom-surface presets                              # what can I ask for?
dicom-surface convert DICOM_DIR -o out.stl         # do the thing
dicom-surface validate mesh.stl                    # is this mesh sound?
dicom-surface repair   mesh.stl -o fixed.stl       # make it watertight
```

Output format follows the extension: `.stl`, `.ply`, `.obj`, `.vtp`.
Coordinates are **LPS** (DICOM patient space), which is what STL consumers expect.

### Presets

| preset | modality | threshold | grid | for |
|---|---|---|---|---|
| `bone` | CT | 300 HU | 0.6 mm | general bone. Denoises without erasing teeth or sutures. |
| `bone-detail` | CT | 300 HU | native | maximum fidelity, very large files |
| `bone-print` | CT | 350 HU | 1.0 mm | smooth low-poly bone for 3D printing |
| `teeth` | CT | 1200 HU | 0.3 mm | enamel and dense dentin |
| `skin` | CT | -300 HU | 1.0 mm | outer skin surface |
| `auto` | any | Otsu | 0.6 mm | MR, CBCT, or any uncalibrated intensity |

Every preset parameter is overridable:

```sh
dicom-surface convert scans/ -o out.stl \
    --series 6 --threshold 250 --closing-mm 3.2 --resample-mm 0.8 --smooth-iters 25
```

### Not every scan is a head CT

- **Series selection.** A study holds scouts, several reconstruction kernels, and
  reformats. `list` shows them all; `--series` takes a number, a UID, or a
  description substring. The default pick prefers small voxels, then the axial
  (acquired) plane, and only then slice count.
- **Modality.** Hounsfield thresholds are meaningful only on CT. Ask for
  `--preset bone` on an MR and the tool refuses rather than emitting a
  confidently wrong mesh. Use `--preset auto` (Otsu) or an explicit
  `--threshold`.
- **Anisotropy.** Filter sizes are given in millimetres and convert to a
  different voxel count on every axis.
- **Multi-part anatomy.** `--all-islands` keeps every bone fragment;
  `--all-components` keeps enclosed cavities (sinuses, marrow space) as separate
  shells.

## Why the output is watertight

Six things silently produce a plausible-looking but wrong mesh. This tool handles
each, and the test suite has a regression test for every one.

**1. Filenames are not slice order.** In the study this was built on, file `1`
held instance 222 and file `231` held instance 226. Sorting by name gives a
scrambled volume that still renders as a convincing blob. Slices are ordered by
`ImagePositionPatient` projected on the slice normal.

**2. Marching cubes does not cap the volume boundary.** When anatomy runs off the
edge of the scan -- a cranial vault truncated by the field of view -- the surface
is left *open* there. The volume is padded with one voxel of background, which
caps it flat. The tool warns when this happens, because that anatomy was never
scanned and no amount of meshing brings it back.

**3. Millimetres are not voxels.** With 0.8 mm slices, `round(1.0mm / 2 / 0.8)`
gives a 3-voxel kernel spanning 2.40 mm -- 2.4x what was asked, quietly erasing
anatomy. Kernels are floored so the realised extent never exceeds the request.
And `2.4 / 0.8` is `2.9999999999999996`, so the flooring needs an epsilon or the
kernel vanishes instead.

**4. An internal cavity is not an island.** Sinuses and trabecular air cells are
*background*, so no labelmap island filter removes them. They survive as separate
closed shells inside the mesh. Only keeping the largest *surface component* drops
them (2,477 of them, on one head CT).

**5. Decimation tears thin walls.** A skull's orbital walls are one voxel thick;
edge collapses weld their opposite faces together. Measured on a 4M-triangle
skull, `vtkQuadricDecimation` leaks defects at *every* reduction -- 6 boundary
edges at 50%, 90 at 85% -- so no backoff converges. `vtkDecimatePro` with
`PreserveTopologyOn` both misses the target and emits non-manifold edges.

Triangle count is therefore controlled by `--resample-mm`, which resamples the
mask onto an isotropic grid before meshing. Marching cubes always returns a
manifold surface, so this cannot introduce a defect. `--target-faces` still
exists, is verified after the fact, and warns if it broke the mesh.

**6. A signed distance field is not free.** It is the textbook way to resample a
mask, but ITK's Maurer transform quantises distance to voxel centres, putting its
zero level half a voxel inside the true boundary -- a 6% volume loss on a 20 mm
sphere. Smoothed fractional occupancy, thresholded at 0.5, keeps the error under
0.5%.

## Sharp kernels

CT reconstruction kernels trade noise against resolution. Sharp ones (Siemens
`Hr68`, `B70f`; anything in the `[BHUY]r?[6-9]x` family) amplify high-frequency
noise by design. Thresholding one at the usual ~200 HU bone value covers the
surface in spurious spikes. `dicom-surface` detects these and warns; the presets
default to 300 HU for exactly this reason.

## Validation

`validate` welds vertices before measuring anything. STL stores a triangle soup --
every triangle carries its own three vertices and no connectivity -- so
"watertight", "manifold" and "hole" are undefined until identical positions are
merged. Without that step every STL on earth reports as 100% boundary edges.

It also never repairs what it measures. `trimesh.split()` builds submeshes with
`repair=True`, which fills holes; a validator using it would report an open mesh
as closed.

Reported: triangle and vertex counts, connected components, watertightness,
winding consistency, boundary edges, non-manifold edges, degenerate faces, genus,
volume, bounding box, and (with the `quality` extra) self-intersecting faces.

Volume is reported as `undefined` on an open mesh rather than printing the
meaningless number the divergence theorem yields.

## Privacy

The tool reads geometry, modality and acquisition parameters. It never reads or
writes patient names, IDs, or dates, and the JSON provenance record contains none.
Meshes derived from a scan are still personal health data -- a skull is
recognisable. Treat outputs accordingly.

## Not a medical device

Output is for visualisation, education and 3D printing. It is not validated for
diagnosis, surgical planning, or any clinical use.

## Development

```sh
uv venv --python 3.12
uv pip install -e '.[dev,repair,quality]'
pytest
```

The suite runs on synthetic volumes -- spheres, hollow spheres, clipped spheres --
and on DICOM files it writes itself. No patient data required.

## License

MIT
