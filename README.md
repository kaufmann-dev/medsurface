# dicom-surface

Extract a DICOM series as a 3D surface mesh, then report whether the result is
watertight and manifold.

> **Safety and privacy:** this is not a medical device and has not been validated
> for diagnosis or surgical planning. A valid mesh can still be anatomically
> wrong, self-intersecting, under-sampled, or incomplete. Derived meshes are
> personal health data; see [Privacy](#privacy) and [Not a medical
> device](#not-a-medical-device).

```console
$ dicom-surface list ~/scans/head-ct
#      MOD  DESCRIPTION                      SLICES VOXEL mm               PLANE     NOTES
------------------------------------------------------------------------------------------------
1      CT   Topogramm 0,60 sag Tr20 MPR           1 -                      -         localizer
2      CT   GS nativ 3,00 ax Hr40 A3 MPR         92 0.345 x 0.345 x 2.000  axial
6      CT   GS nativ 1,00 ax Hr68 A1 MPR        231 0.315 x 0.315 x 0.800  axial     <- default, sharp kernel Hr68
8      CT   GS nativ 1,00 sag Hr68 A1 MPR       242 0.315 x 0.315 x 0.800  sagittal  sharp kernel Hr68
501    CT   Patientenprotokoll                    1 -                      -         only 1 slice(s)
1021.1 CT   mpr ax 3/2                           64 0.289 x 0.289 x 2.000  axial     orientation 1 of 2 in this UID
1021.2 CT   mpr ax 3/2                            1 1.066 x 1.066 x 0.500  -         only 1 slice(s), orientation 2 of 2

$ dicom-surface convert ~/scans/head-ct -o skull.stl --preset bone
...
quality:
  triangles           600,000
  components          1
  watertight          yes
  boundary edges      0
  non-manifold edges  0
  degenerate faces    0
  genus               1225
  volume              385394 mm3
```

That abbreviated quality block is an illustrative result from the development
head CT, not a promised face count, genus, or volume for other data. Unless
`--no-validate` is given, conversion checks the written file in formats the
validator supports. Self-intersections are checked only when
`--self-intersections` is requested.

## Contents

- [Install](#install)
- [Usage](#usage)
- [DICOM compatibility and limits](#dicom-compatibility-and-limits)
- [How closed output is produced and checked](#how-closed-output-is-produced-and-checked)
- [Merging two scans](#merging-two-scans)
- [Printing](#printing)
- [Smoothing](#smoothing)
- [Sharp kernels](#sharp-kernels)
- [Validation](#validation)
- [Evidence and reproducibility](#evidence-and-reproducibility)
- [Privacy](#privacy)
- [Not a medical device](#not-a-medical-device)
- [Development](#development)

## Install

```sh
pip install dicom-surface                 # core
pip install 'dicom-surface[repair]'       # + repair of meshes from elsewhere
```

Python 3.10+. The CI matrix covers CPython 3.10–3.13 on Ubuntu; other operating
systems are not exercised by this repository. Core dependencies include NumPy,
pydicom, SimpleITK, VTK, Trimesh, SciPy, and PyMeshLab. MeshLib is used only by
the `repair` extra. No 3D Slicer installation is needed.

PyMeshLab is declared as `pymeshlab>=2023.12` with no upper pin. Its Linux wheels
load OpenGL-linked code even in a headless workflow. Many desktop distributions
already include the libraries; the Ubuntu CI image and minimal containers need:

```sh
apt install libgl1 libopengl0
```

Without `libgl1`, importing PyMeshLab can fail with `libGL.so.1: cannot open
shared object file`. Repository CI history shows that without `libopengl0` the
meshing plugin can fail to load and decimation later raises PyMeshLab's current,
greppable `Filter does not exists.` message. These package names are
Debian/Ubuntu-specific; PyMeshLab's [upstream installation
page](https://pymeshlab.readthedocs.io/en/latest/installation.html) does not list
distribution packages or a repository-tested headless substitute.

## Usage

```sh
dicom-surface list     DICOM_DIR                   # what's in this folder?
dicom-surface presets                              # what can I ask for?
dicom-surface convert  DICOM_DIR -o out.stl        # do the thing
dicom-surface merge    DIR_A DIR_B -o out.stl      # fuse two scans of one body
dicom-surface validate mesh.stl                    # is this mesh sound?
dicom-surface repair   mesh.stl -o fixed.stl       # make it watertight
```

Output format follows the extension: `.stl`, `.ply`, `.obj`, `.vtp`.
Coordinates are emitted in millimetres in DICOM patient (**LPS**) space. STL stores
no coordinate-system metadata, so downstream tools will treat them as ordinary XYZ.

STL, PLY, and OBJ are covered by write/read validation tests. VTP writing is
implemented with VTK, but the current Trimesh validator cannot read `.vtp` and
raises `NotImplementedError: file_type 'vtp' not supported`; `convert` therefore
needs `--no-validate` for VTP, followed by validation in a VTP-capable tool.
`dicom-surface validate mesh.vtp` has the same limitation.

### Presets

| preset | modality | threshold | target triangles | for |
|---|---|---|---|---|
| `bone` | CT | 300 HU | 600k | general bone. Denoises without erasing teeth or sutures. |
| `bone-detail` | CT | 300 HU | all | maximum fidelity, very large files |
| `teeth` | CT | 1200 HU | 300k | enamel and dense dentin |
| `skin` | CT | -300 HU | 400k | outer skin surface |
| `auto` | any | Otsu | 600k | MR, CBCT, or any uncalibrated intensity |

Every displayed preset value is a starting point. Explicit CLI flags override
the selected preset:

```sh
dicom-surface convert scans/ -o out.stl \
    --series 6 --threshold 250 --closing-mm 3.2 --target-faces 250000 --smooth-iters 25
```

### Not every scan is a head CT

- **Series selection.** A study holds scouts, several reconstruction kernels, and
  reformats. `list` shows them all; `--series` takes an ident (`6`, `1021.2`), a
  UID, or a description substring. The default pick prefers small voxels, then
  the axial (acquired) plane, and only then slice count. Stacks whose geometry
  does not hold up are rejected with a reason, never silently reconstructed.
- **Modality.** Hounsfield thresholds are meaningful only on CT. Ask for
  `--preset bone` on an MR and the tool refuses rather than emitting a
  confidently wrong mesh. Use `--preset auto` (Otsu) or an explicit
  `--threshold`.
- **Anisotropy.** Filter sizes are given in millimetres and convert to a
  different voxel count on every axis.
- **Multi-part anatomy.** `--all-islands` keeps every bone fragment;
  `--all-components` keeps enclosed cavities (sinuses, marrow space) as separate
  shells.

## DICOM compatibility and limits

The loader is intentionally a classic-series loader, not a general DICOM
normalizer. pydicom reads headers for discovery; SimpleITK's GDCM-backed reader
loads the ordered pixel files. Compatibility is therefore bounded by both this
repository's checks and the codecs in the installed SimpleITK build.

| case | current status |
|---|---|
| Classic single-frame CT/MR-style stacks | Supported path and covered by generated, uncompressed `MONOCHROME2` DICOM tests. At least five instances are required for discovery; loading itself requires at least two. |
| Enhanced multi-frame DICOM | Unsupported. A multi-frame object is counted as one file/frame by discovery and is rejected as too few slices. Per-frame functional groups are not parsed. |
| Compressed transfer syntaxes | Pixel decoding is delegated to SimpleITK/GDCM. No pydicom pixel plugin is needed for header discovery, but no compressed syntax is tested here; support depends on the installed SimpleITK build and codec. |
| `RescaleSlope` / `RescaleIntercept` | Relied on through SimpleITK for CT intensity/HU conversion. The repository has no independent cross-codec calibration test. |
| Pixel padding and `MONOCHROME1` | No explicit padding mask or photometric inversion is implemented and neither case is tested. Inspect intensities and use an explicit threshold; treat these as unverified. |
| Consistent oblique acquisitions | Supported: ordering uses the slice normal and the final affine preserves direction cosines. Synthetic oblique coordinate mapping is tested. |
| Gantry tilt, nonparallel slices, varying orientation | No tilt correction is implemented. Orientations within 2° are grouped; larger changes split the UID. Within-group tilt/nonparallel geometry is not validated, so these inputs are unsupported unless independently checked. |
| Irregular inter-slice distance | The median projected distance becomes the regular grid spacing. A spread under 0.01 mm is called uniform; larger variation warns. A stack is rejected only when spread exceeds `max(0.1 mm, 0.5 × median spacing)`. Accepted irregular stacks are regularized by the reader and may shift geometry. |
| Duplicate positions | Not rejected: zero position differences are omitted when estimating spacing, but duplicate files remain in the load list. Treat as unsupported. |
| Missing geometry tags or varying rows/columns/pixel spacing | Discovery has fallbacks for missing orientation/position, but does not establish that the resulting geometry is correct or consistent. Such inputs may fail in SimpleITK or load incorrectly; they are unsupported. |
| Localizers/scouts | Rejected when `ImageType` contains `LOCALIZER`; very short series are also rejected. An unlabelled scout is not recognized by content. |
| Reformats | A consistent classic reformat can load. Axial stacks are preferred after voxel volume because they are commonly acquired, but `ImageType` is not used to prove acquisition versus reformat. |
| Vendor mosaics | Unsupported and untested; mosaic tiles are not unpacked. |

`FrameOfReferenceUID` is neither read nor compared. Matching UIDs do not skip
registration, and differing UIDs do not change the method. Patient coordinates
are used within each volume, but cross-series alignment is never assumed: `merge`
always performs and gates rigid registration.

## How closed output is produced and checked

These are important ways a DICOM directory can produce a plausible-looking but
wrong mesh. The stated synthetic cases have regression tests; historical
development-scan numbers are identified as such.

**1. Filenames are not slice order.** In the study this was built on, file `1`
held instance 222 and file `231` held instance 226. Sorting by name gives a
scrambled volume that still renders as a convincing blob. Slices are ordered by
`ImagePositionPatient` projected on the slice normal.

**1b. One SeriesInstanceUID may hold several orientations.** Nothing in DICOM
forbids it, and reformat series do it routinely: one real study had 64 axial
frames and a single perpendicular frame sharing a UID. Take the slice normal from
whichever instance the filesystem happens to yield first, and every position is
projected onto the wrong axis -- slice spacing collapses from 2.0 mm to
3.9e-07 mm, the ordering scrambles, and nothing raises. Worse, that fake voxel
size then looks like the finest data in the study and wins automatic selection.

Instances are grouped by `(SeriesInstanceUID, orientation)`, with orientations matched
within 2° so float jitter in `ImageOrientationPatient` does not shatter a stack.
Split UIDs get dotted idents (`1021.1`, `1021.2`). Geometry is validated: a stack
whose slice spacing is implausible, or whose spacing varies by more than half its
own median (subject to the 0.1 mm floor above), is rejected with a reason rather
than silently selected. Lesser irregularity is warned and regularized. Directory
and file names are sorted, so discovery does not depend on filesystem iteration
order.

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
them. A historical development-head measurement counted 2,477 such shells; the
scan and benchmark command are not distributed.

**5. The selected decimator is topology-constrained.** PyMeshLab's
`meshing_decimation_quadric_edge_collapse` is called with `preservetopology=True`,
`preserveboundary=True`, `preservenormal=True`, `planarquadric=True`, and
`qualitythr=0.3`. A historical 4M-triangle development skull exposed defects in
both VTK alternatives: `vtkQuadricDecimation` produced 6 boundary edges at 50%
reduction and 90 boundary plus 207 non-manifold edge uses at 85%; a
`vtkDecimatePro` run requested at 600k stopped at 898k and was still
non-manifold. The original scan, commands, VTK/PyMeshLab versions, and hardware
were not retained, so these are scoped observations, not universal VTK results.

`targetfacenum` is a requested budget, not an API guarantee that every constrained
mesh can reach it. Tests cover exact-within-two faces on a synthetic sphere,
genus preservation on a torus, no-op behavior below budget, and closed synthetic
meshes at several reductions. In the audit environment, PyMeshLab 2025.7.post1
hit 600,000 and 250,000 faces exactly on a procedural 786,432-face torus, kept
genus 1 and zero edge defects, and changed volume by -0.0016% and -0.0154%.
Historical skull runs also landed exactly at 600k and 250k while retaining genus
1225 and roughly 385,200 mm³, but those artifacts are unavailable. Conversion
compares boundary/non-manifold edge counts before and after decimation and warns
if they increase; it does not prove a universal face-count or fidelity guarantee.

**6. Resampling can change topology.** The implementation uses
`vtkFlyingEdges3D` at isovalue 0.5 on a padded binary or fractional-occupancy
field. Padding is intended to close field-of-view boundaries, and the tested
spheres remain watertight, but neither Flying Edges nor resampling preserves the
input mask's component count, cavities, tunnels, genus, or anatomical topology.
VTK documents that Flying Edges reuses the Marching Cubes case table and can emit
degenerate triangles; there is no package-level guarantee that every ambiguous
field produces a valid manifold. The written mesh is therefore validated rather
than assumed valid.

A coarser grid can blur a one-voxel membrane below the occupancy threshold and
reopen a pore. In a historical development-head run, `--resample-mm 0.6` added
257 tunnels (genus 1225 → 1482), perforated the brow and orbital walls, and
terraced the vault. The scan, script, command details, and package versions are
not retained. Reproducible tests make the narrower claim that a 0.4 mm synthetic
slab survives a 0.2 mm grid but not a 1.5 mm grid, while a 4 mm slab survives.

The default is therefore the native grid. `--resample-mm` remains available for
speed and memory, and warns when it is coarser than the finest native spacing.

**7. Fractional occupancy is the chosen resampling field.** A historical
20 mm-sphere experiment reported about 6% volume loss for an ITK Maurer signed
distance field and under 0.5% for smoothed occupancy thresholded at 0.5. Its
script and environment were not retained. The current synthetic regression test
only establishes that resampled sphere volume stays within 5% of the analytic
volume over its tested 0.6, 1.0, and 1.5 mm grids; it does not reproduce the
historical sub-0.5% figure.

## Merging two scans

Two studies of one body often cover more together than either does alone: a
facial CT that stops mid-vault, a sinus CT that stops at the maxilla. `merge`
registers one onto the other and fuses them.

```sh
dicom-surface merge facial-ct/ sinus-ct/ --series-a 6 --series-b 2 -o skull.stl
```

`merge` does not inspect `FrameOfReferenceUID` or use it to skip registration.
DICOM patient coordinates are patient-oriented (LPS), but this implementation
does not assume they align across series or studies. One historical development
pair started 847 mm apart with a 10.6° pose difference; that pair is not included
in the repository.

Registration is global first, then local. Both binary masks are nearest-neighbour
resampled onto separate, axis-aligned 2.0 mm world lattices. The global stage
uses full 3D FFT cross-correlation to score every integer translation lag; it
does not search rotation or sub-lattice translation. Point-to-plane ICP then
refines the transform and supplies rotation.

The reproducible chiral-phantom tests recover rotations of 3°, 8°, 25°, and 40°
to under 1° with translation error under 1 mm. They are four samples of one
synthetic shape, not a guaranteed capture range. Failure beyond roughly 60° and
the real-pair comparison (0.12 mm point-to-plane versus point-to-point still at
0.93 mm after 60 iterations) are historical observations without retained
scripts/data.

Surfaces for ICP are extracted from padded masks, smoothed for 10 iterations at
passband 0.1, transformed into patient space, and assigned vertex normals. All
fixed vertices are retained; moving vertices are sampled without replacement to
at most 60,000 with seed 0. Each ICP pass runs up to 80 iterations and rejects
candidate pairs at or beyond 4 mm; it stops early only when fewer than 200
trimmed pairs remain, not on a transform-delta convergence test. The first pass
keeps the closest 45%. The second keeps `clamp(0.9 × matched_fraction, 0.45,
0.95)`. Historical development measurements reported 4.49° error for a fixed
45% trim on a 5° case and 0.44° after widening; no standalone benchmark artifact
is retained.

The reported metrics mean:

| metric | exact implementation |
|---|---|
| point-to-plane RMS / median | Absolute normal-direction residual to the nearest fixed surface vertex, evaluated after fitting on the closest 45% of sampled moving vertices. It is descriptive and not gated. |
| moving→fixed surface overlap | Fraction of transformed sampled moving vertices that lie inside the fixed image field of view and have a fixed surface vertex within 4 mm. Fewer than 100 eligible points yields 0. |
| fixed→moving surface overlap | The reverse calculation after applying the inverse transform; fixed vertices are sampled to at most 60,000 and compared with the sampled moving surface. |
| minimum/symmetric surface overlap | `min(moving→fixed, fixed→moving)`; the gate requires at least 0.60. This is a symmetric acceptance statistic, not a Hausdorff distance. |
| shared-field Dice | Both masks are linearly resampled and thresholded above 0.5 on a 1.5 mm common world grid. Dice is `2|A∩B|/(|A|+|B|)` only where resampled coverage masks say both scans acquired data; the gate requires 0.55. |
| shared field of view | Count of common-coverage voxels on that 1.5 mm grid times 1.5³ mm³; the gate requires 20,000 mm³. |

The union is taken on the **masks**, not the meshes. Boolean-unioning two shells
leaves a seam ridge wherever they disagree by a fraction of a millimetre; fusing
the solids and running marching cubes once gives a single continuous surface. The
fused grid is **isotropic**: resampling a 0.3 mm sinus CT onto a 0.8 mm facial
CT's grid would throw away the finer scan's resolution and stamp the coarser
scan's slice pitch across the whole vault.

An isotropic grid does not, however, *remove* terracing. Each scan's staircase is
baked into its own mask by its own slice pitch, and resampling a 0.8 mm staircase
onto a 0.4 mm grid merely samples it more finely. So the fused surface is smoothed
with the preset's settings, exactly as a single-scan surface is — 20 iterations,
then 25 more after decimation for `bone`. A historical development comparison
found a lighter merge-specific surface stage visibly rougher; this is rationale
for sharing the preset, not a universal visual guarantee.

### Acceptance gates and identity checks

Registration cannot fail on its own. FFT always has a peak, ICP always converges
somewhere. Point two unrelated scans at it and you get a confident transform and a
mesh made of two bodies stuck together. So the answer is checked, not trusted.

The geometric thresholds were empirically selected after observing the following
small development set. The four positive rows are repeated studies or
reconstructions of one skull, not four independent patients; the negative row is
one metal bar phantom assigned that skull's patient identity. Scanner models,
full protocols, commands, software versions, and raw data are not retained, and
there is no record that thresholds were chosen before these cases were seen.
Treat this as exploratory gate selection, not clinical calibration:

| pair | min overlap | dice | rms |
|---|---|---|---|
| 2024 × 2023, cross-study, cross-kernel | 0.922 | 0.855 | 0.123 mm |
| 2024 × 2023 soft kernel | 0.935 | 0.835 | 0.205 mm |
| 2024 × 2024, 3 mm reconstruction | 0.969 | 0.852 | 0.187 mm |
| 2024 × 2024, sagittal reformat | 0.991 | 0.932 | 0.090 mm |
| **bar phantom, unrelated anatomy** | **0.270** | **0.323** | 0.432 mm |
| gate | 0.60 | 0.55 | *not gated* |

The automated suite additionally registers one synthetic unrelated shell/bar
pair and checks that it fails, but it does not add independent patient negatives.
Three implementation choices came from the development table.

**The residual is not evidence.** The impostor's 0.43 mm sits inside any plausible
bound, because ICP drives *some* residual down no matter what it is fitting. It is
reported, not gated on. A gate that has never fired is a false sense of security.

**Surface overlap must be symmetric.** Measured only moving→fixed, the bar scored
96.4%: a small dense object buried in a large one finds a neighbouring surface
almost everywhere. Measured the other way it collapses to 26.6%. The gate takes
the weaker direction, each restricted to the other scan's field of view.

**Geometry cannot establish identity.** A historical 3%-scaled-skull experiment
passed both gates (overlap 0.989, Dice 0.675); the current test injects those
measured metric values into the gate rather than reproducing the scaling
experiment. `merge` therefore also compares selected patient fields.

Patient identity is decided in this order; "both present" means a non-empty value
in both first instances:

| evidence | result without `--force` |
|---|---|
| Either header cannot be read, either side has no patient fields at all, or nothing reaches a rule below | `unknown`: proceed with a warning that identity cannot be verified. |
| `PatientBirthDate`, `PatientSex`, or normalized `PatientName` are both present and conflict | `different`: refuse, even if IDs match. Missing values do not conflict. |
| IDs match and issuers do not both exist and conflict | `same`: accept. Equal pseudonymous IDs therefore corroborate identity unless issuers conflict. |
| Names match and birth dates either match or are not both present | `same`: accept. Sex is not required, though a present sex conflict was already refused. |
| IDs are both present and differ, and no earlier same-person rule matched | `different`: refuse. Matching birth date/sex alone do not override differing IDs. |
| IDs match but both issuers differ, with no qualifying name match | `unknown`: warn and proceed; the shared number is institution-scoped and insufficient. |

Names are case-folded, carets are changed to spaces, and repeated whitespace is
collapsed; other fields are compared as stripped strings. These heuristics are
corroboration, not proof. Fully de-identified inputs warn rather than proceeding
silently. Differently pseudonymized studies are refused unless another same rule
matches or the operator uses `--force`.

`--force` overrides modality mismatch, a `different` patient verdict, and failed
registration gates, recording `forced: true` in provenance. It does not override
selecting the same series twice, discovery/load errors, an excessive fused grid,
or other invalid input.

## Printing

A watertight mesh is not the same as a printable one. Bone is full of pores that
print as fragile holes, orbital walls are thinner than a nozzle, and specks of
bone smaller than a grain of rice cannot be handled. `--print-profile` fixes those,
on both `convert` and `merge`:

```sh
dicom-surface convert scans/ -o skull.stl --print-profile fdm
dicom-surface merge   a/ b/  -o skull.stl --print-profile resin
```

| profile | closing extent floor | island-volume floor | minimum-feature target |
|---|---|---|---|
| `anatomical` (default) | 0 (no profile change) | 0 (no profile change) | 0 (no thickening) |
| `resin` | 3.2 mm | 100 mm³ | 0.6 mm |
| `fdm` | 4.8 mm | 200 mm³ | 1.2 mm |

There is no separate `thicken` parameter. `min_feature_mm` supplies both the
opening/detection target and twice the selective-dilation radius. The profile
values are engineering defaults: closings were tuned on one development head CT,
and the feature targets are typical machine guidance, not printer validation.

`anatomical` has three zero fields, so it adds no print-profile morphology. A
synthetic test asserts that its output mask equals the profile-free mask. Git
history records manual byte-equality checks for a few development conversions,
but serialized bytes, deterministic ordering, metadata, platforms, and every
output format are not covered by an automated byte test. The supported claim is
mask-stage identity, not universal byte-for-byte output.

Without explicit profile-related flags, closing and island floors compose with
the tissue preset using `max()`; for example, a profile cannot lower a larger
preset value. Explicit `--closing-mm` and `--min-island-mm3` are copied to both
the preset and profile before that `max()`, so an explicit flag wins and can
raise or lower the named defaults. `--min-feature-mm` overrides the profile's
feature target. Other explicit flags override only the preset. Triangle budget
stays on `--preset` / `--target-faces`.

### It happens on the mask, not on the mesh

This is the reason it lives here rather than in a tool that post-processes an STL.
Feeding a finished mesh back through voxelisation can cost fidelity before any
printability work happens. The following is a historical measurement on one
600k-triangle development skull with morphology switched off; the input scan,
round-trip script, raster grid, comparison implementation, package versions, and
hardware are unavailable:

| | triangles | volume | mean dihedral | RMS dev | max dev |
|---|---|---|---|---|---|
| `convert --preset bone` | 600,000 | 385,798 mm³ | 11.90° | — | — |
| → STL round trip, no morphology | 600,000 | 388,458 mm³ | 10.22° | 0.064 mm | 0.40 mm |

That development round trip re-rasterised onto a hard binary grid, then ran
surface extraction, smoothing, and decimation a second time. Applying the
profile to the in-memory mask avoids that extra mesh-to-voxel-to-mesh pass; the
table is not a general error bound.

### Selective mask-space thickening

Dilating everything is simpler, and it is what a naive print-prep step does, but it
is dimensionally wrong. A cranial vault is 5 mm of solid bone; inflating it moves
the model's outer surface for no benefit. Only the paper-thin structures — orbital
floor, ethmoid, nasal septum — need material.

The implementation does not use `mask \ opening(mask, r)` as the growth selector,
because voxelized convex corners make that set speckle across otherwise thick
surfaces. It uses the reach of the opened core:

```
r    = min_feature / 2
core = opening(mask, r)         # everything already thick enough
thin = mask \ dilate(core, r)    # beyond the thick material's reach
out  = mask | dilate(thin, r)
```

This happens after thresholding, island filtering, median/opening from the tissue
preset, and the effective closing. If a feature target is active, the mask is
first moved to the selected isotropic printability grid; selective thickening is
then followed by a second island filter. Mesh extraction, initial smoothing,
component filtering, optional decimation, and post-smoothing happen later.

`T` is therefore an intended mask-space feature target, not a measured final wall
thickness. The opening and dilation use integer voxel radii, and a discrete ball
of radius `r` occupies `(2r+1)` voxel centers across an axis even though the code's
reported feature quantity is `2r × spacing`. A 0.3 mm grid with `T=1.2 mm`, for
example, uses `r=2`; a one-voxel controlled slab grew from 0.3 mm to 1.5 mm, not
exactly 1.2 mm. Fractional resampling, the 0.5 isosurface, smoothing, decimation,
and post-smoothing can move the final surface again. No final-mesh local-thickness
measurement exists, so the profile does not guarantee an exact STL or printed
minimum wall thickness.

There is also a known selector blind spot. Material within `r` of the opened thick
core is excluded from `thin`, even when it is itself a short fragile fin or
bridge. At the FDM grid, audit phantoms containing a 0.3 mm short fin and short
bridge selected 0% of those features and left their controlled thickness at
0.3 mm; the global post-profile thin fraction was only about 1%, below the 5%
warning threshold. Long fins/bridges and the isotropic skull-like plate/rim case
were selected and grown. See the [print-profile investigation](docs/print-profile-investigation.md).

Synthetic tests do establish that a sufficiently large solid sphere and cube are
unchanged while a separated thin sheet grows. Historical measurements on one
development skull reported that the FDM selector acted on 2.0% of mask material,
with +18.6% volume and at most +0.38 mm bounding-box change, versus +93.8% and
+2.70 mm for uniform dilation. Those skull measurements are not reproducible from
the repository and are observations, not bounds.

### The scanner's slice pitch must not become the printer's tolerance

Every kernel here is specified in millimetres and built from whole voxels. For the
*floored* kernels that is harmless — they come out a little gentler than asked.
For the thickening radius it is not. That one **ceils**, because undershooting a
thickening leaves a wall too thin to print, and on an anisotropic grid the axes
round independently. The "ball" becomes an ellipsoid:

| scan | spacing (mm) | native-grid dilation radius per axis for `T=1.2 mm` | aspect |
|---|---|---|---|
| sinus CT | 0.45 × 0.45 × 0.30 | 0.90 × 0.90 × 0.60 | 1.5:1 |
| facial CT | 0.32 × 0.32 × 0.80 | 0.63 × 0.63 × 0.80 | 1.3:1 |
| routine head | 0.50 × 0.50 × 2.00 | 1.00 × 1.00 × **2.00** | 2.0:1 |
| survey | 0.90 × 0.90 × 5.00 | 0.90 × 0.90 × **5.00** | 5.6:1 |

On the last row, native-grid dilation would move each z-facing mask boundary by
5 mm. The implementation instead chooses `grid_mm = (T/2)/r` for integer `r`,
starting with the smallest grid no coarser than the finest input spacing. This
makes the code's reported `2r × grid_mm` quantity equal `T`; it does not make the
discrete support or final mesh thickness exactly `T`. A 300-million-voxel budget
can force a coarser grid, and a warning reports the resulting overshoot or loss of
native resolution.

A grid this fine is only built when `min_feature_mm > 0`. `anatomical` does not
request this resampling path.

It is not free. One historical 0.32 × 0.32 × 0.80 mm development head CT put both
profiles on a 0.30 mm, 234-million-voxel grid; reported end-to-end cost rose from
110 s / 3.2 GB to about 270 s / 6.0 GB. The hardware, OS, commands, and package
versions are unavailable, so these are sizing anecdotes, not performance
expectations. The reported probe/dilation quantities became 0.60 and 1.20 mm per
axis; final wall thickness was not measured.

### The thin-feature check

The diagnostic is a mask-space morphological proxy, not an independent
local-thickness measurement:

```
thin_fraction = |mask \ opening(mask, radius = min_feature / 2)| / |mask|
```

This is separate from the selector above. It can mark voxelized corners on thick
objects, and it can dilute a severe but small local miss into a global fraction.
Above 5%, `convert`/`merge` warn; at or below 5%, they are silent. The warning does
not enforce or certify a final-mesh thickness.

The probe is a ball rounded up to whole voxels, so on a coarse grid it tests for
more material than you asked about. A solid sphere at 3 mm voxels reports 3.9% thin
at a 1.2 mm feature size — an artefact of the probe, not of the geometry. When the
realised probe exceeds twice the request, the tool reports *that* instead of a
number.

There is a limit no grid can lift. Through-plane sampling at 0.8 mm does not make
every 0.6 mm anatomical wall impossible: a wall's orientation, reconstruction
kernel, slice sensitivity profile, and partial-volume signal all matter. It does
mean that thickness below the through-plane pitch is not reliably resolved along
that direction, and isotropic upsampling cannot recover information that was not
sampled. A print profile can add mask material based on the reconstructed signal;
it cannot establish the anatomical thickness or guarantee the manufactured part.

### Merging for print

`merge` always registers on **unmodified anatomy**, then re-segments with the
profile for the union. Thickening both scans before registration would inflate
Dice and surface overlap — the very numbers used by the empirically selected gates —
so a print profile would quietly make `merge` easier to fool.

After registration passes, the current order is:

```text
segment+close+select+thicken(A) ─┐
                                ├─ linear resample fractional occupancy ─ max union
segment+close+select+thicken(B) ─┘
```

Only the final dilation primitive distributes over union. Effective closing,
island removal, printability-grid resampling, opened-core selection, and the
complete selective transform do not. The current order preserves each scan's
post-profile fractional boundary during fusion, but it is not mathematically
equivalent to profiling the fused mask.

Focused post-segmentation FDM phantoms (with tissue-preset median/opening set to
no-op so processing order was isolated) found identical results for identical partial fields,
complementary halves, and one fractional-boundary case; a shell with a one-voxel
registration error differed by 0.1% of anatomical volume. Counterexamples were
material: a fused five-layer plate already thick enough at 1.5 mm became 2.7 mm
when its two thin observations were profiled separately (surface p95 0.66 mm),
and a similarly offset plate behaved the same way. Conversely, fusion-first
closing sealed a 3.6 mm inter-scan gap that per-scan closing left open, adding
9.1% volume and joining two components. This makes processing order a fidelity
trade-off, not a safe distributivity identity. No production behavior changes in
this audit; a future implementation task should compare a post-union selector or
full fusion-first profile on representative skull pairs. See the [investigation
report](docs/print-profile-investigation.md).

## Smoothing

Smoothing is the one knob you will feel, and it is a real trade: every iteration
buys a smoother surface by moving it further from the data. Windowed-sinc is used
rather than a plain Laplacian, which would shrink a closed surface toward its
centroid on every pass.

Historical development-head measurement at a 600k target, comparing each final
row with a surface described in the commit as "unsmoothed":

| `--post-smooth-iters` | mean dihedral | creased edges | RMS error | max error |
|---|---|---|---|---|
| 0 | 19.7° | 21.4% | 0.045 mm | 0.23 mm |
| 12 | 15.1° | 13.6% | 0.069 mm | 0.28 mm |
| **25** (default) | **11.7°** | **8.8%** | **0.107 mm** | **0.52 mm** |
| 40 | 10.8° | 7.5% | 0.131 mm | 0.64 mm |

`0` means zero **post-decimation** iterations, not an unsmoothed pipeline: the
`bone` preset still performed its initial 20 windowed-sinc iterations, component
selection, and decimation before that row. This explains its non-zero deviation.
The exact reference export, distance implementation/interpolation, creased-edge
definition, command, packages, and data were not retained, so the table cannot be
reproduced from the repository and should be read only as a historical parameter
choice on one scan. The current synthetic test establishes only that mean
face-adjacency angle decreases monotonically over 0/10/25/40 iterations on a
sphere while volume stays within 2% of the analytic sphere.

The historical rationale chose 25 where the measured curve flattened: 40 bought
0.9° for another 0.024 mm RMS. The associated 0.4 mm stair-step-amplitude and
0.085 mm mean-displacement comparison was specific to that 0.8 mm-pitch scan; it
is not an error bound or proof that smoothing improves every anatomy. Maximum
movement concentrates at sharp crests and thin spicules. `bone-detail` uses fewer
initial iterations and no decimation target, but its historical ~0.02 mm mean
figure is likewise not reproducible from repository artifacts.

## Sharp kernels

CT reconstruction kernels trade noise against resolution. The detector is a
case-insensitive search using this literal regular expression:

```regex
(?:^|[^0-9])(?:[BHUY]r?|BONE|EDGE|LUNG)\s*_?([6-9]\d)
```

It matches `Hr68`, `B70f`, `U90u`, `BONE70`, and related tokens; tests confirm
`Hr68` and `B70f`, while `Hr40`, `B30f`, bare `BONE`, and `I70f` do not match.
There is no manufacturer lookup or token normalization beyond case-insensitivity
and optional whitespace/underscore, so false positives and false negatives are
possible. A match emits a warning that sharp/edge-enhancing reconstruction can
amplify threshold noise and suggests considering ~300 HU rather than ~200 HU or a
smoother series. That is guidance, not automatic threshold adjustment.

## Validation

`validate` loads through Trimesh with `process=True`, which removes non-finite
vertices and merges positions before topology is calculated. STL stores triangle
facets without an explicit shared-vertex graph, so an unwelded reader commonly
sees boundary edges even when facet coordinates coincide.

The weld is approximate, not exact-coordinate-only: Trimesh rounds positions at
digits derived from `trimesh.constants.tol.merge` (`1e-8` in the audited Trimesh
4.12.2 environment, effectively eight decimal places). Because Trimesh is
declared only as `>=4.0`, this implementation detail can vary by installed
version and can change topology for near-coincident vertices. Default processing
does not remove duplicate faces or repair inconsistent winding; those conditions
remain in the measured mesh, with winding reported separately.

It also never repairs what it measures. `trimesh.split()` builds submeshes with
`repair=True`, which fills holes; a validator using it would report an open mesh
as closed.

Reported by default: triangle and vertex counts, connected components,
watertightness, winding consistency, boundary edges, non-manifold edge uses,
degenerate faces, genus when defined, volume when closed, and bounding box.
`--self-intersections` additionally asks PyMeshLab to select intersecting faces
after duplicate-vertex removal. Import/filter errors are reported as text rather
than treated as zero.

Conversion checks boundary and non-manifold edges around decimation, but does not
check self-intersections during morphology, extraction, smoothing, or decimation.
The final CLI validation also omits that expensive check unless the flag is
given. Watertightness, manifold edge counts, consistent winding, and a defined
volume do not by themselves rule out self-intersections.

Volume is reported as `undefined` on an open mesh rather than printing the
meaningless number the divergence theorem yields.

## Evidence and reproducibility

The repository contains no patient DICOM, benchmark directory, or retained
benchmark script. Git history records useful measurements from one development
head CT and related repeat studies, but generally not the source scan, full
command, random seed, comparison tool, package versions, OS, hardware, or whether
I/O time was included. This README keeps those numbers only when labelled
historical and scoped to that dataset.

Reproducible evidence in the repository is synthetic:

| topic | evidence available now |
|---|---|
| Slice ordering/grouping, spacing, selection, sharp-kernel matching | Generated classic single-frame DICOM headers/pixels in `tests/test_series.py`. |
| Capping, coordinate transforms, resampling survival, smoothing trend, decimation topology, file round trips | Synthetic spheres, slabs, clipped volumes, and a procedural torus in `tests/test_surface.py`. |
| Registration transform recovery and refusal | One generated chiral shell over four tested rotations plus an unrelated synthetic bar in `tests/test_merge.py`; sampling seed is 0 in production. |
| Identity decisions and privacy of explicit patient values | Generated DICOM instances in `tests/test_merge.py`. |
| Printability grid, selectors, morphology, warning thresholds, anatomical mask identity | Synthetic masks in `tests/test_printing.py`. |
| Thin-feature and merge-order concerns | Audit-only temporary phantoms and commands summarized in [docs/print-profile-investigation.md](docs/print-profile-investigation.md); no production algorithm or permanent diagnostic test was added. |

The historical real-data claims include VTK defects, 600k/250k skull results,
257 resampling tunnels, signed-distance volume loss, STL round-trip deviation,
registration metrics/gates, print-profile volume/bounds, runtime/memory, smoothing
deviation, and visually observed sharp-kernel behavior. They should not be read
as general performance, accuracy, safety, or compatibility guarantees. The
official [`vtkFlyingEdges3D`
documentation](https://vtk.org/doc/nightly/html/classvtkFlyingEdges3D.html)
describes its case table and degenerate-triangle warning; SimpleITK's [DICOM
series-reader example](https://simpleitk.readthedocs.io/en/master/link_DicomSeriesReader_docs.html)
describes the underlying ordered-file reader.

## Privacy

`list`, `convert`, `validate` and `repair` read geometry, modality and acquisition
parameters only.

`merge` additionally reads a limited set of patient-identifying fields —
`PatientID`, `IssuerOfPatientID`, `PatientName`, `PatientBirthDate`, `PatientSex` —
for the sole purpose of checking whether two studies plausibly belong to the same
person. It has to: geometry alone cannot tell two people apart — see
the [acceptance-gate discussion](#acceptance-gates-and-identity-checks).

Those five patient-field values are compared and discarded. They are never
logged, included in exceptions, or written to JSON provenance — only field names
can appear, as in "`birth_date` differs". Tests cover that property for the
explicit values they create.

That protection does not make all other metadata anonymous. `SeriesInstanceUID`,
series number/ident, `SeriesDescription`, modality, convolution kernel, slice
count, and spacing can appear in console output or provenance; input/output paths
can also appear in messages. UIDs are linkable identifiers, and free-text series
descriptions or paths can contain identifying text. Remove or pseudonymize those
separately when required.

Meshes derived from a scan are still personal health data: a skull is recognisable.
Treat outputs accordingly.

## Not a medical device

Output is for visualisation, education and 3D printing. It is not validated for
diagnosis, surgical planning, or any clinical use.

## Development

```sh
uv venv --python 3.12
uv pip install -e '.[dev,repair]'
pytest
```

The suite runs on synthetic volumes -- spheres, hollow spheres, clipped spheres --
and on DICOM files it writes itself. No patient data required.

## License

MIT
