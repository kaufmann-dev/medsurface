# dicom-surface

Turn a DICOM series into a **watertight** 3D surface mesh.

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

Every conversion is validated before it is handed back. If the mesh is not
watertight, the tool says so rather than letting you find out in a slicer.

## Install

```sh
pip install dicom-surface                 # core
pip install 'dicom-surface[repair]'       # + repair of meshes from elsewhere
```

Python 3.10+. No 3D Slicer.

Desktop Linux already has what you need. On a minimal image — a container, a CI
runner, a headless server — install the OpenGL libraries that pymeshlab's plugins
link against, or decimation fails with pymeshlab's `Filter does not exists.`
(the typo is upstream's, quoted verbatim so it is greppable):

```sh
apt install libgl1 libopengl0
```

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

### Presets

| preset | modality | threshold | triangles | for |
|---|---|---|---|---|
| `bone` | CT | 300 HU | 600k | general bone. Denoises without erasing teeth or sutures. |
| `bone-detail` | CT | 300 HU | all | maximum fidelity, very large files |
| `teeth` | CT | 1200 HU | 300k | enamel and dense dentin |
| `skin` | CT | -300 HU | 400k | outer skin surface |
| `auto` | any | Otsu | 600k | MR, CBCT, or any uncalibrated intensity |

Every preset parameter is overridable:

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

## Why the output is watertight

These are the ways a DICOM directory silently produces a plausible-looking but
wrong mesh. This tool handles each, and the test suite has a regression test for
every one.

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

Instances are grouped by `(SeriesInstanceUID, orientation)`, orientations matched
within 2° so float jitter in `ImageOrientationPatient` does not shatter a stack.
Split UIDs get dotted idents (`1021.1`, `1021.2`). Geometry is validated: a stack
whose slice spacing is implausible, or whose spacing varies by more than half its
own median, is rejected with a reason rather than resampled onto a regular grid it
does not fit. Directory and file names are sorted, so discovery does not depend on
filesystem iteration order.

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

**5. VTK's decimators tear thin walls.** A skull's orbital walls are one voxel
thick; edge collapses weld their opposite faces together. Measured on a
4M-triangle skull, `vtkQuadricDecimation` leaks defects at *every* reduction --
6 boundary edges at 50%, 90 boundary and 207 non-manifold at 85% -- so no backoff
converges. `vtkDecimatePro` with `PreserveTopologyOn` both misses the target
(898k triangles when asked for 600k) and still emits non-manifold edges.

MeshLab's quadric edge collapse with `preservetopology=True` is used instead. It
hits the target exactly and holds genus and volume constant: 1225 and ~385,200
mm³ at 600k *and* at 250k triangles. Decimation is still verified afterwards, and
the tool warns if it ever breaks the mesh.

**6. Resampling erases what the closing sealed.** Resampling the mask onto a
coarser isotropic grid *is* topology-safe -- marching cubes always returns a
manifold surface -- and it is tempting, being 25x faster than decimation. But the
morphological closing seals a pore with a membrane one voxel thick, and a coarser
grid blurs that membrane below the occupancy threshold, so the pore reopens.

On a head CT, `--resample-mm 0.6` added 257 tunnels (genus 1225 → 1482), perforated
the brow and orbital walls, and terraced the vault with interpolation contours.
The default is therefore the native grid. `--resample-mm` remains available for
speed and memory, and warns when it is coarser than the native voxel.

**7. A signed distance field is not free.** If you *do* resample, the textbook
choice is a signed distance field -- but ITK's Maurer transform quantises distance
to voxel centres, putting its zero level half a voxel inside the true boundary: a
6% volume loss on a 20 mm sphere. Smoothed fractional occupancy, thresholded at
0.5, keeps the error under 0.5%.

## Merging two scans

Two studies of one body often cover more together than either does alone: a
facial CT that stops mid-vault, a sinus CT that stops at the maxilla. `merge`
registers one onto the other and fuses them.

```sh
dicom-surface merge facial-ct/ sinus-ct/ --series-a 6 --series-b 2 -o skull.stl
```

Two studies rarely share a usable frame of reference. DICOM patient coordinates are
patient-oriented (LPS), but origin, pose and head tilt are set by the acquisition,
so in one real case the meshes were 847 mm apart with a 10.6° difference in tilt.

Registration is global first, then local. The global stage is an exhaustive
**translation** search by FFT cross-correlation: every integer lag is scored at
once, so the translation needs no initial guess. Rotation is not searched there —
it is left entirely to point-to-plane ICP, and so still depends on that method's
basin of attraction. Measured on a phantom, it recovers 3° to 40° to under 1° and
fails beyond roughly 60°. Point-to-plane converged to 0.12 mm on the real pair,
where point-to-point was still descending at 0.93 mm after 60 iterations.

Correspondences are trimmed, because the scans may overlap only partially — but
the trim then *widens* to the overlap actually measured. A fixed 45% trim breaks
small rotations: the correspondences it discards are the ones furthest from the
rotation axis, which is exactly where the rotational signal lives. A 5°
misalignment of two identical volumes converged 4.49° off. Refitting with the
trim widened brings it to 0.44°.

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
then 25 more after decimation. Smooth it any more lightly and it comes out
visibly rougher than either scan it was built from.

### It refuses input it should refuse

Registration cannot fail on its own. FFT always has a peak, ICP always converges
somewhere. Point two unrelated scans at it and you get a confident transform and a
mesh made of two bodies stuck together. So the answer is checked, not trusted.

The geometric gates are calibrated against four true pairs and one deliberate
impostor — a metal bar phantom, handed the skull's own patient identity:

| pair | min overlap | dice | rms |
|---|---|---|---|
| 2024 × 2023, cross-study, cross-kernel | 0.922 | 0.855 | 0.123 mm |
| 2024 × 2023 soft kernel | 0.935 | 0.835 | 0.205 mm |
| 2024 × 2024, 3 mm reconstruction | 0.969 | 0.852 | 0.187 mm |
| 2024 × 2024, sagittal reformat | 0.991 | 0.932 | 0.090 mm |
| **bar phantom, unrelated anatomy** | **0.270** | **0.323** | 0.432 mm |
| gate | 0.60 | 0.55 | *not gated* |

Three lessons are baked into that table.

**The residual is not evidence.** The impostor's 0.43 mm sits inside any plausible
bound, because ICP drives *some* residual down no matter what it is fitting. It is
reported, not gated on. A gate that has never fired is a false sense of security.

**Surface overlap must be symmetric.** Measured only moving→fixed, the bar scored
96.4%: a small dense object buried in a large one finds a neighbouring surface
almost everywhere. Measured the other way it collapses to 26.6%. The gate takes
the weaker direction, each restricted to the other scan's field of view.

**Geometry cannot tell two people apart.** A skull uniformly scaled by 3% — well
inside person-to-person variation, and something rigid registration cannot absorb
— still passes both gates (overlap 0.989, dice 0.675). So `merge` checks
demographics too. This is not bureaucracy; it is the only check that catches a
second body.

`PatientID` equality is *not* that check. Medical record numbers are scoped to the
issuing institution — which is why DICOM carries `IssuerOfPatientID` — and two
real studies of one skull differed in both `PatientID` (7 vs 15 characters) and
`PatientName` formatting while agreeing exactly on birth date and sex. Conflicting
demographics are treated as evidence of different people; a matching ID, or a
matching name plus birth date, as sufficient corroboration that the studies come
from the same person. Neither is proof — identifiers get re-issued, pseudonymised
and mistyped — which is what `--force` is for. Identifiers are compared, never
logged, raised, or written to the provenance record; only field *names* ever appear.

**De-identified data merges without complaint.** Scans stripped of identifiers
land on "unknown", which warns and proceeds. Studies sharing a pseudonymous ID are
accepted outright. Only two differing IDs with nothing to corroborate them are
refused, and that message names de-identification as the likely cause.

`--force` overrides every gate, and says so in the provenance.

## Printing

A watertight mesh is not the same as a printable one. Bone is full of pores that
print as fragile holes, orbital walls are thinner than a nozzle, and specks of
bone smaller than a grain of rice cannot be handled. `--print-profile` fixes those,
on both `convert` and `merge`:

```sh
dicom-surface convert scans/ -o skull.stl --print-profile fdm
dicom-surface merge   a/ b/  -o skull.stl --print-profile resin
```

| profile | closing | thicken | islands | min feature |
|---|---|---|---|---|
| `anatomical` (default) | — | — | — | — |
| `resin` | ≥ 3.2 mm | 0.4 mm | ≥ 100 mm³ | 0.6 mm |
| `fdm` | ≥ 4.8 mm | 1.2 mm | ≥ 200 mm³ | 1.2 mm |

`anatomical` is the identity: every field is zero and every zero is a no-op, so
the default output is unchanged, byte for byte. There is no "printing disabled"
branch that could drift out of sync with the enabled one.

Profiles compose with any preset using `max()`, never assignment — a profile can
only *add* printability, never relax a preset that already closes harder than the
printer needs (`skin` closes 3.2 mm). Triangle budget stays where it belongs, on
`--preset` and `--target-faces`.

### It happens on the mask, not on the mesh

This is the reason it lives here rather than in a tool that post-processes an STL.
Feeding a finished mesh back through voxelisation costs fidelity before any
printability work happens at all. Measured on a 600k-triangle skull, with the
morphology **switched off entirely**:

| | triangles | volume | mean dihedral | RMS dev | max dev |
|---|---|---|---|---|---|
| `convert --preset bone` | 600,000 | 385,798 mm³ | 11.90° | — | — |
| → STL round trip, no morphology | 600,000 | 388,458 mm³ | 10.22° | 0.064 mm | 0.40 mm |

The round trip re-rasterises onto a hard binary grid — discarding the fractional
occupancy this tool already holds — then runs marching cubes, windowed-sinc
smoothing and quadric decimation a *second* time. Doing the closing and the
dilation on the mask skips all of it.

### Thickening grows only what is too thin

Dilating everything is simpler, and it is what a naive print-prep step does, but it
is dimensionally wrong. A cranial vault is 5 mm of solid bone; inflating it moves
the model's outer surface for no benefit. Only the paper-thin structures — orbital
floor, ethmoid, nasal septum — need material.

Selecting the thin set is subtler than it looks. `mask \ opening(mask, r)` is the
textbook answer and it is **wrong here**. An opening is the union of the balls it
contains, so it cannot reach into a sharp convex corner — and every surface of a
voxelised object is locally sharp. That set is a speckle over the *entire* surface,
and dilating it inflates the whole model: a voxelised sphere grows its bounding box
by the full thickening radius. So thin material is instead material that no
sufficiently thick region can reach:

```
core = opening(mask, feature / 2)         # everything thick enough
thin = mask \ dilate(core, thicken_mm)    # beyond the core's reach
out  = mask | dilate(thin, thicken_mm)
```

Measured: a solid sphere and a solid cube both grow by **+0.0%** with their bounding
boxes unchanged; a one-voxel sheet grows from 1 voxel thick to 5. On the real skull,
the `fdm` profile went from **+93.8% volume and +2.70 mm on the bounding box** with
uniform dilation, to **+23.9% and +0.03 mm** by growing only the 1.6% of material
that is genuinely too thin.

`--thicken-mm` is a *radius*, not a kernel extent. Every other millimetre parameter
here is a kernel extent, floored so the realised kernel never exceeds the request,
because overshooting a filter erases anatomy. Undershooting a thickening leaves a
wall too thin to print — so it ceils, and never rounds a positive request down to
nothing. On 0.8 mm slices a 0.4 mm request realises as 0.8 mm, and the tool says so.

### The thin-feature check

Reported, never enforced. Thin material is exactly the material a ball of the
minimum feature size cannot reach — a morphological opening, no local-thickness
estimator required:

```
thin_fraction = |mask \ opening(mask, radius = min_feature / 2)| / |mask|
```

This is the right measure for *reporting* local thickness and the wrong one for
*selecting* what to grow — see above. It carries a small, curvature-dependent floor
from the voxelised surface: a solid sphere of radius 20 mm reports 0.09% thin at a
3 mm feature size, one of radius 10 mm reports 1.6%. Above 5%, `convert` warns. A printer's minimum feature size is a property of the
printer, not of the anatomy, and the right answer to a thin orbital floor is a
decision, not an automatic edit.

### Merging for print

`merge` always registers on **unmodified anatomy**, then re-segments with the
profile for the union. Thickening both scans before registration would inflate
Dice and surface overlap — the very numbers the impostor gates are calibrated on —
so a print profile would quietly make `merge` easier to fool.

Thickening is safe to apply per scan because dilation distributes over union:
`dilate(A ∪ B) == dilate(A) ∪ dilate(B)`. Closing does not, but closing each scan
separately keeps each one's fractional occupancy field intact, which is worth more
than sealing across the seam.

## Smoothing

Smoothing is the one knob you will feel, and it is a real trade: every iteration
buys a smoother surface by moving it further from the data. Windowed-sinc is used
rather than a plain Laplacian, which would shrink a closed surface toward its
centroid on every pass.

Measured on a head CT (600k triangles, deviation against the unsmoothed surface):

| `--post-smooth-iters` | mean dihedral | creased edges | RMS error | max error |
|---|---|---|---|---|
| 0 | 19.7° | 21.4% | 0.045 mm | 0.23 mm |
| 12 | 15.1° | 13.6% | 0.069 mm | 0.28 mm |
| **25** (default) | **11.7°** | **8.8%** | **0.107 mm** | **0.52 mm** |
| 40 | 10.8° | 7.5% | 0.131 mm | 0.64 mm |

25 is the default because it is where the curve turns. The scan's slices are
0.8 mm apart, so its stair-step artefact has an amplitude around 0.4 mm --
*five times* the 0.085 mm mean displacement the filter introduces to remove it.
Smoothing here erases more error than it creates. Past 25, the returns collapse:
40 iterations buy 0.9° of smoothness for another 0.024 mm of RMS.

The cost is not evenly spread. Mean displacement stays far below one voxel, but
the *maximum* lands on thin spicules and sharp crests -- exactly where a CT is
least trustworthy, and exactly where the geometry is real. If you are measuring
rather than looking, use `bone-detail`, which smooths lightly (~0.02 mm mean) and
lets you see the scanner's stair-steps instead of a filter's opinion of them.

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
volume, bounding box, and self-intersecting faces.

Volume is reported as `undefined` on an open mesh rather than printing the
meaningless number the divergence theorem yields.

## Privacy

`list`, `convert`, `validate` and `repair` read geometry, modality and acquisition
parameters only.

`merge` additionally reads a limited set of patient-identifying fields —
`PatientID`, `IssuerOfPatientID`, `PatientName`, `PatientBirthDate`, `PatientSex` —
for the sole purpose of checking whether two studies plausibly belong to the same
person. It has to: geometry alone cannot tell two people apart — see
"[It refuses input it should refuse](#it-refuses-input-it-should-refuse)".

Those values are compared and discarded. They are never logged, raised in an error,
or written to the JSON provenance record — only *field names* ever appear in
messages, as in "`birth_date` differs". Nothing else in the package reads them.

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
