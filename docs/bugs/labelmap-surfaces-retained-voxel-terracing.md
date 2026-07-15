# Labelmap surfaces retained voxel terracing

Fixed: 2026-07-15 01:40:23 CEST (+0200)

Baseline commit: `81a57388e2966bb62d9e51d48a9cdb52d2ab01ee`

## Symptom

Merging two craniofacial binary labelmaps produced a valid but visibly terraced
surface. The default MeshLib relaxation ran, but topology-safe simplification
was discarded and the output retained all 4,933,828 raw triangles in a 235 MB
STL. Increasing relaxation iterations produced only diminishing improvement.

## Confirmed root cause

Labelmap commands correctly skipped intensity segmentation and morphology, but
they relied exclusively on one-ring mesh relaxation after marching cubes. One
input had 0.8 mm slice spacing, so its discrete boundary already contained
slice terraces. On the 0.4 mm fusion grid, marching cubes represented that
boundary with millions of short edges; topology-based relaxation therefore had
too little physical reach to remove the artifacts. Replaying the former split
20-plus-40 relaxation sequence produced the same surface, confirming that the
single-pass change was not the cause.

## Fix

Labelmap conversion and merge now expose one physical control, `--smooth-mm`,
instead of mesh iteration and force controls. Its default 0.8 mm Gaussian sigma
is applied to the occupancy field before marching cubes; merge applies it after
the registered masks are unioned. A fixed 20-iteration intersection-safe mesh
relaxation removes residual tessellation noise, and `--smooth-mm 0` disables
both smoothing stages. Commands warn that physical smoothing can round, merge,
or erase features near its scale, and JSON provenance records the chosen sigma.

On the reported masks, the default produced a valid, watertight,
self-intersection-free 298,600-triangle STL of 14.2 MB. The fused voxel volume
changed from approximately 731,201 mm³ without field smoothing to 729,975 mm³
with it, about 0.17%.

## Follow-up: recursive Gaussian input contract

Revised: 2026-07-15 01:45:26 CEST (+0200)

Baseline commit: `45016b0e088f10125dc981636742ae388fadd620`

The physical smoothing path now uses `SmoothingRecursiveGaussian` exclusively.
Labelmap conversion and merge reject any input with fewer than four voxels on
an axis before meshing or registration, even when `--smooth-mm 0` is selected.
The error reports the actual dimensions. At this revision, the ordinary
intensity-volume contract remained at least two voxels per axis.

## Follow-up: unified smoothing contract

Revised: 2026-07-15 02:10:32 CEST (+0200)

Baseline commit: `1e50b8de194149546ac7d8a55b4dc7959ed0e4d2`

Normal conversion and merge now use the same recursive-Gaussian occupancy
smoothing and fixed light mesh relaxation as labelmap commands. The normal
iteration and force controls were replaced by `--smooth-mm`, and both merge
paths smooth only after registration and occupancy fusion. The minimum input
size is now four voxels per axis for every volume type, with no smoothing
fallback; this intentionally supersedes the earlier two-voxel ordinary-volume
contract.
