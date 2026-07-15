# Default mask smoothing erased thin bone

Fixed: 2026-07-15 02:47:43 CEST (+0200)

Baseline commit: `06e24ae1f491047c77f1223a57320a5c75daf11f`

## Symptom

A default bone conversion produced numerous artificial perforations through the
frontal bone and orbital region that were absent from the earlier surface. The
mesh still passed structural validation because each perforation was a closed,
consistently wound tunnel rather than an open mesh boundary.

On the reported volume, the unified default produced a 365,398 mm³ surface.
Disabling mask smoothing while keeping the same segmentation produced
387,144 mm³, showing that the pre-mesh operation removed approximately 5.6% of
the segmented volume.

## Confirmed root cause

The labelmap terracing fix had unified every command behind one `--smooth-mm`
control. Normal presets consequently applied a 0.8 mm recursive Gaussian to the
thresholded binary mask before marching cubes. Thin cortical connections and
partial-volume regions fell below the 0.5 occupancy isovalue and disappeared.

Surface relaxation and simplification were not the cause. Relaxation retains
mesh connectivity, and simplification rejects candidates that change the
component/hole/Euler signature. The topology-changing perforations were already
present in the marching-cubes input after Gaussian mask smoothing.

## Fix

All normal and labelmap conversion and merge commands now expose the same two
independent stages. `--mask-smooth-mm` controls the topology-changing Gaussian
before meshing, while `--mesh-smooth-iters` controls fixed-force,
intersection-safe, volume-preserving relaxation afterward. Either value accepts
`0` to disable only its own stage; the former `--smooth-mm` spelling was removed
without an alias.

Normal presets default mask smoothing to `0` and restore topology-preserving
mesh relaxation to 60 iterations for bone and auto, 10 for teeth, and 35 for
skin. External labelmap conversion and merge retain the empirically verified
0.8 mm mask smoothing plus 20 mesh-relaxation iterations. Provenance records
both resolved settings independently, and the CLI warns only when the
topology-changing mask stage is enabled.

Regression coverage verifies both options on all four command surfaces,
independent enable/disable behavior, normal and labelmap defaults, stage order,
validation, warnings, and provenance.

## Follow-up: restore final mesh relaxation

Revised: 2026-07-15

Baseline commit: `954db0af47a6029a4ce9506e4ead03f0bb3cd970`

The restored normal mesh-relaxation totals are again split around
simplification: bone and auto use 20 iterations before and 40 after, skin uses
25 before and 10 after, and teeth remains 10 before and 0 after. Keeping all
iterations before simplification allowed quadric reduction to reintroduce
visible facets. External labelmaps remain at 20 before and 0 after because their
default physical mask smoothing already addresses voxel terracing.
