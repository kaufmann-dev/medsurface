"""Lightweight command defaults shared by the CLI and processing modules."""

#: Fine enough to preserve sub-millimetre bone while keeping a whole-head fused
#: grid within the fusion pipeline's memory limit.
DEFAULT_FUSION_GRID_MM = 0.4

#: Gaussian sigma used to regularize externally segmented labelmap boundaries
#: before marching cubes. This is deliberately a physical distance rather than
#: a mesh-iteration count so its effect does not depend on triangle density.
DEFAULT_LABELMAP_MASK_SMOOTH_MM = 0.8

#: Fixed-force topology-preserving cleanup for external labelmap surfaces.
DEFAULT_LABELMAP_SURFACE_SMOOTH_ITERS = 20

#: Default MeshLib surface-deviation limit for external labelmap surfaces.
DEFAULT_LABELMAP_SIMPLIFY_ERROR_MM = 0.25

#: External labelmaps are already regularized in physical space before meshing,
#: so their default finishing pass does not need additional post-simplification
#: relaxation. Users can enable it explicitly for unusually faceted outputs.
DEFAULT_LABELMAP_POST_SURFACE_SMOOTH_ITERS = 0

#: Recursive Gaussian smoothing requires at least four samples per processed
#: dimension. The shared input contract enforces this from image headers before
#: pixel data are loaded.
MIN_VOLUME_AXIS_VOXELS = 4

#: Mesh relaxation force is fixed so the public surface control needs only one
#: understandable value: its iteration count.
SURFACE_RELAX_FORCE = 0.1

#: Coarse emergency ceiling for source volumes and planned processing grids.
#: This is not a memory guarantee: pixel types and concurrent working images vary.
MAX_VOXELS = 500_000_000

#: Surface formats supported consistently by extraction, validation, and repair.
SUPPORTED_MESH_EXTENSIONS = (".stl", ".ply", ".obj")

#: Atomic single-file volume destinations supported by conversion and fusion.
#: Detached NRRD and MetaImage headers remain input-only because they cannot be
#: published with one atomic filesystem replacement.
SUPPORTED_VOLUME_EXTENSIONS = (".nii", ".nii.gz", ".nrrd", ".mha")

#: Taubin fairing iterations for optional stair-step removal. 600 cleared broad
#: 8-16 mm slice-terrace ripples from a simplified skull dome while a fixed
#: displacement clamp kept the change sub-millimetre.
DEFAULT_DESTEP_ITERS = 600

#: Largest per-vertex displacement, in model mm, that stair-step fairing may
#: introduce. It bounds anatomical change wherever the fairing runs.
DEFAULT_DESTEP_MAX_MM = 1.0

#: Taubin lambda/mu pair. The negative mu step re-inflates what lambda shrinks,
#: so repeated passes remove ripples without steadily shrinking the surface.
DESTEP_LAMBDA = 0.5
DESTEP_MU = -0.53

#: Faces whose normal rotates past this dot product during fairing are treated
#: as folded; their vertices and one surrounding ring return to their input.
DESTEP_UNFOLD_DOT = 0.2

#: Share of vertices at the displacement clamp above which the result warns.
#: Taubin's slight low-frequency gain accumulates with iterations, so many
#: clamped vertices usually mean too many iterations for the triangle size.
DESTEP_CLAMP_WARNING_FRACTION = 0.10

#: Automatic region detection estimates curvature on a fully faired copy of the
#: surface, where slice terraces are gone but anatomy remains. Normals are first
#: averaged over this many neighbour passes to suppress tessellation noise.
DESTEP_AUTO_NORMAL_PASSES = 5

#: Vertices whose faired curvature radius is below this are anatomical detail.
#: On a CT skull this marks teeth, orbital and nasal rims, mandible and
#: zygomatic edges, and processes, while the vault stays clear.
DESTEP_AUTO_DETAIL_RADIUS_MM = 4.0

#: Connected detail smaller than this surface area is an isolated speck. It is
#: faired with its surroundings, because freezing it would pin the ripples
#: around it; on a CT skull these specks covered the vault.
DESTEP_AUTO_SPECK_AREA_MM2 = 20.0

#: Remaining detail is frozen and blends back to full fairing across this
#: geodesic distance.
DESTEP_AUTO_GUARD_FROZEN_MM = 1.0
DESTEP_AUTO_GUARD_FULL_MM = 4.0

#: Detail clusters at least this large, such as the face and dentition, also
#: freeze their surroundings up to the frozen distance and fully fair beyond
#: the full distance. This keeps small smooth patches between facial features
#: unchanged while the open vault is faired, matching a manually banded skull.
DESTEP_AUTO_FEATURE_AREA_MM2 = 100.0
DESTEP_AUTO_FEATURE_FROZEN_MM = 10.0
DESTEP_AUTO_FEATURE_FULL_MM = 20.0
