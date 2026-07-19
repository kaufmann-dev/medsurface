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
