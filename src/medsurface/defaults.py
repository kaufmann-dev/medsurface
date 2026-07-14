"""Lightweight command defaults shared by the CLI and processing modules."""

#: Fine enough to preserve sub-millimetre bone while keeping a whole-head fused
#: grid within the merge pipeline's memory limit.
DEFAULT_MERGE_GRID_MM = 0.4

#: Gaussian sigma used to regularize externally segmented labelmap boundaries
#: before marching cubes. This is deliberately a physical distance rather than
#: a mesh-iteration count so its effect does not depend on triangle density.
DEFAULT_LABELMAP_SMOOTH_MM = 0.8

#: Coarse emergency ceiling for source volumes and planned processing grids.
#: This is not a memory guarantee: pixel types and concurrent working images vary.
MAX_VOXELS = 500_000_000

#: Formats supported consistently by conversion, merge, validation, and repair.
SUPPORTED_MESH_EXTENSIONS = (".stl", ".ply", ".obj")
