"""Lightweight command defaults shared by the CLI and processing modules."""

#: Fine enough to preserve sub-millimetre bone while keeping a whole-head fused
#: grid within the merge pipeline's memory limit.
DEFAULT_MERGE_GRID_MM = 0.4

#: Formats supported consistently by conversion, merge, validation, and repair.
SUPPORTED_MESH_EXTENSIONS = (".stl", ".ply", ".obj")
