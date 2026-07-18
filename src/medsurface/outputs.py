"""Output-format classification shared by CLI and processing layers."""

from __future__ import annotations

import os
from enum import Enum

from .defaults import SUPPORTED_MESH_EXTENSIONS, SUPPORTED_NIFTI_EXTENSIONS


class OutputKind(str, Enum):
    """The two result types produced by merge workflows."""

    MESH = "mesh"
    NIFTI = "nifti"


def extension(path: str | os.PathLike[str]) -> str:
    """Return a normalized extension while preserving compound NIfTI suffixes."""
    normalized = os.fspath(path).casefold()
    if normalized.endswith(".nii.gz"):
        return ".nii.gz"
    return os.path.splitext(normalized)[1]


def merge_output_kind(path: str | os.PathLike[str]) -> OutputKind:
    """Classify a merge destination or reject its unsupported extension."""
    ext = extension(path)
    if ext in SUPPORTED_MESH_EXTENSIONS:
        return OutputKind.MESH
    if ext in SUPPORTED_NIFTI_EXTENSIONS:
        return OutputKind.NIFTI
    raise ValueError(
        "unsupported output extension %r; supported: %s"
        % (
            ext,
            ", ".join((*SUPPORTED_MESH_EXTENSIONS, *SUPPORTED_NIFTI_EXTENSIONS)),
        )
    )


def validate_mesh_output(path: str | os.PathLike[str]) -> None:
    """Require a mesh destination and report compound suffixes accurately."""
    ext = extension(path)
    if ext not in SUPPORTED_MESH_EXTENSIONS:
        raise ValueError(
            "unsupported mesh extension %r; supported: %s"
            % (ext, ", ".join(SUPPORTED_MESH_EXTENSIONS))
        )
