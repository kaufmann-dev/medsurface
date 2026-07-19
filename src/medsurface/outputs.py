"""Output-format classification shared by CLI and processing layers."""

from __future__ import annotations

import os
from dataclasses import dataclass

from .defaults import SUPPORTED_MESH_EXTENSIONS, SUPPORTED_VOLUME_EXTENSIONS


@dataclass(frozen=True)
class VolumeOutput:
    """Suffix-selected storage contract for one atomic volume file."""

    extension: str
    format: str
    compression: str

    @property
    def compressed(self) -> bool:
        return self.compression != "none"


def extension(path: str | os.PathLike[str]) -> str:
    """Return a normalized extension while preserving compound NIfTI suffixes."""
    normalized = os.fspath(path).casefold()
    if normalized.endswith(".nii.gz"):
        return ".nii.gz"
    return os.path.splitext(normalized)[1]


def volume_output(path: str | os.PathLike[str]) -> VolumeOutput:
    """Classify an atomic volume destination or reject its extension."""
    ext = extension(path)
    contracts = {
        ".nii": VolumeOutput(".nii", "NIfTI", "none"),
        ".nii.gz": VolumeOutput(".nii.gz", "NIfTI", "gzip"),
        ".nrrd": VolumeOutput(".nrrd", "NRRD", "gzip"),
        ".mha": VolumeOutput(".mha", "MetaImage", "zlib"),
    }
    try:
        return contracts[ext]
    except KeyError:
        raise ValueError(
            "unsupported volume output extension %r; supported: %s"
            % (ext, ", ".join(SUPPORTED_VOLUME_EXTENSIONS))
        ) from None


def validate_mesh_output(path: str | os.PathLike[str]) -> None:
    """Require a mesh destination and report compound suffixes accurately."""
    ext = extension(path)
    if ext not in SUPPORTED_MESH_EXTENSIONS:
        raise ValueError(
            "unsupported mesh extension %r; supported: %s"
            % (ext, ", ".join(SUPPORTED_MESH_EXTENSIONS))
        )
