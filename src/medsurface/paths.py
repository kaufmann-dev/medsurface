"""Filesystem identity and output-safety helpers."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable


def same_file(left: str | os.PathLike[str], right: str | os.PathLike[str]) -> bool:
    """Compare paths through lexical aliases, symlinks, and hard links."""
    left_path = os.path.realpath(os.path.abspath(left))
    right_path = os.path.realpath(os.path.abspath(right))
    if left_path == right_path:
        return True
    try:
        return os.path.samefile(left_path, right_path)
    except OSError:
        return False


def protect_outputs(
    input_paths: Iterable[Path],
    primary_output: Path,
    json_report: Path | None,
) -> None:
    """Reject output aliases before processing can modify medical image data."""
    if json_report is not None and same_file(primary_output, json_report):
        raise ValueError("output and JSON report must be different files")

    for input_path in input_paths:
        if same_file(primary_output, input_path):
            raise ValueError("output must not overwrite input file %s" % input_path)
        if json_report is not None and same_file(json_report, input_path):
            raise ValueError("JSON report must not overwrite input file %s" % input_path)
