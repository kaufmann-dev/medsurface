"""Read and write label-name tables embedded in NIfTI-1 header extensions.

TotalSegmentator and other tools store ``{label: name}`` tables as a Caret-style
XML extension (``<LabelTable><Label Key="5"><![CDATA[liver]]></Label>``). This
module understands that format without extra dependencies. Every other format
has no standard place for names, so callers fall back to a JSON sidecar.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import struct
import tempfile
import xml.etree.ElementTree as ElementTree
from pathlib import Path
from typing import Mapping

NIFTI1_HEADER_BYTES = 348
_EXTENSION_FLAG_BYTES = 4
_MAX_PREFIX_BYTES = 64 * 1024 * 1024
_SLUG_INVALID = re.compile(r"[^a-z0-9]+")


def _is_gzip(path: Path) -> bool:
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _open(path: Path):
    return gzip.open(path, "rb") if _is_gzip(path) else open(path, "rb")


def _endian(header: bytes) -> str | None:
    if len(header) < NIFTI1_HEADER_BYTES:
        return None
    if struct.unpack("<i", header[:4])[0] == NIFTI1_HEADER_BYTES:
        return "<"
    if struct.unpack(">i", header[:4])[0] == NIFTI1_HEADER_BYTES:
        return ">"
    return None


def _is_single_file_nifti1(header: bytes) -> bool:
    return header[344:348] == b"n+1\x00"


def _vox_offset(header: bytes, endian: str) -> int:
    return int(struct.unpack(endian + "f", header[108:112])[0])


def _parse_label_table(payload: bytes) -> dict[int, str]:
    payload = payload.rstrip(b"\x00")
    if b"<" not in payload:
        return {}
    try:
        root = ElementTree.fromstring(payload)
    except ElementTree.ParseError:
        return {}
    names: dict[int, str] = {}
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1].lower() != "label":
            continue
        key = None
        for attribute, value in element.attrib.items():
            if attribute.lower() in ("key", "index", "value"):
                try:
                    key = int(value)
                except ValueError:
                    key = None
                break
        text = (element.text or "").strip()
        if key is not None and key > 0 and text:
            names[key] = text
    return names


def read_label_names(path: str | os.PathLike[str]) -> dict[int, str]:
    """Return ``{label: name}`` from a NIfTI-1 label table, or ``{}``.

    Non-NIfTI inputs, NIfTI files without extensions, and unparseable
    extensions all return an empty mapping instead of raising.
    """
    source = Path(path)
    name = source.name.lower()
    if not (name.endswith(".nii") or name.endswith(".nii.gz")):
        return {}
    try:
        with _open(source) as handle:
            header = handle.read(NIFTI1_HEADER_BYTES)
            endian = _endian(header)
            if endian is None or not _is_single_file_nifti1(header):
                return {}
            flag = handle.read(_EXTENSION_FLAG_BYTES)
            if len(flag) < 1 or flag[0] == 0:
                return {}
            limit = min(_vox_offset(header, endian), _MAX_PREFIX_BYTES)
            remaining = max(0, limit - NIFTI1_HEADER_BYTES - _EXTENSION_FLAG_BYTES)
            extensions = handle.read(remaining)
    except (OSError, EOFError, struct.error):
        return {}

    names: dict[int, str] = {}
    position = 0
    while position + 8 <= len(extensions):
        size, _code = struct.unpack(endian + "ii", extensions[position : position + 8])
        if size < 8 or size % 16 or position + size > len(extensions):
            break
        names.update(_parse_label_table(extensions[position + 8 : position + size]))
        position += size
    return names


def load_names_json(path: str | os.PathLike[str]) -> dict[int, str]:
    """Load a ``{"label": "name"}`` JSON object with positive integer keys."""
    with open(path, encoding="utf-8") as handle:
        raw = json.load(handle)
    if not isinstance(raw, dict):
        raise ValueError("label-name JSON must be an object mapping label IDs to names")
    names: dict[int, str] = {}
    for key, value in raw.items():
        try:
            label = int(key)
        except (TypeError, ValueError):
            raise ValueError("label-name key %r is not an integer label ID" % key) from None
        if label <= 0:
            raise ValueError("label-name keys must be positive label IDs; got %d" % label)
        if not isinstance(value, str) or not value.strip():
            raise ValueError("label %d needs a non-empty string name" % label)
        names[label] = value.strip()
    return names


def slug(name: str) -> str:
    """File-name-safe lowercase slug of a structure name."""
    value = _SLUG_INVALID.sub("-", name.lower()).strip("-")
    return value[:64].rstrip("-") or "label"


def label_table_xml(names: Mapping[int, str]) -> bytes:
    """Serialize names as the Caret label table TotalSegmentator writes."""
    rows = [
        '<Label Key="%d" Red="1" Green="1" Blue="1" Alpha="1"><![CDATA[%s]]></Label>'
        % (label, names[label].replace("]]>", "]]]]><![CDATA[>"))
        for label in sorted(names)
    ]
    document = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<CaretExtension><VolumeInformation Index="0"><LabelTable>'
        + "".join(rows)
        + "</LabelTable></VolumeInformation></CaretExtension>"
    )
    return document.encode("utf-8")


def embed_label_names(path: str | os.PathLike[str], names: Mapping[int, str]) -> bool:
    """Atomically add a label table extension to a single-file NIfTI-1 volume.

    Existing extensions are replaced. Voxel bytes are copied unchanged. Returns
    ``False`` without touching the file when it is not a single-file NIfTI-1
    image or when ``names`` is empty.
    """
    target = Path(path)
    if not names or not target.name.lower().endswith((".nii", ".nii.gz")):
        return False
    compressed = _is_gzip(target)
    with _open(target) as handle:
        data = handle.read()
    header = bytearray(data[:NIFTI1_HEADER_BYTES])
    endian = _endian(bytes(header))
    if endian is None or not _is_single_file_nifti1(bytes(header)):
        return False
    old_offset = _vox_offset(bytes(header), endian)
    if old_offset < NIFTI1_HEADER_BYTES or old_offset > len(data):
        return False

    payload = label_table_xml(names)
    size = 8 + len(payload)
    size += (-size) % 16
    extension = struct.pack(endian + "ii", size, 0) + payload
    extension += b"\x00" * (size - len(extension))
    new_offset = NIFTI1_HEADER_BYTES + _EXTENSION_FLAG_BYTES + size
    header[108:112] = struct.pack(endian + "f", float(new_offset))
    rebuilt = bytes(header) + b"\x01\x00\x00\x00" + extension + data[old_offset:]

    descriptor, temporary = tempfile.mkstemp(
        prefix=".%s." % target.name, dir=str(target.parent)
    )
    try:
        with os.fdopen(descriptor, "wb") as raw:
            if compressed:
                with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as out:
                    out.write(rebuilt)
            else:
                raw.write(rebuilt)
            raw.flush()
            os.fsync(raw.fileno())
        with _open(Path(temporary)) as check:
            reread = check.read()
        if reread[new_offset:] != data[old_offset:]:
            raise OSError("label table embedding changed voxel bytes")
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return True
