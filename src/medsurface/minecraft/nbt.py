"""Minimal Java Edition NBT encoding and decoding.

Writing needs explicit tag types, so numbers are wrapped in :class:`Byte`,
:class:`Int`, and friends. Decoding returns plain Python values: ``dict`` for
compounds, ``list`` for lists, ``int``/``float``/``str`` for scalars, and numpy
arrays for byte, int, and long arrays.
"""

from __future__ import annotations

import gzip
import os
import struct
from typing import Any

import numpy as np

TAG_END = 0
TAG_BYTE = 1
TAG_SHORT = 2
TAG_INT = 3
TAG_LONG = 4
TAG_FLOAT = 5
TAG_DOUBLE = 6
TAG_BYTE_ARRAY = 7
TAG_STRING = 8
TAG_LIST = 9
TAG_COMPOUND = 10
TAG_INT_ARRAY = 11
TAG_LONG_ARRAY = 12


class Byte(int):
    pass


class Short(int):
    pass


class Int(int):
    pass


class Long(int):
    pass


class Float(float):
    pass


class Double(float):
    pass


class IntArray:
    def __init__(self, values: Any) -> None:
        self.values = np.asarray(values, dtype=np.int32)


class LongArray:
    def __init__(self, values: Any) -> None:
        self.values = np.asarray(values, dtype=np.int64)


class List:
    """A typed NBT list; an empty list is written with item type TAG_End."""

    def __init__(self, item_type: int, items: list[Any] | None = None) -> None:
        self.item_type = item_type
        self.items = list(items or [])


_SCALARS: tuple[tuple[type, int, str], ...] = (
    (Byte, TAG_BYTE, ">b"),
    (Short, TAG_SHORT, ">h"),
    (Int, TAG_INT, ">i"),
    (Long, TAG_LONG, ">q"),
    (Float, TAG_FLOAT, ">f"),
    (Double, TAG_DOUBLE, ">d"),
)


def _tag_type(value: Any) -> int:
    for kind, tag, _ in _SCALARS:
        if isinstance(value, kind):
            return tag
    if isinstance(value, str):
        return TAG_STRING
    if isinstance(value, dict):
        return TAG_COMPOUND
    if isinstance(value, List):
        return TAG_LIST
    if isinstance(value, IntArray):
        return TAG_INT_ARRAY
    if isinstance(value, LongArray):
        return TAG_LONG_ARRAY
    raise TypeError("cannot encode %r as NBT; wrap numbers in an NBT type" % (value,))


def _string(value: str) -> bytes:
    encoded = value.encode("utf-8")
    return struct.pack(">H", len(encoded)) + encoded


def _payload(value: Any, tag: int, out: list[bytes]) -> None:
    for _kind, scalar_tag, fmt in _SCALARS:
        if tag == scalar_tag:
            out.append(struct.pack(fmt, value))
            return
    if tag == TAG_STRING:
        out.append(_string(value))
    elif tag == TAG_COMPOUND:
        for key, item in value.items():
            item_tag = _tag_type(item)
            out.append(bytes((item_tag,)) + _string(key))
            _payload(item, item_tag, out)
        out.append(bytes((TAG_END,)))
    elif tag == TAG_LIST:
        item_type = value.item_type if value.items else TAG_END
        out.append(struct.pack(">bi", item_type, len(value.items)))
        for item in value.items:
            _payload(item, item_type, out)
    elif tag == TAG_INT_ARRAY:
        out.append(struct.pack(">i", len(value.values)))
        out.append(value.values.astype(">i4").tobytes())
    elif tag == TAG_LONG_ARRAY:
        out.append(struct.pack(">i", len(value.values)))
        out.append(value.values.astype(">i8").tobytes())
    else:
        raise TypeError("unsupported NBT tag %d" % tag)


def encode(root: dict[str, Any]) -> bytes:
    """Encode ``root`` as an unnamed root compound."""
    out: list[bytes] = [bytes((TAG_COMPOUND,)) + _string("")]
    _payload(root, TAG_COMPOUND, out)
    return b"".join(out)


def write_gzip(path: str | os.PathLike[str], root: dict[str, Any]) -> None:
    """Write a gzip-compressed NBT file such as ``level.dat``."""
    os.makedirs(os.path.dirname(os.fspath(path)) or ".", exist_ok=True)
    with gzip.open(path, "wb") as handle:
        handle.write(encode(root))


class _Reader:
    def __init__(self, data: bytes) -> None:
        self.data = data
        self.position = 0

    def take(self, count: int) -> bytes:
        end = self.position + count
        if count < 0 or end > len(self.data):
            raise ValueError("truncated NBT data")
        chunk = self.data[self.position : end]
        self.position = end
        return chunk

    def unpack(self, fmt: str) -> Any:
        return struct.unpack(fmt, self.take(struct.calcsize(fmt)))[0]

    def string(self) -> str:
        return self.take(self.unpack(">H")).decode("utf-8", errors="replace")

    def array(self, dtype: str) -> np.ndarray:
        length = self.unpack(">i")
        width = np.dtype(dtype).itemsize
        return np.frombuffer(self.take(length * width), dtype=dtype).astype(dtype[1:])

    def payload(self, tag: int) -> Any:
        if tag == TAG_BYTE:
            return self.unpack(">b")
        if tag == TAG_SHORT:
            return self.unpack(">h")
        if tag == TAG_INT:
            return self.unpack(">i")
        if tag == TAG_LONG:
            return self.unpack(">q")
        if tag == TAG_FLOAT:
            return self.unpack(">f")
        if tag == TAG_DOUBLE:
            return self.unpack(">d")
        if tag == TAG_BYTE_ARRAY:
            return self.array(">i1")
        if tag == TAG_STRING:
            return self.string()
        if tag == TAG_LIST:
            item_type = self.unpack(">b")
            length = self.unpack(">i")
            return [self.payload(item_type) for _ in range(max(length, 0))]
        if tag == TAG_COMPOUND:
            values: dict[str, Any] = {}
            while True:
                item_tag = self.unpack(">B")
                if item_tag == TAG_END:
                    return values
                key = self.string()
                values[key] = self.payload(item_tag)
        if tag == TAG_INT_ARRAY:
            return self.array(">i4")
        if tag == TAG_LONG_ARRAY:
            return self.array(">i8")
        raise ValueError("invalid NBT tag type %d" % tag)


def decode(data: bytes) -> dict[str, Any]:
    """Decode an uncompressed NBT document whose root is a compound."""
    reader = _Reader(data)
    tag = reader.unpack(">B")
    if tag != TAG_COMPOUND:
        raise ValueError("NBT root is not a compound")
    reader.string()
    return reader.payload(TAG_COMPOUND)
