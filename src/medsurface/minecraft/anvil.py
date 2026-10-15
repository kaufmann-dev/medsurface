"""Java Edition Anvil region files, chunk sections, and packed long arrays."""

from __future__ import annotations

import gzip
import os
import struct
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import numpy as np

from . import nbt

MINECRAFT_VERSION = "26.2"
MINECRAFT_DATA_VERSION = 4903
SECTOR_BYTES = 4096
AIR_BLOCKS = frozenset({"minecraft:air", "minecraft:cave_air", "minecraft:void_air"})


@dataclass(frozen=True)
class DimensionBounds:
    """Vertical extent of the generated Overworld dimension."""

    min_y: int
    height: int

    @property
    def max_y(self) -> int:
        return self.min_y + self.height - 1


def region_directory(world: Path) -> Path:
    return world / "dimensions" / "minecraft" / "overworld" / "region"


def region_path(world: Path, region_x: int, region_z: int) -> Path:
    return region_directory(world) / ("r.%d.%d.mca" % (region_x, region_z))


def _bits(count: int) -> int:
    """Bits needed to store values ``0 .. count - 1``."""
    return max(count - 1, 0).bit_length()


def pack(values: np.ndarray, bits: int) -> np.ndarray:
    """Pack unsigned values into longs without spanning, lowest bits first."""
    per_long = 64 // bits
    count = len(values)
    longs = -(-count // per_long)
    padded = np.zeros(longs * per_long, dtype=np.uint64)
    padded[:count] = values
    shifts = np.arange(per_long, dtype=np.uint64) * np.uint64(bits)
    packed = np.bitwise_or.reduce(padded.reshape(longs, per_long) << shifts, axis=1)
    return packed.view(np.int64)


def unpack(longs: np.ndarray, bits: int, count: int) -> np.ndarray:
    """Inverse of :func:`pack`."""
    per_long = 64 // bits
    words = np.asarray(longs, dtype=np.int64).view(np.uint64)
    if len(words) * per_long < count:
        raise ValueError("packed block data is too short")
    shifts = np.arange(per_long, dtype=np.uint64) * np.uint64(bits)
    mask = np.uint64((1 << bits) - 1)
    values = (words[:, None] >> shifts) & mask
    return values.reshape(-1)[:count].astype(np.int64)


def _section(section_y: int, names: Sequence[str], ids: np.ndarray) -> dict[str, Any]:
    present, inverse = np.unique(ids, return_inverse=True)
    if present[0] != 0:
        # Like Minecraft, keep air as the first palette entry.
        present = np.concatenate(([0], present))
        inverse = inverse + 1
    block_states: dict[str, Any] = {
        "palette": nbt.List(
            nbt.TAG_COMPOUND, [{"Name": names[int(value)]} for value in present]
        )
    }
    if len(present) > 1:
        bits = max(4, _bits(len(present)))
        block_states["data"] = nbt.LongArray(pack(inverse.reshape(-1), bits))
    return {
        "Y": nbt.Byte(section_y),
        "block_states": block_states,
        "biomes": {"palette": nbt.List(nbt.TAG_STRING, ["minecraft:the_void"])},
    }


def build_chunk(
    chunk_x: int,
    chunk_z: int,
    column: np.ndarray,
    bounds: DimensionBounds,
    names: Sequence[str],
) -> dict[str, Any]:
    """Chunk NBT for ``column``: global block IDs ``[height, 16 z, 16 x]`` (0 = air)."""
    if column.shape != (bounds.height, 16, 16):
        raise ValueError("chunk column has shape %r" % (column.shape,))
    min_section = bounds.min_y // 16
    max_section = bounds.max_y // 16
    sections = []
    for section_y in range(min_section, max_section + 1):
        start = section_y * 16 - bounds.min_y
        sections.append(_section(section_y, names, column[start : start + 16].reshape(-1)))

    filled = column != 0
    top = bounds.height - np.argmax(filled[::-1], axis=0)
    heights = np.where(filled.any(axis=0), top, 0).reshape(-1)
    packed_heights = pack(heights, _bits(bounds.height + 1))
    heightmaps = {
        key: nbt.LongArray(packed_heights)
        for key in (
            "MOTION_BLOCKING",
            "MOTION_BLOCKING_NO_LEAVES",
            "OCEAN_FLOOR",
            "WORLD_SURFACE",
        )
    }
    return {
        "DataVersion": nbt.Int(MINECRAFT_DATA_VERSION),
        "xPos": nbt.Int(chunk_x),
        "yPos": nbt.Int(min_section),
        "zPos": nbt.Int(chunk_z),
        "Status": "minecraft:full",
        "LastUpdate": nbt.Long(0),
        "InhabitedTime": nbt.Long(0),
        "isLightOn": nbt.Byte(0),
        "sections": nbt.List(nbt.TAG_COMPOUND, sections),
        "Heightmaps": heightmaps,
        "block_entities": nbt.List(nbt.TAG_END),
        "block_ticks": nbt.List(nbt.TAG_END),
        "fluid_ticks": nbt.List(nbt.TAG_END),
        "PostProcessing": nbt.List(
            nbt.TAG_LIST,
            [nbt.List(nbt.TAG_END) for _ in range(min_section, max_section + 1)],
        ),
        "structures": {"starts": {}, "References": {}},
    }


def write_region(
    path: Path, chunks: Iterable[tuple[tuple[int, int], dict[str, Any]]]
) -> None:
    """Write zlib-compressed chunks keyed by absolute ``(chunk_x, chunk_z)``.

    ``chunks`` may be a generator so only one chunk's NBT is held at a time.
    """
    header = bytearray(2 * SECTOR_BYTES)
    body = bytearray()
    for (chunk_x, chunk_z), chunk in chunks:
        compressed = zlib.compress(nbt.encode(chunk))
        payload = struct.pack(">IB", len(compressed) + 1, 2) + compressed
        sectors = -(-len(payload) // SECTOR_BYTES)
        if sectors > 255:
            raise ValueError(
                "chunk (%d, %d) needs %d region sectors; at most 255 are supported"
                % (chunk_x, chunk_z, sectors)
            )
        offset = 2 + len(body) // SECTOR_BYTES
        index = (chunk_x % 32) + (chunk_z % 32) * 32
        header[4 * index : 4 * index + 4] = offset.to_bytes(3, "big") + bytes((sectors,))
        body += payload + bytes(sectors * SECTOR_BYTES - len(payload))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(bytes(header) + bytes(body))


def write_chunks(
    world: Path,
    positions: Iterable[tuple[int, int]],
    make_chunk: Callable[[int, int], dict[str, Any]],
) -> None:
    """Build and write every chunk in ``positions``, one region file at a time."""
    regions: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for chunk_x, chunk_z in positions:
        regions.setdefault((chunk_x // 32, chunk_z // 32), []).append((chunk_x, chunk_z))
    for (region_x, region_z), members in sorted(regions.items()):
        write_region(
            region_path(world, region_x, region_z),
            (((x, z), make_chunk(x, z)) for x, z in sorted(members, key=lambda p: (p[1], p[0]))),
        )


def _decompress(kind: int, data: bytes, where: str) -> bytes:
    try:
        if kind == 1:
            return gzip.decompress(data)
        if kind == 2:
            return zlib.decompress(data)
        if kind == 3:
            return data
    except (OSError, zlib.error) as exc:
        raise ValueError("%s is corrupt: %s" % (where, exc)) from None
    if kind == 4:
        raise ValueError(
            "%s uses LZ4 compression; set region-file-compression=deflate and "
            "re-save the world" % where
        )
    raise ValueError("%s uses unknown compression type %d" % (where, kind))


class RegionReader:
    """Read chunks from a world, loading each region file once."""

    def __init__(self, world: Path) -> None:
        self.world = world
        self._regions: dict[tuple[int, int], bytes | None] = {}

    def _region(self, region_x: int, region_z: int) -> bytes | None:
        key = (region_x, region_z)
        if key not in self._regions:
            path = region_path(self.world, region_x, region_z)
            self._regions[key] = path.read_bytes() if path.is_file() else None
        return self._regions[key]

    def chunk(self, chunk_x: int, chunk_z: int) -> dict[str, Any] | None:
        data = self._region(chunk_x // 32, chunk_z // 32)
        if data is None:
            return None
        where = "chunk (%d, %d)" % (chunk_x, chunk_z)
        if len(data) < 2 * SECTOR_BYTES:
            raise ValueError("region file for %s has a truncated header" % where)
        index = (chunk_x % 32) + (chunk_z % 32) * 32
        offset = int.from_bytes(data[4 * index : 4 * index + 3], "big")
        sectors = data[4 * index + 3]
        if offset == 0 and sectors == 0:
            return None
        start = offset * SECTOR_BYTES
        if start + 5 > len(data):
            raise ValueError("%s points outside its region file" % where)
        length, kind = struct.unpack(">IB", data[start : start + 5])
        if kind & 0x80:
            external = region_directory(self.world) / ("c.%d.%d.mcc" % (chunk_x, chunk_z))
            try:
                payload = external.read_bytes()
            except OSError as exc:
                raise ValueError("cannot read external %s: %s" % (where, exc)) from None
            kind &= 0x7F
        else:
            if length < 1 or start + 4 + length > len(data):
                raise ValueError("%s is truncated" % where)
            payload = data[start + 5 : start + 4 + length]
        try:
            return nbt.decode(_decompress(kind, payload, where))
        except (ValueError, struct.error, UnicodeDecodeError) as exc:
            raise ValueError("%s has invalid NBT: %s" % (where, exc)) from None


def chunk_block_names(
    chunk: dict[str, Any], y_min: int, y_max: int
) -> tuple[list[str], np.ndarray]:
    """Block names in ``y_min..y_max`` as a name table and ``[y, z, x]`` indices.

    Only palette ``Name`` values are read; block-state properties such as the
    ``facing`` of a placed glazed terracotta are ignored. Missing sections are air.
    """
    version = chunk.get("DataVersion")
    if version != MINECRAFT_DATA_VERSION:
        raise ValueError(
            "chunk (%s, %s) has data version %s; expected %d (Java %s)"
            % (
                chunk.get("xPos"),
                chunk.get("zPos"),
                version,
                MINECRAFT_DATA_VERSION,
                MINECRAFT_VERSION,
            )
        )
    names = ["minecraft:air"]
    lookup = {"minecraft:air": 0}
    ids = np.zeros((y_max - y_min + 1, 16, 16), dtype=np.int64)
    for section in chunk.get("sections", []):
        base = int(section.get("Y", 0)) * 16
        low = max(base, y_min)
        high = min(base + 15, y_max)
        states = section.get("block_states")
        if low > high or not states or not states.get("palette"):
            continue
        palette = states["palette"]
        local = []
        for entry in palette:
            name = str(entry.get("Name", "minecraft:air"))
            if name not in lookup:
                lookup[name] = len(names)
                names.append(name)
            local.append(lookup[name])
        if len(palette) == 1:
            values = np.zeros(4096, dtype=np.int64)
        else:
            values = unpack(states.get("data", np.zeros(0)), max(4, _bits(len(palette))), 4096)
            if values.max() >= len(palette):
                raise ValueError("chunk section %d indexes beyond its palette" % (base // 16))
        mapped = np.asarray(local, dtype=np.int64)[values].reshape(16, 16, 16)
        ids[low - y_min : high - y_min + 1] = mapped[low - base : high - base + 1]
    return names, ids


def write_session_lock(world: Path, now_ms: int) -> None:
    with open(os.path.join(world, "session.lock"), "wb") as handle:
        handle.write(struct.pack(">q", now_ms))
