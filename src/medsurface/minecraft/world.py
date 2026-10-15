"""Round-trip discrete labelmaps through editable Minecraft Java worlds.

By default one labelmap voxel becomes one block and label ``0`` is air. The
chosen vertical index axis maps to Minecraft Y; the other two index axes, in
order, map to X and Z. A sea-lantern wireframe, a spawn platform, and a block
legend sit outside the exported volume. A manifest in ``.medsurface/`` records
the grid geometry and label palette needed to export the edited world.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np
import SimpleITK as sitk

from .. import __version__, labelnames
from .. import labelmap as labelmap_mod
from .. import volume as volume_mod
from ..catalog import FileSource, VolumeCandidate
from ..outputs import volume_output
from . import anvil, blockgrid, nbt
from .anvil import MINECRAFT_DATA_VERSION, MINECRAFT_VERSION, DimensionBounds
from .palette import PaletteEntry, PaletteOverride, assign

Logger = Callable[[str], None]

SCHEMA_VERSION = 1
MANIFEST_DIR = ".medsurface"
MANIFEST_FILE = "minecraft.json"
DEFAULT_PACK_ID = "medsurface"
AXES = ("x", "y", "z")
VANILLA_MAX_Y = 319
VANILLA_HEIGHT = 384
MIN_DIMENSION_Y = -2032
MAX_VERTICAL_VOXELS = 4064
_PACK_ID = re.compile(r"[a-z0-9_.-]{1,64}")
_FRAME_BLOCK = "minecraft:sea_lantern"
_PLATFORM_BLOCK = "minecraft:smooth_quartz"


class UnknownBlocksError(ValueError):
    """Blocks inside the export volume that are not in the label palette."""

    def __init__(self, counts: dict[str, int], samples: list[dict[str, Any]]) -> None:
        total = sum(counts.values())
        listed = ", ".join("%s (%d)" % (name, count) for name, count in sorted(counts.items()))
        where = "; ".join(
            "%s at %s" % (sample["block"], " ".join(str(value) for value in sample["position"]))
            for sample in samples[:3]
        )
        super().__init__(
            "found %d blocks inside the export volume that are not in the label "
            "palette: %s; for example %s. Replace them with palette blocks or air "
            "(see 'medsurface minecraft inspect')" % (total, listed, where)
        )
        self.counts = counts
        self.samples = samples


@dataclass(frozen=True)
class Bounds:
    min: tuple[int, int, int]
    max: tuple[int, int, int]

    def as_dict(self) -> dict[str, list[int]]:
        return {"min": list(self.min), "max": list(self.max)}


@dataclass(frozen=True)
class Placement:
    """Where the labelmap grid sits in the world."""

    vertical: int
    horizontal: tuple[int, int]
    dimension: DimensionBounds
    bounds: Bounds

    @classmethod
    def create(cls, size: tuple[int, ...], vertical: int) -> Placement:
        length = int(size[vertical])
        if length > MAX_VERTICAL_VOXELS:
            fitting = [AXES[axis] for axis in range(3) if size[axis] <= MAX_VERTICAL_VOXELS]
            raise ValueError(
                "labelmap %s axis has %d voxels, exceeding Minecraft's 4,064-block "
                "custom-dimension height limit; fitting axes: %s"
                % (AXES[vertical], length, ", ".join(fitting) or "none")
            )
        height = -(-max(length, VANILLA_HEIGHT) // 16) * 16
        dimension = DimensionBounds(max(VANILLA_MAX_Y + 1 - height, MIN_DIMENSION_Y), height)
        first, second = (axis for axis in range(3) if axis != vertical)
        len_x, len_z = int(size[first]), int(size[second])
        min_x = -(len_x // 2)
        min_z = -(len_z // 2)
        min_y = dimension.min_y + (dimension.height - length) // 2
        return cls(
            vertical=vertical,
            horizontal=(first, second),
            dimension=dimension,
            bounds=Bounds(
                (min_x, min_y, min_z),
                (min_x + len_x - 1, min_y + length - 1, min_z + len_z - 1),
            ),
        )

    def axis_mapping(self) -> dict[str, str]:
        names = {self.vertical: "y", self.horizontal[0]: "x", self.horizontal[1]: "z"}
        return {AXES[axis]: names[axis] for axis in range(3)}

    def chunk_bounds(self) -> dict[str, int]:
        return {
            "min_x": self.bounds.min[0] // 16,
            "max_x": self.bounds.max[0] // 16,
            "min_z": self.bounds.min[2] // 16,
            "max_z": self.bounds.max[2] // 16,
        }

    def required_chunks(self) -> list[tuple[int, int]]:
        bounds = self.chunk_bounds()
        return [
            (chunk_x, chunk_z)
            for chunk_z in range(bounds["min_z"], bounds["max_z"] + 1)
            for chunk_x in range(bounds["min_x"], bounds["max_x"] + 1)
        ]

    def to_world_order(self, voxels_ijk: np.ndarray) -> np.ndarray:
        """``[i, j, k]`` array as Minecraft ``[y, z, x]``."""
        return np.transpose(voxels_ijk, (self.vertical, self.horizontal[1], self.horizontal[0]))

    def from_world_order(self, blocks_yzx: np.ndarray) -> np.ndarray:
        """Minecraft ``[y, z, x]`` array as ``[i, j, k]``."""
        order = {self.vertical: 0, self.horizontal[1]: 1, self.horizontal[0]: 2}
        return np.transpose(blocks_yzx, [order[axis] for axis in range(3)])


@dataclass(frozen=True)
class WorldOptions:
    #: Index axis (``x``, ``y``, ``z``) mapped to Minecraft Y; ``None`` means ``z``.
    vertical_axis: str | None = None
    transform: blockgrid.TransformOptions = field(default_factory=blockgrid.TransformOptions)
    #: ``LevelName`` shown in Minecraft's world list (default: folder name).
    world_name: str | None = None
    #: File name recorded in the manifest instead of the real input name.
    source_name: str | None = None
    palette: Mapping[int, PaletteOverride] = field(default_factory=dict)
    #: Namespace and folder of the custom-height dimension data pack.
    pack_id: str = DEFAULT_PACK_ID


@dataclass(frozen=True)
class Manifest:
    """Validated contents of ``.medsurface/minecraft.json``."""

    raw: dict[str, Any]
    placement: Placement
    palette: list[PaletteEntry]
    size: tuple[int, int, int]
    spacing: tuple[float, float, float]
    origin: tuple[float, float, float]
    direction: tuple[float, ...]
    dtype: str


@dataclass
class WorldResult:
    output_path: str
    size: tuple[int, int, int]
    voxel_count: int
    labels: list[PaletteEntry]
    dimension: DimensionBounds
    volume_bounds: Bounds
    transform: blockgrid.TransformRecord | None
    seconds: float

    def as_dict(self) -> dict[str, Any]:
        return {
            "output": self.output_path,
            "minecraft_version": MINECRAFT_VERSION,
            "size": list(self.size),
            "voxel_count": self.voxel_count,
            "nonzero_labels": len(self.labels),
            "palette": [entry.as_dict() for entry in self.labels],
            "dimension_bounds": _dimension_dict(self.dimension),
            "volume_bounds": self.volume_bounds.as_dict(),
            "transform": self.transform.as_dict() if self.transform else None,
            "seconds": self.seconds,
        }


@dataclass
class VolumeResult:
    output_path: str
    format: str
    size: tuple[int, int, int]
    pixel_type: str
    label_voxels: dict[int, int]
    transform: dict[str, Any] | None
    seconds: float

    def as_dict(self) -> dict[str, Any]:
        return {
            "output": self.output_path,
            "format": self.format,
            "size": list(self.size),
            "pixel_type": self.pixel_type,
            "label_voxels": {str(label): count for label, count in self.label_voxels.items()},
            "transform": self.transform,
            "seconds": self.seconds,
        }


@dataclass
class ValidationResult:
    required_chunks: int
    checked_voxels: int
    label_voxels: dict[int, int]

    def as_dict(self) -> dict[str, Any]:
        return {
            "valid": True,
            "minecraft_version": MINECRAFT_VERSION,
            "required_chunks": self.required_chunks,
            "checked_voxels": self.checked_voxels,
            "label_voxels": {str(label): count for label, count in self.label_voxels.items()},
        }


def _dimension_dict(dimension: DimensionBounds) -> dict[str, int]:
    return {"min_y": dimension.min_y, "height": dimension.height}


def _say(log: Logger | None, message: str) -> None:
    if log is not None:
        log(message)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_pack_id(pack_id: str) -> None:
    if not _PACK_ID.fullmatch(pack_id) or pack_id in ("minecraft", "vanilla"):
        raise ValueError(
            "--pack-id must use 1-64 of a-z, 0-9, '_', '-' or '.' and must not be "
            "minecraft or vanilla"
        )


def _vertical_index(options: WorldOptions) -> int:
    anatomical = options.transform.orientation == "anatomical"
    if anatomical and options.vertical_axis not in (None, "z"):
        raise ValueError(
            "--orient anatomical always places superior on Minecraft Y; omit --vertical-axis"
        )
    axis = options.vertical_axis or "z"
    if axis not in AXES:
        raise ValueError("vertical axis must be x, y, or z")
    return AXES.index(axis)


def _load_labelmap(
    candidate: VolumeCandidate, *, allow_large_volume: bool
) -> tuple[sitk.Image, dict[int, str]]:
    volume = labelmap_mod._load_volume(candidate, allow_large_volume=allow_large_volume)
    labelmap_mod._validate_labelmap(volume.image)
    return volume.image, labelmap_mod.embedded_label_names(candidate)


def _label_counts(labels: np.ndarray) -> dict[int, int]:
    values, counts = np.unique(labels, return_counts=True)
    return {int(value): int(count) for value, count in zip(values, counts)}


# ---------------------------------------------------------------- guides
def _guides(
    placement: Placement, palette: list[PaletteEntry]
) -> tuple[dict[tuple[int, int, int], str], tuple[int, int, int]]:
    """Wireframe, platform, and legend blocks, all outside the volume bounds."""
    bounds = placement.bounds
    dimension = placement.dimension
    guides: dict[tuple[int, int, int], str] = {}
    frame_x = (bounds.min[0] - 2, bounds.max[0] + 2)
    frame_z = (bounds.min[2] - 2, bounds.max[2] + 2)
    for y in range(bounds.min[1], bounds.max[1] + 1):
        for x in frame_x:
            for z in frame_z:
                guides[(x, y, z)] = _FRAME_BLOCK
    for y in (bounds.min[1], bounds.max[1]):
        for x in range(frame_x[0], frame_x[1] + 1):
            for z in frame_z:
                guides[(x, y, z)] = _FRAME_BLOCK
        for z in range(frame_z[0], frame_z[1] + 1):
            for x in frame_x:
                guides[(x, y, z)] = _FRAME_BLOCK

    origin_x = bounds.max[0] + 10
    origin_z = bounds.min[2]
    platform_y = min(max(bounds.min[1] + 5, dimension.min_y + 1), dimension.max_y - 3)
    rows = -(-max(len(palette), 1) // 16)
    for x in range(origin_x - 3, origin_x + 19):
        for z in range(origin_z - 5, origin_z + rows + 3):
            guides[(x, platform_y, z)] = _PLATFORM_BLOCK
    for index, entry in enumerate(palette):
        position = (origin_x + index % 16, platform_y + 1, origin_z + index // 16)
        entry.legend_position = position
        guides[position] = entry.block
    spawn = (origin_x, platform_y + 1, origin_z - 3)
    return guides, spawn


# ---------------------------------------------------------- world files
def _data_root(data: dict[str, Any]) -> dict[str, Any]:
    return {"DataVersion": nbt.Int(MINECRAFT_DATA_VERSION), "data": data}


def _dimension_generator(dimension_type: str, biome: str) -> dict[str, Any]:
    return {
        "type": dimension_type,
        "generator": {
            "type": "minecraft:flat",
            "settings": {
                "biome": biome,
                "features": nbt.Byte(0),
                "lakes": nbt.Byte(0),
                "layers": nbt.List(
                    nbt.TAG_COMPOUND, [{"block": "minecraft:air", "height": nbt.Int(1)}]
                ),
                "structure_overrides": nbt.List(nbt.TAG_END),
            },
        },
    }


def _write_dimension_pack(world: Path, dimension: DimensionBounds, pack_id: str) -> None:
    description = (
        "medsurface custom-height Overworld"
        if pack_id == DEFAULT_PACK_ID
        else "Custom-height Overworld"
    )
    pack = {"pack": {"description": description, "min_format": [107, 1], "max_format": 107}}
    dimension_type = {
        "ambient_light": 0.0,
        "attributes": {
            "minecraft:audio/ambient_sounds": {
                "mood": {
                    "block_search_extent": 8,
                    "offset": 2.0,
                    "sound": "minecraft:ambient.cave",
                    "tick_delay": 6000,
                }
            },
            "minecraft:audio/background_music": {
                "creative": {
                    "max_delay": 24000,
                    "min_delay": 12000,
                    "sound": "minecraft:music.creative",
                },
                "default": {
                    "max_delay": 24000,
                    "min_delay": 12000,
                    "sound": "minecraft:music.game",
                },
            },
            "minecraft:gameplay/bed_rule": {
                "can_set_spawn": "always",
                "can_sleep": "when_dark",
                "error_message": {"translate": "block.minecraft.bed.no_sleep"},
            },
            "minecraft:gameplay/nether_portal_spawns_piglin": True,
            "minecraft:gameplay/respawn_anchor_works": False,
            "minecraft:visual/ambient_light_color": "#0a0a0a",
            "minecraft:visual/cloud_color": "#ccffffff",
            "minecraft:visual/cloud_height": 192.33,
            "minecraft:visual/fog_color": "#c0d8ff",
            "minecraft:visual/sky_color": "#78a7ff",
        },
        "coordinate_scale": 1.0,
        "default_clock": "minecraft:overworld",
        "has_ceiling": False,
        "has_ender_dragon_fight": False,
        "has_skylight": True,
        "height": dimension.height,
        "infiniburn": "#minecraft:infiniburn_overworld",
        "logical_height": dimension.height,
        "min_y": dimension.min_y,
        "monster_spawn_block_light_limit": 0,
        "monster_spawn_light_level": {
            "type": "minecraft:uniform",
            "max_inclusive": 7,
            "min_inclusive": 0,
        },
        "timelines": "#minecraft:in_overworld",
    }
    root = world / "datapacks" / pack_id
    directory = root / "data" / pack_id / "dimension_type"
    directory.mkdir(parents=True)
    (root / "pack.mcmeta").write_text(json.dumps(pack, indent=2), encoding="utf-8")
    (directory / "overworld.json").write_text(
        json.dumps(dimension_type, indent=2), encoding="utf-8"
    )


def _write_world_metadata(
    world: Path,
    name: str,
    pack_id: str,
    placement: Placement,
    guides: Mapping[tuple[int, int, int], str],
    spawn: tuple[int, int, int],
) -> None:
    now_ms = int(time.time() * 1000)
    level = {
        "Data": {
            "DataVersion": nbt.Int(MINECRAFT_DATA_VERSION),
            "version": nbt.Int(19133),
            "LevelName": name,
            "LastPlayed": nbt.Long(now_ms),
            "GameType": nbt.Int(1),
            "allowCommands": nbt.Byte(1),
            "initialized": nbt.Byte(1),
            "WasModded": nbt.Byte(0),
            "SpawnX": nbt.Int(spawn[0]),
            "SpawnY": nbt.Int(spawn[1]),
            "SpawnZ": nbt.Int(spawn[2]),
            "SpawnAngle": nbt.Float(0.0),
            # Current Java versions store the world spawn in this compound; the
            # legacy fields above remain for older readers.
            "spawn": {
                "dimension": "minecraft:overworld",
                "pos": nbt.IntArray(spawn),
                "yaw": nbt.Float(0.0),
                "pitch": nbt.Float(0.0),
            },
            "Version": {
                "Id": nbt.Int(MINECRAFT_DATA_VERSION),
                "Name": MINECRAFT_VERSION,
                "Series": "main",
                "Snapshot": nbt.Byte(0),
            },
            "DataPacks": {
                "Enabled": nbt.List(nbt.TAG_STRING, ["vanilla", "file/%s" % pack_id]),
                "Disabled": nbt.List(nbt.TAG_END),
            },
            "enabled_features": nbt.List(nbt.TAG_STRING, ["minecraft:vanilla"]),
            "difficulty_settings": {
                "difficulty": "peaceful",
                "locked": nbt.Byte(1),
                "hardcore": nbt.Byte(0),
            },
            "server_brands": nbt.List(nbt.TAG_END),
        }
    }
    nbt.write_gzip(world / "level.dat", level)
    _write_dimension_pack(world, placement.dimension, pack_id)

    data = world / "data" / "minecraft"
    game_rules = {
        "minecraft:advance_time": nbt.Byte(0),
        "minecraft:advance_weather": nbt.Byte(0),
        "minecraft:spawn_mobs": nbt.Byte(0),
        "minecraft:spawn_patrols": nbt.Byte(0),
        "minecraft:spawn_wandering_traders": nbt.Byte(0),
        "minecraft:spawn_wardens": nbt.Byte(0),
        "minecraft:spawn_phantoms": nbt.Byte(0),
        "minecraft:fire_spread_radius_around_player": nbt.Int(0),
        "minecraft:mob_griefing": nbt.Byte(0),
        "minecraft:random_tick_speed": nbt.Int(0),
        "minecraft:spread_vines": nbt.Byte(0),
        "minecraft:keep_inventory": nbt.Byte(1),
    }
    nbt.write_gzip(data / "game_rules.dat", _data_root(game_rules))
    clocks = {"minecraft:overworld": {"total_ticks": nbt.Long(6000), "paused": nbt.Byte(1)}}
    nbt.write_gzip(data / "world_clocks.dat", _data_root(clocks))
    weather = {
        "clear_weather_time": nbt.Int(1_000_000),
        "rain_time": nbt.Int(0),
        "raining": nbt.Byte(0),
        "thunder_time": nbt.Int(0),
        "thundering": nbt.Byte(0),
    }
    nbt.write_gzip(data / "weather.dat", _data_root(weather))
    world_gen = {
        "bonus_chest": nbt.Byte(0),
        "generate_features": nbt.Byte(0),
        "seed": nbt.Long(0),
        "dimensions": {
            "minecraft:overworld": _dimension_generator(
                "%s:overworld" % pack_id, "minecraft:the_void"
            ),
            "minecraft:the_nether": _dimension_generator(
                "minecraft:the_nether", "minecraft:nether_wastes"
            ),
            "minecraft:the_end": _dimension_generator("minecraft:the_end", "minecraft:the_end"),
        },
    }
    nbt.write_gzip(data / "world_gen_settings.dat", _data_root(world_gen))

    xs = [placement.bounds.min[0], placement.bounds.max[0], *(x for x, _, _ in guides)]
    zs = [placement.bounds.min[2], placement.bounds.max[2], *(z for _, _, z in guides)]
    size = float(max(max(xs) - min(xs) + 129, max(zs) - min(zs) + 129))
    border = {
        "center_x": nbt.Double((min(xs) + max(xs)) / 2.0),
        "center_z": nbt.Double((min(zs) + max(zs)) / 2.0),
        "size": nbt.Double(size),
        "size_lerp_time": nbt.Long(0),
        "size_lerp_target": nbt.Double(size),
        "safe_zone": nbt.Double(5.0),
        "damage_per_block": nbt.Double(0.0),
        "warning_blocks": nbt.Int(5),
        "warning_time": nbt.Int(15),
    }
    nbt.write_gzip(
        world / "dimensions" / "minecraft" / "overworld" / "data" / "world_border.dat",
        _data_root(border),
    )
    anvil.write_session_lock(world, now_ms)


def _write_blocks(
    world: Path,
    placement: Placement,
    labels_ijk: np.ndarray,
    palette: list[PaletteEntry],
    guides: Mapping[tuple[int, int, int], str],
) -> int:
    """Write every volume and guide chunk; return the number of chunks."""
    names = ["minecraft:air"]
    for block in [entry.block for entry in palette] + list(dict.fromkeys(guides.values())):
        if block not in names:
            names.append(block)
    global_id = {name: index for index, name in enumerate(names)}

    ordered = sorted(palette, key=lambda entry: entry.label)
    values = [entry.label for entry in ordered]
    block_ids = np.array([global_id[entry.block] for entry in ordered], dtype=np.uint8)
    if labels_ijk.dtype.kind in "iu" and (not values or values[-1] < 1 << 16):
        lookup = np.zeros(max(values, default=0) + 1, dtype=np.uint8)
        lookup[values] = block_ids
        ids = lookup[labels_ijk]
    else:
        ids = np.zeros(labels_ijk.shape, dtype=np.uint8)
        foreground = labels_ijk != 0
        ids[foreground] = block_ids[np.searchsorted(values, labels_ijk[foreground])]
    volume_yzx = placement.to_world_order(ids)

    bounds = placement.bounds
    positions = np.array(list(guides), dtype=np.int64).reshape(-1, 3)
    low = np.minimum(np.array(bounds.min), positions.min(axis=0))
    high = np.maximum(np.array(bounds.max), positions.max(axis=0))
    blocks = np.zeros((high[1] - low[1] + 1, high[2] - low[2] + 1, high[0] - low[0] + 1), np.uint8)
    blocks[
        bounds.min[1] - low[1] : bounds.max[1] - low[1] + 1,
        bounds.min[2] - low[2] : bounds.max[2] - low[2] + 1,
        bounds.min[0] - low[0] : bounds.max[0] - low[0] + 1,
    ] = volume_yzx
    guide_ids = np.array([global_id[name] for name in guides.values()], dtype=np.uint8)
    blocks[
        positions[:, 1] - low[1], positions[:, 2] - low[2], positions[:, 0] - low[0]
    ] = guide_ids

    chunks = set(placement.required_chunks())
    chunks.update(zip((positions[:, 0] // 16).tolist(), (positions[:, 2] // 16).tolist()))
    dimension = placement.dimension

    def make_chunk(chunk_x: int, chunk_z: int) -> dict[str, Any]:
        column = np.zeros((dimension.height, 16, 16), dtype=np.uint8)
        x0, x1 = max(chunk_x * 16, low[0]), min(chunk_x * 16 + 15, high[0])
        z0, z1 = max(chunk_z * 16, low[2]), min(chunk_z * 16 + 15, high[2])
        if x0 <= x1 and z0 <= z1:
            column[
                low[1] - dimension.min_y : high[1] - dimension.min_y + 1,
                z0 - chunk_z * 16 : z1 - chunk_z * 16 + 1,
                x0 - chunk_x * 16 : x1 - chunk_x * 16 + 1,
            ] = blocks[:, z0 - low[2] : z1 - low[2] + 1, x0 - low[0] : x1 - low[0] + 1]
        return anvil.build_chunk(chunk_x, chunk_z, column, dimension, names)

    anvil.write_chunks(world, chunks, make_chunk)
    return len(chunks)


# ------------------------------------------------------------- manifest
def manifest_path(world: Path) -> Path:
    return world / MANIFEST_DIR / MANIFEST_FILE


def _manifest(
    image: sitk.Image,
    source_name: str,
    source_sha256: str,
    placement: Placement,
    palette: list[PaletteEntry],
    record: blockgrid.TransformRecord | None,
) -> dict[str, Any]:
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "created_by": "medsurface %s" % __version__,
        "minecraft_version": MINECRAFT_VERSION,
        "minecraft_data_version": MINECRAFT_DATA_VERSION,
        "source": {"name": source_name, "sha256": source_sha256},
        "grid": {
            "size": list(image.GetSize()),
            "spacing": list(image.GetSpacing()),
            "origin": list(image.GetOrigin()),
            "direction": list(image.GetDirection()),
            "dtype": sitk.GetArrayViewFromImage(image).dtype.name,
        },
        "vertical_axis": AXES[placement.vertical],
        "axes": placement.axis_mapping(),
        "dimension_bounds": _dimension_dict(placement.dimension),
        "volume_bounds": placement.bounds.as_dict(),
        "required_chunks": placement.chunk_bounds(),
        "palette": [entry.as_dict() for entry in palette],
    }
    if record is not None:
        manifest["transform"] = record.as_dict()
    return manifest


def load_manifest(world: Path) -> Manifest:
    """Read and validate a world's manifest."""
    if not world.is_dir():
        raise ValueError("%s is not a Minecraft world directory" % world)
    path = manifest_path(world)
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValueError(
            "cannot read %s: %s; was this world created by 'medsurface minecraft to-world'?"
            % (path, exc)
        ) from None
    except json.JSONDecodeError as exc:
        raise ValueError("invalid manifest %s: %s" % (path, exc)) from None
    try:
        if raw["schema_version"] != SCHEMA_VERSION:
            raise ValueError(
                "unsupported manifest schema %r; expected %d"
                % (raw["schema_version"], SCHEMA_VERSION)
            )
        if raw["minecraft_data_version"] != MINECRAFT_DATA_VERSION:
            raise ValueError(
                "world targets Minecraft data version %r; expected %d (Java %s)"
                % (raw["minecraft_data_version"], MINECRAFT_DATA_VERSION, MINECRAFT_VERSION)
            )
        grid = raw["grid"]
        size = tuple(int(value) for value in grid["size"])
        if len(size) != 3 or min(size) < 1:
            raise ValueError("manifest grid size is invalid")
        np.dtype(grid["dtype"])
        palette = [
            PaletteEntry(
                label=int(entry["label"]),
                block=str(entry["block"]),
                name=entry.get("name"),
                voxel_count=int(entry["voxel_count"]),
                legend_position=tuple(entry["legend_position"]),
            )
            for entry in raw["palette"]
        ]
        vertical = AXES.index(raw["vertical_axis"])
        placement = Placement.create(size, vertical)
        if (
            raw["volume_bounds"] != placement.bounds.as_dict()
            or raw["dimension_bounds"] != _dimension_dict(placement.dimension)
            or raw["required_chunks"] != placement.chunk_bounds()
        ):
            raise ValueError("manifest axis mapping and volume bounds are inconsistent")
        manifest = Manifest(
            raw=raw,
            placement=placement,
            palette=palette,
            size=size,  # type: ignore[arg-type]
            spacing=tuple(float(value) for value in grid["spacing"]),  # type: ignore[arg-type]
            origin=tuple(float(value) for value in grid["origin"]),  # type: ignore[arg-type]
            direction=tuple(float(value) for value in grid["direction"]),
            dtype=str(grid["dtype"]),
        )
    except (KeyError, TypeError) as exc:
        raise ValueError("invalid manifest %s: missing or malformed field %s" % (path, exc)) from None
    except ValueError as exc:
        raise ValueError("invalid manifest %s: %s" % (path, exc)) from None
    blocks: set[str] = set()
    labels: set[int] = set()
    for entry in palette:
        if entry.block in blocks:
            raise ValueError("manifest maps more than one label to %s" % entry.block)
        if entry.label <= 0 or entry.label in labels:
            raise ValueError("manifest contains invalid or duplicate label %d" % entry.label)
        blocks.add(entry.block)
        labels.add(entry.label)
    return manifest


# ------------------------------------------------------------- decoding
def ensure_unlocked(world: Path) -> None:
    """Refuse to read a world that Minecraft currently holds open."""
    path = world / "session.lock"
    try:
        handle = open(path, "r+b")
    except OSError as exc:
        raise ValueError("cannot open %s: %s" % (path, exc)) from None
    with handle:
        busy = ValueError(
            "the Minecraft world is currently open; close it before validating or exporting"
        )
        if sys.platform == "win32":
            import msvcrt

            try:
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError:
                raise busy from None
            msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl

            # Java's FileChannel.tryLock uses POSIX record locks on Unix.
            try:
                fcntl.lockf(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError:
                raise busy from None
            fcntl.lockf(handle, fcntl.LOCK_UN)


def decode_world(world: Path, manifest: Manifest) -> np.ndarray:
    """Labels of the export volume as ``[i, j, k]``; unknown blocks raise."""
    ensure_unlocked(world)
    if not (world / "level.dat").is_file():
        raise ValueError("world is missing level.dat")
    placement = manifest.placement
    bounds = placement.bounds
    label_of = {entry.block: entry.label for entry in manifest.palette}
    labels_yzx = np.zeros(
        (
            bounds.max[1] - bounds.min[1] + 1,
            bounds.max[2] - bounds.min[2] + 1,
            bounds.max[0] - bounds.min[0] + 1,
        ),
        dtype=np.dtype(manifest.dtype),
    )
    unknown: dict[str, int] = {}
    samples: list[dict[str, Any]] = []
    reader = anvil.RegionReader(world)
    for chunk_x, chunk_z in placement.required_chunks():
        chunk = reader.chunk(chunk_x, chunk_z)
        if chunk is None:
            raise ValueError("required chunk (%d, %d) is missing" % (chunk_x, chunk_z))
        names, ids = anvil.chunk_block_names(chunk, bounds.min[1], bounds.max[1])
        x0, x1 = max(bounds.min[0], chunk_x * 16), min(bounds.max[0], chunk_x * 16 + 15)
        z0, z1 = max(bounds.min[2], chunk_z * 16), min(bounds.max[2], chunk_z * 16 + 15)
        ids = ids[:, z0 - chunk_z * 16 : z1 - chunk_z * 16 + 1, x0 - chunk_x * 16 : x1 - chunk_x * 16 + 1]
        lookup = np.array(
            [0 if name in anvil.AIR_BLOCKS else label_of.get(name, -1) for name in names],
            dtype=np.int64,
        )
        labels = lookup[ids]
        if (labels < 0).any():
            for y, z, x in np.argwhere(labels < 0):
                name = names[ids[y, z, x]]
                unknown[name] = unknown.get(name, 0) + 1
                if len(samples) < 12:
                    samples.append(
                        {
                            "block": name,
                            "position": [int(x0 + x), int(bounds.min[1] + y), int(z0 + z)],
                        }
                    )
            labels = np.maximum(labels, 0)
        labels_yzx[:, z0 - bounds.min[2] : z1 - bounds.min[2] + 1, x0 - bounds.min[0] : x1 - bounds.min[0] + 1] = labels
    if unknown:
        raise UnknownBlocksError(unknown, samples)
    return placement.from_world_order(labels_yzx)


# ------------------------------------------------------------- commands
def _write_volume(
    owner: list[sitk.Image], output: Path, names: Mapping[int, str], present: set[int]
) -> str:
    """Atomically write a verified labelmap and embed names in NIfTI outputs."""
    result = volume_mod.write_verified_volume(owner[0], str(output), release_source=owner.clear)
    kept = {label: name for label, name in names.items() if label in present}
    if result.output.format == "NIfTI" and kept:
        labelnames.embed_label_names(output, kept)
    return result.output.format


def to_world(
    candidate: VolumeCandidate,
    output: Path,
    options: WorldOptions,
    *,
    allow_large_volume: bool = False,
    log: Logger | None = None,
) -> WorldResult:
    """Convert one labelmap into a new Minecraft Java world directory."""
    started = time.time()
    if os.path.lexists(output):
        raise ValueError("output %s already exists; refusing to overwrite it" % output)
    if not output.name:
        raise ValueError("output must name a world directory")
    validate_pack_id(options.pack_id)
    for option, value in (("--world-name", options.world_name), ("--source-name", options.source_name)):
        if value is not None and not value.strip():
            raise ValueError("%s must not be empty" % option)
    blockgrid.validate_options(options.transform)
    vertical = _vertical_index(options)

    _say(log, "Loading labelmap ...")
    image, names = _load_labelmap(candidate, allow_large_volume=allow_large_volume)
    record = None
    if not options.transform.is_identity:
        _say(log, "Resampling and orienting the label grid ...")
        image, record = blockgrid.transform(image, options.transform)
    labels_ijk = sitk.GetArrayViewFromImage(image).transpose(2, 1, 0)
    placement = Placement.create(image.GetSize(), vertical)
    palette = assign(_label_counts(labels_ijk), names, options.palette)
    guides, spawn = _guides(placement, palette)
    source_path = candidate.source.path if isinstance(candidate.source, FileSource) else None
    manifest = _manifest(
        image,
        options.source_name or candidate.source_name,
        _sha256(source_path) if source_path is not None else "",
        placement,
        palette,
        record,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    staging_root = Path(tempfile.mkdtemp(prefix=".%s." % output.name, dir=output.parent))
    try:
        staging = staging_root / output.name
        staging.mkdir()
        _say(log, "Writing Minecraft Java %s chunks ..." % MINECRAFT_VERSION)
        chunk_count = _write_blocks(staging, placement, labels_ijk, palette, guides)
        _say(log, "  %d chunks" % chunk_count)
        _write_world_metadata(
            staging, options.world_name or output.name, options.pack_id, placement, guides, spawn
        )
        (staging / MANIFEST_DIR).mkdir()
        manifest_path(staging).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        os.rename(staging, output)
    finally:
        shutil.rmtree(staging_root, ignore_errors=True)
    size = tuple(int(value) for value in image.GetSize())
    return WorldResult(
        output_path=str(output),
        size=size,  # type: ignore[arg-type]
        voxel_count=int(np.prod(size)),
        labels=palette,
        dimension=placement.dimension,
        volume_bounds=placement.bounds,
        transform=record,
        seconds=time.time() - started,
    )


def resample(
    candidate: VolumeCandidate,
    output: Path,
    options: blockgrid.TransformOptions,
    *,
    allow_large_volume: bool = False,
    log: Logger | None = None,
) -> VolumeResult:
    """Write the derived block grid that ``to_world`` would build."""
    started = time.time()
    volume_output(output)
    _say(log, "Loading labelmap ...")
    image, names = _load_labelmap(candidate, allow_large_volume=allow_large_volume)
    _say(log, "Resampling and orienting the label grid ...")
    derived, record = blockgrid.transform(image, options)
    image = sitk.Image()
    counts = _label_counts(sitk.GetArrayViewFromImage(derived))
    pixel_type = sitk.GetArrayViewFromImage(derived).dtype.name
    size = tuple(int(value) for value in derived.GetSize())
    _say(log, "Writing %s ..." % output)
    output_format = _write_volume([derived], output, names, set(counts))
    del derived
    return VolumeResult(
        output_path=str(output),
        format=output_format,
        size=size,  # type: ignore[arg-type]
        pixel_type=pixel_type,
        label_voxels={label: count for label, count in counts.items() if label != 0},
        transform=record.as_dict(),
        seconds=time.time() - started,
    )


def validate(world: Path, *, log: Logger | None = None) -> ValidationResult:
    """Check the manifest, required chunks, and every block in the export volume."""
    manifest = load_manifest(world)
    _say(log, "Reading Minecraft blocks ...")
    labels = decode_world(world, manifest)
    counts = _label_counts(labels)
    return ValidationResult(
        required_chunks=len(manifest.placement.required_chunks()),
        checked_voxels=int(labels.size),
        label_voxels={label: count for label, count in counts.items() if label != 0},
    )


def to_labelmap(world: Path, output: Path, *, log: Logger | None = None) -> VolumeResult:
    """Export an edited world as a labelmap on the manifest's grid."""
    started = time.time()
    volume_output(output)
    manifest = load_manifest(world)
    _say(log, "Reading Minecraft blocks ...")
    labels = decode_world(world, manifest)
    counts = _label_counts(labels)
    image = sitk.GetImageFromArray(np.ascontiguousarray(labels.transpose(2, 1, 0)))
    del labels
    image.SetSpacing(manifest.spacing)
    image.SetOrigin(manifest.origin)
    image.SetDirection(manifest.direction)
    names = {entry.label: entry.name for entry in manifest.palette if entry.name}
    _say(log, "Writing %s ..." % output)
    output_format = _write_volume([image], output, names, set(counts))
    del image
    return VolumeResult(
        output_path=str(output),
        format=output_format,
        size=manifest.size,
        pixel_type=manifest.dtype,
        label_voxels={label: count for label, count in counts.items() if label != 0},
        transform=manifest.raw.get("transform"),
        seconds=time.time() - started,
    )
