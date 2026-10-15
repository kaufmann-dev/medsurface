"""Inert label blocks, approximate colors, and label-to-block assignment."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Mapping

MAX_LABELS = 127

_COLORS = (
    "red",
    "lime",
    "blue",
    "yellow",
    "magenta",
    "cyan",
    "orange",
    "purple",
    "green",
    "light_blue",
    "pink",
    "brown",
    "black",
    "light_gray",
    "gray",
    "white",
)
_SUFFIXES = ("concrete", "wool", "glazed_terracotta", "terracotta", "stained_glass")
_MATERIALS = (
    "bone_block",
    "calcite",
    "quartz_block",
    "chiseled_quartz_block",
    "quartz_bricks",
    "smooth_quartz",
    "iron_block",
    "gold_block",
    "diamond_block",
    "emerald_block",
    "lapis_block",
    "redstone_block",
    "coal_block",
    "copper_block",
    "amethyst_block",
    "prismarine",
    "stone",
    "cobblestone",
    "mossy_cobblestone",
    "stone_bricks",
    "mossy_stone_bricks",
    "cracked_stone_bricks",
    "chiseled_stone_bricks",
    "granite",
    "polished_granite",
    "diorite",
    "polished_diorite",
    "andesite",
    "polished_andesite",
    "deepslate",
    "polished_deepslate",
    "bricks",
    "oak_planks",
    "spruce_planks",
    "birch_planks",
    "jungle_planks",
    "acacia_planks",
    "dark_oak_planks",
    "mangrove_planks",
    "cherry_planks",
    "bamboo_planks",
    "crimson_planks",
    "warped_planks",
    "sea_lantern",
    "shroomlight",
    "glowstone",
    "crying_obsidian",
)

#: Deterministic, inert, non-falling, non-fluid full cubes without block entities.
BLOCKS: tuple[str, ...] = tuple(
    "minecraft:%s_%s" % (color, suffix) for suffix in _SUFFIXES for color in _COLORS
) + tuple("minecraft:%s" % name for name in _MATERIALS)

# Approximate average texture colors: concrete, wool, glazed terracotta,
# terracotta, stained glass.
_DYE_COLORS: dict[str, tuple[tuple[int, int, int], ...]] = {
    "white": ((207, 213, 214), (233, 236, 236), (188, 212, 202), (209, 178, 161), (255, 255, 255)),
    "orange": ((224, 97, 1), (240, 118, 19), (154, 147, 91), (161, 83, 37), (216, 127, 51)),
    "magenta": ((169, 48, 159), (189, 68, 179), (208, 100, 191), (149, 88, 108), (178, 76, 216)),
    "light_blue": ((36, 137, 199), (58, 175, 217), (94, 164, 208), (113, 108, 137), (102, 153, 216)),
    "yellow": ((241, 175, 21), (248, 198, 39), (234, 192, 88), (186, 133, 35), (229, 229, 51)),
    "lime": ((94, 169, 24), (112, 185, 25), (162, 197, 55), (103, 117, 52), (127, 204, 25)),
    "pink": ((214, 101, 143), (237, 141, 172), (235, 154, 181), (161, 78, 78), (242, 127, 165)),
    "gray": ((55, 58, 62), (62, 68, 71), (83, 90, 93), (57, 42, 35), (76, 76, 76)),
    "light_gray": ((125, 125, 115), (142, 142, 134), (144, 166, 167), (135, 107, 98), (153, 153, 153)),
    "cyan": ((21, 119, 136), (21, 137, 145), (52, 118, 125), (87, 91, 91), (76, 127, 153)),
    "purple": ((100, 32, 156), (121, 42, 172), (109, 48, 152), (118, 70, 86), (127, 63, 178)),
    "blue": ((45, 47, 143), (53, 57, 157), (47, 64, 139), (74, 59, 91), (51, 76, 178)),
    "brown": ((96, 60, 32), (114, 71, 40), (119, 106, 85), (77, 51, 36), (102, 76, 51)),
    "green": ((73, 91, 36), (84, 109, 27), (117, 142, 67), (76, 83, 42), (102, 127, 51)),
    "red": ((142, 33, 33), (161, 39, 34), (181, 59, 53), (143, 61, 46), (153, 51, 51)),
    "black": ((8, 10, 15), (20, 21, 25), (67, 30, 32), (37, 22, 16), (25, 25, 25)),
}
_MATERIAL_COLORS: dict[str, tuple[int, int, int]] = {
    "bone_block": (225, 221, 201),
    "calcite": (223, 224, 220),
    "quartz_block": (236, 230, 223),
    "chiseled_quartz_block": (231, 226, 218),
    "quartz_bricks": (234, 229, 221),
    "smooth_quartz": (236, 230, 223),
    "iron_block": (220, 220, 220),
    "gold_block": (246, 208, 61),
    "diamond_block": (98, 237, 228),
    "emerald_block": (42, 203, 87),
    "lapis_block": (31, 67, 140),
    "redstone_block": (175, 24, 5),
    "coal_block": (16, 16, 16),
    "copper_block": (192, 107, 79),
    "amethyst_block": (133, 97, 191),
    "prismarine": (99, 156, 151),
    "stone": (125, 125, 125),
    "cobblestone": (127, 127, 127),
    "mossy_cobblestone": (110, 118, 94),
    "stone_bricks": (122, 121, 122),
    "mossy_stone_bricks": (115, 121, 105),
    "cracked_stone_bricks": (118, 117, 118),
    "chiseled_stone_bricks": (119, 118, 119),
    "granite": (149, 103, 85),
    "polished_granite": (154, 106, 89),
    "diorite": (188, 188, 188),
    "polished_diorite": (192, 193, 194),
    "andesite": (136, 136, 136),
    "polished_andesite": (132, 134, 133),
    "deepslate": (80, 80, 82),
    "polished_deepslate": (72, 72, 73),
    "bricks": (150, 97, 83),
    "oak_planks": (162, 130, 78),
    "spruce_planks": (114, 84, 48),
    "birch_planks": (192, 175, 121),
    "jungle_planks": (160, 115, 80),
    "acacia_planks": (168, 90, 50),
    "dark_oak_planks": (66, 43, 20),
    "mangrove_planks": (117, 54, 48),
    "cherry_planks": (226, 178, 172),
    "bamboo_planks": (193, 173, 80),
    "crimson_planks": (101, 48, 70),
    "warped_planks": (43, 104, 99),
    "sea_lantern": (172, 199, 190),
    "shroomlight": (240, 146, 70),
    "glowstone": (171, 131, 84),
    "crying_obsidian": (32, 10, 60),
}

# Anatomical keyword families and their preferred blocks, tried in order.
_SEMANTIC: tuple[tuple[tuple[str, ...], tuple[str, ...]], ...] = (
    (
        ("bone", "vertebra", "rib", "sternum", "sacrum", "skull", "femur", "humerus",
         "clavicula", "scapula", "hip"),
        ("bone_block", "calcite", "quartz_block", "smooth_quartz", "light_gray_concrete",
         "white_wool"),
    ),
    (
        ("artery", "aorta", "heart", "myocard", "blood"),
        ("red_concrete", "red_wool", "red_glazed_terracotta", "red_terracotta"),
    ),
    (
        ("vein", "vena_cava", "portal_vein"),
        ("blue_concrete", "blue_wool", "blue_glazed_terracotta", "blue_terracotta"),
    ),
    (
        ("lung", "trachea", "bronch"),
        ("light_blue_concrete", "cyan_concrete", "light_blue_wool", "cyan_wool"),
    ),
    (
        ("brain", "spinal_cord"),
        ("pink_concrete", "magenta_concrete", "pink_wool", "magenta_wool"),
    ),
    (
        ("liver", "spleen", "muscle", "gland", "pancreas"),
        ("brown_concrete", "orange_concrete", "purple_concrete", "brown_wool"),
    ),
    (
        ("kidney", "bladder", "ureter", "prostate", "uterus"),
        ("orange_concrete", "yellow_concrete", "orange_wool", "yellow_wool"),
    ),
)


@dataclass
class PaletteEntry:
    label: int
    block: str
    name: str | None
    voxel_count: int
    legend_position: tuple[int, int, int] = (0, 0, 0)

    def as_dict(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "block": self.block,
            "name": self.name,
            "voxel_count": self.voxel_count,
            "legend_position": list(self.legend_position),
        }


@dataclass(frozen=True)
class PaletteOverride:
    block: str | None = None
    name: str | None = None


def block_color(block: str) -> tuple[int, int, int] | None:
    """Approximate average texture color of a supported label block."""
    block_id = block.removeprefix("minecraft:")
    for index, suffix in enumerate("_" + value for value in _SUFFIXES):
        if block_id.endswith(suffix):
            color = block_id[: -len(suffix)]
            # "light_gray_glazed_terracotta" must not match "_terracotta" first.
            if suffix == "_terracotta" and color.endswith("_glazed"):
                continue
            dyes = _DYE_COLORS.get(color)
            return dyes[index] if dyes else None
    return _MATERIAL_COLORS.get(block_id)


def _normalize_block(block: str) -> str:
    block = block.strip()
    return block if ":" in block else "minecraft:%s" % block


def parse_overrides(text: str) -> dict[int, PaletteOverride]:
    """Parse ``{"5": "bone_block"}`` or ``{"5": {"block": ..., "name": ...}}``."""
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError("palette file is not valid JSON: %s" % exc) from None
    if not isinstance(raw, dict):
        raise ValueError("palette file must be a JSON object keyed by label ID")
    overrides: dict[int, PaletteOverride] = {}
    for key, entry in raw.items():
        try:
            label = int(str(key).strip())
        except ValueError:
            label = 0
        if label <= 0:
            raise ValueError("palette key %r is not a positive label ID" % key)
        if isinstance(entry, str):
            overrides[label] = PaletteOverride(block=_normalize_block(entry))
        elif isinstance(entry, dict):

            def text_field(field: str) -> str | None:
                value = entry.get(field)
                if value is None:
                    return None
                if isinstance(value, str) and value.strip():
                    return value.strip()
                raise ValueError(
                    "palette label %d: %s must be a non-empty string" % (label, field)
                )

            block = text_field("block")
            overrides[label] = PaletteOverride(
                block=_normalize_block(block) if block else None, name=text_field("name")
            )
        else:
            raise ValueError("palette label %d must map to a block ID or an object" % label)
    return overrides


def _semantic_candidates(name: str | None) -> list[str]:
    if name is None:
        return []
    lowered = name.lower()
    for keywords, blocks in _SEMANTIC:
        if any(keyword in lowered for keyword in keywords):
            return ["minecraft:%s" % block for block in blocks]
    return []


def assign(
    counts: Mapping[int, int],
    names: Mapping[int, str],
    overrides: Mapping[int, PaletteOverride] | None = None,
) -> list[PaletteEntry]:
    """Assign one unique block per nonzero label, honoring overrides first."""
    overrides = overrides or {}
    labels = sorted((label, count) for label, count in counts.items() if label != 0)
    if len(labels) > MAX_LABELS:
        raise ValueError(
            "labelmap contains %d nonzero labels; the supported maximum is %d"
            % (len(labels), MAX_LABELS)
        )
    allowed = set(BLOCKS)
    used: set[str] = set()
    for label, override in sorted(overrides.items()):
        if override.block is None:
            continue
        if override.block not in allowed:
            raise ValueError(
                "palette label %d: %s is not one of the %d supported label blocks "
                "(see 'medsurface minecraft blocks')" % (label, override.block, MAX_LABELS)
            )
        if override.block in used:
            raise ValueError("palette maps more than one label to %s" % override.block)
        used.add(override.block)

    single_label = len(labels) == 1
    entries: list[PaletteEntry] = []
    for label, voxel_count in labels:
        chosen = overrides.get(label)
        name = (chosen.name if chosen else None) or names.get(label)
        if chosen is not None and chosen.block is not None:
            block = chosen.block
        else:
            candidates = _semantic_candidates(name)
            if single_label and name is None:
                candidates.insert(0, "minecraft:white_concrete")
            block = next(
                candidate
                for candidate in [*candidates, *BLOCKS]
                if candidate not in used
            )
            used.add(block)
        entries.append(PaletteEntry(label, block, name, int(voxel_count)))
    return entries
