"""Minecraft world round trips, block grids, and the Anvil/NBT codecs."""

from __future__ import annotations

import gzip
import json
import struct
import subprocess
import sys
import zlib
from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk
from typer.testing import CliRunner

from medsurface import catalog, cli, labelnames
from medsurface.minecraft import anvil, blockgrid, nbt, palette, world

runner = CliRunner(env={"_TYPER_FORCE_DISABLE_TERMINAL": "1"})

# i -> Left, j -> Anterior, k -> Superior in LPS, like many DICOM-derived files.
LAI_DIRECTION = (1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0)


def _write(path: Path, labels_kji: np.ndarray, spacing=(1.0, 1.0, 1.0), *, origin=(0.0, 0.0, 0.0), direction=None) -> Path:
    image = sitk.GetImageFromArray(labels_kji)
    image.SetSpacing(tuple(float(value) for value in spacing))
    image.SetOrigin(origin)
    if direction is not None:
        image.SetDirection(direction)
    sitk.WriteImage(image, str(path))
    return path


def _candidate(path: Path):
    return catalog.select(catalog.discover(path), None)


def _labels(path: Path) -> np.ndarray:
    return sitk.GetArrayFromImage(sitk.ReadImage(str(path)))


def _sphere(size_ijk, spacing, radius_mm) -> np.ndarray:
    k, j, i = np.meshgrid(*(np.arange(n) for n in size_ijk[::-1]), indexing="ij")
    center = [(n - 1) * s / 2 for n, s in zip(size_ijk, spacing)]
    distance = (
        (i * spacing[0] - center[0]) ** 2
        + (j * spacing[1] - center[1]) ** 2
        + (k * spacing[2] - center[2]) ** 2
    )
    return (distance <= radius_mm**2).astype(np.uint8)


def _extent(labels_kji: np.ndarray, label: int) -> list[int]:
    k, j, i = np.nonzero(labels_kji == label)
    return [int(axis.max() - axis.min() + 1) for axis in (i, j, k)]


def _centroid(labels_kji: np.ndarray, label: int) -> np.ndarray:
    k, j, i = np.nonzero(labels_kji == label)
    return np.array([i.mean(), j.mean(), k.mean()])


def _replace_in_chunk(world_path: Path, chunk_x: int, chunk_z: int, old: bytes, new: bytes) -> None:
    """Rewrite one chunk's NBT bytes, appending it to its region like Minecraft."""
    path = anvil.region_path(world_path, chunk_x // 32, chunk_z // 32)
    data = bytearray(path.read_bytes())
    index = (chunk_x % 32) + (chunk_z % 32) * 32
    offset = int.from_bytes(data[4 * index : 4 * index + 3], "big") * anvil.SECTOR_BYTES
    length = struct.unpack(">I", data[offset : offset + 4])[0]
    raw = zlib.decompress(bytes(data[offset + 5 : offset + 4 + length]))
    assert old in raw
    compressed = zlib.compress(raw.replace(old, new))
    payload = struct.pack(">IB", len(compressed) + 1, 2) + compressed
    sectors = -(-len(payload) // anvil.SECTOR_BYTES)
    start = len(data) // anvil.SECTOR_BYTES
    data[4 * index : 4 * index + 4] = start.to_bytes(3, "big") + bytes((sectors,))
    data += payload + bytes(sectors * anvil.SECTOR_BYTES - len(payload))
    path.write_bytes(bytes(data))


def _nbt_string(value: str) -> bytes:
    encoded = value.encode()
    return struct.pack(">H", len(encoded)) + encoded


def _name_tag(block: str) -> bytes:
    return b"\x08" + _nbt_string("Name") + _nbt_string(block)


# ------------------------------------------------------------------ codecs
def test_nbt_round_trips_every_tag_type():
    root = {
        "byte": nbt.Byte(-3),
        "short": nbt.Short(1234),
        "int": nbt.Int(-70000),
        "long": nbt.Long(1 << 40),
        "float": nbt.Float(0.5),
        "double": nbt.Double(-2.25),
        "string": "minecraft:stone ü",
        "ints": nbt.IntArray([1, -2, 3]),
        "longs": nbt.LongArray([-1, 1 << 62]),
        "list": nbt.List(nbt.TAG_COMPOUND, [{"Name": "a"}, {"Name": "b"}]),
        "empty": nbt.List(nbt.TAG_STRING),
        "nested": {"inner": {"value": nbt.Byte(1)}},
    }
    decoded = nbt.decode(nbt.encode(root))

    assert decoded["byte"] == -3
    assert decoded["short"] == 1234
    assert decoded["int"] == -70000
    assert decoded["long"] == 1 << 40
    assert decoded["float"] == 0.5
    assert decoded["double"] == -2.25
    assert decoded["string"] == "minecraft:stone ü"
    assert decoded["ints"].tolist() == [1, -2, 3]
    assert decoded["longs"].tolist() == [-1, 1 << 62]
    assert decoded["list"] == [{"Name": "a"}, {"Name": "b"}]
    assert decoded["empty"] == []
    assert decoded["nested"] == {"inner": {"value": 1}}
    with pytest.raises(TypeError):
        nbt.encode({"plain": 1})


@pytest.mark.parametrize("bits", [4, 5, 7, 9, 13])
def test_packed_values_round_trip_without_spanning(bits):
    values = np.random.default_rng(bits).integers(0, 1 << bits, 4096)
    packed = anvil.pack(values, bits)

    assert len(packed) == -(-4096 // (64 // bits))
    assert np.array_equal(anvil.unpack(packed, bits, 4096), values)


def test_chunks_keep_air_first_and_use_at_least_four_bits():
    bounds = anvil.DimensionBounds(-64, 384)
    names = ["minecraft:air", "minecraft:red_concrete", "minecraft:bone_block"]
    column = np.zeros((384, 16, 16), dtype=np.uint8)
    column[-10 + 64, 4, 3] = 1
    column[0:16] = 2  # the bottom section (y -64..-49) is solid bone

    chunk = nbt.decode(nbt.encode(anvil.build_chunk(0, 0, column, bounds, names)))
    section = next(entry for entry in chunk["sections"] if entry["Y"] == -1)
    full = next(entry for entry in chunk["sections"] if entry["Y"] == -4)

    assert [entry["Name"] for entry in section["block_states"]["palette"]] == names[:2]
    assert len(section["block_states"]["data"]) == 256  # 4 bits per block
    assert [entry["Name"] for entry in full["block_states"]["palette"]] == [
        "minecraft:air",
        "minecraft:bone_block",
    ]
    decoded_names, ids = anvil.chunk_block_names(chunk, -64, 319)
    assert decoded_names[ids[-10 + 64, 4, 3]] == "minecraft:red_concrete"
    assert decoded_names[ids[-10 + 64, 4, 2]] == "minecraft:air"
    heights = anvil.unpack(chunk["Heightmaps"]["WORLD_SURFACE"], 9, 256).reshape(16, 16)
    assert heights[4, 3] == -10 + 1 + 64
    assert heights[0, 0] == 16  # top of the bone section: y -49 + 1 - (-64)


def test_regions_use_euclidean_coordinates_for_negative_chunks(tmp_path):
    bounds = anvil.DimensionBounds(-64, 384)
    names = ["minecraft:air", "minecraft:stone"]

    def make(chunk_x, chunk_z):
        column = np.zeros((384, 16, 16), dtype=np.uint8)
        column[100, chunk_z % 16, chunk_x % 16] = 1
        return anvil.build_chunk(chunk_x, chunk_z, column, bounds, names)

    positions = [(-1, -1), (-33, 0), (31, -32), (0, 0)]
    anvil.write_chunks(tmp_path, positions, make)

    assert anvil.region_path(tmp_path, -2, 0).is_file()
    assert anvil.region_path(tmp_path, 0, -1).is_file()
    reader = anvil.RegionReader(tmp_path)
    for chunk_x, chunk_z in positions:
        chunk = reader.chunk(chunk_x, chunk_z)
        assert chunk is not None
        assert (chunk["xPos"], chunk["zPos"]) == (chunk_x, chunk_z)
        names_read, ids = anvil.chunk_block_names(chunk, 36, 36)
        assert names_read[ids[0, chunk_z % 16, chunk_x % 16]] == "minecraft:stone"
    assert reader.chunk(5, 5) is None


# ---------------------------------------------------------------- palette
def test_palette_has_127_unique_colored_blocks():
    assert len(palette.BLOCKS) == palette.MAX_LABELS
    assert len(set(palette.BLOCKS)) == palette.MAX_LABELS
    assert all(palette.block_color(block) for block in palette.BLOCKS)
    assert palette.block_color("minecraft:light_gray_glazed_terracotta") == (144, 166, 167)
    assert palette.block_color("minecraft:light_gray_terracotta") == (135, 107, 98)


def test_overrides_reserve_blocks_and_rename_labels():
    overrides = palette.parse_overrides(
        '{"2": "bone_block", "3": {"name": "Heart", "block": "minecraft:red_wool"}}'
    )
    entries = palette.assign({0: 5, 1: 2, 2: 3, 3: 4}, {1: "femur_left"}, overrides)

    assert [entry.block for entry in entries] == [
        "minecraft:calcite",
        "minecraft:bone_block",
        "minecraft:red_wool",
    ]
    assert entries[2].name == "Heart"


def test_semantic_names_and_single_labels_choose_sensible_blocks():
    entries = palette.assign({0: 5, 1: 2, 2: 3}, {1: "femur_left", 2: "aorta"})
    assert [entry.block for entry in entries] == ["minecraft:bone_block", "minecraft:red_concrete"]
    assert palette.assign({0: 1, 9: 1}, {})[0].block == "minecraft:white_concrete"


@pytest.mark.parametrize(
    "text",
    ['{"1": "minecraft:dirt"}', '{"1": "bone_block", "2": "bone_block"}'],
)
def test_invalid_override_blocks_are_rejected(text):
    with pytest.raises(ValueError):
        palette.assign({1: 2, 2: 3}, {}, palette.parse_overrides(text))


@pytest.mark.parametrize("text", ['{"0": "bone_block"}', "[1]", '{"1": 5}', "{", '{"1": {"name": ""}}'])
def test_malformed_palette_files_are_rejected(text):
    with pytest.raises(ValueError):
        palette.parse_overrides(text)


def test_more_than_127_labels_are_rejected():
    with pytest.raises(ValueError, match="supported maximum is 127"):
        palette.assign({label: 1 for label in range(1, 129)}, {})


# ------------------------------------------------------------- round trip
@pytest.mark.parametrize("suffix", [".nii.gz", ".nrrd", ".mha"])
def test_multilabel_world_round_trips_labels_type_and_geometry(tmp_path, suffix):
    labels = np.zeros((12, 20, 24), dtype=np.uint16)
    labels[2:10, 4:16, 3:20] = 4
    labels[5:8, 8:12, 10:14] = 900
    labels[0, 0, 0] = 2
    source = _write(
        tmp_path / ("labels" + suffix),
        labels,
        (0.8, 0.6, 2.5),
        origin=(10.0, -20.0, 5.0),
        direction=LAI_DIRECTION,
    )
    names = {4: "femur_left", 900: "aorta"}
    if suffix == ".nii.gz":
        assert labelnames.embed_label_names(source, names)

    report = world.to_world(_candidate(source), tmp_path / "world", world.WorldOptions())
    assert report.transform is None
    assert world.validate(tmp_path / "world").label_voxels == {
        2: 1,
        4: int((labels == 4).sum()),
        900: 48,
    }
    output = tmp_path / ("edited" + suffix)
    world.to_labelmap(tmp_path / "world", output)

    stored = sitk.ReadImage(str(output))
    original = sitk.ReadImage(str(source))
    assert np.array_equal(sitk.GetArrayFromImage(stored), labels)
    assert stored.GetPixelIDValue() == original.GetPixelIDValue()
    assert np.allclose(stored.GetSpacing(), original.GetSpacing())
    assert np.allclose(stored.GetOrigin(), original.GetOrigin())
    assert np.allclose(stored.GetDirection(), original.GetDirection())
    manifest = json.loads(world.manifest_path(tmp_path / "world").read_text())
    blocks = {entry["label"]: entry["block"] for entry in manifest["palette"]}
    if suffix == ".nii.gz":
        assert labelnames.read_label_names(output) == names
        assert blocks[4] == "minecraft:bone_block"
        assert blocks[2] == "minecraft:red_concrete"  # unnamed labels are assigned first
        assert blocks[900] == "minecraft:red_wool"  # so the aorta takes the next red
    assert "transform" not in manifest
    assert manifest["grid"]["size"] == [24, 20, 12]


def test_default_world_places_one_voxel_per_block_with_guides(tmp_path):
    labels = np.zeros((6, 5, 8), dtype=np.uint8)
    labels[:, :, 0] = 1  # i = 0 lies on the world's minimum X
    source = _write(tmp_path / "labels.nii.gz", labels)

    world.to_world(_candidate(source), tmp_path / "world", world.WorldOptions())
    manifest = world.load_manifest(tmp_path / "world")
    bounds = manifest.placement.bounds

    assert manifest.raw["axes"] == {"x": "x", "y": "z", "z": "y"}
    assert bounds.min[0] == -4 and bounds.max[0] == 3
    assert bounds.min[2] == -2 and bounds.max[2] == 2
    assert bounds.max[1] - bounds.min[1] == 5
    reader = anvil.RegionReader(tmp_path / "world")
    names, ids = anvil.chunk_block_names(reader.chunk(-1, -1), bounds.min[1], bounds.min[1])
    assert names[ids[0, (bounds.min[2]) % 16, bounds.min[0] % 16]] == "minecraft:white_concrete"
    names, ids = anvil.chunk_block_names(reader.chunk(-1, -1), bounds.min[1], bounds.min[1])
    corner = names[ids[0, (bounds.min[2] - 2) % 16, (bounds.min[0] - 2) % 16]]
    assert corner == "minecraft:sea_lantern"
    legend = manifest.palette[0].legend_position
    assert legend[0] > bounds.max[0]
    level = nbt.decode(gzip.decompress((tmp_path / "world" / "level.dat").read_bytes()))
    assert level["Data"]["DataVersion"] == anvil.MINECRAFT_DATA_VERSION
    assert level["Data"]["LevelName"] == "world"
    assert list(level["Data"]["DataPacks"]["Enabled"]) == ["vanilla", "file/medsurface"]
    assert (tmp_path / "world" / "datapacks" / "medsurface" / "pack.mcmeta").is_file()


def test_existing_world_output_is_never_overwritten(tmp_path):
    source = _write(tmp_path / "labels.nrrd", np.ones((4, 4, 4), dtype=np.uint8))
    (tmp_path / "world").mkdir()
    with pytest.raises(ValueError, match="refusing to overwrite"):
        world.to_world(_candidate(source), tmp_path / "world", world.WorldOptions())


def test_unknown_blocks_abort_validation_and_export(tmp_path):
    labels = np.ones((4, 6, 6), dtype=np.uint8)
    source = _write(tmp_path / "labels.nii.gz", labels)
    world.to_world(_candidate(source), tmp_path / "world", world.WorldOptions())
    _replace_in_chunk(
        tmp_path / "world",
        0,
        0,
        _name_tag("minecraft:white_concrete"),
        _name_tag("minecraft:dirt"),
    )

    with pytest.raises(world.UnknownBlocksError, match="minecraft:dirt") as caught:
        world.validate(tmp_path / "world")
    assert caught.value.counts == {"minecraft:dirt": 4 * 3 * 3}
    output = tmp_path / "edited.nii.gz"
    with pytest.raises(world.UnknownBlocksError):
        world.to_labelmap(tmp_path / "world", output)
    assert not output.exists()


def test_block_state_properties_do_not_change_the_label(tmp_path):
    labels = np.zeros((4, 6, 6), dtype=np.uint8)
    labels[:, 3:, 3:] = 5
    source = _write(tmp_path / "labels.nrrd", labels)
    overrides = palette.parse_overrides('{"5": "orange_glazed_terracotta"}')
    world.to_world(_candidate(source), tmp_path / "world", world.WorldOptions(palette=overrides))
    tag = _name_tag("minecraft:orange_glazed_terracotta")
    properties = b"\x0a" + _nbt_string("Properties") + b"\x08" + _nbt_string("facing") + _nbt_string("north") + b"\x00"
    _replace_in_chunk(tmp_path / "world", 0, 0, tag + b"\x00", tag + properties + b"\x00")

    world.to_labelmap(tmp_path / "world", tmp_path / "edited.nrrd")
    assert np.array_equal(_labels(tmp_path / "edited.nrrd"), labels)


def test_tall_volume_uses_a_custom_dimension_height(tmp_path):
    labels = np.zeros((400, 4, 4), dtype=np.uint8)
    labels[::7, 1, 2] = 3
    source = _write(tmp_path / "tall.mha", labels)

    world.to_world(_candidate(source), tmp_path / "world", world.WorldOptions())
    manifest = json.loads(world.manifest_path(tmp_path / "world").read_text())
    dimension_type = json.loads(
        (
            tmp_path / "world" / "datapacks" / "medsurface" / "data" / "medsurface"
            / "dimension_type" / "overworld.json"
        ).read_text()
    )

    assert manifest["dimension_bounds"] == {"min_y": -80, "height": 400}
    assert (dimension_type["min_y"], dimension_type["height"]) == (-80, 400)
    world.to_labelmap(tmp_path / "world", tmp_path / "edited.mha")
    assert np.array_equal(_labels(tmp_path / "edited.mha"), labels)


def test_vertical_axis_above_the_dimension_limit_is_rejected(tmp_path):
    source = _write(tmp_path / "tall.nrrd", np.ones((4065, 4, 4), dtype=np.uint8))
    with pytest.raises(ValueError, match="fitting axes: x, y"):
        world.to_world(_candidate(source), tmp_path / "world", world.WorldOptions())
    world.to_world(
        _candidate(source), tmp_path / "sideways", world.WorldOptions(vertical_axis="x")
    )


def test_pack_id_renames_the_dimension_data_pack(tmp_path):
    source = _write(tmp_path / "labels.nrrd", np.ones((4, 4, 4), dtype=np.uint8))
    world.to_world(
        _candidate(source), tmp_path / "world", world.WorldOptions(pack_id="clinic.demo")
    )
    level = nbt.decode(gzip.decompress((tmp_path / "world" / "level.dat").read_bytes()))
    generator = nbt.decode(
        gzip.decompress((tmp_path / "world/data/minecraft/world_gen_settings.dat").read_bytes())
    )

    assert list(level["Data"]["DataPacks"]["Enabled"]) == ["vanilla", "file/clinic.demo"]
    assert (
        generator["data"]["dimensions"]["minecraft:overworld"]["type"] == "clinic.demo:overworld"
    )
    pack = json.loads((tmp_path / "world/datapacks/clinic.demo/pack.mcmeta").read_text())
    assert pack["pack"]["description"] == "Custom-height Overworld"
    for invalid in ("", "Upper", "minecraft", "vanilla", "a/b"):
        with pytest.raises(ValueError):
            world.validate_pack_id(invalid)


def test_open_worlds_are_refused(tmp_path):
    source = _write(tmp_path / "labels.nrrd", np.ones((4, 4, 4), dtype=np.uint8))
    world.to_world(_candidate(source), tmp_path / "world", world.WorldOptions())
    holder = subprocess.Popen(
        [
            sys.executable,
            "-c",
            "import fcntl, sys\n"
            "handle = open(sys.argv[1], 'r+b')\n"
            "fcntl.lockf(handle, fcntl.LOCK_EX)\n"
            "print('locked', flush=True)\n"
            "sys.stdin.read()\n",
            str(tmp_path / "world" / "session.lock"),
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
    )
    try:
        assert holder.stdout is not None and holder.stdout.readline().strip() == "locked"
        with pytest.raises(ValueError, match="currently open"):
            world.validate(tmp_path / "world")
    finally:
        holder.communicate("")
    assert world.validate(tmp_path / "world").checked_voxels == 64


# ------------------------------------------------------------- block grid
def test_block_size_parses_millimetres_and_auto():
    assert blockgrid.parse_block_mm("2") == 2.0
    assert blockgrid.parse_block_mm("AUTO") == "auto"
    for invalid in ("0", "-1", "x", "nan"):
        with pytest.raises(ValueError):
            blockgrid.parse_block_mm(invalid)


def test_transform_options_are_validated():
    for options in (
        blockgrid.TransformOptions(crop=True),
        blockgrid.TransformOptions(block_mm=1.0, fill_threshold=0.0),
        blockgrid.TransformOptions(block_mm=1.0, fill_threshold=1.5),
    ):
        with pytest.raises(ValueError):
            blockgrid.validate_options(options)


def test_block_resampling_makes_anisotropic_spheres_round_and_keeps_volume(tmp_path):
    spacing = (0.5, 0.5, 2.5)
    labels = _sphere((60, 60, 12), spacing, 12.0)
    source = _write(tmp_path / "sphere.nrrd", labels, spacing)
    raw = _extent(labels, 1)
    assert raw[0] > 4 * raw[2]

    output = tmp_path / "sphere-1mm.nii.gz"
    report = world.resample(
        _candidate(source), output, blockgrid.TransformOptions(block_mm=1.0)
    )
    derived = sitk.ReadImage(str(output))
    blocks = sitk.GetArrayFromImage(derived)

    assert np.allclose(derived.GetSpacing(), (1.0, 1.0, 1.0))
    assert report.size == (30, 30, 30)
    assert all(abs(value - 24) <= 2 for value in _extent(blocks, 1))
    assert abs((blocks == 1).sum() / ((labels == 1).sum() * 0.625) - 1.0) < 0.03


def test_thin_structures_survive_downsampling_by_coverage(tmp_path):
    labels = np.zeros((20, 20, 20), dtype=np.uint8)
    labels[:, :, 8:10] = 7  # a 1 mm slab of two 0.5 mm voxels
    source = _write(tmp_path / "slab.nrrd", labels, (0.5, 0.5, 0.5))

    world.resample(
        _candidate(source), tmp_path / "slab.nrrd.mha", blockgrid.TransformOptions(block_mm=1.0)
    )
    assert _extent(_labels(tmp_path / "slab.nrrd.mha"), 7) == [1, 10, 10]


def test_ties_prefer_the_smaller_label_and_coverage_controls_filling(tmp_path):
    labels = np.zeros((4, 4, 4), dtype=np.uint8)
    labels[:, :, 0] = 9
    labels[:, :, 1] = 3
    tie = sitk.GetImageFromArray(labels)
    half_covered = sitk.GetImageFromArray(np.where(labels == 9, 9, 0).astype(np.uint8))

    tied, _ = blockgrid.transform(tie, blockgrid.TransformOptions(block_mm=2.0))
    filled, _ = blockgrid.transform(half_covered, blockgrid.TransformOptions(block_mm=2.0))
    strict, _ = blockgrid.transform(
        half_covered, blockgrid.TransformOptions(block_mm=2.0, fill_threshold=0.75)
    )
    assert sitk.GetArrayFromImage(tied)[:, :, 0].tolist() == [[3, 3], [3, 3]]
    assert sitk.GetArrayFromImage(filled)[:, :, 0].tolist() == [[9, 9], [9, 9]]
    assert not sitk.GetArrayFromImage(strict).any()


def _oriented_fixture(path: Path) -> Path:
    labels = np.zeros((14, 10, 12), dtype=np.uint8)
    labels[4:10, 3:7, 4:8] = 1
    labels[5:8, 4:6, 1] = 2  # patient right: small i, because +i points left
    labels[13, 5, 5:7] = 3  # superior: large k
    labels[6, 0, 5:7] = 4  # posterior: small j, because +j points anterior
    return _write(path, labels, (1.0, 1.0, 2.0), origin=(30.0, -10.0, 5.0), direction=LAI_DIRECTION)


def _assert_anatomical_axes(labels_kji: np.ndarray) -> None:
    body = _centroid(labels_kji, 1)
    assert _centroid(labels_kji, 2)[0] > body[0] + 2  # +i right
    assert _centroid(labels_kji, 3)[2] > body[2] + 2  # +k superior
    assert _centroid(labels_kji, 4)[1] > body[1] + 1  # +j posterior


def test_anatomical_orientation_permutes_without_resampling(tmp_path):
    source = _oriented_fixture(tmp_path / "oriented.nii.gz")
    output = tmp_path / "anatomical.nii.gz"

    report = world.resample(
        _candidate(source), output, blockgrid.TransformOptions(orientation="anatomical")
    )
    derived = sitk.ReadImage(str(output))
    original = sitk.ReadImage(str(source))

    assert report.size == (12, 10, 14)
    assert np.allclose(derived.GetSpacing(), (1.0, 1.0, 2.0))
    _assert_anatomical_axes(sitk.GetArrayFromImage(derived))
    values, counts = np.unique(sitk.GetArrayFromImage(original), return_counts=True)
    assert report.label_voxels == {int(v): int(c) for v, c in zip(values, counts) if v}
    # Every voxel keeps its physical position.
    index = (5, 4, 6)
    point = derived.TransformIndexToPhysicalPoint(index)
    assert derived[index] == original[original.TransformPhysicalPointToIndex(point)]


def test_anatomical_resampling_builds_an_upright_unmirrored_world(tmp_path):
    source = _oriented_fixture(tmp_path / "oriented.nrrd")
    options = blockgrid.TransformOptions(block_mm=1.0, orientation="anatomical")

    report = world.to_world(
        _candidate(source), tmp_path / "world", world.WorldOptions(transform=options)
    )
    manifest = world.load_manifest(tmp_path / "world")
    labels_ijk = world.decode_world(tmp_path / "world", manifest)
    # Minecraft: +X east = patient right, +Y up = superior, +Z south = posterior.
    blocks_yzx = manifest.placement.to_world_order(labels_ijk)
    position = {label: np.argwhere(blocks_yzx == label).mean(axis=0) for label in (1, 2, 3, 4)}

    assert report.transform is not None
    assert manifest.raw["axes"] == {"x": "x", "y": "z", "z": "y"}
    assert position[2][2] > position[1][2] + 2  # X
    assert position[3][0] > position[1][0] + 2  # Y
    assert position[4][1] > position[1][1] + 1  # Z
    with pytest.raises(ValueError, match="omit --vertical-axis"):
        world.to_world(
            _candidate(source),
            tmp_path / "other",
            world.WorldOptions(vertical_axis="x", transform=options),
        )


def test_scaled_world_round_trips_on_the_block_grid(tmp_path):
    spacing = (0.5, 0.5, 2.5)
    labels = _sphere((40, 40, 8), spacing, 8.0)
    labels[:, :20][labels[:, :20] == 1] = 3
    source = _write(tmp_path / "labels.nii.gz", labels, spacing)
    labelnames.embed_label_names(source, {1: "liver", 3: "spleen"})
    options = blockgrid.TransformOptions(block_mm=1.0, crop=True)

    world.resample(_candidate(source), tmp_path / "derived.nii.gz", options)
    report = world.to_world(
        _candidate(source), tmp_path / "world", world.WorldOptions(transform=options)
    )
    exported = world.to_labelmap(tmp_path / "world", tmp_path / "edited.nii.gz")

    assert report.transform is not None and report.transform.cropped
    assert exported.transform is not None and exported.transform["block_mm"] == 1.0
    assert np.array_equal(_labels(tmp_path / "edited.nii.gz"), _labels(tmp_path / "derived.nii.gz"))
    edited = sitk.ReadImage(str(tmp_path / "edited.nii.gz"))
    assert np.allclose(edited.GetSpacing(), (1.0, 1.0, 1.0))
    assert labelnames.read_label_names(tmp_path / "edited.nii.gz") == {1: "liver", 3: "spleen"}


# --------------------------------------------------------------------- CLI
def test_cli_round_trip_with_json_reports(tmp_path):
    labels = np.zeros((6, 8, 8), dtype=np.uint8)
    labels[1:5, 2:6, 2:6] = 2
    source = _write(tmp_path / "labels.nii.gz", labels)
    world_path = tmp_path / "my world"

    created = runner.invoke(
        cli.app,
        ["minecraft", "to-world", str(source), "-o", str(world_path), "--json", str(tmp_path / "w.json")],
        prog_name="medsurface",
    )
    assert created.exit_code == 0, created.output
    assert "minecraft to-labelmap" in created.stdout
    assert json.loads((tmp_path / "w.json").read_text())["nonzero_labels"] == 1

    inspected = runner.invoke(cli.app, ["minecraft", "inspect", str(world_path)], prog_name="medsurface")
    assert inspected.exit_code == 0
    assert "minecraft:white_concrete" in inspected.stdout

    validated = runner.invoke(
        cli.app,
        ["minecraft", "validate", str(world_path), "--progress", "json"],
        prog_name="medsurface",
    )
    events = [json.loads(line) for line in validated.stderr.splitlines()]
    assert validated.exit_code == 0
    assert events[-1]["event"] == "result"
    assert events[-1]["result"]["label_voxels"] == {"2": 64}

    exported = runner.invoke(
        cli.app,
        ["minecraft", "to-labelmap", str(world_path), "-o", str(tmp_path / "edited.nrrd")],
        prog_name="medsurface",
    )
    assert exported.exit_code == 0, exported.output
    assert np.array_equal(_labels(tmp_path / "edited.nrrd"), labels)


@pytest.mark.parametrize(
    "extra, message",
    [
        (["--crop"], "--crop requires --block-mm"),
        (["--fill-threshold", "0.3"], "--fill-threshold requires --block-mm"),
        (["--block-mm", "0"], "positive number"),
        (["--block-mm", "1", "--fill-threshold", "2"], "at most 1"),
        (["--orient", "anatomical", "--vertical-axis", "x"], "omit --vertical-axis"),
        (["--pack-id", "minecraft"], "--pack-id"),
        (["--world-name", " "], "--world-name must not be empty"),
    ],
)
def test_cli_rejects_invalid_world_options(tmp_path, extra, message):
    source = _write(tmp_path / "labels.nrrd", np.ones((4, 4, 4), dtype=np.uint8))
    result = runner.invoke(
        cli.app,
        ["minecraft", "to-world", str(source), "-o", str(tmp_path / "world"), *extra],
        prog_name="medsurface",
    )
    assert result.exit_code == 2
    assert message in " ".join(result.stderr.split())
    assert not (tmp_path / "world").exists()


def test_cli_blocks_lists_every_block_as_json(tmp_path):
    result = runner.invoke(
        cli.app,
        ["minecraft", "blocks", "--json", str(tmp_path / "blocks.json")],
        prog_name="medsurface",
    )
    blocks = json.loads((tmp_path / "blocks.json").read_text())["blocks"]

    assert result.exit_code == 0
    assert len(blocks) == palette.MAX_LABELS
    assert blocks[0] == {"block": "minecraft:red_concrete", "rgb": [142, 33, 33], "hex": "#8e2121"}
