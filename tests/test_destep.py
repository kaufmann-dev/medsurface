"""Masked stair-step fairing on a rippled surface with thin detail."""

from __future__ import annotations

import json

import numpy as np
import pytest
import SimpleITK as sitk
from typer.testing import CliRunner

from medsurface import cli, destep, pipeline, segment, surface, validate
from medsurface.presets import DestepSettings, validate_destep
from tests.mesh_helpers import sphere

RADIUS = 30.0
RIPPLE_MM = 10.0
RIPPLE_AMPLITUDE = 0.4
SPACING = 0.6
HALF = RADIUS + 10.0
FIN_X = (-12.0, -4.0, 4.0, 12.0)


def _rippled_labelmap() -> sitk.Image:
    """A sphere with broad z ripples and four thin, tall fins below it."""
    axis = np.arange(int(2 * HALF / SPACING)) * SPACING - HALF
    z, y, x = np.meshgrid(axis, axis, axis, indexing="ij")
    radius = RADIUS + RIPPLE_AMPLITUDE * np.sin(2 * np.pi * z / RIPPLE_MM)
    mask = np.sqrt(x**2 + y**2 + z**2) <= radius
    for center in FIN_X:
        mask |= (
            (np.abs(x - center) < 0.8)
            & (np.abs(y) < 6.0)
            & (z < -RADIUS + 15.0)
            & (z > -RADIUS - 7.0)
        )
    image = sitk.GetImageFromArray(mask.astype(np.uint8))
    image.SetSpacing((SPACING,) * 3)
    image.SetOrigin((-HALF,) * 3)
    return image


@pytest.fixture(scope="module")
def rippled():
    grid = segment.smooth_occupancy(segment.pad(_rippled_labelmap(), 1), 0.5)
    mesh = surface.transform(
        surface.marching_cubes(grid, segment.ISO_OCCUPANCY),
        surface.index_to_physical(grid),
    )
    finished = pipeline.finish_surface(
        mesh,
        surface_smooth_iters=20,
        simplify_error_mm=0.25,
        post_surface_smooth_iters=40,
        keep_largest_component=True,
        step=lambda _message, function: function(),
        log=lambda _message: None,
    )
    vertices = surface.to_arrays(finished.poly)[0]
    radius = np.linalg.norm(vertices, axis=1)
    dome = vertices[:, 2] > 10.0
    fins = (vertices[:, 2] < -RADIUS + 15.0) & (np.abs(vertices[:, 1]) < 6.5)
    fins &= radius > RADIUS + 1.0
    assert dome.sum() > 500 and fins.sum() > 50
    return finished.poly, vertices, dome, fins


def _ripple_amplitude(vertices: np.ndarray, selected: np.ndarray) -> float:
    """Least-squares amplitude of the known radial ripple."""
    radius = np.linalg.norm(vertices[selected], axis=1)
    phase = 2 * np.pi * vertices[selected, 2] / RIPPLE_MM
    design = np.c_[np.ones_like(phase), np.sin(phase), np.cos(phase)]
    coefficients = np.linalg.lstsq(design, radius, rcond=None)[0]
    return float(np.hypot(coefficients[1], coefficients[2]))


def _run(mesh, settings):
    result, stats = destep.destep_safely(mesh, settings)
    vertices, faces = surface.to_arrays(result)
    report = validate.validate_arrays(vertices, faces)
    assert report["valid"], report["problems"]
    return vertices, stats


def test_band_fairs_its_full_side_and_never_moves_the_frozen_side(rippled):
    mesh, before, dome, fins = rippled
    after, stats = _run(mesh, DestepSettings("band", full_mm=10.0, frozen_mm=-10.0))

    frozen = before[:, 2] <= -10.0
    assert np.array_equal(after[frozen], before[frozen])
    assert _ripple_amplitude(after, dome) < 0.8 * _ripple_amplitude(before, dome)
    assert stats.safeguard.accepted
    assert 0 < stats.faired_vertex_fraction < 1
    assert 0 < stats.frozen_vertex_fraction < 1


def test_reversed_band_fairs_the_low_side_instead(rippled):
    mesh, before, dome, fins = rippled
    after, _stats = _run(mesh, DestepSettings("band", full_mm=-10.0, frozen_mm=10.0))

    assert np.array_equal(after[dome], before[dome])
    assert np.abs(after[fins] - before[fins]).max() > 0


def test_auto_freezes_thin_detail_and_fairs_broad_surfaces(rippled):
    mesh, before, dome, fins = rippled
    after, stats = _run(mesh, DestepSettings("auto"))

    assert np.array_equal(after[fins], before[fins])
    assert _ripple_amplitude(after, dome) < 0.8 * _ripple_amplitude(before, dome)
    assert stats.faired_vertex_fraction > 0.5
    assert stats.frozen_vertex_fraction > 0


def test_all_fairs_every_region_within_the_displacement_limit(rippled):
    mesh, before, dome, fins = rippled
    after, stats = _run(mesh, DestepSettings("all", max_displacement_mm=0.3))

    moved = np.linalg.norm(after - before, axis=1)
    # MeshLib stores float32 coordinates.
    assert moved.max() <= 0.3 + 1e-5
    assert moved[fins].max() > 0
    assert stats.faired_vertex_fraction == 1.0
    assert stats.frozen_vertex_fraction == 0.0
    assert stats.clamped_vertices > 0
    assert stats.clamped_vertex_fraction == pytest.approx(
        stats.clamped_vertices / len(before)
    )


def test_band_that_misses_the_surface_names_its_axis_range(rippled):
    mesh = rippled[0]
    with pytest.raises(ValueError, match=r"whose x range is -\d"):
        destep.destep_safely(
            mesh, DestepSettings("band", axis="x", full_mm=200.0, frozen_mm=150.0)
        )


def test_revert_folds_restores_a_flipped_neighbourhood():
    original, faces = sphere(2, radius=10.0)
    folded = original.copy()
    folded[0] = -0.5 * original[0]
    adjacent = destep.adjacency(len(original), faces)

    restored, reverted = destep.revert_folds(folded, original, faces, adjacent)

    assert reverted >= 1
    assert np.array_equal(restored[0], original[0])


@pytest.mark.parametrize(
    ("settings", "message"),
    [
        (DestepSettings("everything"), "region"),
        (DestepSettings("all", iterations=0), "iterations"),
        (DestepSettings("all", max_displacement_mm=0.0), "maximum displacement"),
        (DestepSettings("all", max_displacement_mm=float("inf")), "maximum"),
        (DestepSettings("band", axis="w", full_mm=1, frozen_mm=0), "axis"),
        (DestepSettings("band", full_mm=1.0), "both full and frozen"),
        (DestepSettings("band", full_mm=1.0, frozen_mm=1.0), "must differ"),
        (DestepSettings("band", full_mm=float("nan"), frozen_mm=1.0), "finite"),
        (DestepSettings("auto", full_mm=1.0), "require the band region"),
    ],
)
def test_invalid_destep_settings_are_rejected(settings, message):
    with pytest.raises(ValueError, match=message):
        validate_destep(settings)


def test_finish_surface_records_destep_provenance_and_warnings(rippled):
    mesh = rippled[0]
    finished = pipeline.finish_surface(
        mesh,
        surface_smooth_iters=0,
        simplify_error_mm=0.0,
        post_surface_smooth_iters=0,
        keep_largest_component=False,
        step=lambda _message, function: function(),
        log=lambda _message: None,
        destep=DestepSettings("all", iterations=2400, max_displacement_mm=0.05),
    )

    record = finished.provenance["destep"]
    assert record["region"] == "all"
    assert record["iterations"] == 2400
    assert record["safeguard"]["accepted"]
    assert any("displacement limit" in message for message in finished.warnings)


def test_labelmap_extract_command_publishes_a_valid_destepped_mesh(tmp_path):
    source = tmp_path / "rippled.nii.gz"
    sitk.WriteImage(_rippled_labelmap(), str(source))
    output = tmp_path / "rippled.stl"
    report = tmp_path / "rippled.json"

    result = CliRunner().invoke(
        cli.app,
        [
            "labelmap",
            "extract",
            str(source),
            "-o",
            str(output),
            "--destep",
            "auto",
            "--json",
            str(report),
        ],
        prog_name="medsurface",
    )

    assert result.exit_code == 0, result.output
    assert "destep fairing" in result.stderr
    payload = json.loads(report.read_text())
    record = payload["provenance"]["surface_finishing"]["destep"]
    assert record["region"] == "auto"
    assert record["safeguard"]["accepted"]
    assert payload["provenance"]["surface"]["destep"]["region"] == "auto"
    assert payload["quality"]["valid"]
    assert validate.validate(str(output))["valid"]
