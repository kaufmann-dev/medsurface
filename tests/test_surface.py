"""Surface extraction: capping, cavity removal, coordinate mapping, volume fidelity."""

import math
import os

import numpy as np
import pytest
import SimpleITK as sitk

from dicom_surface import segment, surface, validate
from tests.mesh_helpers import (
    box,
    concatenate,
    mean_adjacency_angle,
    sphere,
    torus,
    transformed,
    volume,
    write,
)


def _mesh(image, cap=True, smooth_iters=0, largest=False):
    if cap:
        image = segment.pad(image, 1)
    affine = surface.index_to_physical(image)
    poly = surface.marching_cubes(image)
    if smooth_iters:
        poly = surface.smooth(poly, smooth_iters, 0.1)
    if largest:
        poly, _n = surface.largest_component(poly)
    return surface.transform(poly, affine)


def _write(poly, tmp_path, name="m.stl"):
    p = os.path.join(str(tmp_path), name)
    surface.write(poly, p)
    return p


def test_solid_sphere_is_watertight_with_correct_volume(solid_sphere, tmp_path):
    image, radius = solid_sphere
    poly = _mesh(image, smooth_iters=10)
    report = validate.validate(_write(poly, tmp_path))

    assert report["watertight"]
    assert report["valid"]
    assert report["problems"] == []
    assert report["components"] == 1
    assert report["boundary_edges"] == 0
    assert report["genus"] == 0

    expected = 4.0 / 3.0 * math.pi * radius**3
    assert report["volume_mm3"] == pytest.approx(expected, rel=0.02)


def test_uncapped_boundary_leaves_the_mesh_open(clipped_sphere, tmp_path):
    """Regression: marching cubes does not close the surface at the image edge.

    Anatomy truncated by the scanner's field of view therefore yields an OPEN
    mesh unless the volume is padded with background first.
    """
    assert segment.FOREGROUND == 1
    from dicom_surface.volume import touches_boundary

    assert touches_boundary(clipped_sphere), "fixture must actually touch the edge"

    open_report = validate.validate(
        _write(_mesh(clipped_sphere, cap=False), tmp_path, "open.stl"),
    )
    assert not open_report["watertight"]
    assert not open_report["valid"]
    assert open_report["boundary_edges"] > 0

    capped_report = validate.validate(
        _write(_mesh(clipped_sphere, cap=True), tmp_path, "capped.stl"),
    )
    assert capped_report["watertight"]
    assert capped_report["boundary_edges"] == 0


def test_largest_component_removes_internal_cavities(hollow_sphere, tmp_path):
    """Regression: an enclosed cavity is background, so no labelmap island filter
    can remove it. It survives as a second, disconnected surface shell."""
    image, outer_r, inner_r = hollow_sphere

    # The labelmap has exactly ONE foreground component ...
    cc = sitk.ConnectedComponent(image)
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(cc)
    assert len(stats.GetLabels()) == 1

    # ... yet the surface has two shells.
    all_shells = _mesh(image)
    report_all = validate.validate(_write(all_shells, tmp_path, "shells.stl"))
    assert report_all["components"] == 2

    kept = _mesh(image, largest=True)
    report_one = validate.validate(_write(kept, tmp_path, "outer.stl"))
    assert report_one["components"] == 1
    assert report_one["watertight"]

    # Keeping the outer shell should give the solid sphere's volume, not the
    # hollow shell's volume.
    solid = 4.0 / 3.0 * math.pi * outer_r**3
    hollow = solid - 4.0 / 3.0 * math.pi * inner_r**3
    assert report_one["volume_mm3"] == pytest.approx(solid, rel=0.03)
    assert abs(report_one["volume_mm3"] - hollow) > 0.05 * solid


def test_index_to_physical_matches_simpleitk():
    """The affine must agree with SimpleITK for oblique orientations too."""
    arr = np.zeros((4, 5, 6), dtype=np.uint8)
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((0.3, 0.7, 2.0))
    img.SetOrigin((-13.0, 22.0, 5.5))
    theta = 0.37
    d = [
        math.cos(theta), -math.sin(theta), 0.0,
        math.sin(theta), math.cos(theta), 0.0,
        0.0, 0.0, 1.0,
    ]
    img.SetDirection(d)

    m = surface.index_to_physical(img)
    for idx in [(0, 0, 0), (5, 4, 3), (2, 1, 0)]:
        expected = img.TransformIndexToPhysicalPoint(idx)
        got = (m[:3, :3] @ np.asarray(idx) + m[:3, 3]).tolist()
        assert got == pytest.approx(list(expected), abs=1e-9)


def test_padding_preserves_physical_coordinates(solid_sphere):
    image, _ = solid_sphere
    padded = segment.pad(image, 1)
    # index (0,0,0) of the original is index (1,1,1) of the padded image
    assert padded.TransformIndexToPhysicalPoint((1, 1, 1)) == pytest.approx(
        image.TransformIndexToPhysicalPoint((0, 0, 0))
    )
    assert padded.GetSize() == tuple(s + 2 for s in image.GetSize())


def test_looser_simplification_error_produces_no_more_faces(solid_sphere):
    image, _ = solid_sphere
    poly = _mesh(image, smooth_iters=5)
    before = poly.topology.numValidFaces()
    strict = surface.decimate(poly, 0.05)
    loose = surface.decimate(poly, 0.25)
    assert strict.topology.numValidFaces() < before
    assert loose.topology.numValidFaces() <= strict.topology.numValidFaces()


def test_smoothing_trades_roughness_for_displacement(solid_sphere):
    """More iterations means a smoother surface that sits further from the data.

    Guards the smoothing knob users actually feel. Volume-preserving relaxation
    must not contract the surface toward its centroid.
    """
    image, radius = solid_sphere
    padded = segment.pad(image, 1)
    affine = surface.index_to_physical(padded)
    raw = surface.transform(surface.marching_cubes(padded), affine)

    exact = 4.0 / 3.0 * math.pi * radius**3

    roughness = []
    for iterations in (0, 10, 25, 40):
        poly = surface.smooth(raw, iterations, 0.1) if iterations else raw
        mesh = surface.to_arrays(poly)
        roughness.append(mean_adjacency_angle(mesh))
        # no shrinkage: volume stays within 2% of the analytic sphere
        assert volume(mesh) == pytest.approx(exact, rel=0.02), iterations

    assert roughness == sorted(roughness, reverse=True), roughness
    assert roughness[0] > roughness[-1]


def test_safe_smoothing_leaves_a_clean_result_exactly_unchanged(solid_sphere):
    image, _ = solid_sphere
    original = _mesh(image)
    requested = surface.smooth(original, 10, 0.1)

    guarded, stats = surface.protect_smoothed_surface(original, requested, 10)

    requested_verts, requested_faces = surface.to_arrays(requested)
    guarded_verts, guarded_faces = surface.to_arrays(guarded)
    assert np.array_equal(guarded_faces, requested_faces)
    assert np.array_equal(guarded_verts, requested_verts)
    assert stats.initial_self_intersecting_faces == 0
    assert stats.protected_vertices == 0


def test_smoothing_safeguard_changes_only_a_collision_neighborhood():
    left = sphere(subdivisions=2, radius=1.0)
    right = transformed(sphere(subdivisions=2, radius=1.0), translation=(2.2, 0.0, 0.0))
    original_verts, original_faces = concatenate(left, right)

    moved = original_verts.copy()
    moved[len(left[0]):, 0] -= 0.4
    original = surface.from_arrays(original_verts, original_faces)
    intersecting = surface.from_arrays(moved, original_faces)
    assert surface.selected_self_intersecting_faces(intersecting).any()

    guarded, stats = surface.protect_smoothed_surface(original, intersecting, 25)

    assert not surface.selected_self_intersecting_faces(guarded).any()
    assert 0 < stats.protected_vertices < len(original_verts)
    assert surface.count_defects(guarded) == surface.count_defects(original)


def test_meshlib_marks_both_faces_in_an_intersecting_pair():
    import meshlib.mrmeshnumpy as mrmeshnumpy
    import meshlib.mrmeshpy as mrmeshpy
    mesh = concatenate(
        sphere(subdivisions=2, radius=1.0),
        transformed(sphere(subdivisions=2, radius=1.0), translation=(1.6, 0.0, 0.0)),
    )
    poly = surface.from_arrays(*mesh)

    selected = surface.selected_self_intersecting_faces(poly)
    mr_mesh = mrmeshnumpy.meshFromFacesVerts(mesh[1].astype(np.int32), mesh[0])
    pairs = mrmeshpy.findSelfCollidingTriangles(mrmeshpy.MeshPart(mr_mesh))

    assert pairs
    pair = pairs[0]
    assert selected[int(pair.aFace)]
    assert selected[int(pair.bFace)]


def test_count_defects_matches_validate(clipped_sphere, solid_sphere):
    closed = _mesh(solid_sphere[0])
    assert surface.count_defects(closed) == (0, 0)

    opened = _mesh(clipped_sphere, cap=False)
    boundary, holes = surface.count_defects(opened)
    assert boundary > 0
    assert holes > 0


@pytest.mark.parametrize("mm", [0.6, 1.0, 1.5])
def test_isotropic_resampling_stays_watertight(solid_sphere, tmp_path, mm):
    """The tested padded sphere stays closed after occupancy resampling."""
    image, radius = solid_sphere
    padded = segment.pad(image, 1)
    grid = segment.resample_isotropic(padded, mm)

    poly = surface.marching_cubes(grid, segment.ISO_OCCUPANCY)
    poly = surface.smooth(poly, 10, 0.1)
    poly = surface.transform(poly, surface.index_to_physical(grid))

    assert surface.count_defects(poly) == (0, 0)
    report = validate.validate(_write(poly, tmp_path, "iso%s.stl" % mm))
    assert report["watertight"]
    assert report["components"] == 1
    expected = 4.0 / 3.0 * math.pi * radius**3
    assert report["volume_mm3"] == pytest.approx(expected, rel=0.05)


def test_coarser_resampling_yields_fewer_triangles(solid_sphere):
    image, _ = solid_sphere
    padded = segment.pad(image, 1)
    counts = []
    for mm in (0.5, 1.0, 2.0):
        grid = segment.resample_isotropic(padded, mm)
        poly = surface.marching_cubes(grid, segment.ISO_OCCUPANCY)
        counts.append(poly.topology.numValidFaces())
    assert counts[0] > counts[1] > counts[2]


def test_resampling_preserves_physical_placement(solid_sphere):
    image, _ = solid_sphere
    padded = segment.pad(image, 1)
    grid = segment.resample_isotropic(padded, 1.0)
    assert grid.GetSpacing() == pytest.approx((1.0, 1.0, 1.0))
    # the resampled grid must fully contain the original extent
    for axis in range(3):
        assert grid.GetSize()[axis] * 1.0 >= padded.GetSize()[axis] * padded.GetSpacing()[axis]


@pytest.mark.parametrize("error_mm", [0.05, 0.1, 0.25])
def test_simplification_never_opens_a_closed_mesh(solid_sphere, tmp_path, error_mm):
    """MeshLib simplification must keep a closed surface watertight."""
    image, radius = solid_sphere
    poly = _mesh(image, smooth_iters=10)
    assert surface.count_defects(poly) == (0, 0)

    smaller = surface.decimate(poly, error_mm)

    assert surface.count_defects(smaller) == (0, 0)
    report = validate.validate(_write(smaller, tmp_path, "dec_%s.stl" % error_mm))
    assert report["watertight"], report
    assert report["valid"], report["problems"]
    assert report["components"] == 1
    assert report["genus"] == 0
    expected = 4.0 / 3.0 * math.pi * radius**3
    assert report["volume_mm3"] == pytest.approx(expected, rel=0.05)


def test_zero_simplification_error_is_a_noop(solid_sphere):
    image, _ = solid_sphere
    poly = _mesh(image, smooth_iters=5)
    assert surface.decimate(poly, 0).topology.numValidFaces() == poly.topology.numValidFaces()


def test_decimation_preserves_genus(tmp_path):
    """A skull has genus >1000. Decimation must not close its tunnels."""
    p = os.path.join(str(tmp_path), "torus.stl")
    torus_mesh = torus(10.0, 3.0, 192, 96)
    write(p, torus_mesh)
    assert validate.validate(p)["genus"] == 1

    poly = surface.from_arrays(*torus_mesh)

    smaller = surface.decimate(poly, 0.25)
    after = validate.validate(_write(smaller, tmp_path, "torus_dec.stl"))
    assert after["watertight"]
    assert after["genus"] == 1, "decimation closed the tunnel"


def test_unsafe_decimation_is_discarded(monkeypatch):
    left = sphere(subdivisions=2, radius=1.0)
    right = transformed(sphere(subdivisions=2, radius=1.0), translation=(2.2, 0.0, 0.0))
    original_verts, original_faces = concatenate(left, right)
    original = surface.from_arrays(original_verts, original_faces)

    moved = original_verts.copy()
    moved[len(left[0]):, 0] -= 0.4
    intersecting = surface.from_arrays(moved, original_faces)
    signature = (2, 0, 4)
    monkeypatch.setattr(
        surface,
        "_decimate_candidate",
        lambda *_args: (intersecting, signature, signature, 0.25),
    )

    result, stats = surface.decimate_safely(original, 0.25)

    assert result is original
    assert not stats.accepted
    assert stats.remaining_self_intersecting_faces > 0
    assert "self-intersecting" in stats.rejection_reason


def test_decimation_protects_source_patches_until_the_candidate_is_clean(monkeypatch):
    source_mesh = concatenate(
        sphere(subdivisions=3, radius=1.0),
        transformed(sphere(subdivisions=3, radius=1.0), translation=(2.2, 0.0, 0.0)),
    )
    original = surface.from_arrays(*source_mesh)

    left = sphere(subdivisions=2, radius=1.0)
    valid_mesh = concatenate(
        left,
        transformed(sphere(subdivisions=2, radius=1.0), translation=(2.2, 0.0, 0.0)),
    )
    valid = surface.from_arrays(*valid_mesh)
    moved = valid_mesh[0].copy()
    moved[len(left[0]):, 0] -= 0.4
    intersecting = surface.from_arrays(moved, valid_mesh[1])
    signature = (2, 0, 4)

    def candidate(_poly, _error, protected):
        result = valid if protected.any() else intersecting
        return result, signature, signature, 0.25

    monkeypatch.setattr(surface, "_decimate_candidate", candidate)

    result, stats = surface.decimate_safely(original, 0.25)

    assert result is valid
    assert stats.accepted
    assert stats.repair_attempts == 1
    assert stats.initial_self_intersecting_faces > 0
    assert stats.remaining_self_intersecting_faces == 0
    assert stats.protected_input_faces > 0


def _slab_survives(thickness_voxels, mm, spacing=0.2):
    """Does a slab of the given thickness still register as foreground after
    resampling onto an `mm` isotropic grid?"""
    import numpy as np
    import SimpleITK as sitk

    arr = np.zeros((40, 60, 40), dtype=np.uint8)
    y0 = 30
    arr[5:35, y0:y0 + thickness_voxels, 5:35] = 1
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((spacing, spacing, spacing))

    grid = segment.resample_isotropic(img, mm)
    a = sitk.GetArrayViewFromImage(grid)
    return bool((a > segment.ISO_OCCUPANCY).any())


def test_resampling_erases_structures_thinner_than_the_target_voxel():
    """Resampling onto a coarser grid destroys sub-voxel structure.

    The morphological closing seals a pore with a membrane one voxel thick;
    resampling to 0.6 mm blurs that membrane below the 0.5 occupancy level and the
    pore reopens. On a head CT that adds 257 tunnels (genus 1225 -> 1482) and
    terraces the vault. Hence the native grid by default.
    """
    # 2 voxels = 0.4 mm of bone.
    assert _slab_survives(2, 0.2), "a membrane must survive its own resolution"
    assert not _slab_survives(2, 1.5), "0.4 mm cannot survive a 1.5 mm grid"

    # Thick cortical bone is unaffected: 20 voxels = 4 mm.
    assert _slab_survives(20, 1.5), "4 mm of bone must survive a 1.5 mm grid"


def test_simplification_with_zero_error_is_a_noop(solid_sphere):
    image, _ = solid_sphere
    poly = _mesh(image)
    n = poly.topology.numValidFaces()
    assert surface.decimate(poly, 0).topology.numValidFaces() == n


def test_unsupported_extension_is_rejected(solid_sphere, tmp_path):
    image, _ = solid_sphere
    poly = _mesh(image)
    with pytest.raises(ValueError, match="unsupported output extension"):
        surface.write(poly, os.path.join(str(tmp_path), "mesh.xyz"))


def test_invalid_serialized_output_does_not_replace_destination(solid_sphere, tmp_path, monkeypatch):
    image, _ = solid_sphere
    destination = tmp_path / "mesh.stl"
    destination.write_bytes(b"known valid destination")
    monkeypatch.setattr(
        validate,
        "validate",
        lambda _path: {"valid": False, "problems": ["test serialization failure"]},
    )

    with pytest.raises(ValueError, match="serialized output mesh is invalid"):
        surface.write_validated(_mesh(image), str(destination))

    assert destination.read_bytes() == b"known valid destination"
    assert not list(tmp_path.glob(".mesh.stl.*"))


@pytest.mark.parametrize("ext", [".stl", ".ply", ".obj"])
def test_roundtrip_formats(solid_sphere, tmp_path, ext):
    image, _ = solid_sphere
    poly = _mesh(image, smooth_iters=5)
    path = _write(poly, tmp_path, "m" + ext)
    assert os.path.getsize(path) > 0
    report = validate.validate(path)
    assert report["triangles"] > 0
    assert report["valid"], report["problems"]


def test_vtp_is_not_supported(solid_sphere, tmp_path):
    image, _ = solid_sphere
    with pytest.raises(ValueError, match="unsupported output extension"):
        _write(_mesh(image, smooth_iters=5), tmp_path, "m.vtp")

    input_path = tmp_path / "input.vtp"
    input_path.write_text("<VTKFile />")
    with pytest.raises(ValueError, match="unsupported mesh extension"):
        validate.validate(str(input_path))


def test_disjoint_closed_shells_are_valid(tmp_path):
    """Multiple closed printable parts do not make a mesh invalid."""
    path = str(tmp_path / "two-shells.ply")
    write(path, concatenate(box(), transformed(box(), translation=(3.0, 0.0, 0.0))))

    report = validate.validate(path)
    assert report["components"] == 2
    assert report["valid"], report["problems"]


def test_transversely_overlapping_closed_shells_are_invalid(tmp_path):
    """Watertightness alone does not rule out intersecting surfaces."""
    path = str(tmp_path / "intersecting-shells.ply")
    crossing = transformed(
        box((2.0, 2.0, 2.0)),
        translation=(0.35, 0.1, 0.2),
        rotation_z_deg=30.0,
    )
    write(path, concatenate(box((2.0, 2.0, 2.0)), crossing))

    report = validate.validate(path)
    assert report["watertight"]
    assert report["self_intersecting_faces"] > 0
    assert not report["valid"]
    assert any("self-intersecting" in problem for problem in report["problems"])


def test_inconsistent_winding_fails_validation(tmp_path):
    vertices, faces = box()
    faces = faces.copy()
    faces[0] = faces[0][::-1]
    path = str(tmp_path / "flipped-face.ply")
    write(path, (vertices, faces))

    report = validate.validate(path)
    assert report["disoriented_faces"] > 0
    assert not report["winding_consistent"]
    assert not report["valid"]


def test_meshlib_native_report_omits_raw_geometry_counters(tmp_path):
    path = write(tmp_path / "box.ply", box())
    report = validate.validate(path)
    assert "degenerate_faces" not in report
    assert "nonmanifold_edge_uses" not in report
    assert report["holes"] == 0
    assert report["disoriented_faces"] == 0


def test_failed_intersection_measurement_fails_validation(tmp_path, monkeypatch):
    path = write(tmp_path / "box.stl", box())
    monkeypatch.setattr(validate, "_self_intersections", lambda _mesh: "error: Test")

    report = validate.validate(path)
    assert not report["valid"]
    assert "did not complete" in " ".join(report["problems"])
    assert "valid               NO" in validate.summarise(report)
