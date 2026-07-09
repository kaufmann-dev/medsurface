"""Surface extraction: capping, cavity removal, coordinate mapping, volume fidelity."""

import math
import os

import numpy as np
import pytest
import SimpleITK as sitk

from dicom_surface import segment, surface, validate


def _mesh(image, cap=True, smooth_iters=0, largest=False):
    if cap:
        image = segment.pad(image, 1)
    affine = surface.index_to_physical(image)
    poly = surface.marching_cubes(surface.to_vtk_image(image))
    if smooth_iters:
        poly = surface.smooth(poly, smooth_iters, 0.1)
    if largest:
        poly, _n = surface.largest_component(poly)
    return surface.transform(poly, affine)


def _write(poly, tmp_path, name="m.stl"):
    p = os.path.join(str(tmp_path), name)
    surface.write(surface.compute_normals(poly), p)
    return p


def test_solid_sphere_is_watertight_with_correct_volume(solid_sphere, tmp_path):
    image, radius = solid_sphere
    poly = _mesh(image, smooth_iters=10)
    report = validate.validate(_write(poly, tmp_path), self_intersections=False)

    assert report["watertight"]
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
        self_intersections=False,
    )
    assert not open_report["watertight"]
    assert open_report["boundary_edges"] > 0

    capped_report = validate.validate(
        _write(_mesh(clipped_sphere, cap=True), tmp_path, "capped.stl"),
        self_intersections=False,
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
    report_all = validate.validate(_write(all_shells, tmp_path, "shells.stl"),
                                   self_intersections=False)
    assert report_all["components"] == 2

    kept = _mesh(image, largest=True)
    report_one = validate.validate(_write(kept, tmp_path, "outer.stl"),
                                   self_intersections=False)
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
        got = [
            sum(m.GetElement(r, c) * idx[c] for c in range(3)) + m.GetElement(r, 3)
            for r in range(3)
        ]
        assert got == pytest.approx(list(expected), abs=1e-9)


def test_padding_preserves_physical_coordinates(solid_sphere):
    image, _ = solid_sphere
    padded = segment.pad(image, 1)
    # index (0,0,0) of the original is index (1,1,1) of the padded image
    assert padded.TransformIndexToPhysicalPoint((1, 1, 1)) == pytest.approx(
        image.TransformIndexToPhysicalPoint((0, 0, 0))
    )
    assert padded.GetSize() == tuple(s + 2 for s in image.GetSize())


def test_decimation_hits_the_triangle_budget(solid_sphere, tmp_path):
    image, _ = solid_sphere
    poly = _mesh(image, smooth_iters=5)
    before = poly.GetNumberOfPolys()
    target = before // 4
    smaller = surface.decimate(poly, target)
    assert smaller.GetNumberOfPolys() < before
    assert smaller.GetNumberOfPolys() == pytest.approx(target, rel=0.25)


def test_smoothing_trades_roughness_for_displacement(solid_sphere):
    """More iterations means a smoother surface that sits further from the data.

    Guards the one knob users actually feel. Windowed-sinc must not shrink the
    surface either -- a plain Laplacian would contract it toward the centroid.
    """
    import numpy as np
    import trimesh

    image, radius = solid_sphere
    padded = segment.pad(image, 1)
    affine = surface.index_to_physical(padded)
    raw = surface.transform(surface.marching_cubes(surface.to_vtk_image(padded)), affine)

    ref = trimesh.Trimesh(*surface.to_arrays(raw), process=True)
    exact = 4.0 / 3.0 * math.pi * radius**3

    roughness = []
    for iterations in (0, 10, 25, 40):
        poly = surface.smooth(raw, iterations, 0.1) if iterations else raw
        m = trimesh.Trimesh(*surface.to_arrays(poly), process=True)
        roughness.append(float(np.degrees(m.face_adjacency_angles).mean()))
        # no shrinkage: volume stays within 2% of the analytic sphere
        assert m.volume == pytest.approx(exact, rel=0.02), iterations

    assert roughness == sorted(roughness, reverse=True), roughness
    assert roughness[0] > roughness[-1]


def test_count_defects_matches_validate(clipped_sphere, solid_sphere):
    closed = _mesh(solid_sphere[0])
    assert surface.count_defects(closed) == (0, 0)

    opened = _mesh(clipped_sphere, cap=False)
    boundary, nonmanifold = surface.count_defects(opened)
    assert boundary > 0
    assert nonmanifold == 0


@pytest.mark.parametrize("mm", [0.6, 1.0, 1.5])
def test_isotropic_resampling_stays_watertight(solid_sphere, tmp_path, mm):
    """Resampling is the topology-safe way to shed triangles: marching cubes
    always returns a manifold surface, so no defect can be introduced."""
    image, radius = solid_sphere
    padded = segment.pad(image, 1)
    grid = segment.resample_isotropic(padded, mm)

    poly = surface.marching_cubes(surface.to_vtk_image(grid), segment.ISO_OCCUPANCY)
    poly = surface.smooth(poly, 10, 0.1)
    poly = surface.transform(poly, surface.index_to_physical(grid))

    assert surface.count_defects(poly) == (0, 0)
    report = validate.validate(_write(poly, tmp_path, "iso%s.stl" % mm),
                               self_intersections=False)
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
        poly = surface.marching_cubes(surface.to_vtk_image(grid), segment.ISO_OCCUPANCY)
        counts.append(poly.GetNumberOfPolys())
    assert counts[0] > counts[1] > counts[2]


def test_resampling_preserves_physical_placement(solid_sphere):
    image, _ = solid_sphere
    padded = segment.pad(image, 1)
    grid = segment.resample_isotropic(padded, 1.0)
    assert grid.GetSpacing() == pytest.approx((1.0, 1.0, 1.0))
    # the resampled grid must fully contain the original extent
    for axis in range(3):
        assert grid.GetSize()[axis] * 1.0 >= padded.GetSize()[axis] * padded.GetSpacing()[axis]


@pytest.mark.parametrize("keep", [0.5, 0.2, 0.05])
def test_decimation_never_opens_a_closed_mesh(solid_sphere, tmp_path, keep):
    """Regression: both of VTK's decimators tear a closed surface.

    vtkQuadricDecimation leaks boundary edges at every reduction (6 at 50%, 90 at
    85% on a real skull); vtkDecimatePro with PreserveTopologyOn misses the target
    and still emits non-manifold edges. The MeshLab decimator holds topology, so a
    watertight mesh stays watertight and its volume barely moves.
    """
    image, radius = solid_sphere
    poly = _mesh(image, smooth_iters=10)
    assert surface.count_defects(poly) == (0, 0)

    target = max(64, int(poly.GetNumberOfPolys() * keep))
    smaller = surface.decimate(poly, target)

    assert surface.count_defects(smaller) == (0, 0)
    report = validate.validate(_write(smaller, tmp_path, "dec_%s.stl" % keep),
                               self_intersections=False)
    assert report["watertight"], report
    assert report["components"] == 1
    assert report["genus"] == 0
    expected = 4.0 / 3.0 * math.pi * radius**3
    assert report["volume_mm3"] == pytest.approx(expected, rel=0.05)


def test_decimation_hits_the_target(solid_sphere):
    """Within a triangle or two: a closed surface has an even face count, so an
    odd target cannot be hit exactly."""
    image, _ = solid_sphere
    poly = _mesh(image, smooth_iters=5)
    target = poly.GetNumberOfPolys() // 3
    assert surface.decimate(poly, target).GetNumberOfPolys() == pytest.approx(target, abs=2)


def test_decimation_preserves_genus(tmp_path):
    """A skull has genus >1000. Decimation must not close its tunnels."""
    import trimesh
    import vtk as _vtk

    torus = trimesh.creation.torus(major_radius=10.0, minor_radius=3.0,
                                   major_sections=192, minor_sections=96)
    p = os.path.join(str(tmp_path), "torus.stl")
    torus.export(p)
    assert validate.validate(p, self_intersections=False)["genus"] == 1

    reader = _vtk.vtkSTLReader()
    reader.SetFileName(p)
    reader.Update()
    poly = reader.GetOutput()

    smaller = surface.decimate(poly, max(64, poly.GetNumberOfPolys() // 6))
    after = validate.validate(_write(smaller, tmp_path, "torus_dec.stl"),
                              self_intersections=False)
    assert after["watertight"]
    assert after["genus"] == 1, "decimation closed the tunnel"


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
    """Regression for the quality bug that shipped.

    Resampling onto a coarser grid destroys sub-voxel structure. The morphological
    closing seals a pore with a membrane one voxel thick; resampling to 0.6 mm
    blurs that membrane below the 0.5 occupancy level and the pore reopens. On a
    head CT this added 257 tunnels (genus 1225 -> 1482) and terraced the vault.
    """
    # 2 voxels = 0.4 mm of bone.
    assert _slab_survives(2, 0.2), "a membrane must survive its own resolution"
    assert not _slab_survives(2, 1.5), "0.4 mm cannot survive a 1.5 mm grid"

    # Thick cortical bone is unaffected: 20 voxels = 4 mm.
    assert _slab_survives(20, 1.5), "4 mm of bone must survive a 1.5 mm grid"


def test_decimation_is_a_noop_when_already_under_budget(solid_sphere):
    image, _ = solid_sphere
    poly = _mesh(image)
    n = poly.GetNumberOfPolys()
    assert surface.decimate(poly, n * 10).GetNumberOfPolys() == n
    assert surface.decimate(poly, 0).GetNumberOfPolys() == n


def test_unsupported_extension_is_rejected(solid_sphere, tmp_path):
    image, _ = solid_sphere
    poly = _mesh(image)
    with pytest.raises(ValueError, match="unsupported output extension"):
        surface.write(poly, os.path.join(str(tmp_path), "mesh.xyz"))


@pytest.mark.parametrize("ext", [".stl", ".ply", ".obj"])
def test_roundtrip_formats(solid_sphere, tmp_path, ext):
    image, _ = solid_sphere
    poly = _mesh(image, smooth_iters=5)
    path = _write(poly, tmp_path, "m" + ext)
    assert os.path.getsize(path) > 0
    report = validate.validate(path, self_intersections=False)
    assert report["triangles"] > 0
