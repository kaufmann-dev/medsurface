"""Kernel sizing. Every assertion here guards a silent data-quality bug."""

import math

import pytest

from dicom_surface.geometry import (
    is_anisotropic,
    kernel_extent_mm,
    kernel_radius_voxels,
    mm3_to_voxels,
    voxel_volume_mm3,
)

# The real spacing of the head CT this tool was developed against.
HEAD_CT = (0.315, 0.315, 0.8)


def test_realised_kernel_never_exceeds_request():
    """The core contract: for any active axis, (2r+1)*spacing <= requested mm.

    A radius of 0 is the identity filter on that axis and is exempt -- a kernel
    can never be narrower than a single voxel.
    """
    for mm in (0.5, 0.8, 1.0, 1.5, 2.0, 2.4, 3.2, 5.0):
        for spacing in (HEAD_CT, (1.0, 1.0, 1.0), (0.4, 0.4, 3.0), (0.7, 0.7, 0.7)):
            radii = kernel_radius_voxels(mm, spacing)
            for r, s, extent in zip(radii, spacing, kernel_extent_mm(radii, spacing)):
                if r > 0:
                    assert extent <= mm + 1e-9, (mm, spacing, radii, extent)
                else:
                    assert extent == pytest.approx(s)


def test_exact_fit_survives_binary_float_representation():
    """Regression: 2.4 / 0.8 == 2.9999999999999996, so a naive floor yields r=0
    and silently disables the filter along z."""
    assert 2.4 / 0.8 != 3.0  # the trap itself
    assert kernel_radius_voxels(2.4, (0.8, 0.8, 0.8)) == [1, 1, 1]
    assert kernel_radius_voxels(0.9, (0.3, 0.3, 0.3)) == [1, 1, 1]
    assert kernel_radius_voxels(1.5, (0.5, 0.5, 0.5)) == [1, 1, 1]


def test_epsilon_does_not_round_up_a_genuine_near_miss():
    """The guard must not turn a kernel that genuinely does not fit into one
    that does: 2.39 mm cannot hold three 0.8 mm voxels."""
    assert kernel_radius_voxels(2.39, (0.8, 0.8, 0.8)) == [0, 0, 0]


def test_regression_round_would_overshoot_in_z():
    """A 1 mm median against 0.8 mm slices.

    round(0.5/0.8) == 1 gives a 3-voxel kernel spanning 2.40 mm -- 2.4x the
    request -- which erases real anatomy while appearing to honour the parameter.
    Flooring the *extent* keeps z to a single voxel.
    """
    radii = kernel_radius_voxels(1.0, HEAD_CT)
    assert radii == [1, 1, 0]

    naive_round = [max(0, int(round((1.0 / 2.0) / s))) for s in HEAD_CT]
    assert naive_round == [2, 2, 1]  # what the buggy version produced
    assert radii != naive_round

    extents = kernel_extent_mm(radii, HEAD_CT)
    assert extents == pytest.approx([0.945, 0.945, 0.8])
    assert max(extents) <= 1.0


def test_closing_kernel_matches_validated_pipeline():
    """2.4 mm closing on the head CT reproduces the 7x7x3 kernel that was verified."""
    assert kernel_radius_voxels(2.4, HEAD_CT) == [3, 3, 1]
    extents = kernel_extent_mm([3, 3, 1], HEAD_CT)
    assert extents == pytest.approx([2.205, 2.205, 2.4])


def test_zero_and_negative_mm_disable_the_kernel():
    assert kernel_radius_voxels(0.0, HEAD_CT) == [0, 0, 0]
    assert kernel_radius_voxels(-1.0, HEAD_CT) == [0, 0, 0]


def test_kernel_collapses_on_thick_slices():
    """A 1 mm kernel cannot reach across 3 mm slices; z radius must be 0, not 1."""
    assert kernel_radius_voxels(1.0, (0.4, 0.4, 3.0)) == [0, 0, 0]


def test_radius_grows_monotonically_with_request():
    prev = -1
    for mm in [x / 10 for x in range(1, 100)]:
        r = kernel_radius_voxels(mm, (0.5, 0.5, 0.5))[0]
        assert r >= prev
        prev = r


def test_rejects_nonpositive_spacing():
    with pytest.raises(ValueError):
        kernel_radius_voxels(1.0, (0.5, 0.0, 0.5))


def test_voxel_volume_and_island_conversion():
    assert voxel_volume_mm3(HEAD_CT) == pytest.approx(0.315 * 0.315 * 0.8)
    # 50 mm3 of bone at head-CT resolution
    assert mm3_to_voxels(50.0, HEAD_CT) == pytest.approx(
        round(50.0 / voxel_volume_mm3(HEAD_CT)), abs=1
    )
    assert mm3_to_voxels(0.0, HEAD_CT) == 0
    # never rounds a positive request down to "no filtering at all"
    assert mm3_to_voxels(1e-6, HEAD_CT) == 1


def test_anisotropy_detection():
    assert is_anisotropic(HEAD_CT)
    assert not is_anisotropic((0.5, 0.5, 0.5))
    assert not math.isnan(voxel_volume_mm3(HEAD_CT))
