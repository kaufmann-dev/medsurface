import numpy as np
import pytest
import SimpleITK as sitk


def sphere_array(shape, center, radius, spacing):
    """Binary sphere on an anisotropic grid, indexed (z, y, x)."""
    zz, yy, xx = np.meshgrid(
        np.arange(shape[0]) * spacing[2],
        np.arange(shape[1]) * spacing[1],
        np.arange(shape[2]) * spacing[0],
        indexing="ij",
    )
    cz, cy, cx = center
    d2 = (zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2
    return (d2 <= radius**2).astype(np.uint8)


def as_image(arr, spacing):
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(tuple(float(s) for s in spacing))
    return img


@pytest.fixture
def spacing():
    return (0.5, 0.5, 1.0)


@pytest.fixture
def solid_sphere(spacing):
    """A 20 mm-radius sphere with comfortable background margin."""
    shape = (60, 120, 120)  # z, y, x  -> 60mm x 60mm x 60mm
    arr = sphere_array(shape, center=(30.0, 30.0, 30.0), radius=20.0, spacing=spacing)
    return as_image(arr, spacing), 20.0


@pytest.fixture
def hollow_sphere(spacing):
    """20 mm sphere with a 10 mm concentric cavity: two disconnected shells."""
    shape = (60, 120, 120)
    outer = sphere_array(shape, (30.0, 30.0, 30.0), 20.0, spacing)
    inner = sphere_array(shape, (30.0, 30.0, 30.0), 10.0, spacing)
    return as_image((outer & ~inner).astype(np.uint8), spacing), 20.0, 10.0


@pytest.fixture
def clipped_sphere(spacing):
    """A sphere cut off by the volume boundary, like FOV-truncated anatomy."""
    shape = (30, 120, 120)  # only 30 mm of z: the sphere runs off the top
    arr = sphere_array(shape, (20.0, 30.0, 30.0), 20.0, spacing)
    img = as_image(arr, spacing)
    return img
