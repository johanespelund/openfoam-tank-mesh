import numpy as np
import pytest

from openfoam_tank_mesh.Profile import KSiteProfile, SphereProfile

INCH = 0.0254
KSITE_A = 0.5 * 73 * INCH  # cap height
KSITE_B = 0.5 * 87.6 * INCH  # cylinder radius
KSITE_C = 1.5 * INCH  # cylinder height


@pytest.fixture
def ksite_profile():
    return KSiteProfile(
        fill_level=0.5,
        outlet_radius=0.0127,
        bulk_cell_size=25e-3,
        wall_tan_cell_size=5.0e-3,
        wall_cell_size=5.0e-3,
        r_BL=1.05,
    )


@pytest.fixture
def sphere_profile():
    return SphereProfile(
        radius=1.0,
        fill_level=0.5,
        outlet_radius=0.05,
        bulk_cell_size=0.1,
        wall_tan_cell_size=0.05,
        wall_cell_size=0.01,
        r_BL=1.1,
    )


def test_ksite_profile_dimensions(ksite_profile):
    assert ksite_profile.cylinder_radius == pytest.approx(KSITE_B)
    assert ksite_profile.cylinder_height == pytest.approx(KSITE_C)
    assert ksite_profile.cap_height == pytest.approx(KSITE_A)


def test_ksite_profile_radius_at_bottom(ksite_profile):
    """The radius at the bottom pole should be zero."""
    assert ksite_profile.get_radius(ksite_profile.y_start) == pytest.approx(0.0, abs=1e-9)


def test_ksite_profile_interface_inside_tank(ksite_profile):
    """At 50% fill, the interface must lie inside the tank."""
    assert ksite_profile.y_start < ksite_profile.y_interface < ksite_profile.y_end


def test_ksite_profile_fill_level(ksite_profile):
    """Liquid volume / total volume should match the specified fill level."""
    assert ksite_profile.volume_liquid / ksite_profile.volume == pytest.approx(0.5, rel=1e-2)


def test_ksite_profile_volume_conservation(ksite_profile):
    """Liquid + gas volumes must sum to the total volume."""
    assert ksite_profile.volume_liquid + ksite_profile.volume_gas == pytest.approx(ksite_profile.volume, rel=1e-3)


def test_sphere_profile_dimensions(sphere_profile):
    assert sphere_profile.cylinder_radius == pytest.approx(1.0)
    assert sphere_profile.cylinder_height == pytest.approx(0.0)
    assert sphere_profile.cap_height == pytest.approx(1.0)


def test_sphere_profile_radius_at_poles(sphere_profile):
    """The radius at the bottom pole of a sphere should be zero."""
    assert sphere_profile.get_radius(sphere_profile.y_start) == pytest.approx(0.0, abs=1e-9)


def test_sphere_profile_fill_level(sphere_profile):
    """Liquid volume / total volume should match the specified fill level."""
    assert sphere_profile.volume_liquid / sphere_profile.volume == pytest.approx(0.5, rel=1e-2)


def test_sphere_profile_interface_at_equator(sphere_profile):
    """At 50% fill, interface should be at the equator (y = radius = 1.0)."""
    assert sphere_profile.y_interface == pytest.approx(1.0, abs=0.02)


def test_sphere_profile_volume_conservation(sphere_profile):
    """Liquid + gas volumes must sum to the total volume."""
    assert np.isclose(sphere_profile.volume_liquid + sphere_profile.volume_gas, sphere_profile.volume, rtol=1e-3)
