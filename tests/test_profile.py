import numpy as np
import pytest

from openfoam_tank_mesh.Profile import CylinderCapsTankProfile, CylinderTankProfile, KSiteProfile, SphereProfile

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


# ---------------------------------------------------------------------------
# CylinderCapsTankProfile
# ---------------------------------------------------------------------------

CYL_RADIUS = 0.5
CYL_HEIGHT = 1.0
CAP_HEIGHT = 0.4


@pytest.fixture
def cylinder_caps_profile():
    return CylinderCapsTankProfile(
        cylinder_radius=CYL_RADIUS,
        cylinder_height=CYL_HEIGHT,
        cap_height=CAP_HEIGHT,
        fill_level=0.5,
        outlet_radius=0.03,
        bulk_cell_size=0.05,
        wall_tan_cell_size=0.01,
        wall_cell_size=0.005,
        r_BL=1.1,
    )


@pytest.fixture
def cylinder_caps_profile_via_diameter():
    """Same geometry as cylinder_caps_profile but specified with diameter."""
    return CylinderCapsTankProfile(
        cylinder_radius=0.0,  # will be overridden by cylinder_diameter
        cylinder_height=CYL_HEIGHT,
        cap_height=CAP_HEIGHT,
        fill_level=0.5,
        outlet_radius=0.03,
        bulk_cell_size=0.05,
        wall_tan_cell_size=0.01,
        wall_cell_size=0.005,
        r_BL=1.1,
        cylinder_diameter=CYL_RADIUS * 2,
    )


@pytest.fixture
def cylinder_caps_profile_zero_cylinder():
    """Degenerate case: cylinder_height = 0 gives a purely spheroidal tank."""
    return CylinderCapsTankProfile(
        cylinder_radius=0.6,
        cylinder_height=0.0,
        cap_height=0.6,
        fill_level=0.5,
        outlet_radius=0.03,
        bulk_cell_size=0.05,
        wall_tan_cell_size=0.01,
        wall_cell_size=0.005,
        r_BL=1.1,
    )


def test_cylinder_caps_profile_dimensions(cylinder_caps_profile):
    assert cylinder_caps_profile.cylinder_radius == pytest.approx(CYL_RADIUS)
    assert cylinder_caps_profile.cylinder_height == pytest.approx(CYL_HEIGHT)
    assert cylinder_caps_profile.cap_height == pytest.approx(CAP_HEIGHT)


def test_cylinder_caps_profile_total_height(cylinder_caps_profile):
    """Total height (y_end attribute) should be 2 * cap_height + cylinder_height."""
    expected_height = 2 * CAP_HEIGHT + CYL_HEIGHT
    assert cylinder_caps_profile.y_end == pytest.approx(expected_height)


def test_cylinder_caps_profile_radius_at_poles(cylinder_caps_profile):
    """Radius at the bottom pole must be zero."""
    assert cylinder_caps_profile.get_radius(cylinder_caps_profile.y_start) == pytest.approx(0.0, abs=1e-9)


def test_cylinder_caps_profile_radius_at_equator(cylinder_caps_profile):
    """Radius at cap-cylinder junction must equal cylinder_radius."""
    y_junction = CAP_HEIGHT
    assert cylinder_caps_profile.get_radius(y_junction) == pytest.approx(CYL_RADIUS, rel=1e-6)


def test_cylinder_caps_profile_fill_level(cylinder_caps_profile):
    """Liquid volume / total volume should match the specified fill level."""
    assert cylinder_caps_profile.volume_liquid / cylinder_caps_profile.volume == pytest.approx(0.5, rel=1e-2)


def test_cylinder_caps_profile_volume_conservation(cylinder_caps_profile):
    """Liquid + gas volumes must sum to the total volume."""
    assert cylinder_caps_profile.volume_liquid + cylinder_caps_profile.volume_gas == pytest.approx(
        cylinder_caps_profile.volume, rel=1e-3
    )


def test_cylinder_caps_profile_interface_inside_tank(cylinder_caps_profile):
    """Interface must lie inside the tank."""
    assert cylinder_caps_profile.y_start < cylinder_caps_profile.y_interface < cylinder_caps_profile.y_end


def test_cylinder_caps_profile_diameter_equivalent(cylinder_caps_profile, cylinder_caps_profile_via_diameter):
    """Specifying cylinder_diameter should give the same geometry as cylinder_radius."""
    assert cylinder_caps_profile_via_diameter.cylinder_radius == pytest.approx(CYL_RADIUS)
    assert cylinder_caps_profile_via_diameter.volume == pytest.approx(cylinder_caps_profile.volume, rel=1e-6)


def test_cylinder_caps_profile_zero_cylinder_height(cylinder_caps_profile_zero_cylinder):
    """A zero cylinder_height tank should still build and conserve volume."""
    p = cylinder_caps_profile_zero_cylinder
    assert p.cylinder_height == pytest.approx(0.0)
    assert p.volume_liquid + p.volume_gas == pytest.approx(p.volume, rel=1e-3)


def test_cylinder_caps_profile_wall_thickness_parameter():
    """Wall-point offset should match the configured wall thickness."""
    wall_thickness = 4.0e-3
    p = CylinderCapsTankProfile(
        cylinder_radius=CYL_RADIUS,
        cylinder_height=CYL_HEIGHT,
        cap_height=CAP_HEIGHT,
        fill_level=0.5,
        outlet_radius=0.03,
        bulk_cell_size=0.05,
        wall_tan_cell_size=0.01,
        wall_cell_size=0.005,
        r_BL=1.1,
        wall_thickness=wall_thickness,
    )
    pts = p.get_mesh_points()
    offset = np.linalg.norm(pts.outer_points[0] - pts.wall_points[0])
    assert offset == pytest.approx(wall_thickness, rel=1e-6)


@pytest.fixture
def cylinder_profile():
    return CylinderTankProfile(
        cylinder_radius=0.5,
        cylinder_height=1.2,
        fill_level=0.4,
        outlet_radius=0.05,
        bulk_cell_size=0.05,
        wall_tan_cell_size=0.01,
        wall_cell_size=0.005,
        r_BL=1.1,
    )


def test_cylinder_profile_dimensions(cylinder_profile):
    assert cylinder_profile.cylinder_radius == pytest.approx(0.5)
    assert cylinder_profile.cylinder_height == pytest.approx(1.2)
    assert cylinder_profile.cap_height == pytest.approx(0.0)


def test_cylinder_profile_flat_endpoints_and_sidewall(cylinder_profile):
    eps = 1e-9
    outer_points, _ = cylinder_profile.get_profile_points()
    assert cylinder_profile.get_radius(cylinder_profile.y_start) == pytest.approx(0.0, abs=1e-9)
    assert outer_points[-1][0] == pytest.approx(cylinder_profile.outlet_radius)
    assert cylinder_profile.get_radius(cylinder_profile.y_end - eps) == pytest.approx(cylinder_profile.cylinder_radius)
    assert cylinder_profile.get_radius(0.5 * cylinder_profile.y_end) == pytest.approx(cylinder_profile.cylinder_radius)


def test_cylinder_profile_volume(cylinder_profile):
    expected = np.pi * cylinder_profile.cylinder_radius**2 * cylinder_profile.cylinder_height
    assert cylinder_profile.volume == pytest.approx(expected, rel=1e-6)


def test_cylinder_profile_mesh_points_include_wall(cylinder_profile):
    pts = cylinder_profile.get_mesh_points()
    assert len(pts.outer_points) >= 4
    assert len(pts.wall_points) == len(pts.outer_points)
    assert np.isfinite(np.array(pts.outer_points)).all()
