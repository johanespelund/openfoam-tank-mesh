"""
Tests for CylinderCapsMesh.

The geometry-comparison tests run in CI (no OpenFOAM needed).
The full generate() test is skipped in CI because it requires an OpenFOAM
installation and a valid case directory.
"""

import shutil

import pytest

from openfoam_tank_mesh.Profile import CylinderCapsTankProfile, KSiteProfile

HAS_OPENFOAM = shutil.which("blockMesh") is not None

# K-Site dimensions (Stochl & Knoll, 1991)
INCH = 0.0254
KSITE_A = 0.5 * 73 * INCH  # cap height
KSITE_B = 0.5 * 87.6 * INCH  # cylinder radius
KSITE_C = 1.5 * INCH  # cylinder height

# Common mesh parameters shared by both tests below.
_MESH_PARAMS = {
    "fill_level": 0.49,
    "wall_cell_size": 5.0e-3,
    "wall_tan_cell_size": 5.0e-3,
    "bulk_cell_size": 25e-3,
    "r_BL": 1.05,
    "tri_bulk": False,
    "outlet_radius": 0.0127,
    "internal_outlet": 0.0127 * 4,
    "debug": False,
    "revolve": 0,
    "n_revolve": 0,
    "n_wall_layers": 6,
    "VoF": False,
}

# CylinderCapsMesh additionally needs the geometry parameters.
CAPS_KSITE_PARAMS = dict(
    **_MESH_PARAMS,
    cylinder_radius=KSITE_B,
    cylinder_height=KSITE_C,
    cap_height=KSITE_A,
)


# ---------------------------------------------------------------------------
# Geometry-comparison tests (run in CI, no OpenFOAM required)
# ---------------------------------------------------------------------------


@pytest.fixture
def ksite_profile():
    return KSiteProfile(
        fill_level=_MESH_PARAMS["fill_level"],
        outlet_radius=_MESH_PARAMS["outlet_radius"],
        bulk_cell_size=_MESH_PARAMS["bulk_cell_size"],
        wall_tan_cell_size=_MESH_PARAMS["wall_tan_cell_size"],
        wall_cell_size=_MESH_PARAMS["wall_cell_size"],
        r_BL=_MESH_PARAMS["r_BL"],
        internal_outlet=_MESH_PARAMS["internal_outlet"],
    )


@pytest.fixture
def caps_ksite_profile():
    return CylinderCapsTankProfile(
        cylinder_radius=KSITE_B,
        cylinder_height=KSITE_C,
        cap_height=KSITE_A,
        fill_level=_MESH_PARAMS["fill_level"],
        outlet_radius=_MESH_PARAMS["outlet_radius"],
        bulk_cell_size=_MESH_PARAMS["bulk_cell_size"],
        wall_tan_cell_size=_MESH_PARAMS["wall_tan_cell_size"],
        wall_cell_size=_MESH_PARAMS["wall_cell_size"],
        r_BL=_MESH_PARAMS["r_BL"],
        internal_outlet=_MESH_PARAMS["internal_outlet"],
    )


def test_cylinder_caps_ksite_dimensions(caps_ksite_profile):
    """CylinderCapsTankProfile with K-Site dims stores the correct dimensions."""
    assert caps_ksite_profile.cylinder_radius == pytest.approx(KSITE_B)
    assert caps_ksite_profile.cylinder_height == pytest.approx(KSITE_C)
    assert caps_ksite_profile.cap_height == pytest.approx(KSITE_A)


def test_cylinder_caps_matches_ksite_volume(caps_ksite_profile, ksite_profile):
    """CylinderCapsTankProfile with K-Site dims should have the same volume as KSiteProfile."""
    assert caps_ksite_profile.volume == pytest.approx(ksite_profile.volume, rel=1e-4)


def test_cylinder_caps_matches_ksite_interface(caps_ksite_profile, ksite_profile):
    """Interface height should match between CylinderCaps and KSite profiles."""
    assert caps_ksite_profile.y_interface == pytest.approx(ksite_profile.y_interface, rel=1e-3)


def test_cylinder_caps_matches_ksite_liquid_volume(caps_ksite_profile, ksite_profile):
    """Liquid volume should match between CylinderCaps and KSite profiles."""
    assert caps_ksite_profile.volume_liquid == pytest.approx(ksite_profile.volume_liquid, rel=1e-3)


def test_cylinder_caps_matches_ksite_cylinder_radius(caps_ksite_profile, ksite_profile):
    """cylinder_radius attribute must match."""
    assert caps_ksite_profile.cylinder_radius == pytest.approx(ksite_profile.cylinder_radius, rel=1e-6)


@pytest.mark.skipif(not HAS_OPENFOAM, reason="OpenFOAM is not available")
def test_cylinder_caps_mesh_missing_radius_raises():
    """CylinderCapsMesh should raise ValueError when neither cylinder_radius nor cylinder_diameter is given."""
    from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh

    bad_params = dict(CAPS_KSITE_PARAMS)
    del bad_params["cylinder_radius"]

    with pytest.raises(ValueError, match="cylinder_radius.*cylinder_diameter"):
        CylinderCapsMesh(input_parameters=bad_params)


@pytest.mark.skipif(not HAS_OPENFOAM, reason="OpenFOAM is not available")
def test_mirror_without_empty_2d_raises():
    """mirror=True without empty_2d=True should raise MirrorRequiresEmpty2D (a ValueError)."""
    from openfoam_tank_mesh.exceptions import MirrorRequiresEmpty2D
    from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh

    params = dict(CAPS_KSITE_PARAMS, mirror=True)  # empty_2d defaults to False

    with pytest.raises(MirrorRequiresEmpty2D):
        CylinderCapsMesh(input_parameters=params)


@pytest.mark.skipif(not HAS_OPENFOAM, reason="OpenFOAM is not available")
def test_extrude_cylinder_without_empty_2d_raises():
    """extrude_cylinder > 0 without empty_2d=True should raise ExtrudeCylinderRequiresEmpty2D."""
    from openfoam_tank_mesh.exceptions import ExtrudeCylinderRequiresEmpty2D
    from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh

    params = dict(CAPS_KSITE_PARAMS, extrude_cylinder=0.5)  # empty_2d defaults to False

    with pytest.raises(ExtrudeCylinderRequiresEmpty2D):
        CylinderCapsMesh(input_parameters=params)


@pytest.mark.skipif(not HAS_OPENFOAM, reason="OpenFOAM is not available")
def test_extrude_cylinder_n_layers():
    """extrude_cylinder with empty_2d=True should compute n_layers = round(flow/bulk_cell_size)."""
    from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh

    bulk = _MESH_PARAMS["bulk_cell_size"]  # 0.025 m
    flow = bulk * 10  # 10 layers expected
    params = dict(CAPS_KSITE_PARAMS, empty_2d=True, extrude_cylinder=flow)

    mesh = CylinderCapsMesh(input_parameters=params)
    n_layers = max(1, round(mesh.extrude_cylinder / mesh.bulk_cell_size))
    assert n_layers == 10


# ---------------------------------------------------------------------------
# Full mesh-generation tests (skipped in CI - require OpenFOAM)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_OPENFOAM, reason="OpenFOAM is not available")
def test_cylinder_caps_mesh():
    """CylinderCapsMesh with K-Site dimensions should generate a valid mesh."""
    from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh

    mesh = CylinderCapsMesh(input_parameters=dict(CAPS_KSITE_PARAMS))
    mesh.generate()
    assert mesh.tank.get_radius(mesh.tank.y_start) == pytest.approx(0.0, abs=1e-9)
    assert mesh.tank.cylinder_radius == pytest.approx(KSITE_B)


@pytest.mark.skipif(not HAS_OPENFOAM, reason="OpenFOAM is not available")
def test_cylinder_caps_mesh_matches_ksite_mesh():
    """CylinderCapsMesh with K-Site dims should produce the same mesh as KSiteMesh."""
    from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh, KSiteMesh

    ksite_mesh = KSiteMesh(input_parameters=dict(_MESH_PARAMS))
    caps_mesh = CylinderCapsMesh(input_parameters=dict(CAPS_KSITE_PARAMS))

    # Geometry must be identical.
    assert caps_mesh.tank.cylinder_radius == pytest.approx(ksite_mesh.tank.cylinder_radius, rel=1e-6)
    assert caps_mesh.tank.cylinder_height == pytest.approx(ksite_mesh.tank.cylinder_height, rel=1e-6)
    assert caps_mesh.tank.cap_height == pytest.approx(ksite_mesh.tank.cap_height, rel=1e-6)
    assert caps_mesh.tank.volume == pytest.approx(ksite_mesh.tank.volume, rel=1e-4)
    assert caps_mesh.tank.y_interface == pytest.approx(ksite_mesh.tank.y_interface, rel=1e-3)

    # Both meshes should generate without errors.
    ksite_mesh.generate()
    caps_mesh.generate()


if __name__ == "__main__":
    test_cylinder_caps_mesh()
    print("test_cylinder_caps_mesh passed")
