import shutil

import pytest

from openfoam_tank_mesh.Profile import CylinderTankProfile

HAS_OPENFOAM = shutil.which("blockMesh") is not None


def test_cylinder_profile_diameter_equivalent():
    by_radius = CylinderTankProfile(
        cylinder_radius=0.4,
        cylinder_height=1.0,
        fill_level=0.5,
        outlet_radius=0.04,
        bulk_cell_size=0.04,
        wall_tan_cell_size=0.01,
        wall_cell_size=0.005,
        r_BL=1.05,
    )
    by_diameter = CylinderTankProfile(
        cylinder_radius=0.0,
        cylinder_height=1.0,
        fill_level=0.5,
        outlet_radius=0.04,
        bulk_cell_size=0.04,
        wall_tan_cell_size=0.01,
        wall_cell_size=0.005,
        r_BL=1.05,
        cylinder_diameter=0.8,
    )
    assert by_diameter.cylinder_radius == pytest.approx(by_radius.cylinder_radius)
    assert by_diameter.volume == pytest.approx(by_radius.volume, rel=1e-9)


def test_cylinder_profile_outlet_validation():
    with pytest.raises(ValueError, match="outlet_radius"):
        CylinderTankProfile(
            cylinder_radius=0.4,
            cylinder_height=1.0,
            fill_level=0.5,
            outlet_radius=0.4,
            bulk_cell_size=0.04,
            wall_tan_cell_size=0.01,
            wall_cell_size=0.005,
            r_BL=1.05,
        )


@pytest.mark.skipif(not HAS_OPENFOAM, reason="OpenFOAM is not available")
def test_cylinder_mesh_missing_radius_raises():
    from openfoam_tank_mesh.mesh_builders import CylinderMesh

    params = {
        "cylinder_height": 1.0,
        "fill_level": 0.5,
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
    with pytest.raises(ValueError, match="cylinder_radius.*cylinder_diameter"):
        CylinderMesh(input_parameters=params)
