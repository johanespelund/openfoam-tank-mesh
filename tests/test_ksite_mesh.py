import os

import pytest

from openfoam_tank_mesh.KSiteMesh import KSiteMesh


@pytest.mark.skipif("CI" in os.environ, reason="OpenFOAM is not available in CI")
def test_ksite_mesh():
    mesh = KSiteMesh(
        input_parameters={
            "fill_level": 0.49,
            "wall_cell_size": 3e-3,
            "bulk_cell_size": 9e-3,
            "outlet_radius": 0.01,
            "debug": False,
            "revolve": 15,
        }
    )
    assert mesh.tank.fill_level == 0.49
    assert mesh.tank.outlet_radius == 0.01
    assert mesh.tank.get_radius(0) == mesh.tank.cylinder_radius
    assert mesh.tank.get_radius(mesh.tank.y2) == 0
    assert mesh.tank.get_radius(mesh.tank.y1) == 0
    return None


if __name__ == "__main__":
    test_ksite_mesh()
    print("test_ksite_mesh passed")
    pass
