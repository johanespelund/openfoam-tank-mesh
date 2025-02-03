import os

import pytest

from openfoam_tank_mesh.NASA1m3Mesh import NASA1m3Mesh


@pytest.mark.skipif("CI" in os.environ, reason="OpenFOAM is not available in CI")
def test_mesh():
    mesh = NASA1m3Mesh(
        input_parameters={
            "fill_level": 0.81,
            "wall_cell_size": 4.0e-3,
            "bulk_cell_size": 10e-3,
            "outlet_radius": 0.02,
            "debug": True,
            "revolve": 30,
        }
    )
    mesh.generate()
    assert mesh.tank.get_radius(0) == mesh.tank.cylinder_radius
    assert mesh.tank.get_radius(mesh.tank.y2) == 0
    assert mesh.tank.get_radius(mesh.tank.y1) == 0
    return None


if __name__ == "__main__":
    test_mesh()
    print("test_ksite_mesh passed")
    pass
