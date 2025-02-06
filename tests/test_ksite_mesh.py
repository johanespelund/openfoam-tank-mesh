import os

import pytest

from openfoam_tank_mesh.KSiteMesh import KSiteMesh


@pytest.mark.skipif("CI" in os.environ, reason="OpenFOAM is not available in CI")
def test_ksite_mesh():
    mesh = KSiteMesh(
        input_parameters={
            "fill_level": 0.49,
            "wall_cell_size": 24e-3,
            "bulk_cell_size": 48e-3,
            "outlet_radius": 0.024,
            "debug": False,
            "revolve": 90,
        }
    )
    # mesh.cfMesh(nLayers=3)
    # mesh.add_wall(wall_thickness=0.002, n_layers=10)
    assert mesh.tank.get_radius(0) == mesh.tank.cylinder_radius
    assert mesh.tank.get_radius(mesh.tank.y2) == 0
    assert mesh.tank.get_radius(mesh.tank.y1) == 0
    return None


if __name__ == "__main__":
    test_ksite_mesh()
    print("test_ksite_mesh passed")
    pass
