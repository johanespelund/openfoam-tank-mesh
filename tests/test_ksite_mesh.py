import os

import pytest

from openfoam_tank_mesh.TwoPhaseMesh import KSiteMesh #, MHTBMesh


@pytest.mark.skipif("CI" in os.environ, reason="OpenFOAM is not available in CI")
def test_ksite_mesh():
    mesh = KSiteMesh(
    # mesh = MHTBMesh(
        input_parameters={
            "fill_level": 0.49,
            "wall_cell_size": 5.0e-3,
            "wall_tan_cell_size": 5.0e-3,
            "bulk_cell_size": 25e-3,
            "r_BL": 1.05,
            "tri_bulk": False,
            "outlet_radius": 0.0127,
            # "outlet_radius": 0.0127*2,
            "internal_outlet": 0.0127*4,
            "debug": True,
            "revolve": 0,
            "n_revolve": 0,
            "n_wall_layers": 6
        }
    )
    mesh.generate()
    # mesh.tank.write_sample_lines([15, 30, 45, 60, 75], 50e-3)
    # mesh.cfMesh(nLayers=0)
    # mesh.add_wall(wall_thickness=0.002, n_layers=10)
    # assert mesh.tank.get_radius(0) == mesh.tank.cylinder_radius
    # assert mesh.tank.get_radius(mesh.tank.y2) == 0
    # assert mesh.tank.get_radius(mesh.tank.y1) == 0
    return None


if __name__ == "__main__":
    test_ksite_mesh()
    print("test_ksite_mesh passed")
    pass
