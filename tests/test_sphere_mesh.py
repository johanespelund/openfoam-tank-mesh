import os

import pytest

from openfoam_tank_mesh.TwoPhaseMesh import SphereMesh


@pytest.mark.skipif("CI" in os.environ, reason="OpenFOAM is not available in CI")
def test_ksite_mesh():
    mesh = SphereMesh(
        input_parameters={
            "radius": 1.4,
            "fill_level": 0.1,
            "wall_cell_size": 5.0e-3,
            "wall_tan_cell_size": 25.0e-3,
            "bulk_cell_size": 25.0e-3,
            "r_BL": 1.2,
            "tri_bulk": False,
            "outlet_radius": 0.07,
            "internal_outlet": 0,
            "debug": False,
            "revolve": 0,
            "n_revolve": 0,
            "n_wall_layers": 6,
        }
    )
    mesh.generate()
    mesh.remove_wall()
    return None


if __name__ == "__main__":
    test_ksite_mesh()
    print("test_ksite_mesh passed")
    pass
