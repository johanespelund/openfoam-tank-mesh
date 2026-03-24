import os

import pytest


@pytest.mark.skipif("CI" in os.environ, reason="OpenFOAM is not available in CI")
def test_ksite_mesh():
    from openfoam_tank_mesh.TwoPhaseMesh import KSiteMesh

    mesh = KSiteMesh(
        input_parameters={
            "fill_level": 0.49,
            "wall_cell_size": 5.0e-3,
            "wall_tan_cell_size": 5.0e-3,
            "bulk_cell_size": 25e-3,
            "r_BL": 1.05,
            "tri_bulk": False,
            "outlet_radius": 0.0127,
            "internal_outlet": 0.0127 * 4,
            "debug": True,
            "revolve": 0,
            "n_revolve": 0,
            "n_wall_layers": 6,
            "VoF": False
        }
    )
    mesh.generate()
    assert mesh.tank.get_radius(mesh.tank.y_start) == pytest.approx(0.0, abs=1e-9)
    assert mesh.tank.cylinder_radius > 0


if __name__ == "__main__":
    test_ksite_mesh()
    print("test_ksite_mesh passed")
