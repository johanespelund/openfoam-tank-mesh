import os

import pytest

from openfoam_tank_mesh.NASA1m3Mesh import NASA1m3Mesh


@pytest.mark.skipif("CI" in os.environ, reason="OpenFOAM is not available in CI")
def test_mesh():
    mesh = NASA1m3Mesh(
        input_parameters={
            "fill_level": 0.82,
            "wall_cell_size": 0.25e-3,
            "wall_tan_cell_size": 3.0e-3,
            "bulk_cell_size": 4e-3,
            "r_BL": 1.1,
            "tri_bulk": False,
            "outlet_radius": 0.024,
            "debug": True,
            "revolve": 0,
            "insulation_type": "bubbles",
            "cargo": "LH2",
        }
    )
    mesh.generate()
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # mesh.tank.plot_tank(ax)
    # plt.show()
    # mesh.tank.write_sample_lines([15, 30, 45, 60, 75], 50e-3)

    # mesh.remove_wall()
    # mesh.cfMesh(5)
    assert mesh.tank.get_radius(0) == mesh.tank.cylinder_radius
    assert mesh.tank.get_radius(mesh.tank.y2) == 0
    assert mesh.tank.get_radius(mesh.tank.y1) == 0
    return None


if __name__ == "__main__":
    test_mesh()
    print("test_ksite_mesh passed")
    pass
