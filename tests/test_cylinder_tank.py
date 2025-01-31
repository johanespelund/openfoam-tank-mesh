import numpy as np

from openfoam_tank_mesh import CylinderTank as ct


def test_spherical_tank():
    tank = ct.CylinderTank(
        name="test",
        fill_level=0.5,
        outlet_radius=0.01,
        cylinder_radius=1,
        cylinder_height=0.02,
        cap_height=0.75,
    )

    # import matplotlib.pyplot as plt

    # fig, ax = plt.subplots()

    # tank.plot_tank(ax)

    # plt.show()

    assert tank.get_radius(0) == 1
    # assert tank.get_radius(0.6) == 0
    # assert tank.get_radius(-0.6) == 0
    assert np.isclose(tank.calculate_interface_position(), 0)


# if __name__ == "__main__":
#     test_spherical_tank()
