import numpy as np

from openfoam_tank_mesh import SphericalTank as st


def test_spherical_tank():
    tank = st.SphericalTank("test", 0.5, 0.01, 1)
    assert tank.get_radius(0) == 1
    assert tank.get_radius(1) == 0
    assert tank.get_radius(-1) == 0
    assert tank.get_radius_derivative(0) == 0
    assert tank.get_volume() == np.pi * 4 / 3 * 1**3
    assert tank.get_partial_volume(-1, 1) == np.pi * 4 / 3 * 1**3
    assert tank.get_partial_volume(-1, 0) == np.pi * 4 / 3 * 1**3 / 2
    assert np.isclose(tank.calculate_interface_position(), 0)
