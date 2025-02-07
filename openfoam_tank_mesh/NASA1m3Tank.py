import numpy as np

from openfoam_tank_mesh.CylinderTank import CylinderTank

# The radius of the NASA 1m^3 tank,
# using the volume of a sphere, V = 4/3 * pi * r^3
RADIUS = (3 / (4 * np.pi)) ** (1 / 3)


class NASA1m3Tank(CylinderTank):
    """
    A class to represent the
    """

    def __init__(self, fill_level: float, outlet_radius: float, insulation_type: str, cargo: str) -> None:
        self.insulation_type = insulation_type
        self.cargo = cargo
        super().__init__(
            name="NASA1m3",
            fill_level=fill_level,
            outlet_radius=outlet_radius,
            cylinder_radius=RADIUS,
            cylinder_height=0,
            cap_height=RADIUS,
        )
        return None
