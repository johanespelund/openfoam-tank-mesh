import numpy as np

from openfoam_tank_mesh.Tank import Tank


class CylinderTank(Tank):
    """
    A tank with a cylindrical body and two spheroidal caps.
    """

    def __init__(
        self,
        name: str,
        fill_level: float,
        outlet_radius: float,
        cylinder_radius: float,
        cylinder_height: float,
        cap_height: float,
    ) -> None:
        self.cylinder_radius = cylinder_radius
        self.cylinder_height = cylinder_height
        self.cap_height = cap_height
        super().__init__(name, fill_level, outlet_radius)
        return None

    @property
    def height(self) -> float:
        return self.cylinder_height + 2 * self.cap_height

    def get_radius(self, y: float) -> float:
        self.validate_y_range(y)

        A = self.cylinder_radius
        B = self.cap_height
        C = self.cylinder_height

        if y > C / 2:
            return float(A * np.sqrt(1 - (y - C / 2) ** 2 / B**2))
        elif y > -C / 2:
            return A
        else:
            return float(A * np.sqrt(1 - (y + C / 2) ** 2 / B**2))

    def get_radius_derivative(self, y: float) -> float:
        self.validate_y_range(y)

        A = self.cylinder_radius
        B = self.cap_height
        C = self.cylinder_height

        if y > C / 2:
            return -float(A * (y - C / 2) / (B * np.sqrt(B**2 - (y - C / 2) ** 2)))
        elif y > -C / 2:
            return 0
        else:
            return -float(A * (y + C / 2) / (B * np.sqrt(B**2 - (y + C / 2) ** 2)))
