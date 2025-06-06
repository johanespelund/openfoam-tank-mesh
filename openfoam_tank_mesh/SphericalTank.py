import numpy as np

from openfoam_tank_mesh.Tank import Tank


class SphericalTank(Tank):
    def __init__(self, name: str, fill_level: float, outlet_radius: float, radius: float) -> None:
        self.radius = radius
        super().__init__(name, fill_level, outlet_radius)
        return None

    @property
    def height(self) -> float:
        return 2 * self.radius

    def get_radius(self, y: float) -> float:
        self.validate_y_range(y)
        return float(np.sqrt(self.radius**2 - y**2))

    def get_radius_derivative(self, y: float) -> float:
        self.validate_y_range(y)
        return -float(y / np.sqrt(self.radius**2 - y**2))
