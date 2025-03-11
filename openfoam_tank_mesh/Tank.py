import pathlib
from abc import ABC, abstractmethod

import numpy as np
import scipy.integrate as spi  # type: ignore[import-untyped]
import scipy.optimize as spo  # type: ignore[import-untyped]
from matplotlib.axes import Axes

from openfoam_tank_mesh.exceptions import OutOfRange


class Tank(ABC):
    """
    Base class for tank geometries.
    """

    def __init__(self, name: str, fill_level: float, outlet_radius: float) -> None:
        self.name = name
        self.fill_level = fill_level
        self.outlet_radius = outlet_radius

        self.y1 = -self.height / 2
        self.y2 = self.height / 2

        self.volume = self.get_volume()
        self.y_interface = self.calculate_interface_position()
        self.interface_radius = self.get_radius(self.y_interface)
        self.y_outlet = self.calculate_outlet_position()
        self.gas_height = self.y_outlet - self.y_interface
        self.volume_liquid = self.get_partial_volume(self.y1, self.y_interface)
        self.volume_gas = self.get_partial_volume(self.y_interface, self.y2)
        self.area_liquid = self.get_partial_area(self.y1, self.y_interface)
        self.area_gas = self.get_partial_area(self.y_interface, self.y2)
        super().__init__()

    @property
    @abstractmethod
    def height(self) -> float:
        """
        The height of the tank.
        """
        pass

    @property
    @abstractmethod
    def cylinder_radius(self) -> float:
        """
        The radius of the cylindrical part of the tank, if applicable.
        """
        pass

    @property
    @abstractmethod
    def cylinder_height(self) -> float:
        """
        The height of the cylindrical part of the tank, if applicable.
        """
        pass

    @property
    @abstractmethod
    def cap_height(self) -> float:
        """
        The height of the cap of the tank, if applicable.
        """
        pass

    def get_volume(self) -> float:
        return self.get_partial_volume(self.y1, self.y2)

    @abstractmethod
    def get_radius(self, y: float) -> float:
        """
        Get the radius at a given height y, where y is in the range [-height/2, height/2].
        """
        pass

    @abstractmethod
    def get_radius_derivative(self, y: float) -> float:
        """
        Get the derivative of the radius at a given height y, where y is in the range [-height/2, height/2].
        """
        pass

    def get_tangent(self, y: float) -> np.ndarray:
        """
        Return the normalized tangent to the curve r(y),
        using the derivative of the curve.
        """
        self.validate_y_range(y)
        dy_dx = self.get_radius_derivative(y)
        norm = np.sqrt(dy_dx**2 + 1)
        return np.array([dy_dx / norm, 1 / norm])

    def get_normal(self, y: float) -> np.ndarray:
        """
        Return the normalized normal to the curve r(y),
        using the derivative of the curve.
        Points towards the center of the tank.
        """
        self.validate_y_range(y)
        dy_dx = self.get_radius_derivative(y)
        norm = np.sqrt(dy_dx**2 + 1)
        return np.array([-1 / norm, dy_dx / norm])

    def get_partial_volume(self, y1: float, y2: float) -> float:
        """
        Get the volume between y1 and y2, where y1 < y2
        and y1, y2 are in the range [-height/2, height/2].
        """
        r = lambda y: self.get_radius(y)
        integrand = lambda y: np.pi * r(y) ** 2
        return float(spi.quad(integrand, y1, y2)[0])

    def get_partial_area(self, y1: float, y2: float) -> float:
        """
        Get the area between y1 and y2, where y1 < y2
        and y1, y2 are in the range [-height/2, height/2].
        """

        r = lambda y: self.get_radius(y)
        drdy = lambda y: self.get_radius_derivative(y)
        integrand = lambda y: 2 * np.pi * r(y) * np.sqrt(1 + drdy(y) ** 2)
        return float(spi.quad(integrand, y1, y2)[0])

    def calculate_interface_position(self) -> float:
        """
        Calculate the position of the interface between the liquid and the gas.
        """

        def objective(y: float) -> float:
            current_volume = self.get_partial_volume(self.y1, y)
            return current_volume - self.fill_level * self.volume

        return float(spo.fsolve(objective, 0)[0])

    def calculate_outlet_position(self) -> float:
        """
        Calculate the position of the outlet.
        """

        def objective(y: float) -> float:
            current_radius = self.get_radius(y)
            return current_radius - self.outlet_radius

        return float(spo.least_squares(objective, 0.95 * self.y2, bounds=(0, self.y2)).x)

    def plot_tank(self, ax: Axes) -> None:
        """
        Plot the tank on the given axis.
        """

        y = np.linspace(self.y1, self.y_outlet, 1000)
        r = np.array([self.get_radius(yi) for yi in y])

        ax.plot(r, y, label=self.name)
        ax.fill_betweenx(y, r, 0, where=y < self.y_interface, alpha=0.5)  # type: ignore[arg-type]

        normal = self.get_normal(self.y_interface)
        n = normal / np.linalg.norm(normal)
        tangent = self.get_tangent(self.y_interface)
        t = tangent / np.linalg.norm(tangent)

        # Plot normal and tangent
        ax.arrow(self.interface_radius, self.y_interface, n[0], n[1], color="r", head_width=0.05)
        ax.arrow(self.interface_radius, self.y_interface, t[0], t[1], color="g", head_width=0.05)

        ax.set_aspect("equal")

    def validate_y_range(self, y: float) -> None:
        if y < self.y1 or y > self.y2:
            raise OutOfRange(y)

    def write_sample_lines(self, angles: list[float], length: float) -> None:
        """
        Write sample lines that are normal to the wall of the tank.
        angle = 0 is top of the tank (y-axis)
        """

        output_file = pathlib.Path("system/lineDefinitions/normal_lines")
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, "w") as f:
            for angle in angles:
                theta = np.deg2rad(angle)
                y = self.y2 * np.cos(theta)
                x = self.get_radius(y)
                n = self.get_normal(y)
                dx, dy = n * length
                x_end, y_end = x + dx, y + dy
                # Move the start point slightly outside the tank (1% of length)
                x += 1e-3 * dx
                y += 1e-3 * dy

                f.write(f"normal_{angle}_deg\n{{\n")
                f.write("\ttype lineCellFace;\n")
                f.write("\taxis distance;\n")
                start = f"\t( {x:.6f} {y:.6f} 0 )"
                end = f"\t( {x_end:.6f} {y_end:.6f} 0 )"
                f.write(f"\tstart{start};\n")
                f.write(f"\tend{end};\n")
                f.write("}\n\n")
