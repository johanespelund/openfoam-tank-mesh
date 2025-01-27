# Base class for tank
from abc import ABC, abstractmethod

import numpy as np
import scipy.integrate as spi
import scipy.optimize as spo


class OutOfRange(Exception):
    def __init__(self, y):
        super().__init__(f"y = {y} is out of range.")


class Tank(ABC):
    def __init__(self, name, fill_level, outlet_radius):
        self.name = name
        self.fill_level = fill_level
        self.outlet_radius = outlet_radius

        self.y1 = -self.height / 2
        self.y2 = self.height / 2

        self.volume = self.get_volume()
        self.y_interface = self.calculate_interface_position()
        self.y_outlet = self.calculate_outlet_position()
        self.gas_height = self.y_outlet - self.y_interface
        self.volume_liquid = self.get_partial_volume(self.y1, self.y_interface)
        self.volume_gas = self.get_partial_volume(self.y_interface, self.y2)
        self.area_liquid = self.get_partial_area(self.y1, self.y_interface)
        self.area_gas = self.get_partial_area(self.y_interface, self.y2)

    @property
    @abstractmethod
    def height(self) -> float:
        """
        The height of the tank.
        """
        pass

    def get_volume(self):
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

    def get_partial_volume(self, y1: float, y2: float) -> float:
        """
        Get the volume between y1 and y2, where y1 < y2
        and y1, y2 are in the range [-height/2, height/2].
        """
        r = lambda y: self.get_radius(y)
        integrand = lambda y: np.pi * r(y) ** 2
        return spi.quad(integrand, y1, y2)[0]

    def get_partial_area(self, y1: float, y2: float) -> float:
        """
        Get the area between y1 and y2, where y1 < y2
        and y1, y2 are in the range [-height/2, height/2].
        """

        r = lambda y: self.get_radius(y)
        drdy = lambda y: self.get_radius_derivative(y)
        integrand = lambda y: 2 * np.pi * r(y) * np.sqrt(1 + drdy(y) ** 2)
        return spi.quad(integrand, y1, y2)[0]

    def calculate_interface_position(self) -> float:
        """
        Calculate the position of the interface between the liquid and the gas.
        """

        def objective(y):
            current_volume = self.get_partial_volume(self.y1, y)
            return current_volume - self.fill_level * self.volume

        return np.around(spo.fsolve(objective, 0)[0])

    def calculate_outlet_position(self) -> float:
        """
        Calculate the position of the outlet.
        """

        def objective(y):
            current_radius = self.get_radius(y)
            return current_radius - self.outlet_radius

        return np.around(spo.least_squares(objective, 0.95 * self.y2, bounds=(0, self.y2)).x)

    def plot_tank(self, ax) -> None:
        """
        Plot the tank on the given axis.
        """

        y = np.linspace(self.y1, self.y_outlet, 1000)
        r = np.array([self.get_radius(yi) for yi in y])

        ax.plot(r, y, label=self.name)
        ax.fill_betweenx(y, r, 0, where=y < self.y_interface, alpha=0.5)
        ax.set_aspect("equal")

    def validate_y_range(self, y) -> None:
        if y < self.y1 or y > self.y2:
            raise OutOfRange(y)
