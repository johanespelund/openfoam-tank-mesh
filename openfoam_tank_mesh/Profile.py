import numpy as np
from abc import ABC, abstractmethod
import scipy.integrate as spi  # type: ignore[import-untyped]
import scipy.optimize as spo  # type: ignore[import-untyped]
from openfoam_tank_mesh.exceptions import OutOfRange
import matplotlib.pyplot as plt
import copy
from typing import Optional
from dataclasses import dataclass


@dataclass
class PointCoords:
    # Class to store coordinates and indices for use in gmsh
    points: dict[str, np.ndarray]
    inner_points: list[np.ndarray]
    outer_points: list[np.ndarray]
    wall_points: list[np.ndarray]
    outlet_points: list[np.ndarray]
    internal_outlet_points: list[np.ndarray]
    axis_points: list[np.ndarray]
    i_bl_lower: int
    i_bl: int
    i_bl_upper: int
    y_int_outlet: float


def calculate_boundary_layer(
    r_BL: float,
    wall_cell_size: float,
    wall_tan_cell_size: float,
) -> tuple[int, float, float]:
    """
    Calculate boundary layer parameters.
    returns:
        n_bl: int = number of boundary layer cells
        t_bl: float = thickness of boundary layer
        e_bl: float = boundary layer expansion
    """
    if r_BL == 1 or wall_cell_size == wall_tan_cell_size:
        r_BL = 1.0
        return 4, wall_cell_size * 4, 1
    n = 1
    x = wall_cell_size
    t = wall_cell_size
    while x <= wall_tan_cell_size / r_BL:
        n += 1
        x *= r_BL
        t += x
    return n, t, x / wall_cell_size


def closest_odd(n: float) -> int:
    """
    Return the closest odd integer greater than or equal to n.
    """
    return max(3, int(n) // 2 * 2 + 1)


def closest_even(n: float) -> int:
    return max(2, int(n) // 2 * 2 + 2)


class Segment(ABC):
    def __init__(self, name: str, y_start: float) -> None:
        self.name: str = name
        self.y_start: float = y_start
        self.y_end: float = 0.0 
        self.r_start: float = 0.0
        self.r_end: float = 0.0
        self.length: float = 0.0
        self.upperNeighbor: Segment
        self.lowerNeighbor: Segment
        self.length_scale: float = 0.0
        self.N: int = 0
        self.r: float = 1.0

    @abstractmethod
    def get_radius(self, y: float) -> float:
        pass

    @abstractmethod
    def get_radius_derivative(self, y: float) -> float:
        pass

    @abstractmethod
    def get_length(self) -> float:
        pass

    def __str__(self) -> str:
        return f"{self.name}:\n  y: {self.y_start} - {self.y_end}\n  r: {self.r_start} - {self.r_end}\n  N: {self.N}, r: {self.r}"

    def get_tangent(self, y: float) -> np.ndarray:
        dy_dx = self.get_radius_derivative(y)
        norm = np.sqrt(dy_dx**2 + 1)
        return np.array([dy_dx / norm, 1 / norm])

    def get_normal(self, y: float) -> np.ndarray:
        dy_dx = self.get_radius_derivative(y)
        norm = np.sqrt(dy_dx**2 + 1)
        return np.array([-1 / norm, dy_dx / norm])


class LineSegment(Segment):
    def __init__(
        self,
        name: str,
        y_start: float,
        y_end: float,
        r_start: float,
        r_end: float,
        length_scale: float = 0.0,
    ) -> None:
        super().__init__(name, y_start)
        self.y_end = y_end
        self.r_start = r_start
        self.r_end = r_end
        self.length_scale = length_scale
        self.N = closest_even(self.get_length() / length_scale)

    def get_radius(self, y: float) -> float:
        if self.y_start is None or self.y_end is None or self.r_start is None or self.r_end is None:
            raise RuntimeError("Segment not fully initialized.")
        if y < self.y_start or y > self.y_end:
            raise OutOfRange(y)
        return self.r_start + (self.r_end - self.r_start) * (y - self.y_start) / (self.y_end - self.y_start)

    def get_radius_derivative(self, y: float) -> float:
        if self.r_start is None or self.r_end is None or self.y_start is None or self.y_end is None:
            raise RuntimeError("Segment not fully initialized.")
        return (self.r_end - self.r_start) / (self.y_end - self.y_start)

    def get_length(self) -> float:
        if self.r_start is None or self.r_end is None or self.y_start is None or self.y_end is None:
            raise RuntimeError("Segment not fully initialized.")
        return float(np.sqrt((self.r_end - self.r_start) ** 2 + (self.y_end - self.y_start) ** 2))


class EllipseArc(Segment):
    def __init__(
        self,
        name: str,
        y_start: float,
        y_end: float,
        axis_major: float,
        axis_minor: float,
        y_offset: float = 0.0,
        length_scale: float = 0.0,
    ) -> None:
        super().__init__(name, y_start)
        self.axis_major: float = axis_major
        self.axis_minor: float = axis_minor
        self.y_offset = y_offset
        self.y_start += y_offset
        self.y_end: float = y_end + y_offset
        self.r_start: float = self.get_radius(self.y_start)
        self.r_end: float = self.get_radius(self.y_end)
        self.length_scale: float = length_scale
        self.N: int = closest_even(self.get_length() / length_scale)
        self.major_point: np.typing.NDArray = np.array([self.axis_major, self.axis_minor + self.y_offset])
        self.origo: np.typing.NDArray = np.array([0, self.axis_minor + self.y_offset])

    def get_radius(self, y: float) -> float:
        A, B = self.axis_major, self.axis_minor
        _y = y - B - self.y_offset
        if y < self.y_start or y > self.y_end:
            raise OutOfRange(y)
        return A * float(np.sqrt(max(0.0, 1.0 - (_y / B) ** 2)))

    def get_radius_derivative(self, y: float) -> float:
        A, B = self.axis_major, self.axis_minor
        _y = y - B - self.y_offset
        if y < self.y_start or y > self.y_end:
            raise OutOfRange(y)
        denom = np.sqrt(max(0.0, 1.0 - (_y / B) ** 2))
        return -float(A * _y / (B**2 * denom) if denom != 0 else 0.0)

    def get_length(self) -> float:
        y_vals = np.linspace(self.y_start, self.y_end, 1000)
        drdy = np.array([self.get_radius_derivative(yi) for yi in y_vals])
        return float(spi.trapezoid(np.sqrt(drdy**2 + 1), y_vals))


class CircleArc(EllipseArc):
    def __init__(
        self,
        name: str,
        y_start: float,
        y_end: float,
        radius: float,
        y_offset: float = 0.0,
        length_scale: float = 0.0,
    ) -> None:
        super().__init__(name, y_start, y_end, radius, radius, y_offset, length_scale)


class Profile:
    def __init__(self, segments: list[Segment]) -> None:
        self.segments: list[Segment] = segments
        self.sort_segments()
        self.check_segment_connectivity()

        self.volume: float = self.get_partial_volume(
            self.segments[0].y_start, self.segments[-1].y_end
        )
        self.area: float = self.get_partial_area(
            self.segments[0].y_start, self.segments[-1].y_end
        )
        self.y_start: float = self.segments[0].y_start
        self.y_end: float = self.segments[-1].y_end
        self.i_interface: int = 0
        self.n_upper_bl_segments: int = 0
        self.n_lower_bl_segments: int = 0

    def sort_segments(self) -> None:
        """
        Sort the segments by their y_start value.
        """
        self.segments.sort(key=lambda x: x.y_start)
        for i in range(len(self.segments) - 1):
            self.segments[i].upperNeighbor = self.segments[i + 1]
            self.segments[i + 1].lowerNeighbor = self.segments[i]

    def check_segment_connectivity(self) -> None:
        """
        Check the connectivity of the segments.
        """
        for i in range(len(self.segments) - 1):
            if self.segments[i].y_end != self.segments[i + 1].y_start:
                print(self.segments[i].y_end, self.segments[i + 1].y_start)

                raise ValueError(f"Segments {self.segments[i].name} and {self.segments[i + 1].name} are not connected.")
            if self.segments[i].r_end != self.segments[i + 1].r_start:
                print(self.segments[i])
                print(self.segments[i + 1])

                raise ValueError(f"Segments {self.segments[i].name} and {self.segments[i + 1].name} are not connected.")

    def get_radius(self, y: float) -> float:
        """
        Get the radius at a given height y, where y is in the range [y_start, y_end].
        """
        for segment in self.segments:
            if segment.y_start <= y <= segment.y_end:
                return segment.get_radius(y)
        raise OutOfRange(y)

    def get_radius_derivative(self, y: float) -> float:
        """
        Get the derivative of the radius at a given height y, where y is in the range [y_start, y_end].
        """
        for segment in self.segments:
            if segment.y_start <= y <= segment.y_end:
                return segment.get_radius_derivative(y)
        raise OutOfRange(y)

    def get_normal(self, y: float) -> np.ndarray:
        """
        Return the normalized normal to the curve r(y),
        using the derivative of the curve.
        Points towards the center of the tank.
        """
        for segment in self.segments:
            if segment.y_start <= y <= segment.y_end:
                return segment.get_normal(y)
        raise OutOfRange(y)

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



class TankProfile(Profile):
    """
    Tank profile class. It takes  fill_level and outlet_radius as argument,
    and splits the profile into two parts: liquid and gas.
    """

    def __init__(
        self,
        segments: list[Segment],
        fill_level: float,
        outlet_radius: float,
        internal_outlet: float = 0,
    ) -> None:
        super().__init__(segments=segments)
        self.fill_level: float = fill_level
        self.outlet_radius: float = outlet_radius
        self.internal_outlet: float = internal_outlet
        self.t_BL: float = 0
        self.N: int = 0
        self.y_interface: float = self.calculate_interface_position()
        self.y_outlet: float = self.calculate_outlet_position()
        self.interface_radius: float = self.get_radius(self.y_interface)
        self.area_liquid: float = self.get_partial_area(self.y_start, self.y_interface)
        self.area_gas: float = self.get_partial_area(self.y_interface, self.y_end)
        self.volume_liquid: float = self.get_partial_volume(self.y_start, self.y_interface)
        self.volume_gas: float = self.get_partial_volume(self.y_interface, self.y_end)

        self.split_profile(self.y_outlet, tol=0.000)
        self.segments.pop()

    def calculate_interface_position(self) -> float:
        """
        Calculate the position of the interface between the liquid and the gas.
        """

        def objective(y: float) -> float:
            current_volume = self.get_partial_volume(0, y)
            return current_volume - self.fill_level * self.volume

        guess = 0.5 * (self.y_start + self.y_end)
        return float(spo.fsolve(objective, guess)[0])

    def calculate_outlet_position(self) -> float:
        """
        Calculate the position of the outlet.
        """

        def objective(y: float) -> float:
            current_radius = self.get_radius(y)
            return current_radius - self.outlet_radius

        result = float(spo.least_squares(objective, 0.9 * self.y_end, bounds=(0.6 * self.y_end, self.y_end)).x)
        return result

    def merge_segments(self, segment1: Segment, segment2: Segment) -> Segment:
        """
        Replace two segments with a new LineSegment.
        """
        # Find the lower and upper segments
        if segment1.y_start < segment2.y_start:
            lower_segment = segment1
            upper_segment = segment2
        else:
            lower_segment = segment2
            upper_segment = segment1
        # Create a new LineSegment
        new_segment = LineSegment(
            name=f"{lower_segment.name}_{upper_segment.name}",
            y_start=lower_segment.y_start,
            y_end=upper_segment.y_end,
            r_start=lower_segment.r_start,
            r_end=upper_segment.r_end,
            length_scale=lower_segment.length_scale,
        )

        # Remove the old segments and add the new one
        self.segments.remove(lower_segment)
        self.segments.remove(upper_segment)
        self.segments.append(new_segment)
        self.sort_segments()
        self.check_segment_connectivity()
        return new_segment

    def get_profile_points(self) ->tuple[list[np.ndarray], int]:
        points = []

        for segment in self.segments:
            points.append(np.array([segment.r_start, segment.y_start]))
        points.append(np.array([self.segments[-1].r_end, self.segments[-1].y_end]))

        self.i_interface = 0
        for i, p in enumerate(points):
            if p[1] == self.y_interface:
                self.i_interface = i
                break

        return points, self.i_interface

    def get_profile_normals(self) -> list[np.ndarray]:
        normals = []

        for segment in self.segments:
            normals.append(segment.get_normal(segment.y_start))
        normals.append(self.segments[-1].get_normal(self.segments[-1].y_end))
        return normals

    def ymin(self) -> float:
        """
        Get the minimum y value of the profile.
        """
        return self.segments[0].y_start

    def ymax(self) -> float:
        """
        Get the maximum y value of the profile.
        """
        return self.segments[-1].y_end

    def split_profile(self, y_split: float, tol: float = 10e-3) -> tuple[Segment,Segment]:
        """
        Go through the segments and split the one where the interface is located.
        """
        for segment in self.segments:
            if segment.y_start <= y_split <= segment.y_end:
                # Replace this segment with two segments of same type,
                # meeting at the interface position.
                y_start = segment.y_start
                y_end = segment.y_end
                r_split = segment.get_radius(y_split)

                segment.y_start = y_start
                segment.y_end = y_split
                segment.r_start = segment.get_radius(y_start)
                segment.r_end = segment.get_radius(y_split)
                lower_segment = copy.deepcopy(segment)
                lower_segment.name = f"{segment.name}_lower"

                # If lower segment is shorter than tol, we need to
                # extend it downwards, and shorten the segment below.
                if abs(segment.get_length()) < tol:
                    ln = segment.lowerNeighbor
                    ln.y_end = y_start - tol
                    ln.r_end = ln.get_radius(ln.y_end)
                    # Create a LineSegment to connect the two segments
                    lower_segment = LineSegment(
                        name=f"{segment.name}_lower_extendDown",
                        y_start=ln.y_end,
                        y_end=y_split,
                        r_start=ln.r_end,
                        r_end=r_split,
                        length_scale=segment.length_scale,
                    )

                upper_segment = copy.copy(segment)
                upper_segment.name = f"{segment.name}_upper"
                upper_segment.y_start = y_split
                upper_segment.y_end = y_end
                upper_segment.r_start = upper_segment.get_radius(y_split)
                upper_segment.r_end = upper_segment.get_radius(y_end)

                if abs(upper_segment.get_length()) < tol:
                    # print(f"Segment {upper_segment.name} is too short, extending it upwards.")
                    un = upper_segment.upperNeighbor
                    un.y_start = y_end + tol
                    un.r_start = un.get_radius(un.y_start)

                    # Create a LineSegment to connect the two segments
                    upper_segment = LineSegment(
                        name=f"{segment.name}_upper_extendUp",
                        y_start=y_split,
                        y_end=un.y_start,
                        r_start=r_split,
                        r_end=un.r_start,
                        length_scale=segment.length_scale,
                    )

                break

        upper_segment.N = closest_even(upper_segment.get_length() / upper_segment.length_scale)
        lower_segment.N = closest_even(lower_segment.get_length() / lower_segment.length_scale)

        # Remove old segment and add new ones
        self.segments.remove(segment)
        self.segments.append(lower_segment)
        self.segments.append(upper_segment)
        self.sort_segments()

        self.check_segment_connectivity()

        return lower_segment, upper_segment

    def add_boundary_layers(self, x_wall: float = 0.5e-3, r_BL: float = 1.1) -> None:
        """
        Add boundary layer around the interface (on each side).
        """
        x_bulk = self.segments[-1].length_scale
        n, t, _ = calculate_boundary_layer(r_BL, x_wall, x_bulk)
        self.t_BL = t
        self.N = n

        tol = min(5 * x_wall, 2 * x_bulk)
        for offset in [t, -t, 0]:
            self.split_profile(self.y_interface + offset, tol=tol)

        def get_bl_segments(y_start: float, y_end: float, reverse: bool=False) -> list[Segment]:
            segments = self.segments[::-1] if reverse else self.segments
            return [seg for seg in segments if y_start <= seg.y_start and seg.y_end <= y_end]

        def distribute_cells(bl_segments: list[Segment], r_sign: int) -> None:
            t_local, n_local = 0, 0
            for segment in bl_segments:
                L = segment.get_length()
                N = 0
                t_accum = 0.0
                while t_accum < L / r_BL or (N % 2 == 1 and N > 3):
                    t_accum += x_wall * r_BL**n_local
                    n_local += 1
                    N += 1
                segment.N = N
                segment.r = r_sign * r_BL

            if n_local != self.N or segment.N <= 1:
                segment.N += self.N - n_local
                if segment.N <= 1:
                    neighbor = segment.lowerNeighbor if r_sign > 0 else segment.upperNeighbor
                    new_N = segment.N + neighbor.N
                    new_segment = self.merge_segments(segment, neighbor)
                    new_segment.N = new_N
                    new_segment.r = r_sign * r_BL

        upper_bl_segments = get_bl_segments(self.y_interface, self.y_interface + self.t_BL)
        lower_bl_segments = get_bl_segments(self.y_interface - self.t_BL, self.y_interface, reverse=True)

        self.n_upper_bl_segments = len(upper_bl_segments)
        self.n_lower_bl_segments = len(lower_bl_segments)

        distribute_cells(upper_bl_segments, +1)
        distribute_cells(lower_bl_segments, -1)

    def plot(self) -> None:
        fig, ax = plt.subplots()

        for segment in self.segments:
            y = np.linspace(segment.y_start, segment.y_end, 1000)
            r = np.array([segment.get_radius(yi) for yi in y])
            if len(y) == 0:
                y = np.array([segment.y_start, segment.y_end])
                r = np.array([segment.r_start, segment.r_end])

            ax.plot(r, y, label=segment.name)
            # ax.fill_betweenx(y, r, 0, where=y < segment.y_end, alpha=0.5)

        ax.hlines(
            self.y_interface,
            0,
            self.get_radius(self.y_interface) + 0.1,
            color="k",
            label="Interface",
            ls="--",
        )

        # points = self.get_mesh_points()
        # normals = self.get_profile_normals()
        # for i, (key, item) in enumerate(points.items()):
        #     try:
        #         int(key)
        #         ax.plot(item[0], item[1], "ro")
        #         ax.text(item[0], item[1], f"{key}", fontsize=8, ha="right")
        #         x1, y1 = item[0], item[1]
        #         n = normals[i]
        #         x2 = x1 + n[0]
        #         y2 = y1 + n[1]
        #         ax.plot([x1, x2], [y1, y2])
        # except:
        #      pass

        ax.set_aspect("equal")
        plt.legend()
        plt.show()

    def get_mesh_points(self) -> PointCoords:
        """
        Get the mesh points for use in gmsh.
        """
        points = []
        self.sort_segments()
        self.check_segment_connectivity()
        s = self.segments

        profile_points, i_interface = self.get_profile_points()
        profile_points.append(np.array([0, s[-1].y_end]))

        i_bl_lower = 0
        i_bl = 0
        i_bl_upper = 0
        for i, point in enumerate(profile_points):
            if point[1] == self.y_interface - self.t_BL:
                i_bl_lower = i
            if point[1] == self.y_interface + self.t_BL:
                i_bl_upper = i
            if point[1] == self.y_interface:
                i_bl = i

        profile_normals = self.get_profile_normals()
        # Make sure the normals at start and outlet are vertical
        profile_normals.append(np.array([0, -1]))
        profile_normals[0] = np.array([0, 1])
        profile_normals[-2] = np.array([0, -1])

        inner_points = [p + n * self.t_BL for p, n in zip(profile_points, profile_normals)]

        # Interface needs to be horizontal and held the same t_BL
        b = self.t_BL / profile_normals[i_interface][0]
        inner_points[i_interface] = profile_points[i_interface] + b * np.array([1, 0])

        # Need to adjust the neighbotr points as well, but with optional relaxation
        r = 0.0

        for j in range(i_bl_lower, i_bl_upper + 1):
            b = self.t_BL / profile_normals[j][0]
            norm = self.t_BL * profile_normals[j]
            hor = b * np.array([1, 0])
            inner_points[j] = profile_points[j] + r * norm + (1 - r) * hor

        points += profile_points
        points += inner_points

        tw = 2.08e-3
        wall_points = []  # Points for the optional wall region, outer - n*tw
        for i, point in enumerate(profile_points[:-1]):
            wall_points.append(point - tw * profile_normals[i])

        points.append(np.array([0, self.y_interface - self.t_BL]))
        points.append(np.array([0, self.y_interface]))
        points.append(np.array([0, self.y_interface + self.t_BL]))
        axis_points = [inner_points[0]] + points[-3:]

        # Add the outlet points

        y_int_outlet = self.y_outlet - max(self.t_BL + 2 * s[-1].length_scale, self.internal_outlet)

        points.append(np.array([0, y_int_outlet]))
        points.append(np.array([self.outlet_radius, y_int_outlet]))
        internal_outlet_points = inner_points[-2:] + points[-2:]
        axis_points.append(points[-2])
        points += wall_points

        return PointCoords(
            points={str(i): p for i, p in enumerate(points)},
            inner_points=inner_points,
            outer_points=profile_points,
            wall_points=wall_points,
            outlet_points=inner_points[-2:] + profile_points[-2:][::-1],
            internal_outlet_points=internal_outlet_points,
            axis_points=axis_points,
            i_bl_lower=i_bl_lower,
            i_bl=i_bl,
            i_bl_upper=i_bl_upper,
            y_int_outlet=y_int_outlet,
        )

    def get_curve_groups(self) -> dict[str, list[Segment]]:
        """
        Returns three groups of curves:
        liquid, interface and gas.
        """
        groups: dict = {
            "liquid": [],
            "interface_liquid": [],
            "interface_gas": [],
            "gas": [],
        }
        for seg in self.segments:
            if seg.y_end <= self.y_interface - self.t_BL:
                groups["liquid"].append(seg)
            # elif self.y_interface - self.t_BL < seg.y_end <= self.y_interface + self.t_BL:
            #     groups["interface"].append(seg)
            elif self.y_interface - self.t_BL <= seg.y_start < self.y_interface:
                groups["interface_liquid"].append(seg)
            elif self.y_interface <= seg.y_start < self.y_interface + self.t_BL:
                groups["interface_gas"].append(seg)
            else:
                groups["gas"].append(seg)

        return groups

    def get_partial_area(self, y1: float, y2: float) -> float:
        """
        Get the area between y1 and y2, where y1 < y2
        and y1, y2 are in the range [-height/2, height/2].
        """

        r = lambda y: self.get_radius(y)
        drdy = lambda y: self.get_radius_derivative(y)
        integrand = lambda y: 2 * np.pi * r(y) * np.sqrt(1 + drdy(y) ** 2)
        return float(spi.quad(integrand, y1, y2)[0])

    def get_partial_volume(self, y1: float, y2: float) -> float:
        """
        Get the volume between y1 and y2, where y1 < y2
        and y1, y2 are in the range [-height/2, height/2].
        """
        r = lambda y: self.get_radius(y)
        integrand = lambda y: np.pi * r(y) ** 2
        return float(spi.quad(integrand, y1, y2)[0])


INCH = 0.0254
A = 0.5 * 73 * INCH
B = 0.5 * 87.6 * INCH
C = 1.5 * INCH


class KSiteProfile(TankProfile):
    def __init__(
        self,
        fill_level: float,
        outlet_radius: float,
        bulk_cell_size: float,
        wall_tan_cell_size: float,
        wall_cell_size: float,
        r_BL: float = 1.2,
        internal_outlet: float = 0,
    ) -> None:
        super().__init__(
            segments=[
                EllipseArc("ellipse1", 0, A, B, A, length_scale=wall_tan_cell_size),
                LineSegment("line1", A, A + C, B, B, length_scale=wall_tan_cell_size),
                EllipseArc(
                    "ellipse2",
                    A,
                    2 * A,
                    B,
                    A,
                    y_offset=C,
                    length_scale=wall_tan_cell_size,
                ),
            ],
            fill_level=fill_level,
            outlet_radius=outlet_radius,
            internal_outlet=internal_outlet,
        )
        self.name = "KSite"
        self.add_boundary_layers(x_wall=wall_cell_size, r_BL=r_BL)
        self.cap_height = A
        self.cylinder_radius = B
        self.cylinder_height = C


# FEET = 0.3048
# A = 3.05 / 4 #10 * FEET / 4
# B = 3.05 / 2 #10 * FEET / 2
# C = 3.05 / 2 #5 * FEET


# class MHTB(TankProfile):
#     def __init__(
#         self,
#         fill_level: float,
#         outlet_radius: float,
#         bulk_cell_size: float,
#         wall_tan_cell_size: float,
#         wall_cell_size: float,
#         r_BL: float = 1.2,
#         internal_outlet: float = 0,
#     ) -> None:
#         super().__init__(
#             segments=[
#                 EllipseArc("ellipse1", 0, A, B, A, length_scale=wall_tan_cell_size),
#                 LineSegment("line1", A, A + C, B, B, length_scale=wall_tan_cell_size),
#                 EllipseArc("ellipse2", A, 2 * A, B, A, y_offset=C, length_scale=wall_tan_cell_size),
#             ],
#             fill_level=fill_level,
#             outlet_radius=outlet_radius,
#             internal_outlet=internal_outlet,
#         )
#         self.name = "MHTB"
#         # self.plot()
#         self.add_boundary_layers(x_wall=wall_cell_size, r_BL=r_BL)
#         self.cap_height = A
#         self.cylinder_radius = B
#         self.cylinder_height = C

# if __name__ == "__main__":
#     tank = MHTB(
#         fill_level=0.5,
#         outlet_radius=0.01,
#         bulk_cell_size=0.01,
#         wall_tan_cell_size=0.0001,
#         wall_cell_size=0.005,
#         r_BL=1.2,
#     )
#     print(tank.volume)
#     print(tank.area)
