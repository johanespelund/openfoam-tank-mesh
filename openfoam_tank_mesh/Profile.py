import copy
import importlib.util
import logging
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass
from math import sqrt

# ---------------------------------------------------------------------------
# mpl_toolkits workaround
# The system-installed mpl_toolkits (/usr/lib/python3/dist-packages) is a
# namespace package that shadows the pip-installed version under ~/.local and
# is incompatible with the pip matplotlib.  Load the pip copy explicitly
# before matplotlib.pyplot is imported so that the '3d' projection registers.
# ---------------------------------------------------------------------------
_pip_mplot3d = "/home/johan/.local/lib/python3.10/site-packages/mpl_toolkits/mplot3d"
for _mod_name, _rel in [
    ("mpl_toolkits.mplot3d", "__init__.py"),
    ("mpl_toolkits.mplot3d.proj3d", "proj3d.py"),
    ("mpl_toolkits.mplot3d.art3d", "art3d.py"),
    ("mpl_toolkits.mplot3d.axis3d", "axis3d.py"),
    ("mpl_toolkits.mplot3d.axes3d", "axes3d.py"),
]:
    if _mod_name not in sys.modules:
        _spec = importlib.util.spec_from_file_location(
            _mod_name,
            f"{_pip_mplot3d}/{_rel}",
            submodule_search_locations=[_pip_mplot3d],
        )
        _m = importlib.util.module_from_spec(_spec)
        sys.modules[_mod_name] = _m
        _spec.loader.exec_module(_m)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy.integrate as spi  # type: ignore[import-untyped]  # noqa: E402
import scipy.optimize as spo  # type: ignore[import-untyped]  # noqa: E402

from openfoam_tank_mesh.exceptions import (  # noqa: E402
    BoundaryLayerTooThick,
    OutOfRange,
    SegmentNotInitialized,
    SegmentsNotConnected,
)

logger = logging.getLogger(__name__)

MM_PER_INCH = 25.4
ELSEVIER_5P_PAGEWIDTH_MM = 190.0
ONE_HALF_PAGEWIDTH_INCH = ELSEVIER_5P_PAGEWIDTH_MM / MM_PER_INCH / 2.0


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

    @abstractmethod
    def get_origo(self) -> tuple[float, float]:
        pass

    @abstractmethod
    def get_rmax(self) -> float:
        pass

    def __str__(self) -> str:
        return (
            f"{self.name}:\n  y: {self.y_start} - {self.y_end}"
            f"\n  r: {self.r_start} - {self.r_end}\n  N: {self.N}, r: {self.r}"
        )

    def get_tangent(self, y: float) -> np.ndarray:
        dy_dx = self.get_radius_derivative(y)
        norm = sqrt(dy_dx**2 + 1)
        return np.array([dy_dx / norm, 1 / norm])

    def get_normal(self, y: float) -> np.ndarray:
        dy_dx = self.get_radius_derivative(y)
        norm = sqrt(dy_dx**2 + 1)
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

    def get_origo(self) -> tuple[float, float]:
        return (0.0, (self.y_start + self.y_end) / 2)

    def get_rmax(self) -> float:
        return max(self.r_start, self.r_end)

    def get_radius(self, y: float) -> float:
        if self.y_start is None or self.y_end is None or self.r_start is None or self.r_end is None:
            raise SegmentNotInitialized()
        if y < self.y_start or y > self.y_end:
            raise OutOfRange(y)
        return self.r_start + (self.r_end - self.r_start) * (y - self.y_start) / (self.y_end - self.y_start)

    def get_radius_derivative(self, y: float) -> float:
        if self.r_start is None or self.r_end is None or self.y_start is None or self.y_end is None:
            raise SegmentNotInitialized()
        return (self.r_end - self.r_start) / (self.y_end - self.y_start)

    def get_length(self) -> float:
        if self.r_start is None or self.r_end is None or self.y_start is None or self.y_end is None:
            raise SegmentNotInitialized()
        return float(sqrt((self.r_end - self.r_start) ** 2 + (self.y_end - self.y_start) ** 2))


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

    def get_origo(self) -> tuple[float, float]:
        return (0.0, self.axis_minor + self.y_offset)

    def get_rmax(self) -> float:
        return self.axis_major

    def get_major_point(self) -> tuple[float, float]:
        return (0, self.axis_major)

    def get_radius(self, y: float) -> float:
        A, B = self.axis_major, self.axis_minor
        _y = y - B - self.y_offset
        if y < self.y_start or y > self.y_end:
            raise OutOfRange(y)
        return A * float(sqrt(max(0.0, 1.0 - (_y / B) ** 2)))

    def get_radius_derivative(self, y: float) -> float:
        A, B = self.axis_major, self.axis_minor
        _y = y - B - self.y_offset
        if y < self.y_start or y > self.y_end:
            raise OutOfRange(y)
        denom = sqrt(max(0.0, 1.0 - (_y / B) ** 2))
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
    def __init__(
        self,
        segments: list[Segment],
        extrude: bool = False,
        extrude_thickness: float = 0.0,
    ) -> None:
        self.segments: list[Segment] = segments
        self.extrude: bool = extrude
        self.extrude_thickness: float = extrude_thickness
        self.sort_segments()
        self.check_segment_connectivity()

        self.volume: float = self.get_partial_volume(self.segments[0].y_start, self.segments[-1].y_end)
        self.area: float = self.get_partial_area(self.segments[0].y_start, self.segments[-1].y_end)
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
                logger.debug(
                    "Segment connectivity mismatch: %s y_end=%s, %s y_start=%s",
                    self.segments[i].name,
                    self.segments[i].y_end,
                    self.segments[i + 1].name,
                    self.segments[i + 1].y_start,
                )
                raise SegmentsNotConnected(self.segments[i].name, self.segments[i + 1].name)
            if self.segments[i].r_end != self.segments[i + 1].r_start:
                logger.debug("Segment r mismatch: %s", self.segments[i])
                logger.debug("Next segment: %s", self.segments[i + 1])
                raise SegmentsNotConnected(self.segments[i].name, self.segments[i + 1].name)

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

        For a revolved (axisymmetric) geometry the standard solid-of-revolution
        formula is used: V = ∫ π r(y)² dy.

        For an extruded (empty_2d) geometry the cross-section is a rectangle of
        width 2·r(y) and depth ``extrude_thickness``:
        V = ∫ 2 · r(y) · t dy
        """
        r = lambda y: self.get_radius(y)
        if self.extrude:
            t = self.extrude_thickness
            integrand = lambda y: 2 * r(y) * t
        else:
            integrand = lambda y: np.pi * r(y) ** 2
        return float(spi.quad(integrand, y1, y2)[0])

    def get_partial_area(self, y1: float, y2: float) -> float:
        """
        Get the lateral wall area between y1 and y2, where y1 < y2
        and y1, y2 are in the range [-height/2, height/2].

        For a revolved (axisymmetric) geometry the surface-of-revolution formula
        is used: A = ∫ 2π r(y) √(1 + (dr/dy)²) dy.

        For an extruded (empty_2d) geometry the lateral wall consists of two
        flat side walls (front and back) of thickness ``extrude_thickness``:
        A = ∫ 2 · t · √(1 + (dr/dy)²) dy
        The flat end caps (the empty_pos / empty_neg faces) are excluded.
        """
        r = lambda y: self.get_radius(y)
        drdy = lambda y: self.get_radius_derivative(y)
        if self.extrude:
            t = self.extrude_thickness
            integrand = lambda y: 2 * t * np.sqrt(1 + drdy(y) ** 2)
        else:
            integrand = lambda y: 2 * np.pi * r(y) * np.sqrt(1 + drdy(y) ** 2)
        return float(spi.quad(integrand, y1, y2)[0])

    def get_interface_area(self, y: float) -> float:
        """
        Get cross-sectional area at height ``y``.

        For a revolved (axisymmetric) geometry this is a disk area,
        ``pi * r(y)^2``.

        For an extruded (empty_2d) geometry this is a rectangle area,
        ``2 * r(y) * extrude_thickness``.
        """
        r = self.get_radius(y)
        if self.extrude:
            return float(2 * r * self.extrude_thickness)
        return float(np.pi * r**2)

    def get_y(self, radius: float, ymin: float, ymax: float) -> float:
        def objective(y: float) -> float:
            current_radius = self.get_radius(y)
            return current_radius - radius

        return float(spo.brentq(objective, ymin, ymax))


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
        wall_thickness: float = 2.08e-3,
        tolerance: float = 10e-3,
        extrude: bool = False,
        extrude_thickness: float = 0.0,
    ) -> None:
        super().__init__(segments=segments, extrude=extrude, extrude_thickness=extrude_thickness)
        self.name: str = ""
        self.fill_level: float = fill_level
        self.outlet_radius: float = outlet_radius
        self.internal_outlet: float = internal_outlet
        self.wall_thickness: float = wall_thickness
        self.tolerance: float = tolerance
        self.t_BL: float = 0
        self.N: int = 0
        self.cylinder_radius: float = 0.0
        self.cylinder_height: float = 0.0
        self.cap_height: float = 0.0
        self.y_interface: float = self.calculate_interface_position()
        self.y_outlet: float = self.calculate_outlet_position()
        self.interface_radius: float = self.get_radius(self.y_interface)
        self.area_wall_liquid: float = self.get_area_wall_liquid()
        self.area_wall_gas: float = self.get_area_wall_gas()
        self.area_interface: float = self.get_area_interface()
        self.volume_liquid: float = self.get_partial_volume(self.y_start, self.y_interface)
        self.volume_gas: float = self.get_partial_volume(self.y_interface, self.y_end)

        # Save full segments (before outlet truncation) for visualization.
        self._viz_segments: list[Segment] = copy.deepcopy(self.segments)

        self.split_profile(self.y_outlet, tol=0.000)
        self.segments.pop()

    def calculate_interface_position(self) -> float:
        """
        Calculate the position of the interface between the liquid and the gas.
        """

        def objective(y: float) -> float:
            current_volume = self.get_partial_volume(0, float(y))
            return current_volume - self.fill_level * self.volume

        return float(spo.brentq(objective, self.y_start, self.y_end))

    def get_area_wall_liquid(self) -> float:
        """Return wall area in the liquid region (from ``y_start`` to ``y_interface``)."""
        if hasattr(self, "area_wall_liquid"):
            return float(self.area_wall_liquid)
        return self.get_partial_area(self.y_start, self.y_interface)

    def get_area_wall_gas(self) -> float:
        """Return wall area in the gas region (from ``y_interface`` to ``y_end``)."""
        if hasattr(self, "area_wall_gas"):
            return float(self.area_wall_gas)
        return self.get_partial_area(self.y_interface, self.y_end)

    def get_area_interface(self) -> float:
        """Return liquid-gas interface area at ``y_interface``."""
        if hasattr(self, "area_interface"):
            return float(self.area_interface)
        return self.get_interface_area(self.y_interface)

    @property
    def area_liquid(self) -> float:
        """Backward-compatible alias for ``area_wall_liquid``."""
        return self.area_wall_liquid

    @property
    def area_gas(self) -> float:
        """Backward-compatible alias for ``area_wall_gas``."""
        return self.area_wall_gas

    def calculate_outlet_position(self) -> float:
        """
        Calculate the position of the outlet.
        """

        def objective(y: float) -> float:
            current_radius = self.get_radius(y)
            return current_radius - self.outlet_radius

        return float(spo.brentq(objective, 0.6 * self.y_end, self.y_end))

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

    def get_profile_points(self) -> tuple[list[np.ndarray], int]:
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

    def split_profile(self, y_split: float, tol: float = 10e-3) -> tuple[Segment, Segment | None]:
        """
        Go through the segments and split the one where the interface is located.
        """

        for segment in self.segments:
            if segment.y_start <= y_split <= segment.y_end:
                # Replace this segment with two segments of same type,
                # meeting at the interface position.
                y_start = segment.y_start
                y_end = segment.y_end
                r_start = segment.r_start
                r_end = segment.r_end
                r_split = segment.get_radius(y_split)

                segment.y_start = y_start
                segment.y_end = y_split
                segment.r_start = segment.get_radius(y_start)
                segment.r_end = r_split
                lower_segment = copy.deepcopy(segment)
                lower_segment.name = f"{segment.name}_lower"

                # If lower segment is shorter than tol, we need to
                # extend it downwards, and shorten the segment below.
                if abs(lower_segment.get_length()) < tol:
                    ln = segment.lowerNeighbor
                    new_ln = ln
                    if ln.get_length() <= 2.5 * tol:
                        y_start = ln.y_start
                        r_start = ln.r_start
                        new_ln = ln.lowerNeighbor
                        self.segments.remove(ln)
                    else:
                        y_start = y_start - 1.0 * tol
                        r_start = ln.get_radius(y_start)
                        ln.y_end = y_start
                        ln.r_end = r_start
                    # Create a LineSegment to connect the two segments
                    lower_segment = LineSegment(
                        name=f"{segment.name}_lower_extendDown",
                        y_start=y_start,
                        y_end=y_split,
                        r_start=r_start,
                        r_end=r_split,
                        length_scale=segment.length_scale,
                    )
                    lower_segment.lowerNeighbor = new_ln

                upper_segment = copy.copy(segment)
                upper_segment.name = f"{segment.name}_upper"
                upper_segment.y_start = y_split
                upper_segment.y_end = y_end
                upper_segment.r_start = r_split
                upper_segment.r_end = r_end

                if abs(upper_segment.get_length()) < 1.0 * tol:
                    un = upper_segment.upperNeighbor
                    if un.get_length() <= 2.5 * tol:
                        y_end = un.y_end
                        r_end = un.r_end
                        new_un = un.upperNeighbor
                        self.segments.remove(un)
                    else:
                        y_end = y_end + 1.0 * tol
                        r_end = un.get_radius(y_end)
                        un.y_start = y_end
                        un.r_start = r_end
                        new_un = un

                    # Create a LineSegment to connect the two segments
                    upper_segment = LineSegment(
                        name=f"{segment.name}_upper_extendUp",
                        y_start=y_split,
                        y_end=y_end,
                        r_start=r_split,
                        r_end=r_end,
                        length_scale=segment.length_scale,
                    )
                    upper_segment.upperNeighbor = new_un

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

        # Check that the boundary layer is not too large.
        if 2 * self.t_BL > min(self.y_interface, self.y_outlet - self.y_interface):
            raise BoundaryLayerTooThick()

        tol = min(5 * x_wall, 2 * x_bulk)
        for offset in [t, -t, 0]:
            self.split_profile(self.y_interface + offset, tol=tol)

        def get_bl_segments(y_start: float, y_end: float, reverse: bool = False) -> list[Segment]:
            segments = self.segments[::-1] if reverse else self.segments
            return [seg for seg in segments if y_start <= seg.y_start and seg.y_end <= y_end]

        def distribute_cells(bl_segments: list[Segment], r_sign: int) -> None:
            n_local = 0
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

            segment = bl_segments[-1]
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
        _fig, ax = plt.subplots()

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

        # point_coords = self.get_mesh_points()
        # normals = self.get_profile_normals()
        # points = point_coords.points
        # for _i, (key, item) in enumerate(points.items()):
        #     ax.plot(item[0], item[1], "ro")
        #     ax.text(item[0], item[1], f"{key}", fontsize=8, ha="right")
        # x1, y1 = item[0], item[1]
        # n = normals[i]
        # x2 = x1 + n[0]
        # y2 = y1 + n[1]
        # ax.plot([x1, x2], [y1, y2])

        ax.set_aspect("equal")
        plt.legend()
        plt.show()

    def plot_3d(  # noqa: C901
        self,
        wedge_angle: float = 20.0,
        save_path: str | None = None,
        dpi: int = 600,
        show: bool = True,
        sim_domain_plane_only: bool = False,
    ) -> None:
        """
        Plot the tank geometry in 3D.

        The full tank wall is shown as a surface of revolution around the
        vertical (y) axis.  Liquid and gas regions are coloured differently.
        The simulation wedge domain (±wedge_angle/2 degrees) is overlaid as
        wireframe lines.

        Parameters
        ----------
        wedge_angle:
            Total wedge angle in degrees (default 30°).  The wedge is centred
            on the xz-plane (phi = 0) and spans ±wedge_angle/2.
        save_path:
            Optional output path for saving the figure as an image.
        dpi:
            Resolution used when saving the figure.
        show:
            If True, display the figure interactively.
        sim_domain_plane_only:
            If True, draw the simulation domain as a single symmetry plane
            (phi=0) instead of a wedge volume.
        """
        N_phi = 180  # azimuthal resolution for the surface
        N_y = 500  # axial resolution for the surface

        fig = plt.figure()  # figsize=(ONE_HALF_PAGEWIDTH_INCH, ONE_HALF_PAGEWIDTH_INCH * (10.0 / 8.0)))
        ax = fig.add_subplot(111, projection="3d")

        # Coordinate convention: the symmetry axis (y in Profile) maps to the
        # vertical Z axis of the 3D plot so that "up" is visually correct.
        # Radial directions map to X and Y of the 3D plot.
        # Profile y  →  plot Z
        # r·cos(φ)   →  plot X
        # r·sin(φ)   →  plot Y

        # ------------------------------------------------------------------
        # Use full pre-outlet segments for a smooth surface to the axis
        # ------------------------------------------------------------------
        viz_segs = getattr(self, "_viz_segments", self.segments)
        y_viz_start = viz_segs[0].y_start
        y_viz_end = viz_segs[-1].y_end

        def _viz_radius(y: float) -> float:
            for seg in viz_segs:
                if seg.y_start <= y <= seg.y_end:
                    return seg.get_radius(y)
            raise OutOfRange(y)

        y_all = np.linspace(y_viz_start, y_viz_end, N_y)
        r_all = np.array([_viz_radius(yi) for yi in y_all])

        half = np.radians(wedge_angle / 2.0)
        phi_full = np.linspace(0, 2 * np.pi, N_phi)
        n_wedge = max(24, int(N_phi * wedge_angle / 360.0))
        phi_wedge = np.linspace(-half, half, n_wedge)

        def _surface(y_arr: np.ndarray, r_arr: np.ndarray, phi_arr: np.ndarray, color: str, alpha: float) -> None:
            Z_s, P = np.meshgrid(y_arr, phi_arr)
            R = np.tile(r_arr, (len(phi_arr), 1))
            X_s = R * np.cos(P)
            Y_s = R * np.sin(P)
            ax.plot_surface(
                X_s,
                Y_s,
                Z_s,
                color=color,
                alpha=alpha,
                linewidth=0,
                antialiased=True,
                shade=False,
                rcount=Z_s.shape[0],
                ccount=Z_s.shape[1],
            )

        # Liquid region (y_start → y_interface)
        y_liq = np.linspace(y_viz_start, self.y_interface, N_y)
        r_liq = np.array([_viz_radius(yi) for yi in y_liq])
        _surface(y_liq, r_liq, phi_full, color="#4499cc", alpha=0.3)
        if not sim_domain_plane_only:
            _surface(y_liq, r_liq, phi_wedge, color="#4499cc", alpha=0.6)

        # Gas region (y_interface → y_end)
        y_gas = np.linspace(self.y_interface, y_viz_end, N_y)
        r_gas = np.array([_viz_radius(yi) for yi in y_gas])
        _surface(y_gas, r_gas, phi_full, color="#dddddd", alpha=0.10)
        if not sim_domain_plane_only:
            _surface(y_gas, r_gas, phi_wedge, color="#dddddd", alpha=0.38)

        # Liquid-gas interface disk
        r_int = _viz_radius(self.y_interface)
        r_disk = np.linspace(0, r_int, 40)
        disk_sets: list[tuple[np.ndarray, float]] = [(phi_full, 0.14)]
        if not sim_domain_plane_only:
            disk_sets.append((phi_wedge, 0.50))
        for phi_disk, alpha_disk in disk_sets:
            P_disk, R_disk = np.meshgrid(phi_disk, r_disk)
            X_disk = R_disk * np.cos(P_disk)
            Y_disk = R_disk * np.sin(P_disk)
            Z_disk = np.full_like(X_disk, self.y_interface)
            ax.plot_surface(
                X_disk,
                Y_disk,
                Z_disk,
                color="#4499cc",
                alpha=alpha_disk,
                linewidth=0,
                antialiased=True,
                shade=True,
                rcount=Z_disk.shape[0],
                ccount=Z_disk.shape[1],
            )

        if sim_domain_plane_only:
            s_plane = np.linspace(0.0, 1.0, 52)

            def _sim_plane(y_arr: np.ndarray, r_arr: np.ndarray, color: str, alpha: float) -> None:
                S, Z_p = np.meshgrid(s_plane, y_arr, indexing="ij")
                X_p = S * r_arr[None, :]
                Y_p = np.zeros_like(X_p)
                ax.plot_surface(
                    X_p,
                    Y_p,
                    Z_p,
                    color=color,
                    alpha=alpha,
                    linewidth=0,
                    antialiased=True,
                    shade=False,
                    rcount=Z_p.shape[0],
                    ccount=Z_p.shape[1],
                )

            _sim_plane(y_liq, r_liq, color="#4499cc", alpha=0.48)
            _sim_plane(y_gas, r_gas, color="#dddddd", alpha=0.30)
        else:
            # Filled wedge side faces (phi = ±half), split by phase at y_interface.
            s_face = np.linspace(0.0, 1.0, 28)

            def _wedge_side_face(
                y_arr: np.ndarray, r_arr: np.ndarray, phi_face: float, color: str, alpha: float
            ) -> None:
                S, Z_f = np.meshgrid(s_face, y_arr, indexing="ij")
                R_f = S * r_arr[None, :]
                X_f = R_f * np.cos(phi_face)
                Y_f = R_f * np.sin(phi_face)
                ax.plot_surface(
                    X_f,
                    Y_f,
                    Z_f,
                    color=color,
                    alpha=alpha,
                    linewidth=0,
                    antialiased=True,
                    shade=True,
                    rcount=Z_f.shape[0],
                    ccount=Z_f.shape[1],
                )

            for phi_w in [half, -half]:
                _wedge_side_face(y_liq, r_liq, phi_w, color="#4499cc", alpha=0.40)
                _wedge_side_face(y_gas, r_gas, phi_w, color="#dddddd", alpha=0.24)

        # ------------------------------------------------------------------
        # Simulation domain wireframe
        # ------------------------------------------------------------------
        lw, col = 1.0, "k"
        r_if = _viz_radius(self.y_interface)
        if sim_domain_plane_only:
            ax.plot(r_all, np.zeros_like(r_all), y_all, color=col, lw=lw, zorder=100)
            ax.plot([0, r_all[0]], [0, 0], [y_viz_start, y_viz_start], color=col, lw=lw, zorder=100)
            ax.plot([0, r_all[-1]], [0, 0], [y_viz_end, y_viz_end], color=col, lw=lw, zorder=100)
            ax.plot([0, r_if], [0, 0], [self.y_interface, self.y_interface], color=col, lw=lw, ls="--", zorder=100)
        else:
            wedge_phis = [half, -half]
            for phi_w in wedge_phis:
                xw = r_all * np.cos(phi_w)
                yw = r_all * np.sin(phi_w)

                # Profile outline on this wedge face
                ax.plot(xw, yw, y_all, color=col, lw=lw, zorder=100)

                # Bottom radial edge (axis → wall)
                ax.plot([0, xw[0]], [0, yw[0]], [y_viz_start, y_viz_start], color=col, lw=lw, zorder=100)

                # Top radial edge (axis → wall)
                ax.plot([0, xw[-1]], [0, yw[-1]], [y_viz_end, y_viz_end], color=col, lw=lw, zorder=100)

                # Interface line (axis → wall at y_interface)
                ax.plot(
                    [0, r_if * np.cos(phi_w)],
                    [0, r_if * np.sin(phi_w)],
                    [self.y_interface, self.y_interface],
                    color=col,
                    lw=lw,
                    ls="--",
                    zorder=100,
                )

            # Arcs connecting the two wedge faces at y_start, y_interface, y_end
            for y_arc, r_arc in [
                (y_viz_start, r_all[0]),
                (self.y_interface, r_if),
                (y_viz_end, r_all[-1]),
            ]:
                phi_arc = np.linspace(-half, half, 60)
                ax.plot(
                    r_arc * np.cos(phi_arc),
                    r_arc * np.sin(phi_arc),
                    np.full_like(phi_arc, y_arc),
                    color=col,
                    lw=lw,
                    zorder=100,
                )

        # Axis line — extends 10 % of tank height above and below
        tank_height = y_viz_end - y_viz_start
        ax_ext = 0.20 * tank_height
        ax.plot(
            [0, 0],
            [0, 0],
            [y_viz_start - ax_ext, y_viz_end + ax_ext],
            color="k",
            lw=1.5,
            ls="dashdot",
            alpha=0.6,
            zorder=-5,
        )

        # Axis line on each wedge face (vertical edge at r=0)
        if not sim_domain_plane_only:
            wedge_phis = [half, -half]
            for _phi_w in wedge_phis:
                ax.plot([0, 0], [0, 0], [y_viz_start, y_viz_end], color=col, lw=lw)

        # ------------------------------------------------------------------
        # Formatting
        # ------------------------------------------------------------------
        ax.set_axis_off()
        r_max = np.max(r_all)
        tank_height = y_viz_end - y_viz_start
        z_min = y_viz_start - ax_ext
        z_max = y_viz_end + ax_ext

        # Keep limits tight so the tank uses the full canvas.
        ax.set_xlim(-r_max, r_max)
        ax.set_ylim(-r_max, r_max)
        ax.set_zlim(z_min, z_max)
        ax.set_box_aspect([2 * r_max, 2 * r_max, z_max - z_min])

        # 3D + axis-off can still leave large default subplot margins.
        ax.set_position([0.0, 0.0, 1.0, 1.0])
        fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight", pad_inches=0)
        if show:
            plt.show()
        else:
            plt.close(fig)

    def plot_3d_horizontal(  # noqa: C901
        self,
        save_path: str | None = None,
        dpi: int = 600,
        show: bool = True,
        sim_domain_plane_only: bool = False,
    ) -> None:
        """
        Plot the tank as a horizontal cylinder for the empty_2d (extruded) case.

        The tank cross-section (the r(y) profile) is in the X-Z plane with Z
        vertical (profile y-axis → plot Z).  The cylinder is extruded along Y
        (into the page) to a compact visualization length of 2·diameter.  A rectangular
        simulation slab of depth 0.5·diameter is shown as wireframe lines
        centred at y=0.  The cut ends are softened with a small fillet.
        Liquid and gas regions are shown in different colours.

        Coordinate mapping
        ------------------
        Profile x (radius, ±r(y))  →  plot X
        Profile y (height)         →  plot Z   (vertical)
        Extrusion direction        →  plot Y   (into the page)

        Parameters
        ----------
        save_path:
            Optional output path for saving the figure as an image.
        dpi:
            Resolution used when saving the figure.
        show:
            If True, display the figure interactively.
        sim_domain_plane_only:
            If True, draw the simulation domain as a single symmetry plane
            (y=0) instead of an extrusion slab.
        """
        N_profile = 500  # points along the profile outline

        fig = plt.figure(figsize=(1.5 * ONE_HALF_PAGEWIDTH_INCH, 1.5 * ONE_HALF_PAGEWIDTH_INCH * (9.0 / 11.0)))
        ax = fig.add_subplot(111, projection="3d")

        # ------------------------------------------------------------------
        # Profile geometry (use full pre-outlet segments)
        # ------------------------------------------------------------------
        viz_segs = getattr(self, "_viz_segments", self.segments)
        y_viz_start = viz_segs[0].y_start
        y_viz_end = viz_segs[-1].y_end

        def _viz_radius(y: float) -> float:
            for seg in viz_segs:
                if seg.y_start <= y <= seg.y_end:
                    return seg.get_radius(y)
            raise OutOfRange(y)

        # Sample the profile
        z_prof = np.linspace(y_viz_start, y_viz_end, N_profile)  # plot Z
        r_prof = np.array([_viz_radius(zi) for zi in z_prof])

        diameter = 2.0 * np.max(r_prof)
        half_len = 1.5 * diameter  # half of total 2d extrusion length (plot Y)
        slab_depth = 0.20 * diameter  # simulation domain: symmetry plane to one side

        # ------------------------------------------------------------------
        # Ruled wall surface helper:
        # surface spans x_out/z_out outline at two y (extrusion) positions
        # ------------------------------------------------------------------
        def _wall_surface(x_out: np.ndarray, z_out: np.ndarray, y0: float, y1: float, color: str, alpha: float) -> None:
            n = len(x_out)
            X_s = np.array([x_out, x_out])
            Z_s = np.array([z_out, z_out])
            Y_s = np.array([[y0] * n, [y1] * n])
            ax.plot_surface(
                X_s,
                Y_s,
                Z_s,
                color=color,
                alpha=alpha,
                linewidth=0,
                antialiased=True,
                shade=False,
                rcount=Z_s.shape[0],
                ccount=Z_s.shape[1],
            )

        # ------------------------------------------------------------------
        # Split outline at interface height
        i_if = np.searchsorted(z_prof, self.y_interface)
        r_if = _viz_radius(self.y_interface)

        x_liq = np.concatenate([r_prof[: i_if + 1], -r_prof[: i_if + 1][::-1], [r_prof[0]]])
        z_liq = np.concatenate([z_prof[: i_if + 1], z_prof[: i_if + 1][::-1], [z_prof[0]]])

        x_gas = np.concatenate([r_prof[i_if:], -r_prof[i_if:][::-1], [r_prof[i_if]]])
        z_gas = np.concatenate([z_prof[i_if:], z_prof[i_if:][::-1], [z_prof[i_if]]])

        # Full geometry as context (more transparent).
        _wall_surface(x_liq, z_liq, -half_len, half_len, "#4499cc", 0.20)
        _wall_surface(x_gas, z_gas, -half_len, half_len, "#cccccc", 0.12)

        # Simulation slab highlighted (less transparent), only on sim-side half.
        x_liq_side = -r_prof[: i_if + 1]
        z_liq_side = z_prof[: i_if + 1]
        x_gas_side = -r_prof[i_if:]
        z_gas_side = z_prof[i_if:]
        if not sim_domain_plane_only:
            _wall_surface(x_liq_side, z_liq_side, 0.0, slab_depth, "#4499cc", 0.58)
            _wall_surface(x_gas_side, z_gas_side, 0.0, slab_depth, "#cccccc", 0.36)

        # ------------------------------------------------------------------
        # End cap disks at y = ±half_len, split liquid/gas
        # ------------------------------------------------------------------
        def _end_disk(y_val: float, alpha_scale: float = 1.0, half_only: bool = False) -> None:
            z_d = np.linspace(y_viz_start, self.y_interface, 200)
            r_d = np.array([_viz_radius(zi) for zi in z_d])
            X_d = np.array([-r_d, np.zeros_like(r_d)]) if half_only else np.array([-r_d, r_d])
            Z_d = np.array([z_d, z_d])
            Y_d = np.full_like(X_d, y_val)
            ax.plot_surface(
                X_d,
                Y_d,
                Z_d,
                color="#4499cc",
                alpha=0.45 * alpha_scale,
                linewidth=0,
                antialiased=True,
                shade=False,
                rcount=Z_d.shape[0],
                ccount=Z_d.shape[1],
            )
            z_g = np.linspace(self.y_interface, y_viz_end, 200)
            r_g = np.array([_viz_radius(zi) for zi in z_g])
            X_g = np.array([-r_g, np.zeros_like(r_g)]) if half_only else np.array([-r_g, r_g])
            Z_g = np.array([z_g, z_g])
            Y_g = np.full_like(X_g, y_val)
            ax.plot_surface(
                X_g,
                Y_g,
                Z_g,
                color="#cccccc",
                alpha=0.25 * alpha_scale,
                linewidth=0,
                antialiased=True,
                shade=True,
                rcount=Z_g.shape[0],
                ccount=Z_g.shape[1],
            )

        _end_disk(-half_len, alpha_scale=0.45)
        _end_disk(half_len, alpha_scale=0.45)
        _end_disk(0.0, alpha_scale=0.95, half_only=True)
        if not sim_domain_plane_only:
            _end_disk(slab_depth, alpha_scale=0.95, half_only=True)

        # ------------------------------------------------------------------
        # Interface plane (rectangle at z = y_interface, full length in Y)
        # ------------------------------------------------------------------
        y_range = np.array([-half_len, half_len])
        x_range = np.array([-r_if, r_if])
        X_if, Y_if = np.meshgrid(x_range, y_range)
        Z_if = np.full_like(X_if, self.y_interface)
        ax.plot_surface(
            X_if,
            Y_if,
            Z_if,
            color="#88bbdd",
            alpha=0.12,
            linewidth=0,
            antialiased=True,
            shade=False,
            rcount=Z_if.shape[0],
            ccount=Z_if.shape[1],
        )

        # Interface slab highlight (simulation domain only).
        y_range_slab = np.array([0.0, 0.0]) if sim_domain_plane_only else np.array([0.0, slab_depth])
        x_range_slab = np.array([-r_if, 0.0])
        X_if_slab, Y_if_slab = np.meshgrid(x_range_slab, y_range_slab)
        Z_if_slab = np.full_like(X_if_slab, self.y_interface)
        ax.plot_surface(
            X_if_slab,
            Y_if_slab,
            Z_if_slab,
            color="#88bbdd",
            alpha=0.34,
            linewidth=0,
            antialiased=True,
            shade=False,
            rcount=Z_if_slab.shape[0],
            ccount=Z_if_slab.shape[1],
        )

        # ------------------------------------------------------------------
        # Simulation slab wireframe
        # Domain is a quarter-slice: x in [0, r(y)], y in [0, slab_depth].
        # x=0 is the tank symmetry axis, y=0 is the extrusion symmetry plane.
        # Show: right-side profile outline on both y-faces + interface lines.
        # ------------------------------------------------------------------
        lw, col = 1.0, "k"

        # Left-side half of the profile outline (x <= 0)
        x_right = -r_prof  # x = -r(z)
        z_right = z_prof  # z = profile height

        sim_faces = [0.0] if sim_domain_plane_only else [0.0, slab_depth]
        for y_face in sim_faces:
            # Right-side profile outline (wall)
            ax.plot(x_right, np.full_like(x_right, y_face), z_right, color=col, lw=lw, zorder=100)
            # # Axis line on this face (x=0, full height)
            # ax.plot([0, 0], [y_face, y_face], [y_viz_start, y_viz_end], color=col, lw=lw, zorder=100)
            # Interface line on this face (axis to wall)
            ax.plot(
                [0, -r_if],
                [y_face, y_face],
                [self.y_interface, self.y_interface],
                color=col,
                lw=lw,
                ls="--",
                zorder=100,
            )

        if not sim_domain_plane_only:
            # Connecting interface edge between the two y-faces (at x=-r_if)
            ax.plot([-r_if, -r_if], [0.0, slab_depth], [self.y_interface, self.y_interface], color=col, lw=lw, ls="--")
            # Connecting interface edge at the axis (x=0)
            ax.plot([0, 0], [0.0, slab_depth], [self.y_interface, self.y_interface], color=col, lw=lw, ls="--")
            # Connecting axis edge between the two y-faces (at x=0, top and bottom)
            for z_edge in [y_viz_start, y_viz_end]:
                ax.plot([0, 0], [0.0, slab_depth], [z_edge, z_edge], color=col, lw=lw)
            # Connecting wall edge between the two y-faces (at x=-r, top and bottom)
            for z_edge, r_edge in [(y_viz_start, r_prof[0]), (y_viz_end, r_prof[-1])]:
                ax.plot([-r_edge, -r_edge], [0.0, slab_depth], [z_edge, z_edge], color=col, lw=lw)

        # Axis line — extends 10 % of tank height above and below
        tank_height = y_viz_end - y_viz_start
        ax_ext = 0.20 * tank_height
        ax.plot(
            [0, 0],
            [0, 0],
            [y_viz_start - ax_ext, y_viz_end + ax_ext],
            color="k",
            lw=1.5,
            ls="dashdot",
            alpha=0.6,
            zorder=-5,
        )

        # ------------------------------------------------------------------
        # Formatting
        # ------------------------------------------------------------------
        ax.set_axis_off()

        # Keep limits tight so the geometry fills the figure area.
        ax.set_xlim(-0.25 * diameter, 0.25 * diameter)
        ax.set_ylim(-half_len, half_len)
        ax.set_zlim(y_viz_start, y_viz_end)
        # ax.set_box_aspect([diameter, 2 * half_len, y_viz_end - y_viz_start])
        ax.set_aspect("equal")
        ax.view_init(elev=18, azim=116)

        # 3D + axis-off can still leave large default subplot margins.
        ax.set_position([0.0, 0.0, 1.0, 1.0])
        fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
        if save_path is not None:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight", pad_inches=0)
        if show:
            plt.show()
        else:
            plt.close(fig)

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
            if np.isclose(point[1], self.y_interface - self.t_BL):
                i_bl_lower = i
            if np.isclose(point[1], self.y_interface + self.t_BL):
                i_bl_upper = i
            if np.isclose(point[1], self.y_interface):
                i_bl = i

        profile_normals = self.get_profile_normals()
        # Make sure the normals at start and outlet are vertical
        profile_normals.append(np.array([0, -1]))
        profile_normals[0] = np.array([0, 1])
        profile_normals[-2] = np.array([0, -1])

        inner_points = [p + n * self.t_BL for p, n in zip(profile_points, profile_normals, strict=True)]

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

        wall_points = []  # Points for the optional wall region, outer - n*tw
        for i, point in enumerate(profile_points):
            wall_points.append(point - self.wall_thickness * profile_normals[i])

        points.append(np.array([0, self.y_interface - self.t_BL]))
        points.append(np.array([0, self.y_interface]))
        points.append(np.array([0, self.y_interface + self.t_BL]))
        axis_points = [inner_points[0], *points[-3:]]

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
        wall_thickness: float = 2.08e-3,
        extrude: bool = False,
        extrude_thickness: float = 0.0,
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
            wall_thickness=wall_thickness,
            extrude=extrude,
            extrude_thickness=extrude_thickness,
        )
        self.name = "KSite"
        self.add_boundary_layers(x_wall=wall_cell_size, r_BL=r_BL)
        self.cap_height = A
        self.cylinder_radius = B
        self.cylinder_height = C

        self.r_lid = 0.302
        self.y_lid = self.get_y(self.r_lid, 0.5 * self.ymax(), self.ymax())
        self.split_profile(self.y_lid)


class SphereProfile(TankProfile):
    def __init__(
        self,
        radius: float,
        fill_level: float,
        outlet_radius: float,
        bulk_cell_size: float,
        wall_tan_cell_size: float,
        wall_cell_size: float,
        r_BL: float = 1.2,
        internal_outlet: float = 0,
        wall_thickness: float = 2.08e-3,
        extrude: bool = False,
        extrude_thickness: float = 0.0,
    ) -> None:
        super().__init__(
            segments=[
                EllipseArc("ellipse1", 0, radius, radius, radius, length_scale=wall_tan_cell_size),
                EllipseArc(
                    "ellipse2",
                    radius,
                    radius * 2,
                    radius,
                    radius,
                    length_scale=wall_tan_cell_size,
                ),
            ],
            fill_level=fill_level,
            outlet_radius=outlet_radius,
            internal_outlet=internal_outlet,
            wall_thickness=wall_thickness,
            extrude=extrude,
            extrude_thickness=extrude_thickness,
        )
        self.name = "Sphere"
        self.add_boundary_layers(x_wall=wall_cell_size, r_BL=r_BL)
        self.cap_height = radius
        self.cylinder_radius = radius
        self.cylinder_height = 0


class CylinderCapsTankProfile(TankProfile):
    """
    General tank profile with a cylindrical midsection and ellipsoidal caps.

    The tank geometry consists of three sections stacked along the y-axis:

    1. **Bottom cap** - an ellipsoidal cap of height ``cap_height`` and
       equatorial radius ``cylinder_radius`` (semi-minor axis = ``cap_height``,
       semi-major axis = ``cylinder_radius``).
    2. **Cylinder** - a straight cylindrical section of height
       ``cylinder_height`` and radius ``cylinder_radius``.
    3. **Top cap** - a mirror of the bottom cap.

    Total tank height = ``2 * cap_height + cylinder_height``.

    This is a generalisation of :class:`KSiteProfile` that accepts
    user-supplied dimensions instead of the hard-coded K-Site values.

    Parameters
    ----------
    cylinder_radius:
        Radius of the cylindrical section (and the equatorial radius of each
        ellipsoidal cap).  Alternatively supply ``cylinder_diameter`` and this
        value will be computed as ``cylinder_diameter / 2``.
    cylinder_height:
        Height (axial length) of the cylindrical section.  May be zero for a
        purely spheroidal tank.
    cap_height:
        Height of each ellipsoidal cap (semi-minor axis in the axial
        direction).  When ``cap_height == cylinder_radius`` the caps are
        hemispheres.
    fill_level:
        Fractional liquid fill level in [0, 1].
    outlet_radius:
        Radius of the outlet opening at the top of the tank.
    bulk_cell_size:
        Target cell size in the bulk mesh regions.
    wall_tan_cell_size:
        Target cell size along the wall (tangential direction).
    wall_cell_size:
        Target first-cell size normal to the wall (used for boundary-layer
        mesh generation).
    r_BL:
        Boundary-layer cell growth ratio.  Default is 1.2.
    internal_outlet:
        Depth of an internal outlet pipe extending into the tank.
        Set to 0 (default) for a flush outlet.
    cylinder_diameter:
        Optional alternative to ``cylinder_radius``.  If provided,
        ``cylinder_radius`` is set to ``cylinder_diameter / 2`` and any
        explicit ``cylinder_radius`` argument is ignored.
    """

    def __init__(
        self,
        cylinder_radius: float,
        cylinder_height: float,
        cap_height: float,
        fill_level: float,
        outlet_radius: float,
        bulk_cell_size: float,
        wall_tan_cell_size: float,
        wall_cell_size: float,
        r_BL: float = 1.2,
        internal_outlet: float = 0,
        cylinder_diameter: float | None = None,
        wall_thickness: float = 2.08e-3,
        extrude: bool = False,
        extrude_thickness: float = 0.0,
    ) -> None:
        if cylinder_diameter is not None:
            cylinder_radius = cylinder_diameter / 2.0

        segments: list[Segment] = [
            EllipseArc(
                "ellipse1",
                0,
                cap_height,
                cylinder_radius,
                cap_height,
                length_scale=wall_tan_cell_size,
            ),
        ]
        if cylinder_height > 0:
            segments.append(
                LineSegment(
                    "line1",
                    cap_height,
                    cap_height + cylinder_height,
                    cylinder_radius,
                    cylinder_radius,
                    length_scale=wall_tan_cell_size,
                )
            )
        segments.append(
            EllipseArc(
                "ellipse2",
                cap_height,
                2 * cap_height,
                cylinder_radius,
                cap_height,
                y_offset=cylinder_height,
                length_scale=wall_tan_cell_size,
            )
        )

        super().__init__(
            segments=segments,
            fill_level=fill_level,
            outlet_radius=outlet_radius,
            internal_outlet=internal_outlet,
            wall_thickness=wall_thickness,
            extrude=extrude,
            extrude_thickness=extrude_thickness,
        )
        self.name = "CylinderCaps"
        self.add_boundary_layers(x_wall=wall_cell_size, r_BL=r_BL)
        self.cap_height = cap_height
        self.cylinder_radius = cylinder_radius
        self.cylinder_height = cylinder_height
        # self.plot()


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

if __name__ == "__main__":
    tank = SphereProfile(
        radius=1,
        fill_level=0.5,
        outlet_radius=0.01,
        bulk_cell_size=0.025,
        wall_tan_cell_size=0.005,
        wall_cell_size=0.025,
        r_BL=1.2,
    )
    tank.plot()
    print(tank.volume)
    print(tank.area)
