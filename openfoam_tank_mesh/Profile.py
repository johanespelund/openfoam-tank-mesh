import numpy as np
from abc import ABC, abstractmethod
import scipy.integrate as spi  # type: ignore[import-untyped]
import scipy.optimize as spo  # type: ignore[import-untyped]
from openfoam_tank_mesh.exceptions import OutOfRange
import matplotlib.pyplot as plt
import copy

def calculate_boundary_layer(r_BL, wall_cell_size, wall_tan_cell_size):
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
    """
    Return the closest even integer greater than or equal to n.
    """
    return max(2, int(n) // 2 * 2 + 2)

class Segment(ABC):
    """
    Base class for tank segments. Will represent
    either lines, arc or ellipse segments.
    """

    def __init__(self, name: str, y_start: float) -> None:
        self.name = name
        self.y_start = y_start
        self.y_end = None
        self.r_start = None
        self.r_end = None
        self.length = None
        self.upperNeighbor = None
        self.lowerNeighbor = None
        self.length_scale = 0
        self.N = 0  # Number of cells in the segment
        self.r = 1  # Growth rate of the cells
        super().__init__()

    @abstractmethod
    def get_radius(self, y: float) -> float:
        """
        Get the radius at a given height y, where y is in the range [-height/2, height/2].
        """
        pass

    @abstractmethod
    def get_radius_derivative(self, y: float) -> float:
        """
        Get the derivative of the radius at a given height y, where y is in the range [y_start, y_end].
        """
        pass

    @abstractmethod
    def get_length(self) -> float:
        """
        Get the length of the segment.
        """
        pass

    def __str__(self) -> str:
        return f"{self.name}:\n  y: {self.y_start} - {self.y_end}\n  r: {self.r_start} - {self.r_end}\n  N: {self.N}, r: {self.r}"

    def get_tangent(self, y: float) -> np.ndarray:
        """
        Return the normalized tangent to the curve r(y),
        using the derivative of the curve.
        """
        dy_dx = self.get_radius_derivative(y)
        norm = np.sqrt(dy_dx**2 + 1)
        return np.array([dy_dx / norm, 1 / norm])

    def get_normal(self, y: float) -> np.ndarray:
        """
        Return the normalized normal to the curve r(y),
        using the derivative of the curve.
        Points towards the center of the tank.
        """
        dy_dx = self.get_radius_derivative(y)
        norm = np.sqrt(dy_dx**2 + 1)
        return np.array([-1 / norm, dy_dx / norm])


class LineSegment(Segment):
    """
    Line segment class.
    """

    def __init__(
        self,
        name: str,
        y_start: float,
        y_end: float,
        r_start: float,
        r_end: float,
        length_scale: float = 0,
    ) -> None:
        super().__init__(name, y_start)
        self.y_end = y_end
        self.r_start = r_start
        self.r_end = r_end
        self.length_scale = length_scale
        self.N = closest_even(self.get_length() / length_scale)

    def get_radius(self, y: float) -> float:
        """
        Get the radius at a given height y, where y is in the range [y_start, y_end].
        """
        if y < self.y_start or y > self.y_end:
            raise OutOfRange(y)
        return self.r_start + (self.r_end - self.r_start) * (y - self.y_start) / (self.y_end - self.y_start)

    def get_radius_derivative(self, y: float) -> float:
        return (self.r_end - self.r_start) / (self.y_end - self.y_start)

    def get_length(self) -> float:
        return np.sqrt((self.r_end - self.r_start) ** 2 + (self.y_end - self.y_start) ** 2)


class EllipseArc(Segment):
    def __init__(
        self,
        name: str,
        y_start: float,
        y_end: float,
        axis_major: float,  # major axis (horizontal)
        axis_minor: float,  # minor axis (vertical)
        y_offset: float = 0,
        length_scale: float = 0,
    ) -> None:
        super().__init__(name, y_start)
        self.y_end = y_end
        self.axis_major = axis_major
        self.axis_minor = axis_minor
        self.y_offset = y_offset
        self.y_start += y_offset
        self.y_end += y_offset
        self.r_start = self.get_radius(self.y_start)
        self.r_end = self.get_radius(self.y_end)
        self.length_scale = length_scale
        self.N = closest_even(self.get_length() / length_scale)

        # Needed for gmsh
        self.major_point = np.array([self.axis_major, self.axis_minor + self.y_offset])
        self.origo = np.array([0, self.axis_minor + self.y_offset])

    def get_radius(self, y: float) -> float:
        """
        Get the radius at a given height y, where y is in the range [y_start, y_end].
        """
        if y < self.y_start or y > self.y_end:
            print(f"y: {y}, y_start: {self.y_start}, y_end: {self.y_end}")
            raise OutOfRange(y)

        A, B = self.axis_major, self.axis_minor
        _y = y - B - self.y_offset
        return A * np.sqrt(1 - (_y / B) ** 2)

    def get_radius_derivative(self, y: float) -> float:
        if y < self.y_start or y > self.y_end:
            print(f"y: {y}, y_start: {self.y_start}, y_end: {self.y_end}")
            raise OutOfRange(y)

        A, B = self.axis_major, self.axis_minor
        _y = y - B - self.y_offset

        denom = np.sqrt(1 - (_y / B) ** 2)
        return -A * _y / (B**2 * denom) if denom != 0 else 0

    def get_length(self) -> float:
        # Numerically calculate the length of the ellipse arc
        # using the trapezoidal rule
        y = np.linspace(self.y_start, self.y_end, 100)
        r = np.array([self.get_radius(yi) for yi in y])
        drdy = np.array([self.get_radius_derivative(yi) for yi in y])
        length = np.trapz(np.sqrt(drdy**2 + 1), y)
        return length


class CircleArc(EllipseArc):
    def __init__(
        self,
        name: str,
        y_start: float,
        y_end: float,
        radius: float,  # major axis (horizontal)
        y_offset: float = 0,
        length_scale: float = 0,
    ) -> None:
        super().__init__(
            name,
            y_start,
            y_end,
            radius,
            radius,
            y_offset=y_offset,
            length_scale=length_scale,
        )


class Profile:
    def __init__(self, segments: list[Segment]) -> None:
        self.segments = segments
        self.sort_segments()
        self.check_segment_connectivity()
        self.volume = self.get_partial_volume(self.segments[0].y_start, self.segments[-1].y_end)
        self.y_start = self.segments[0].y_start
        self.y_end = self.segments[-1].y_end
        self.i_interface = 0
        self.n_upper_bl_segments = 0
        self.n_lower_bl_segments = 0

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

    def get_profile_points(self) -> list[np.ndarray]:
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
        self.fill_level = fill_level
        self.outlet_radius = outlet_radius
        self.internal_outlet = internal_outlet
        self.t_BL = 0
        self.N = 0
        self.y_interface = self.calculate_interface_position()
        self.y_outlet = self.calculate_outlet_position()
        self.interface_radius = self.get_radius(self.y_interface)
        self.area_liquid = self.get_partial_volume(self.y_start, self.y_interface)
        self.area_gas = self.get_partial_volume(self.y_interface, self.y_end)
        self.volume_liquid = self.get_partial_volume(self.y_start, self.y_interface)
        self.volume_gas = self.get_partial_volume(self.y_interface, self.y_end)

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

        result = float(spo.least_squares(objective, 0.7 * self.y_end, bounds=(0, self.y_end)).x)
        return result

    def merge_segments(self, segment1: Segment, segment2: Segment) -> None:
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

    def ymin(self):
        """
        Get the minimum y value of the profile.
        """
        return self.segments[0].y_start

    def ymax(self):
        """
        Get the maximum y value of the profile.
        """
        return self.segments[-1].y_end

    def split_profile(self, y_split: float, tol: float = 10e-3) -> None:
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
                    print(f"Segment {segment.name} is too short, extending it downwards.")
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
                    print(f"Segment {upper_segment.name} is too short, extending it upwards.")
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
        ### Upper segment
        x_bulk = self.segments[-1].length_scale
        n, t, _ = calculate_boundary_layer(r_BL, x_wall, x_bulk)

        # Split the profile
        tol = min(5 * x_wall, 2 * x_bulk)
        bl_segment, _ = self.split_profile(self.y_interface + t, tol=tol)

        self.sort_segments()
        self.check_segment_connectivity()

        ## Lower segment
        x_bulk = self.segments[0].length_scale
        n, t, _ = calculate_boundary_layer(r_BL, x_wall, x_bulk)
        self.t_BL = t
        self.N = n

        tol = min(5 * x_wall, 2 * x_bulk)
        _, bl_segment = self.split_profile(self.y_interface - self.t_BL, tol=tol)

        tol = min(5 * x_wall, 2 * x_bulk)
        _, bl_segment = self.split_profile(self.y_interface, tol=tol)

        bl_segments = []
        for segment in self.segments:
            if segment.y_start >= self.y_interface and segment.y_end <= self.y_interface + self.t_BL:
                bl_segments.append(segment)
        self.n_upper_bl_segments = len(bl_segments)
        # Now we need to distribute the number of cells between the segments

        t, n, N = 0, 0, 0
        L = 0
        for segment in bl_segments:
            # Get the length of the segment
            length = segment.get_length()
            L += length
            N = 0
            while t < length / r_BL or N % 2 == 1 and N > 3:
                t += x_wall * r_BL**n
                n += 1
                N += 1
            segment.N = N
            segment.r = r_BL
        if n != self.N:
            segment.N += self.N - n
            if segment.N <= 1:
                new_N = segment.N + segment.lowerNeighbor.N
                new_segment = self.merge_segments(segment, segment.lowerNeighbor)
                new_segment.N = new_N
                new_segment.r = r_BL
        total_cells_in_bl = sum([seg.N for seg in bl_segments])

        # Do the same for the boundary layer below the interface
        bl_segments = []
        for segment in self.segments[::-1]:
            if segment.y_start >= self.y_interface - self.t_BL and segment.y_end <= self.y_interface:
                bl_segments.append(segment)
        self.n_lower_bl_segments = len(bl_segments)
        # Now we need to distribute the number of cells between the segments
        t, n, N = 0, 0, 0
        L = 0
        for segment in bl_segments:
            # Get the length of the segment
            length = segment.get_length()
            L += length
            N = 0
            while t < length / r_BL or N % 2 == 1 and N > 3:
                t += x_wall * r_BL**n
                n += 1
                N += 1
            segment.N = N
            segment.r = -r_BL
        if n != self.N:
            segment.N += self.N - n
            if segment.N <= 1:
                new_N = segment.N + segment.upperNeighbor.N
                new_segment = self.merge_segments(segment, segment.upperNeighbor)
                new_segment.N = new_N
                new_segment.r = -r_BL


    def insert_interface(self, tol: float = 10e-3, x_wall: float = 1e-3) -> None:
        self.add_boundary_layers(x_wall=x_wall, r_BL=1.2)

    def plot(self):
        fig, ax = plt.subplots()

        for segment in self.segments:
            y = np.linspace(segment.y_start, segment.y_end, 1000)
            r = np.array([segment.get_radius(yi) for yi in y])
            if len(y) == 0:
                y = np.array([segment.y_start, segment.y_end])
                r = np.array([segment.r_start, segment.r_end])

            ax.plot(r, y, label=segment.name)
            ax.fill_betweenx(y, r, 0, where=y < segment.y_end, alpha=0.5)

        ax.hlines(
            self.y_interface,
            0,
            self.get_radius(self.y_interface) + 0.1,
            color="k",
            label="Interface",
            ls="--",
        )

        points = self.get_mesh_points()
        for key, item in points.items():
            try:
                int(key)
                ax.plot(item[0], item[1], "ro")
                ax.text(item[0], item[1], f"{key}", fontsize=8, ha="right")
            except:
                pass

        ax.set_aspect("equal")
        plt.legend()
        plt.show()

    def get_mesh_points(self):
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

        # Interface needs to be horizontal and hayield the same t_BL
        b = self.t_BL / profile_normals[i_interface][0]
        inner_points[i_interface] = profile_points[i_interface] + b * np.array([1, 0])

        # Need to adjust the neighbotr points as well, but with optional relaxation
        r = 0.0

        for k in [-1, 1]:
            j = i_interface - k
            b = self.t_BL / profile_normals[j][0]
            norm = self.t_BL * profile_normals[j]
            hor = b * np.array([1, 0])
            inner_points[j] = profile_points[j] + r * norm + (1 - r) * hor



        # Here we will define some points to use for the optional
        # TransfiniteTri strategy in gmsh, which is used for very high/low
        # fill levels.
        # ABANDONDED

#         transTriPoints = []
#         def func(y):
#             r = self.get_radius(y)
#             n = self.get_normal(y)

#             p = np.array([r, y])
#             d = inner_points[i_interface+1] - p
#             prod = np.dot(d/np.linalg.norm(d), n)
#             return 1 - abs(prod)

        # Find root of func:
        # result = spo.fminbound(func, self.y_interface, self.y_interface + 4*self.t_BL)
        # profile_points[i_interface + 1] = np.array([self.get_radius(result), result])
        # inner_points[i_interface][0] = inner_points[i_interface + 1][0]


        points += profile_points
        points += inner_points



        tw = 2.08e-3
        wall_points = [] # Points for the optional wall region, outer - n*tw
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

        data = {str(i): p for i, p in enumerate(points)}
        data.update({
            "inner_points": inner_points,
            "outer_points": profile_points,
            "wall_points": wall_points,
            "outlet_points": inner_points[-2:] + profile_points[-2:][::-1],
            "internal_outlet_points": internal_outlet_points,
            "axis_points": axis_points,
            "i_bl_lower": i_bl_lower,
            "i_bl": i_bl,
            "i_bl_upper": i_bl_upper,
            "y_int_outlet": y_int_outlet,
        })
        return data

    def get_curve_groups(self):
        """
        Returns three groups of curves:
        liquid, interface and gas.
        """
        groups = {
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
                EllipseArc("ellipse2", A, 2 * A, B, A, y_offset=C, length_scale=wall_tan_cell_size),
            ],
            fill_level=fill_level,
            outlet_radius=outlet_radius,
            internal_outlet=internal_outlet,
        )
        self.name = "KSite"
        self.insert_interface(tol=2*wall_tan_cell_size, x_wall=wall_cell_size)
        self.cap_height = A
        self.cylinder_radius = B
        self.cylinder_height = C
