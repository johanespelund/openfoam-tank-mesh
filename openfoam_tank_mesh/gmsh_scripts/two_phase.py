from __future__ import annotations

import gmsh  # type: ignore[import-untyped]
import numpy as np
import matplotlib.pyplot as plt

from openfoam_tank_mesh.Profile import EllipseArc, LineSegment, TankProfile, Profile

from openfoam_tank_mesh import KSiteMesh
from openfoam_tank_mesh.gmsh_scripts.utilities import (
    add_curve_loop,
    add_ellipse,
    add_line,
    add_physical_surface,
    add_point,
    add_surface,
    closest_odd,
    get_N_outlet,
    gmsh_setup,
    print_debug,
)

def get_coords(pointID: int) -> tuple[float, float]:
    """
    Get the coordinates of a point in the gmsh model.
    """
    x, y, z = gmsh.model.getValue(0, pointID, [])
    return x, y


def find_point(coords, tol=1e-6):
    """
    Loop through all points in the gmsh model,
    and return the point label for the point which is within
    a tolerance distance from coords.
    """

    # Get all points
    points = gmsh.model.getEntities(dim=0)
    coords = (coords[0], coords[1], 0)
    for point in points[::-1]:
        # Get the coordinates of the point
        x, y, z = gmsh.model.getValue(0, point[1], [])
        # Check if the point is within the tolerance distance from coords
        if np.allclose([x, y, z], coords, atol=tol):
            return point[1]

    return -1


def find_line(start, end, tol=1e-6):
    """
    Loop through all lines in the gmsh model,
    and return the line label for the line which is within
    a tolerance distance from start and end.
    """

    # Get all lines
    lines = gmsh.model.getEntities(dim=1)
    start = (start[0], start[1], 0)
    end = (end[0], end[1], 0)
    for line in lines:
        # Get the start and end points of the line
        # p1, p2 = gmsh.model.getValue(1, line[1], [])
        # print(f"Line {line[1]}: ({p1}, {p2})")
        result = gmsh.model.getBoundary([line], oriented=True)
        assert len(result) == 2, "Line should have two points"
        i1 = result[0][1]
        i2 = result[1][1]

        p1 = gmsh.model.getValue(0, i1, [])
        p2 = gmsh.model.getValue(0, i2, [])

        # Need to account for the fact that the line may be reversed,
        # if so, return negative label value.

        if np.allclose([p1, p2], [start, end], atol=tol):
            return line[1]
        elif np.allclose([p2, p1], [start, end], atol=tol):
            return -line[1]
    print(f"Line not found: ({start}, {end})")
    return -1
    # raise ValueError("Line not found")


def sort_xy(points):
    x = np.array([_point[0] for _point in points])
    y = np.array([_point[1] for _point in points])

    x0 = np.mean(x)
    y0 = np.mean(y)

    r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
    angles = np.where((y - y0) > 0, np.arccos((x - x0) / r), 2 * np.pi - np.arccos((x - x0) / r))

    mask = np.argsort(angles)

    x_sorted = x[mask]
    y_sorted = y[mask]

    # Find the index of the point with lowest y, breaking ties with largest x
    min_y = np.min(y_sorted)
    candidates = np.where(y_sorted == min_y)[0]
    i0 = candidates[np.argmax(x_sorted[candidates])]

    # Rotate arrays
    _x = np.concatenate((x_sorted[i0:], x_sorted[:i0]))
    _y = np.concatenate((y_sorted[i0:], y_sorted[:i0]))

    _points = [(x, y, 0) for x, y in zip(_x, _y)]

    return _points


def run(mesh: KSiteMesh.KSiteMesh) -> None:
    tank = mesh.tank
    y_outlet = tank.y_outlet
    y_interface = tank.y_interface
    wedge_angle = mesh.wedge_angle
    revolve = mesh.revolve
    wall_cell_size = mesh.wall_cell_size
    lc = mesh.wall_tan_cell_size
    n_BL = mesh.n_BL + 1
    r_BL = mesh.r_BL

    debug = mesh.debug

    nw = 10

    gmsh_setup()

    p, lines = generate_points_and_lines(mesh)


    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write("mesh.msh")
    # log = gmsh.logger.get()
    # with open("log.gmsh", "w") as f:
    #     f.writelines([line + "\n" for line in log])

    # Run the GUI
    if debug:
        gmsh.fltk.run()
        gmsh.finalize()


def generate_points_and_lines(
    mesh: KSiteMesh.KSiteMesh,
) -> tuple[dict[str, int], dict[str, int]]:
    """
    Generate points and lines.
    """
    tank = mesh.tank
    r_outlet = tank.outlet_radius
    y_outlet = tank.y_outlet
    y_interface = tank.y_interface
    revolve = mesh.revolve
    wedge_angle = mesh.wedge_angle
    bulk_cell_size = mesh.bulk_cell_size
    lc = mesh.wall_tan_cell_size
    t_BL = mesh.t_BL

    a, c = tank.cylinder_radius, tank.cylinder_height
    b = tank.cap_height
    y_cylinder = c / 2

    y_bl = y_interface + t_BL
    z0 = 0

    tw = 2e-3

    p = {}
    lines = {}

    origo = add_point(0, b + y_cylinder, z0, lc)
    major_point = add_point(a, b + y_cylinder, z0, lc)
    y_cylinder_liq = add_point(a, b - y_cylinder, z0, lc)
    major_point_wall = add_point(a + tw, b + y_cylinder, z0, lc)
    p["origo"] = origo
    p["major_point"] = major_point
    p["y_cylinder_liq"] = y_cylinder_liq
    majorPoint = {
        "liquid": add_point(a, b, z0, lc),
        "interface_liquid": add_point(a, b, z0, lc),
        "interface_gas": add_point(a, b + y_cylinder, z0, lc),
        "gas": add_point(a, b + y_cylinder, z0, lc),
    }
    origoGroup = {
        "liquid": add_point(0, b, z0, lc),
        "interface_liquid": add_point(0, b, z0, lc),
        "interface_gas": add_point(0, b + y_cylinder, z0, lc),
        "gas": add_point(0, b + y_cylinder, z0, lc),
    }

    # Liquid phase points:
    tank_profile = mesh.tank #create_tank_profile(mesh)
    n_segments = len(tank_profile.segments)

    profile_points = tank_profile.get_mesh_points()

    outer_points = profile_points["outer_points"]
    inner_points = profile_points["inner_points"]
    wall_points = profile_points["wall_points"]
    outlet_points = profile_points["outlet_points"]
    internal_outlet_points = profile_points["internal_outlet_points"]
    axis_points = profile_points["axis_points"]
    i_bl = profile_points["i_bl"]

    for k, v in profile_points.items():
        try:
            int(k)
            p[k] = add_point(v[0], v[1], z0, lc)
        except Exception:
            pass

    for point in p:
        gmsh.model.setEntityName(0, p[point], point)

    gmsh.model.geo.synchronize()

    curve_groups = tank_profile.get_curve_groups()
    lines = {}
    line_groups = {group: [] for group in list(curve_groups.keys()) + ["outlet", "internal_outlet"]}
    normal_lines = []
    wall_lines = []

    i = 0
    # for _ in range(2):
    #     for group in curve_groups:
    #         for seg in curve_groups[group]:
    #             if isinstance(seg, EllipseArc):
    #                 lines[f"{i}-{i + 1}"] = add_ellipse(p[str(i)], origoGroup[group], majorPoint[group], p[str(i + 1)])
    #             elif isinstance(seg, LineSegment):
    #                 lines[f"{i}-{i + 1}"] = add_line(p[str(i)], p[str(i + 1)])
    #             line_groups[group].append(lines[f"{i}-{i + 1}"])
    #             i += 1
    #     lines[f"{i}-{i + 1}"] = add_line(p[str(i)], p[str(i + 1)])
    #     line_groups["outlet"].append(lines[f"{i}-{i + 1}"])
    #     i += 2

    # line_groups["internal_outlet"].append(line_groups["outlet"][-1])

    i = 0
    for point_group in [outer_points, inner_points, wall_points]:
        pts = point_group
        for group in curve_groups:
            for seg in curve_groups[group]:
                p1 = find_point(pts[i])
                p2 = find_point(pts[i + 1])
                if isinstance(seg, EllipseArc):
                    add_ellipse(
                        p1, origoGroup[group], majorPoint[group], p2
                    )
                elif isinstance(seg, LineSegment):
                    add_line(p1, p2)
                i += 1
        i = 0


    # Create lines between inner and outer points at start and end of groups
    i = 0
    for group in curve_groups:
        p1 = find_point(outer_points[i])
        p2 = find_point(inner_points[i])
        add_line(p1, p2)
        n_seg = len(curve_groups[group])
        i += n_seg

    for _ in range(2):
        p1 = find_point(outer_points[i])
        p2 = find_point(inner_points[i])
        add_line(p1, p2)
        i += 1
    p1 = find_point(outer_points[-2])
    p2 = find_point(outer_points[-1])
    add_line(p1, p2)


    # Crreate lines between outer and wall points for only first and last point
    wall_normal_curves = []
    for i in (0, len(wall_points) - 1):
        p1 = find_point(outer_points[i])
        p2 = find_point(wall_points[i])
        wall_normal_curves.append(add_line(p1, p2))

    # Create outlet lines
    for i in range(len(outlet_points)):
        i1 = i
        i2 = (i + 1) % len(outlet_points)

        # p1 = find_point(outlet_points[i1])
        # p2 = find_point(outlet_points[i2])
        # add_line(p1, p2)
    # Create outlet line_

        p1 = find_point(internal_outlet_points[i1])
        p2 = find_point(internal_outlet_points[i2])
        add_line(p1, p2)

    for i in range(len(axis_points) - 1):
        p1 = find_point(axis_points[i])
        p2 = find_point(axis_points[i + 1])
        add_line(p1, p2)

    for i in [-1, 0, 1]:
        p1 = find_point(axis_points[2 + i])
        p2 = find_point(inner_points[i_bl + i])
        add_line(p1, p2)


    # for _ in range(2):
    #     lines[f"{i}-{i + n_segments + 2}"] = add_line(p[str(i)], p[str(i + n_segments + 2)])
    #     normal_lines.append(lines[f"{i}-{i + n_segments + 2}"])
    #     i += 1

    # p1 = n_segments + 2
    # p2 = 2 * (n_segments + 2)
    # lines[f"{p1}-{p2}"] = add_line(p[str(p1)], p[str(p2)])

    # for _ in range(3):
    #     p1 = p2
    #     p2 = p1 + 1
    #     lines[f"{p1}-{p2}"] = add_line(p[str(p1)], p[str(p2)])

    # p1 = p2
    # p2 = p1 + 1
    # lines[f"{p1}-{p2}"] = add_line(p[str(p1)], p[str(p2)])
    # line_groups["internal_outlet"].append(lines[f"{p1}-{p2}"])

    # p1 = p2
    # p2 = 2 * (n_segments + 1)
    # lines[f"{p1}-{p2}"] = add_line(p[str(p1)], p[str(p2)])
    # line_groups["internal_outlet"].append(lines[f"{p1}-{p2}"])

    # p1 -= 1
    # p2 = 2 * (n_segments + 1) + 1
    # lines[f"{p1}-{p2}"] = add_line(p[str(p1)], p[str(p2)])
    # line_groups["internal_outlet"].append(lines[f"{p1}-{p2}"])

    gmsh.model.geo.synchronize()

    for l in wall_normal_curves:
        gmsh.model.geo.mesh.setTransfiniteCurve(l, 7)

    outlet_line = find_line(outer_points[-2], outer_points[-1])
    N_outlet = get_N_outlet(mesh)
    gmsh.model.geo.mesh.setTransfiniteCurve(outlet_line, N_outlet)

    for l in line_groups["internal_outlet"]:
        result = gmsh.model.getBoundary([[1, l]], oriented=True)
        i1 = result[0][1]
        i2 = result[1][1]
        p1 = gmsh.model.getValue(0, i1, [])
        p2 = gmsh.model.getValue(0, i2, [])

        d = np.linalg.norm(np.array(p1) - np.array(p2))

        if p1[1] == p2[1]:
            N = N_outlet
        else:
            N = closest_odd(d / lc)
        gmsh.model.geo.mesh.setTransfiniteCurve(l, N, "Progression", 1)


    for sgn in [-1, 0, 1]:
        i1 = find_point((0, tank_profile.y_interface + sgn * tank_profile.t_BL))
        p1 = list(p.keys())[list(p.values()).index(i1)]
        if sgn == -1:
            key = "i_bl_lower"
        elif sgn == 0:
            key = "i_bl"
        else:
            key = "i_bl_upper"
        (
            x,
            y,
        ) = inner_points[profile_points[key]]
        i2 = find_point((x, y))
        p2 = list(p.keys())[list(p.values()).index(i2)]

        lines[f"{p1}-{p2}"] = add_line(p[str(p1)], p[str(p2)])
        gmsh.model.geo.synchronize()

    gmsh.model.geo.synchronize()

    ## LIQUID REGION
    _points = []
    for i in range(profile_points["i_bl_lower"] + 1):
        _points.append(inner_points[i])
    _points.append((0, tank_profile.y_interface - tank_profile.t_BL))

    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clLiquid = gmsh.model.geo.addCurveLoops(_lines)
    sLiquid = gmsh.model.geo.addPlaneSurface(clLiquid)
    line_groups["liquid"] = [abs(l) for l in _lines]

    ## GAS REGION
    _points = []
    for i in range(profile_points["i_bl_upper"], len(inner_points) - 1):
        _points.append(inner_points[i])
    # Add the inner_outlet_point
    _points.append((r_outlet, profile_points["y_int_outlet"]))
    _points.append((0, profile_points["y_int_outlet"]))
    _points.append((0, tank_profile.y_interface + tank_profile.t_BL))

    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clGas = int(gmsh.model.geo.addCurveLoops(_lines))
    sGas = add_surface(clGas)
    line_groups["gas"] = [abs(l) for l in _lines]

    ## OUTER LIQUID BOUNDARY LAYER
    _points = []
    for i in range(profile_points["i_bl_lower"], profile_points["i_bl"] + 1):
        _points.append(outer_points[i])
        _points.append(inner_points[i])

    _points = sort_xy(_points)
    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clOuterLiquidBL = gmsh.model.geo.addCurveLoops(_lines)
    sOuterLiquidBL = gmsh.model.geo.addPlaneSurface(clOuterLiquidBL)

    cornersOuterLiquidBL = [
        find_point(p)
        for p in [
            outer_points[profile_points["i_bl_lower"]],
            outer_points[profile_points["i_bl"]],
            inner_points[profile_points["i_bl"]],
            inner_points[profile_points["i_bl_lower"]],
        ]
    ]

    # INNER LIQUID BOUNDARY LAYER
    _points = []
    for i in range(profile_points["i_bl_lower"], profile_points["i_bl"] + 1):
        _points.append(inner_points[i])
    _points.append((0, tank_profile.y_interface))
    _points.append((0, tank_profile.y_interface - tank_profile.t_BL))

    _points = sort_xy(_points)
    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clInnerLiquidBL = gmsh.model.geo.addCurveLoops(_lines)
    sInnerLiquidBL = gmsh.model.geo.addPlaneSurface(clInnerLiquidBL)

    cornersInnerLiquidBL = [
        find_point(p)
        for p in [
            inner_points[profile_points["i_bl_lower"]],
            inner_points[profile_points["i_bl"]],
            (0, tank_profile.y_interface),
            (0, tank_profile.y_interface - tank_profile.t_BL),
        ]
    ]

    ## OUTER GAS BOUNDARY LAYER
    _points = []
    for i in range(profile_points["i_bl"], profile_points["i_bl_upper"] + 1):
        _points.append(outer_points[i])
        _points.append(inner_points[i])

    _points = sort_xy(_points)
    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clOuterGasBL = gmsh.model.geo.addCurveLoops(_lines)
    sOuterGasBL = gmsh.model.geo.addPlaneSurface(clOuterGasBL)

    cornersOuterGasBL = [
        find_point(p)
        for p in [
            outer_points[profile_points["i_bl"]],
            outer_points[profile_points["i_bl_upper"]],
            inner_points[profile_points["i_bl_upper"]],
            inner_points[profile_points["i_bl"]],
        ]
    ]

    # INNER GAS BOUNDARY LAYER
    _points = []
    for i in range(profile_points["i_bl"], profile_points["i_bl_upper"] + 1):
        _points.append(inner_points[i])
    _points.append((0, tank_profile.y_interface))
    _points.append((0, tank_profile.y_interface + tank_profile.t_BL))

    _points = sort_xy(_points)
    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clInnerGasBL = gmsh.model.geo.addCurveLoops(_lines)
    sInnerGasBL = gmsh.model.geo.addPlaneSurface(clInnerGasBL)

    cornersInnerGasBL = [
        find_point(p)
        for p in [
            inner_points[profile_points["i_bl"]],
            inner_points[profile_points["i_bl_upper"]],
            (0, tank_profile.y_interface + tank_profile.t_BL),
            (0, tank_profile.y_interface),
        ]
    ]

    ## GAS WALL
    _points = []
    for i in range(profile_points["i_bl_upper"], len(inner_points) - 1):
        _points.append(outer_points[i])
    for i in reversed(range(profile_points["i_bl_upper"], len(inner_points) - 1)):
        _points.append(inner_points[i])
    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clGasWall = gmsh.model.geo.addCurveLoops(_lines)
    sGasWall = gmsh.model.geo.addPlaneSurface(clGasWall)
    cornersGasWall = [
        find_point(p)
        for p in [
            outer_points[profile_points["i_bl_upper"]],
            outer_points[len(outer_points) - 2],
            inner_points[len(inner_points) - 2],
            inner_points[profile_points["i_bl_upper"]],
        ]
    ]

    ## LIQUID WALL
    _points = []
    for i in range(0, profile_points["i_bl_lower"] + 1):
        _points.append(outer_points[i])
    for i in reversed(range(0, profile_points["i_bl_lower"] + 1)):
        _points.append(inner_points[i])

    # TODO: Sort function not working when mass center is outside surface!

    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clLiquidWall = gmsh.model.geo.addCurveLoops(_lines)
    sLiquidWall = gmsh.model.geo.addPlaneSurface(clLiquidWall)
    cornersLiquidWall = [
        find_point(p)
        for p in [
            outer_points[0],
            outer_points[profile_points["i_bl_lower"]],
            inner_points[profile_points["i_bl_lower"]],
            inner_points[0],
        ]
    ]

    ## INTERNAL OUTLET
    _points = [
        (r_outlet, profile_points["y_int_outlet"]),
        inner_points[-2],
        inner_points[-1],
        (0, profile_points["y_int_outlet"]),
    ]
    _points = sort_xy(_points)
    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clInternalOutlet = gmsh.model.geo.addCurveLoops(_lines)
    sInternalOutlet = gmsh.model.geo.addPlaneSurface(clInternalOutlet)

    ## OUTLET
    _points = [
        outer_points[-2],
        outer_points[-1],
        inner_points[-1],
        inner_points[-2],
    ]
    _points = sort_xy(_points)
    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    clOutlet = gmsh.model.geo.addCurveLoops(_lines)
    sOutlet = gmsh.model.geo.addPlaneSurface(clOutlet)

    # Add wall normal transfinite curves
    for i in range(len(profile_points["outer_points"])):
        p1 = outer_points[i]
        p2 = inner_points[i]
        l = find_line(p1, p2)
        sign = l/abs(l)  # Get the sign of the line, to determine direction
        print(l, sign)
        gmsh.model.geo.mesh.setTransfiniteCurve(l, tank_profile.N + 1, "Progression", sign*mesh.r_BL)
        # TODO: All of these lines are not defined actually!

    # Add wall tangential transfinite curves
    for i, seg in enumerate(tank_profile.segments):
        inner_line = find_line(
            inner_points[i],
            inner_points[i + 1]
        )
        outer_line = find_line(
            outer_points[i],
            outer_points[i + 1],
        )
        wall_line = find_line(
            wall_points[i],
            wall_points[i + 1],
        )
        N = seg.N
        r_BL = seg.r
        gmsh.model.geo.mesh.setTransfiniteCurve(inner_line, N + 1, "Progression", r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(outer_line, N + 1, "Progression", r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(wall_line, N + 1, "Progression", r_BL)

    # Finally, set N for horizontal interface (+BL) lines:
    N_hor = closest_odd(inner_points[profile_points["i_bl"]][0] / lc)
    bl_up = find_line(
        (0, tank_profile.y_interface + tank_profile.t_BL),
        inner_points[profile_points["i_bl_upper"]]
    )
    interface_line = find_line((0, tank_profile.y_interface), inner_points[profile_points["i_bl"]])
    bl_down = find_line(
        (0, tank_profile.y_interface - tank_profile.t_BL),
        inner_points[profile_points["i_bl_lower"]],
    )
    gmsh.model.geo.mesh.setTransfiniteCurve(bl_up, N_hor)
    gmsh.model.geo.mesh.setTransfiniteCurve(interface_line, N_hor)
    gmsh.model.geo.mesh.setTransfiniteCurve(bl_down, N_hor)

    # Need to set for boundary layer lines on y axis (above and below interface)
    above = find_line(
        (0, tank_profile.y_interface),
        (0, tank_profile.y_interface + tank_profile.t_BL),
    )
    below = find_line(
        (0, tank_profile.y_interface),
        (0, tank_profile.y_interface - tank_profile.t_BL),
    )

    N_above = tank_profile.N + 1  # - tank_profile.n_upper_bl_segments
    gmsh.model.geo.mesh.setTransfiniteCurve(above, N_above, "Progression", mesh.r_BL)

    N_below = tank_profile.N + 1  # - tank_profile.n_lower_bl_segments
    gmsh.model.geo.mesh.setTransfiniteCurve(below, N_below, "Progression", -mesh.r_BL)


    ## WALL REGION
    _points = []

    for i in range(len(wall_points)):
        _points.append(wall_points[i])
    for i in reversed(range(len(wall_points))):
        _points.append(outer_points[i])

    # _points = sort_xy(_points)
    print(_points)
    _lines = [find_line(_points[i], _points[(i + 1) % len(_points)]) for i in range(len(_points))]
    print(_lines)
    clWall = gmsh.model.geo.addCurveLoops(_lines)
    sWall = gmsh.model.geo.addPlaneSurface(clWall)
    cornersWall = [find_point(p) for p in [
        outer_points[0],
        wall_points[0],
        wall_points[-1],
        outer_points[-2]
    ]]





    gmsh.model.geo.synchronize()

    gmsh.model.geo.mesh.setTransfiniteSurface(sOuterLiquidBL, "Left", cornersOuterLiquidBL)
    gmsh.model.geo.mesh.setTransfiniteSurface(sInnerLiquidBL, "Left", cornersInnerLiquidBL)
    gmsh.model.geo.mesh.setTransfiniteSurface(sOuterGasBL, "Left", cornersOuterGasBL)
    gmsh.model.geo.mesh.setTransfiniteSurface(sInnerGasBL, "Left", cornersInnerGasBL)
    gmsh.model.geo.mesh.setTransfiniteSurface(sGasWall, "Left", cornersGasWall)
    gmsh.model.geo.mesh.setTransfiniteSurface(sLiquidWall, "Left", cornersLiquidWall)
    gmsh.model.geo.mesh.setTransfiniteSurface(sInternalOutlet, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(sOutlet, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(sWall, "Left", cornersWall)

    for s in [
        sOuterLiquidBL,
        sOuterGasBL,
        sInnerLiquidBL,
        sInnerGasBL,
        sLiquidWall,
        sGasWall,
        sInternalOutlet,
        sOutlet,
        sWall
    ] + [sGas, sLiquid]:
        # sInternalOutlet] + [sGas, sLiquid]:
        gmsh.model.geo.mesh.setRecombine(2, s)
        # gmsh.model.geo.mesh.setRecombine(2, s4)
    # gmsh.option.setNumber("Mesh.Algorithm", 6)  # 5 or 6
    # else:
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # 5 or 6
    # gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)  # 2 or 3
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)  # 2 or 3
    # gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.model.geo.synchronize()

    bc = mesh.bulk_cell_size
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [])
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", line_groups["gas"] + line_groups["liquid"])
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [])
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", bc)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 2 * r_outlet)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 3 * bc)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.geo.synchronize()
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    # We now want to extrude the wall mesh outwards to generate a new
    # 2D mesh, which we call wall:
    # TODO: Probably easier to just mesh it manually!

    # surfaces: list[tuple[int, int]] = gmsh.model.getEntities(dim=2)
    angle = 2 * np.pi * revolve / 360 if revolve else wedge_angle * np.pi / 180
    n_angle = closest_odd(2 * np.pi * revolve / (360 * lc)) if revolve else 1

    regionSurfaces = {
        "gas": [sGas, sOuterGasBL, sInnerGasBL, sGasWall, sOutlet, sInternalOutlet],
        "liquid": [sLiquid, sOuterLiquidBL, sInnerLiquidBL, sLiquidWall],
        "metal": [sWall],
    }

    regionVolumes = {}

    for region in regionSurfaces:
        surfaces = [(2, s) for s in regionSurfaces[region]]
        result = gmsh.model.geo.revolve(
                surfaces,
                0, 0, 0, 0, 1, 0, angle, numElements=[n_angle], recombine=True)
        gmsh.model.geo.synchronize()
        regionVolumes[region] = [res[1] for res in result if res[0] == 3]

    # gas_volumes = []
    # liquid_volumes = []
    # for v in volumes:
    #     _, ymin, _, _, ymax, _ = gmsh.model.getBoundingBox(*v)
    #     y = (ymin + ymax) / 2
    #     if y > y_interface:
    #         gas_volumes.append(v[1])
    #     else:
    #         liquid_volumes.append(v[1])

    # gas = gmsh.model.addPhysicalGroup(3, gas_volumes)
    # liquid = gmsh.model.addPhysicalGroup(3, liquid_volumes)
    # gmsh.model.setPhysicalName(3, gas, "gas")
    # gmsh.model.setPhysicalName(3, liquid, "liquid")

    gas = gmsh.model.addPhysicalGroup(3, regionVolumes["gas"])
    liquid = gmsh.model.addPhysicalGroup(3, regionVolumes["liquid"])
    walls = gmsh.model.addPhysicalGroup(3, regionVolumes["metal"])
    gmsh.model.setPhysicalName(3, gas, "gas")
    gmsh.model.setPhysicalName(3, liquid, "liquid")
    gmsh.model.setPhysicalName(3, walls, "metal")

    gmsh.model.mesh.generate(3)
    gmsh.model.geo.synchronize()
    # gmsh.model.mesh.recombine()
    gmsh.model.mesh.optimize()

    return p, lines

