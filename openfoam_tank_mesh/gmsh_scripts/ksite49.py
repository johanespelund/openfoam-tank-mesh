from __future__ import annotations

import gmsh  # type: ignore[import-untyped]
import numpy as np
import matplotlib.pyplot as plt

from openfoam_tank_mesh.Tank import EllipseArc, LineSegment, TankProfile, Profile

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
    _, y4, _ = gmsh.model.getValue(0, p["4"], [])
    _, y5, _ = gmsh.model.getValue(0, p["5"], [])
    x7 = tank.get_radius(y4)
    cl0 = add_curve_loop(
        [
            lines["2_3"], -lines["8_3"], lines["8_9"], lines["9_2"]
        ]
    )

    if y4 < y5:
        print_debug(mesh, "y4<y5 / y_BL < y_cylinder")
        cl1 = add_curve_loop(
            [
                # lines["2_3"],
                lines["3_5"],
                -lines["4_5"],
                -lines["10_4"],
                -lines["5B_10"],
                -lines["8_5B"],
                lines["8_3"],
                # lines["8_9"],
                # lines["9_2"],
            ]
        )
        cl2 = add_curve_loop([lines["10_4"], lines["4_6"], lines["6_7"], lines["7_10"]])
        cl3 = add_curve_loop(
            [-lines["7_10"], lines["7_1"], lines["1_11"], -lines["10_11"]]
        )
        cl4 = add_curve_loop(
            [
                lines["8_5B"],
                lines["5B_10"],
                lines["10_11"],
                lines["11_12"],
                lines["12_13"],
                -lines["8_13"],
            ]
        )
        cl5 = add_curve_loop(
            [-lines["8_13"], lines["8_9"], -lines["12_9"], lines["12_13"]]
        )
        # cl_wall = add_curve_loop([
        #     lines["3_3w"],
        #     lines["3_5w"],
        #     -lines["4_5w"],
        #     lines["4_6w"],
        #     -lines["6_6w"],
        #     -lines["4_6"],
        #     lines["4_5"],
        #     -lines["3_5"],
        # ])
    else:
        cl1 = add_curve_loop(
            [
                # lines["2_3"],
                lines["3_4"],
                -lines["10_4"],
                -lines["8_10"],
                lines["8_3"]
                # lines["8_9"],
                # lines["9_2"],
            ]
        )
        cl2 = add_curve_loop(
            [
                lines["10_4"],
                lines["4_5"],
                lines["5_6"],
                lines["6_7"],
                -lines["5B_7"],
                lines["5B_10"],
            ]
        )
        cl3 = add_curve_loop(
            [
                -lines["5B_10"],
                lines["5B_7"],
                lines["7_1"],
                lines["1_11"],
                -lines["10_11"],
            ]
        )
        # cl4 = add_curve_loop([-lines["8_9"], lines["8_10"], lines["10_11"], lines["11_9"]])
        cl4 = add_curve_loop(
            [
                lines["8_10"],
                lines["10_11"],
                lines["11_12"],
                lines["12_13"],
                -lines["8_13"],
            ]
        )
        cl5 = add_curve_loop(
            [-lines["8_13"], lines["8_9"], -lines["12_9"], lines["12_13"]]
        )
        # cl_wall = add_curve_loop([
        #     lines["3_3w"],
        #     lines["3_4w"],
        #     lines["4_5w"],
        #     lines["5_6w"],
        #     -lines["6_6w"],
        #     -lines["5_6"],
        #     -lines["4_5"],
        #     -lines["3_4"],
        # ])
    # lines["1_21"] = add_line(p["1"], p["21"])
    # lines["21_22"] = add_line(p["21"], p["22"])
    # lines["22_23"] = add_line(p["22"], p["23"])
    # lines["23_24"] = add_ellipse(p["23"], origo, y_cylinder_liq, p["24"])
    # lines["24_25"] = add_line(p["24"], p["25"])
    # lines["25_6"] = add_line(p["25"], p["6"])
    # lines["7_26"] = add_line(p["7"], p["26"])
    # lines["26_24"] = add_line(p["26"], p["24"])
    # lines["26_22"] = add_ellipse(p["26"], origo, y_cylinder_liq, p["22"])
    # lines["26_21"] = add_line(p["26"], p["21"])

    cl10 = add_curve_loop(
        [
            lines["1_21"], -lines["26_21"], -lines["7_26"], lines["7_1"]
        ]
    )
    cl11 = add_curve_loop(
        [
            lines["26_24"], lines["24_25"], lines["25_6"], lines["6_7"],
            lines["7_26"]
        ]
    )
    cl12 = add_curve_loop(
        [
            lines["26_22"], lines["22_23"], lines["23_24"], -lines["26_24"]
        ]
    )
    cl13 = add_curve_loop(
        [
            lines["21_22"], -lines["26_22"], lines["26_21"]
        ]
    )

    s0 = add_surface(cl0)
    s1 = add_surface(cl1)
    s2 = add_surface(cl2)
    s3 = add_surface(cl3)
    s4 = add_surface(cl4)
    s5 = add_surface(cl5)

    s10 = add_surface(cl10)
    s11 = add_surface(cl11)
    s12 = add_surface(cl12)
    s13 = add_surface(cl13)

    # swall = add_surface(cl_wall)

    N_2_3 = get_N_outlet(mesh) + 0 
    print_debug(mesh, f"N_2_3 = {N_2_3}")
    N_1_7 = closest_odd(x7 / lc) + 2
    print_debug(mesh, f"N_1_7 = {N_1_7}")

    L_9_12 = (
        max(mesh.t_BL + 2 * mesh.wall_tan_cell_size, mesh.internal_outlet) - mesh.t_BL
    )
    N_9_12 = closest_odd(L_9_12 / lc)
    print_debug(mesh, f"N_9_12 = {N_9_12}")

    gmsh.model.geo.mesh.setTransfiniteCurve(lines["9_2"], n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["2_3"], N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["8_9"], N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["8_3"], n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["12_13"], N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["6_7"], n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["10_4"], n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["1_11"], n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["7_1"], N_1_7)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["10_11"], N_1_7)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["1_21"], n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["7_26"], n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["22_23"], n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["26_24"], n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["26_21"], N_1_7)

    # gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_3w"], nw)
    # gmsh.model.geo.mesh.setTransfiniteCurve(lines["6_6w"], nw)

    if y4 > y5:
        y = np.linspace(y_outlet, y4)
        x = np.array([tank.get_radius(yi) for yi in y])
        L_3_4 = np.sum(
            [
                np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2)
                for i in range(len(y) - 1)
            ]
        )
        N_3_4 = closest_odd(L_3_4 / lc) + 2
        print_debug(mesh, f"L_3_4 = {L_3_4}")
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_4"], N_3_4)
        # gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_4w"], N_3_4)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["8_10"], N_3_4)
        gmsh.model.geo.mesh.setTransfiniteSurface(
            s1, "Left", [p["3"], p["4"], p["10"], p["8"]]
        )
        N_6_5 = 2
        t = wall_cell_size
        i = 1
        while t < y5 - y_interface and i < n_BL - 2:
            t += wall_cell_size * r_BL**i
            i += 1
            N_6_5 += 1
        N_6_5 -= 0
        N_5_4 = n_BL - N_6_5 + 1
        print_debug(mesh, f"N_6_5 = {N_6_5}")
        gmsh.model.geo.mesh.setTransfiniteCurve(
            lines["5B_10"], N_5_4, "Progression", r_BL
        )
        gmsh.model.geo.mesh.setTransfiniteCurve(
            lines["4_5"], N_5_4, "Progression", -r_BL
        )
        # gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_5w"], N_5_4, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(
            lines["5_6"], N_6_5, "Progression", -r_BL
        )
        # gmsh.model.geo.mesh.setTransfiniteCurve(lines["5_6w"], N_6_5, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(
            lines["5B_7"], N_6_5, "Progression", -r_BL
        )
        gmsh.model.geo.mesh.setTransfiniteSurface(
            s2, "Left", [p["6"], p["7"], p["10"], p["4"]]
        )
        gmsh.model.geo.mesh.setTransfiniteSurface(
            s3, "Left", [p["11"], p["10"], p["7"], p["1"]]
        )
    else:
        y = np.linspace(y_outlet, y5)
        x = np.array([tank.get_radius(yi) for yi in y])
        L_3_5 = np.sum(
            [
                np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2)
                for i in range(len(y) - 1)
            ]
        )
        N_3_5 = closest_odd(max(1, int(np.floor(np.around(L_3_5 / lc))) + 1))
        L_4_5 = y5 - y4
        N_4_5 = closest_odd(max(1, int(np.floor(np.around(L_4_5 / lc))) + 1))

        gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_5"], N_3_5)
        # gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_5w"], N_3_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_5"], N_4_5)
        # gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_5w"], N_4_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["8_5B"], N_3_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["5B_10"], N_4_5)
        # gmsh.model.geo.mesh.setTransfiniteCurve(l["10_4"], n_BL, "Progression", -r_BL)
        # gmsh.model.geo.mesh.setTransfiniteCurve(l["6_7"], n_BL, "Progression", r_BL)
        gmsh.model.geo.mesh.setTransfiniteSurface(
            s1, "Left", [p["3"], p["4"], p["10"], p["8"]]
        )
        gmsh.model.geo.mesh.setTransfiniteCurve(
            lines["4_6"], n_BL, "Progression", -r_BL
        )
        # gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_6w"], n_BL, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(
            lines["7_10"], n_BL, "Progression", r_BL
        )
        gmsh.model.geo.mesh.setTransfiniteSurface(
            s2, "Left", [p["6"], p["7"], p["10"], p["4"]]
        )
        gmsh.model.geo.mesh.setTransfiniteSurface(
            s3, "Left", [p["11"], p["10"], p["7"], p["1"]]
        )

    gmsh.model.geo.mesh.setTransfiniteCurve(lines["12_9"], N_9_12)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["8_13"], N_9_12)
    gmsh.model.geo.mesh.setTransfiniteSurface(s5, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(s0, "Left")
    # gmsh.model.geo.mesh.setTransfiniteSurface(swall, "Left", [p["3"], p["w3"], p["w6"], p["6"]])
    gmsh.model.geo.synchronize()
    # gmsh.model.mesh.generate(2)
    N_23_24 = N_3_4

    gmsh.model.geo.mesh.setTransfiniteCurve(lines["23_24"], N_23_24)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["26_22"], N_23_24)
    gmsh.model.geo.mesh.setTransfiniteSurface(s10, "Left")
    gmsh.model.geo.mesh.setTransfiniteSurface(s12, "Left")

    gmsh.model.geo.mesh.setTransfiniteCurve(lines["25_6"], 6, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["24_25"], n_BL - 5, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteSurface(s11, "Left", [p["6"], p["7"], p["26"], p["24"]])

    if y4 < y5:
        curves_list = [
            lines[c]
            for c in ["8_9", "8_5B", "5B_10", "10_11", "11_12", "12_13", "8_13"]
        ]
    else:
        # curves_list = [lines[c] for c in ["8_9", "8_10", "10_11", "11_9"]]
        curves_list = [lines[c] for c in ["8_10", "10_11", "11_12", "12_13", "8_13"]]

    curves_list += [
        lines[c] for c in ["21_22", "26_22", "26_21"]]

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [])
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", curves_list)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [])
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", mesh.wall_tan_cell_size)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", mesh.bulk_cell_size)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 1 * mesh.tank.outlet_radius)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 4 * mesh.bulk_cell_size)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.geo.synchronize()

    # for s in [s1, s2, s3, s5, swall]:
    for s in [s0, s1, s2, s3, s5, s10, s11, s12, s13]:
        gmsh.model.geo.mesh.setRecombine(2, s)
    if not mesh.tri_bulk:
        gmsh.model.geo.mesh.setRecombine(2, s4)
        gmsh.option.setNumber("Mesh.Algorithm", 8)  # 5 or 6
    else:
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # 5 or 6
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)  # 2 or 3
    gmsh.model.geo.synchronize()

    # gmsh.model.mesh.generate(2)
    gmsh.model.geo.synchronize()
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.recombine()
    gmsh.model.geo.synchronize()
    angle = 2 * np.pi * revolve / 360 if revolve else wedge_angle * np.pi / 180
    n_angle = closest_odd(2 * np.pi * revolve / (360 * lc)) if revolve else 1
    _ = gmsh.model.geo.revolve(
        [(2, s1), (2, s2), (2, s3), (2, s4), (2, s5),
         (2, s10), (2, s11), (2, s12), (2, s13)],
        0,
        0,
        0,  # Point on the axis of revolution
        0,
        1,
        0,  # Direction of the axis of revolution
        angle,  # Angle of revolution
        numElements=[n_angle],
        recombine=True,
    )

    gmsh.model.geo.synchronize()

    volumes = gmsh.model.getEntities(dim=3)
    gas = gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes[:]])
    # metal = gmsh.model.addPhysicalGroup(3, [volumes[-1][1]])
    gmsh.model.setPhysicalName(3, gas, "gas")
    # gmsh.model.setPhysicalName(3, metal, "metal")

    # # Generate the 3D mesh
    gmsh.model.mesh.generate(3)
    gmsh.model.geo.synchronize()
    # gmsh.model.mesh.recombine()
    gmsh.model.mesh.optimize()

    surfaces: list[tuple[int, int]] = gmsh.model.getEntities(dim=2)

    # if y4 > y5:
    #     add_physical_surface([0, 1, 2, 3, 4, 5], "cyclic_pos_gmsh", surfaces)
    #     add_physical_surface([11, 17, 20, 23, 24, 30], "cyclic_neg_gmsh", surfaces)
    #     add_physical_surface([26, 27, 28], "walls_gmsh", surfaces)
    #     add_physical_surface([6], "outlet", surfaces)
    #     add_physical_surface([25], "metal_outlet", surfaces)
    #     add_physical_surface([14, 18, 29], "bottom_gmsh", surfaces)
    # else:
    #     add_physical_surface([0, 1, 2, 3, 4, 5], "cyclic_pos_gmsh", surfaces)
    #     add_physical_surface([13, 17, 20, 23, 24, 30], "cyclic_neg_gmsh", surfaces)
    #     add_physical_surface([26, 27, 28], "walls_gmsh", surfaces)
    #     add_physical_surface([6], "outlet", surfaces)
    #     add_physical_surface([25], "metal_outlet", surfaces)
    #     add_physical_surface([15, 18, 29], "bottom_gmsh", surfaces)

    gmsh.model.geo.synchronize()
    gmsh.model.geo.synchronize()

    # # Set the MshFileVersion option to 2.2 (for MSH 2 format)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write("KSite49.msh")
    log = gmsh.logger.get()
    with open("log.gmsh", "w") as f:
        f.writelines([line + "\n" for line in log])

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
    bulk_cell_size = mesh.bulk_cell_size
    lc = bulk_cell_size
    t_BL = mesh.t_BL

    a, c = tank.cylinder_radius, tank.cylinder_height
    y_cylinder = c / 2
    y_bl = y_interface + t_BL
    z0 = 0

    tw = 2e-3

    p = {}
    lines = {}

    origo = add_point(0, y_cylinder, z0, lc)
    major_point = add_point(a, y_cylinder, z0, lc)
    y_cylinder_liq = add_point(a, -y_cylinder, z0, lc)
    major_point_wall = add_point(a + tw, y_cylinder, z0, lc)
    p["origo"] = origo
    p["major_point"] = major_point
    p["y_cylinder_liq"] = y_cylinder_liq

    p["1"] = add_point(0, y_interface, z0, lc)
    p["2"] = add_point(0, y_outlet, z0, lc)
    x3, y3 = r_outlet, y_outlet
    p["3"] = add_point(r_outlet, y_outlet, z0, lc)

    if abs(y_cylinder - y_bl) < mesh.wall_tan_cell_size * 2:
        y_cylinder = y_bl + mesh.wall_tan_cell_size * 2

    y4 = y_bl
    x4 = tank.get_radius(y4)
    p["4"] = add_point(x4, y4, z0, lc)

    x5, y5 = a, y_cylinder
    p["5"] = add_point(x5, y5, z0, lc)
    p["5B"] = add_point(x5 - t_BL, y5, z0, lc)

    x6, y6 = x5, y_interface
    p["6"] = add_point(x6, y6, z0, lc)

    x7, y7 = x5 - t_BL, y_interface
    p["7"] = add_point(x7, y7, z0, lc)

    x8, y8 = x3, y3 - t_BL
    p["8"] = add_point(x8, y8, z0, lc)

    x9, y9 = 0, y_outlet - t_BL
    p["9"] = add_point(x9, y9, z0, lc)

    x10, y10 = x7, y_bl
    p["10"] = add_point(x10, y10, z0, lc)

    x11, y11 = 0, y_bl
    p["11"] = add_point(x11, y11, z0, lc)

    x12, y12 = 0, y_outlet - max(
        t_BL + 2 * mesh.wall_tan_cell_size, mesh.internal_outlet
    )
    p["12"] = add_point(x12, y12, z0, lc)

    x13, y13 = r_outlet, y_outlet - max(
        t_BL + 2 * mesh.wall_tan_cell_size, mesh.internal_outlet
    )
    p["13"] = add_point(x13, y13, z0, lc)


    # Liquid phase points:
    tank_profile = create_tank_profile(mesh)

    p["21"] = add_point(0, y_interface - t_BL, z0, lc)
    p["22"] = add_point(0, mesh.tank.y1 + t_BL, z0, lc)
    p["23"] = add_point(0, mesh.tank.y1, z0, lc)
    # p["24"] = add_point(a, -y_cylinder, z0, lc)

    y_interface = tank.y_interface
    y_BL = tank.y_interface - t_BL
    y_cylinder = -tank.cylinder_height / 2 - 1e-3

    print_debug(mesh, f"{y_interface=}, {y_BL=}, {y_cylinder=}")

    sorting_points = [y_interface, y_BL, y_cylinder]
    sorting_points.sort()
    x24, y24 = tank.get_radius(sorting_points[0]), sorting_points[0]
    x25, y25 = tank.get_radius(sorting_points[1]), sorting_points[1]
    print_debug(mesh, f"{x24=}, {y24=}, {x25=}, {y25=}")
    p["24"] = add_point(mesh.tank.get_radius(sorting_points[0]), sorting_points[0], z0, lc)
    p["25"] = add_point(mesh.tank.get_radius(sorting_points[1]), sorting_points[1], z0, lc)

    # gmsh.model.setEntityName(0, p["24"], "p24")
    # gmsh.model.setEntityName(0, p["25"], "p25")

    p["26"] = add_point(mesh.tank.get_radius(y_BL) - t_BL, y_BL, z0, lc)
    

    for point in p:
        gmsh.model.setEntityName(0, p[point], point)


    # if t_BL > abs(y_interface + y_cylinder):
    #     _x = tank.get_radius(y_interface - t_BL)
    #     p["25"] = add_point(_x, y_interface - t_BL, z0, lc)
    # else:
    #     _x = tank.get_radius(y_interface - t_BL)
    #     p["25"] = add_point(_x, y_interface - t_BL, z0, lc)
    #     # p["25"] = add_point(a, y_interface - t_BL, z0, lc)
    # p["26"] = add_point(a - t_BL, y_interface - t_BL, z0, lc)

    # Add lines
    lines["11_12"] = add_line(p["11"], p["12"])
    lines["12_9"] = add_line(p["12"], p["9"])
    lines["8_13"] = add_line(p["8"], p["13"])
    lines["8_3"] = add_line(p["8"], p["3"])
    lines["12_13"] = add_line(p["12"], p["13"])
    lines["9_2"] = add_line(p["9"], p["2"])
    lines["2_3"] = add_line(p["2"], p["3"])
    lines["8_9"] = add_line(p["8"], p["9"])
    lines["10_4"] = add_line(p["10"], p["4"])

    if y4 < y5:
        lines["5B_10"] = add_line(p["5B"], p["10"])
        lines["4_5"] = add_line(p["4"], p["5"])
    else:
        lines["5B_10"] = add_line(p["5B"], p["10"])
        lines["4_5"] = add_ellipse(p["4"], origo, major_point, p["5"])

    lines["6_7"] = add_line(p["6"], p["7"])
    lines["7_1"] = add_line(p["7"], p["1"])
    lines["1_11"] = add_line(p["1"], p["11"])
    lines["10_11"] = add_line(p["10"], p["11"])

    if y4 < y5:
        lines["3_5"] = add_ellipse(p["3"], origo, major_point, p["5"])
        lines["8_5B"] = add_ellipse(p["8"], origo, major_point, p["5B"])
        lines["4_6"] = add_line(p["4"], p["6"])
        lines["7_10"] = add_line(p["7"], p["10"])
    else:
        lines["3_4"] = add_ellipse(p["3"], origo, major_point, p["4"])
        lines["8_10"] = add_ellipse(p["8"], origo, major_point, p["10"])
        lines["5_6"] = add_line(p["5"], p["6"])
        lines["5B_7"] = add_line(p["5B"], p["7"])


    # Liquid phase lines

    lines["1_21"] = add_line(p["1"], p["21"])
    lines["21_22"] = add_line(p["21"], p["22"])
    lines["22_23"] = add_line(p["22"], p["23"])
    lines["23_24"] = add_ellipse(p["23"], origo, y_cylinder_liq, p["24"])
    # lines["24_25"] = add_line(p["24"], p["25"])
    lines["24_25"] = add_ellipse(p["24"], origo, y_cylinder_liq, p["25"])
    lines["25_6"] = add_line(p["25"], p["6"])
    lines["7_26"] = add_line(p["7"], p["26"])
    lines["26_24"] = add_line(p["26"], p["24"])
    lines["26_22"] = add_ellipse(p["26"], origo, y_cylinder_liq, p["22"])
    lines["26_21"] = add_line(p["26"], p["21"])


    return p, lines


def create_tank_profile(mesh: KSiteMesh.KSiteMesh) -> TP.TankProfile:
    INCH = 0.0254
    A = 0.5 * 73 * INCH
    B = 0.5 * 87.6 * INCH
    C = 1.5 * INCH
    X_BULK = mesh.wall_tan_cell_size

    ellipse1 = EllipseArc(
        name="ellipse1",
        y_start=0,
        y_end=A,
        axis_major=B,
        axis_minor=A,
        length_scale=mesh.wall_tan_cell_size,
    )

    line1 = LineSegment(name="line1", y_start=A, y_end=A + C, r_start=B, r_end=B, length_scale=X_BULK)

    ellipse2 = EllipseArc(
        name="ellipse2",
        y_start=A,
        y_end=2 * A - 0.01,
        axis_major=B,
        axis_minor=A,
        y_offset=C,
        length_scale=X_BULK
    )


    tank_profile = TankProfile(
        segments=[ellipse1, line1, ellipse2],
        fill_level=0.49,
        outlet_radius=0.5,
    )

    tank_profile.insert_interface(tol=0.02, x_wall=mesh.wall_cell_size)
