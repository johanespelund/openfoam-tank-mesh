from __future__ import annotations

import gmsh  # type: ignore[import-untyped]
import numpy as np

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
    bulk_cell_size = mesh.bulk_cell_size
    wall_cell_size = mesh.wall_cell_size
    lc = bulk_cell_size
    t_BL = mesh.t_BL
    n_BL = mesh.n_BL + 1
    r_BL = mesh.r_BL

    debug = mesh.debug

    c = tank.cylinder_height
    y_cylinder = c / 2
    y_bl = y_interface + t_BL

    nw = 10

    gmsh_setup()

    p, lines = generate_points_and_lines(mesh)
    y4, y5 = tank.y_interface + mesh.t_BL, tank.cylinder_height / 2
    x7 = tank.get_radius(y4)

    if y4 < y5:
        print_debug(mesh, "y4<y5 / y_BL < y_cylinder")
        cl1 = add_curve_loop([
            lines["2_3"],
            lines["3_5"],
            -lines["4_5"],
            -lines["10_4"],
            -lines["5B_10"],
            -lines["8_5B"],
            lines["8_9"],
            lines["9_2"],
        ])
        cl2 = add_curve_loop([lines["10_4"], lines["4_6"], lines["6_7"], lines["7_10"]])
        cl3 = add_curve_loop([-lines["7_10"], lines["7_1"], lines["1_11"], -lines["10_11"]])
        cl4 = add_curve_loop([-lines["8_9"], lines["8_5B"], lines["5B_10"], lines["10_11"], lines["11_9"]])
        cl_wall = add_curve_loop([
            lines["3_3w"],
            lines["3_5w"],
            -lines["4_5w"],
            lines["4_6w"],
            -lines["6_6w"],
            -lines["4_6"],
            lines["4_5"],
            -lines["3_5"],
        ])
    else:
        cl1 = add_curve_loop([lines["2_3"], lines["3_4"], -lines["10_4"], -lines["8_10"], lines["8_9"], lines["9_2"]])
        cl2 = add_curve_loop([lines["10_4"], lines["4_5"], lines["5_6"], lines["6_7"], -lines["5B_7"], lines["5B_10"]])
        cl3 = add_curve_loop([-lines["5B_10"], lines["5B_7"], lines["7_1"], lines["1_11"], -lines["10_11"]])
        cl4 = add_curve_loop([-lines["8_9"], lines["8_10"], lines["10_11"], lines["11_9"]])
        cl_wall = add_curve_loop([
            lines["3_3w"],
            lines["3_4w"],
            lines["4_5w"],
            lines["5_6w"],
            -lines["6_6w"],
            -lines["5_6"],
            -lines["4_5"],
            -lines["3_4"],
        ])

    s1 = add_surface(cl1)
    s2 = add_surface(cl2)
    s3 = add_surface(cl3)
    s4 = add_surface(cl4)
    swall = add_surface(cl_wall)

    N_2_3 = get_N_outlet(mesh)
    N_1_7 = closest_odd(x7 / lc) + 2

    gmsh.model.geo.mesh.setTransfiniteCurve(lines["9_2"], n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["2_3"], N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["8_9"], N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["6_7"], n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["10_4"], n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["1_11"], n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["7_1"], N_1_7)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["10_11"], N_1_7)

    gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_3w"], nw)
    gmsh.model.geo.mesh.setTransfiniteCurve(lines["6_6w"], nw)

    if y4 > y5:
        y = np.linspace(y_outlet, y4)
        x = np.array([tank.get_radius(yi) for yi in y])
        L_3_4 = np.sum([np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2) for i in range(len(y) - 1)])
        N_3_4 = closest_odd(L_3_4 / lc) + 2
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_4"], N_3_4)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_4w"], N_3_4)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["8_10"], N_3_4)
        gmsh.model.geo.mesh.setTransfiniteSurface(s1, "Left", [p["2"], p["4"], p["10"], p["9"]])
        N_6_5 = 2
        t = wall_cell_size
        i = 1
        while t < y5 - y_interface and i < n_BL - 2:
            t += wall_cell_size * r_BL**i
            i += 1
            N_6_5 += 1
        N_6_5 -= 0
        N_5_4 = n_BL - N_6_5 + 1
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["5B_10"], N_5_4, "Progression", r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_5"], N_5_4, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_5w"], N_5_4, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["5_6"], N_6_5, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["5_6w"], N_6_5, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["5B_7"], N_6_5, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteSurface(s2, "Left", [p["6"], p["7"], p["10"], p["4"]])
        gmsh.model.geo.mesh.setTransfiniteSurface(s3, "Left", [p["11"], p["10"], p["7"], p["1"]])
    else:
        y = np.linspace(y_outlet, y5)
        x = np.array([tank.get_radius(yi) for yi in y])
        L_3_5 = np.sum([np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2) for i in range(len(y) - 1)])
        N_3_5 = closest_odd(max(1, int(np.floor(np.around(L_3_5 / lc))) + 1))
        L_4_5 = y5 - y4
        N_4_5 = closest_odd(max(1, int(np.floor(np.around(L_4_5 / lc))) + 1))

        gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_5"], N_3_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["3_5w"], N_3_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_5"], N_4_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_5w"], N_4_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["8_5B"], N_3_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["5B_10"], N_4_5)
        # gmsh.model.geo.mesh.setTransfiniteCurve(l["10_4"], n_BL, "Progression", -r_BL)
        # gmsh.model.geo.mesh.setTransfiniteCurve(l["6_7"], n_BL, "Progression", r_BL)
        gmsh.model.geo.mesh.setTransfiniteSurface(s1, "Left", [p["2"], p["4"], p["10"], p["9"]])
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_6"], n_BL, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["4_6w"], n_BL, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(lines["7_10"], n_BL, "Progression", r_BL)
        gmsh.model.geo.mesh.setTransfiniteSurface(s2, "Left", [p["6"], p["7"], p["10"], p["4"]])
        gmsh.model.geo.mesh.setTransfiniteSurface(s3, "Left", [p["11"], p["10"], p["7"], p["1"]])

    gmsh.model.geo.mesh.setTransfiniteSurface(swall, "Left", [p["3"], p["w3"], p["w6"], p["6"]])
    gmsh.model.geo.synchronize()
    # gmsh.model.mesh.generate(2)

    gmsh.model.geo.synchronize()

    # # points = gmsh.model.getEntities(dim=0)
    # # lines = gmsh.model.getEntities(dim=1)
    # # surfaces = gmsh.model.getEntities(dim=2)
    # # volumes = gmsh.model.getEntities(dim=3)
    # # gmsh.model.geo.translate(points, 0, y_offset, 0)
    # # gmsh.model.geo.translate(lines, 0, y_offset, 0)
    # # gmsh.model.geo.translate(surfaces, 0, y_offset, 0)
    # # gmsh.model.geo.translate(volumes, 0, y_offset, 0)

    # # Set meshing algorithm to Frontal-Delaunay for quads
    # # # gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1) # 2 or 3

    # # # Synchronize to finalize geometry operations

    # # # Generate the 2D mesh
    # # # gmsh.model.mesh.optimize('Laplace2D')
    # # # gmsh.model.mesh.optimize('Laplace2D')

    # # gmsh.model.geo.mesh.setRecombine(2, s1)
    # # gmsh.model.geo.mesh.setRecombine(2, s2)
    # # gmsh.model.geo.mesh.setRecombine(2, s3)
    # # gmsh.model.geo.mesh.setRecombine(2, s4)

    for s in [s1, s2, s3, swall]:
        gmsh.model.geo.mesh.setRecombine(2, s)
    gmsh.model.geo.mesh.setRecombine(2, s4)
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # 5 or 6
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
        [(2, s1), (2, s2), (2, s3), (2, s4), (2, swall)],
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
    gas = gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes[:-1]])
    metal = gmsh.model.addPhysicalGroup(3, [volumes[-1][1]])
    gmsh.model.setPhysicalName(3, gas, "gas")
    gmsh.model.setPhysicalName(3, metal, "metal")

    # # Recombine algorithm for quads/hexes
    # # gmsh.option.setNumber("Mesh.RecombineAll", 1)

    # # Generate the 3D mesh
    gmsh.model.mesh.generate(3)
    gmsh.model.geo.synchronize()
    # gmsh.model.mesh.recombine()
    gmsh.model.mesh.optimize()

    # for i in range(len(s)):
    #     # Format an int string with leading zeros
    #     ind = i
    #     int_string = f"s_{ind:02d}"
    #     print(f"Adding physical surface {int_string}")
    #     add_physical_surface([ind], int_string)

    surfaces: list[tuple[int, int]] = gmsh.model.getEntities(dim=2)

    if y_bl > y_cylinder:
        add_physical_surface([0, 1, 2, 3, 4], "cyclic_pos_gmsh", surfaces)
        add_physical_surface([10, 16, 19, 20, 26], "cyclic_neg_gmsh", surfaces)
        add_physical_surface([22, 23, 24], "walls_gmsh", surfaces)
        add_physical_surface([5], "outlet", surfaces)
        add_physical_surface([21], "metal_outlet", surfaces)
        add_physical_surface([13, 17, 25], "bottom_gmsh", surfaces)
    else:
        add_physical_surface([0, 1, 2, 3, 4], "cyclic_pos_gmsh", surfaces)
        add_physical_surface([12, 16, 19, 20, 26], "cyclic_neg_gmsh", surfaces)
        add_physical_surface([22, 23, 24], "walls_gmsh", surfaces)
        add_physical_surface([5], "outlet", surfaces)
        add_physical_surface([21], "metal_outlet", surfaces)
        add_physical_surface([14, 17, 25], "bottom_gmsh", surfaces)

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


def generate_points_and_lines(mesh: KSiteMesh.KSiteMesh) -> tuple[dict[str, int], dict[str, int]]:
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
    major_point_wall = add_point(a + tw, y_cylinder, z0, lc)
    p["origo"] = origo
    p["major_point"] = major_point

    p["1"] = add_point(0, y_interface, z0, lc)
    p["2"] = add_point(0, y_outlet, z0, lc)
    x3, y3 = r_outlet, y_outlet
    p["3"] = add_point(r_outlet, y_outlet, z0, lc)
    p["w3"] = add_point(r_outlet, y_outlet + tw, z0, lc)

    y4 = y_bl
    x4 = tank.get_radius(y4)
    p["4"] = add_point(x4, y4, z0, lc)
    p["w4"] = add_point(x4 + tw, y4, z0, lc)

    x5, y5 = a, y_cylinder
    p["5"] = add_point(x5, y5, z0, lc)
    p["w5"] = add_point(x5 + tw, y5, z0, lc)
    p["5B"] = add_point(x5 - t_BL, y5, z0, lc)

    x6, y6 = x5, y_interface
    p["6"] = add_point(x6, y6, z0, lc)
    p["w6"] = add_point(x6 + tw, y6, z0, lc)

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

    # Add lines

    lines["11_9"] = add_line(p["11"], p["9"])
    lines["9_2"] = add_line(p["9"], p["2"])
    lines["2_3"] = add_line(p["2"], p["3"])
    lines["8_9"] = add_line(p["8"], p["9"])
    lines["10_4"] = add_line(p["10"], p["4"])

    if y4 < y5:
        lines["5B_10"] = add_line(p["5B"], p["10"])
        lines["4_5"] = add_line(p["4"], p["5"])
        lines["4_5w"] = add_line(p["w4"], p["w5"])
    else:
        lines["5B_10"] = add_line(p["5B"], p["10"])
        # l["4_5"] = add_line(p["4"], p["5"])
        lines["4_5"] = add_ellipse(p["4"], origo, major_point, p["5"])
        lines["4_5w"] = add_line(p["w4"], p["w5"])

    lines["6_7"] = add_line(p["6"], p["7"])
    lines["7_1"] = add_line(p["7"], p["1"])
    lines["1_11"] = add_line(p["1"], p["11"])
    lines["10_11"] = add_line(p["10"], p["11"])
    lines["3_3w"] = add_line(p["3"], p["w3"])
    lines["6_6w"] = add_line(p["6"], p["w6"])

    if y4 < y5:
        lines["3_5"] = add_ellipse(p["3"], origo, major_point, p["5"])
        lines["8_5B"] = add_ellipse(p["8"], origo, major_point, p["5B"])
        lines["4_6"] = add_line(p["4"], p["6"])
        lines["7_10"] = add_line(p["7"], p["10"])
        lines["3_5w"] = add_ellipse(p["w3"], origo, major_point_wall, p["w5"])
        lines["4_6w"] = add_line(p["w4"], p["w6"])
    else:
        lines["3_4"] = add_ellipse(p["3"], origo, major_point, p["4"])
        lines["3_4w"] = add_ellipse(p["w3"], origo, major_point_wall, p["w4"])
        lines["8_10"] = add_ellipse(p["8"], origo, major_point, p["10"])
        lines["5_6"] = add_line(p["5"], p["6"])
        lines["5_6w"] = add_line(p["w5"], p["w6"])
        lines["5B_7"] = add_line(p["5B"], p["7"])

    return p, lines
