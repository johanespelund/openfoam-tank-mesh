import gmsh  # type: ignore[import-untyped]
import numpy as np

from openfoam_tank_mesh import KSiteMesh

# from typing import TYPE_CHECKING
# if TYPE_CHECKING:
#     from openfoam_tank_mesh.KSiteMesh import KSiteMesh, KSiteTank


def run(mesh: "KSiteMesh.KSiteMesh") -> None:
    tank = mesh.tank
    r_outlet = tank.outlet_radius
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

    a, c = tank.cylinder_radius, tank.cylinder_height
    y_cylinder = c / 2
    y_bl = y_interface + t_BL
    z0 = 0

    tw = 2e-3
    nw = 10

    gmsh_setup()

    origo = add_point(0, y_cylinder, z0, lc)
    major_point = add_point(a, y_cylinder, z0, lc)
    major_point_wall = add_point(a + tw, y_cylinder, z0, lc)

    p1 = add_point(0, y_interface, z0, lc)
    p2 = add_point(0, y_outlet, z0, lc)
    x3, y3 = r_outlet, y_outlet
    p3 = add_point(r_outlet, y_outlet, z0, lc)
    pw3 = add_point(r_outlet, y_outlet + tw, z0, lc)

    y4 = y_bl
    x4 = tank.get_radius(y4)
    p4 = add_point(x4, y4, z0, lc)
    pw4 = add_point(x4 + tw, y4, z0, lc)

    x5, y5 = a, y_cylinder
    p5 = add_point(x5, y5, z0, lc)
    pw5 = add_point(x5 + tw, y5, z0, lc)
    p5B = add_point(x5 - t_BL, y5, z0, lc)

    x6, y6 = x5, y_interface
    p6 = add_point(x6, y6, z0, lc)
    pw6 = add_point(x6 + tw, y6, z0, lc)

    x7, y7 = x5 - t_BL, y_interface
    p7 = add_point(x7, y7, z0, lc)

    x8, y8 = x3, y3 - t_BL
    p8 = add_point(x8, y8, z0, lc)

    x9, y9 = 0, y_outlet - t_BL
    p9 = add_point(x9, y9, z0, lc)

    x10, y10 = x7, y_bl
    p10 = add_point(x10, y10, z0, lc)

    x11, y11 = 0, y_bl
    p11 = add_point(x11, y11, z0, lc)

    # Add lines

    l_11_9 = add_line(p11, p9)
    l_9_2 = add_line(p9, p2)
    l_2_3 = add_line(p2, p3)
    l_8_9 = add_line(p8, p9)
    l_10_4 = add_line(p10, p4)

    if y4 < y5:
        l_5B_10 = add_line(p5B, p10)
        l_4_5 = add_line(p4, p5)
        l_4_5w = add_line(pw4, pw5)
    else:
        l_5B_10 = add_line(p5B, p10)
        l_4_5 = add_line(p4, p5)
        l_4_5 = add_ellipse(p4, origo, major_point, p5)
        l_4_5w = add_line(pw4, pw5)

    l_6_7 = add_line(p6, p7)
    l_7_1 = add_line(p7, p1)
    l_1_11 = add_line(p1, p11)
    l_10_11 = add_line(p10, p11)
    l_7_10 = add_line(p7, p10)
    l_3_3w = add_line(p3, pw3)
    l6_6w = add_line(p6, pw6)

    if y4 < y5:
        l_3_5 = add_ellipse(p3, origo, major_point, p5)
        l_8_5B = add_ellipse(p8, origo, major_point, p5B)
        l_4_6 = add_line(p4, p6)
        l_7_10 = add_line(p7, p10)
        l_3_5w = add_ellipse(pw3, origo, major_point_wall, pw5)
        l_4_6w = add_line(pw4, pw6)
        cl1 = add_curve_loop([l_2_3, l_3_5, -l_4_5, -l_10_4, -l_5B_10, -l_8_5B, l_8_9, l_9_2])
        cl2 = add_curve_loop([l_10_4, l_4_6, l_6_7, l_7_10])
        cl3 = add_curve_loop([-l_7_10, l_7_1, l_1_11, -l_10_11])
        cl4 = add_curve_loop([-l_8_9, l_8_5B, l_5B_10, l_10_11, l_11_9])
        cl_wall = add_curve_loop([l_3_3w, l_3_5w, -l_4_5w, l_4_6w, -l6_6w, -l_4_6, l_4_5, -l_3_5])
    else:
        l_3_4 = add_ellipse(p3, origo, major_point, p4)
        l_3_4w = add_ellipse(pw3, origo, major_point_wall, pw4)
        l_8_10 = add_ellipse(p8, origo, major_point, p10)
        l_5_6 = add_line(p5, p6)
        l_5_6w = add_line(pw5, pw6)
        l_5B_7 = add_line(p5B, p7)
        cl1 = add_curve_loop([l_2_3, l_3_4, -l_10_4, -l_8_10, l_8_9, l_9_2])
        cl2 = add_curve_loop([l_10_4, l_4_5, l_5_6, l_6_7, -l_5B_7, l_5B_10])
        cl3 = add_curve_loop([-l_5B_10, l_5B_7, l_7_1, l_1_11, -l_10_11])
        cl4 = add_curve_loop([-l_8_9, l_8_10, l_10_11, l_11_9])
        cl_wall = add_curve_loop([l_3_3w, l_3_4w, l_4_5w, l_5_6w, -l6_6w, -l_5_6, -l_4_5, -l_3_4])

    s1 = add_surface(cl1)
    s2 = add_surface(cl2)
    s3 = add_surface(cl3)
    s4 = add_surface(cl4)
    swall = add_surface(cl_wall)

    N_2_3 = get_N_outlet(mesh)
    N_1_7 = closest_odd(x7 / lc) + 2

    gmsh.model.geo.mesh.setTransfiniteCurve(l_9_2, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_2_3, N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_8_9, N_2_3)
    # gmsh.model.geo.mesh.setTransfiniteCurve(l_7_10, n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_6_7, n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_10_4, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_1_11, n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_7_1, N_1_7)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_10_11, N_1_7)

    gmsh.model.geo.mesh.setTransfiniteCurve(l_3_3w, nw)
    gmsh.model.geo.mesh.setTransfiniteCurve(l6_6w, nw)

    if y4 > y5:
        y = np.linspace(y_outlet, y4)
        x = np.array([tank.get_radius(yi) for yi in y])
        L_3_4 = np.sum([np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2) for i in range(len(y) - 1)])
        N_3_4 = closest_odd(L_3_4 / lc) + 2
        gmsh.model.geo.mesh.setTransfiniteCurve(l_3_4, N_3_4)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_3_4w, N_3_4)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_8_10, N_3_4)
        gmsh.model.geo.mesh.setTransfiniteSurface(s1, "Left", [p2, p4, p10, p9])
        N_6_5 = 2
        t = wall_cell_size
        i = 1
        while t < y5 - y_interface and i < n_BL - 2:
            t += wall_cell_size * r_BL**i
            i += 1
            N_6_5 += 1
        N_6_5 -= 0
        N_5_4 = n_BL - N_6_5 + 1
        gmsh.model.geo.mesh.setTransfiniteCurve(l_5B_10, N_5_4, "Progression", r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_4_5, N_5_4, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_4_5w, N_5_4, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_5_6, N_6_5, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_5_6w, N_6_5, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_5B_7, N_6_5, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteSurface(s2, "Left", [p6, p7, p10, p4])
        gmsh.model.geo.mesh.setTransfiniteSurface(s3, "Left", [p11, p10, p7, p1])
    else:
        y = np.linspace(y_outlet, y5)
        x = np.array([tank.get_radius(yi) for yi in y])
        L_3_5 = np.sum([np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2) for i in range(len(y) - 1)])
        N_3_5 = closest_odd(max(1, int(np.floor(np.around(L_3_5 / lc))) + 1))
        L_4_5 = y5 - y4
        N_4_5 = closest_odd(max(1, int(np.floor(np.around(L_4_5 / lc))) + 1))
        gmsh.model.geo.mesh.setTransfiniteCurve(l_3_5, N_3_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_3_5w, N_3_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_4_5, N_4_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_4_5w, N_4_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_8_5B, N_3_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_5B_10, N_4_5)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_10_4, n_BL, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteSurface(s1, "Left", [p2, p4, p10, p9])
        gmsh.model.geo.mesh.setTransfiniteCurve(l_4_6, n_BL, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_4_6w, n_BL, "Progression", -r_BL)
        gmsh.model.geo.mesh.setTransfiniteCurve(l_7_10, n_BL, "Progression", r_BL)
        gmsh.model.geo.mesh.setTransfiniteSurface(s2, "Left", [p6, p7, p10, p4])
        gmsh.model.geo.mesh.setTransfiniteSurface(s3, "Left", [p11, p10, p7, p1])

    gmsh.model.geo.mesh.setTransfiniteSurface(swall, "Left", [p3, pw3, pw6, p6])
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

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

    surfaces: list[tuple[int, int]] = gmsh.model.getEntities(dim=2)

    def add_physical_surface(ind: list[int], name: str) -> None:
        # = [surfaces[i[1]] for i in indices]
        m = gmsh.model.addPhysicalGroup(2, [surfaces[i][1] for i in ind])
        gmsh.model.setPhysicalName(2, m, name)

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

    if y_bl > y_cylinder:
        add_physical_surface([0, 1, 2, 3, 4], "cyclic_pos_gmsh")
        add_physical_surface([10, 16, 19, 20, 26], "cyclic_neg_gmsh")
        add_physical_surface([22, 23, 24], "walls_gmsh")
        add_physical_surface([5], "outlet")
        add_physical_surface([21], "metal_outlet")
        add_physical_surface([13, 17, 25], "bottom_gmsh")
    else:
        add_physical_surface([0, 1, 2, 3, 4], "cyclic_pos_gmsh")
        add_physical_surface([12, 16, 19, 20, 26], "cyclic_neg_gmsh")
        add_physical_surface([22, 23, 24], "walls_gmsh")
        add_physical_surface([5], "outlet")
        add_physical_surface([21], "metal_outlet")
        add_physical_surface([14, 17, 25], "bottom_gmsh")

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


def closest_odd(n: float) -> int:
    return int(n) // 2 * 2 + 1


def get_N_outlet(mesh: "KSiteMesh.KSiteMesh") -> int:
    if mesh.bulk_cell_size >= mesh.outlet_radius:
        return 2
    else:
        return closest_odd(np.ceil(mesh.outlet_radius / mesh.bulk_cell_size))


def gmsh_setup() -> None:
    gmsh.initialize()
    gmsh.logger.start()
    gmsh.model.add("KSite49")
    gmsh.option.setNumber("Geometry.Tolerance", 1e-9)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)


def add_point(x: float, y: float, z: float, lc: float) -> int:
    return int(gmsh.model.geo.addPoint(x, y, z, lc))


def add_line(p1: int, p2: int) -> int:
    return int(gmsh.model.geo.addLine(p1, p2))


def add_ellipse(start: int, origo: int, majorPoint: int, end: int) -> int:
    return int(gmsh.model.geo.addEllipseArc(start, origo, majorPoint, end))


def add_curve_loop(lines: list[int]) -> int:
    return int(gmsh.model.geo.addCurveLoop(lines))


def add_surface(loop: int) -> int:
    return int(gmsh.model.geo.addPlaneSurface([loop]))
