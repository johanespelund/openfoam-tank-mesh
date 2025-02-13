# ruff: noqa
import gmsh  # type: ignore[import-untyped]
import numpy as np

from openfoam_tank_mesh import TankMesh


def run(mesh: "TankMesh.TankMesh") -> None:
    tank = mesh.tank
    r_outlet = tank.outlet_radius
    y_outlet = tank.y_outlet
    y_interface = tank.y_interface
    wedge_angle = mesh.wedge_angle
    revolve = mesh.revolve
    bulk_cell_size = mesh.bulk_cell_size
    wall_cell_size = mesh.wall_cell_size
    lc = mesh.wall_tan_cell_size
    t_BL = mesh.t_BL
    n_BL = mesh.n_BL + 1
    r_BL = mesh.r_BL

    # mesh.bulk_cell_size = bulk_cell_size/4
    # lc = mesh.bulk_cell_size
    # mesh.t_BL = mesh.t_BL/4
    # t_BL = mesh.t_BL
    # n_BL = n_BL

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
    major_point_bl = add_point(a - t_BL, y_cylinder, z0, lc)
    major_point_wall = add_point(a + tw, y_cylinder, z0, lc)

    p1 = add_point(0, y_interface, z0, lc)
    p2 = add_point(0, y_outlet, z0, lc)
    x3, y3 = r_outlet, y_outlet
    p3 = add_point(r_outlet, y_outlet, z0, lc)
    pw3 = add_point(r_outlet, y_outlet + tw, z0, lc)

    x6, y6 = tank.get_radius(y_interface), y_interface
    p6 = add_point(x6, y6, z0, lc)
    pw6 = add_point(x6 + tw, y6, z0, lc)

    # normal = tank.get_normal(y4)
    # n = normal/np.linalg.norm(normal)
    # print(n)
    # x10, y10 = np.array([x4, y4]) + n*t_BL

    n = tank.get_normal(y_bl)
    t = tank.get_tangent(y_bl)
    d = n[0] + t[0]
    # x10, y10 = tank.get_radius(y_bl) + t_BL*d, y_bl
    x10, y10 = get_corner_coords(mesh)
    p10 = add_point(x10, y10, z0, lc)

    normal = tank.get_normal(y10)
    for _ in range(10):
        y4 = y10 - normal[1] * t_BL * 0.5
        x4 = tank.get_radius(y4)
        normal = tank.get_normal(y4)
    # y4 = 1.1*y_bl
    # x4 = tank.get_radius(y4)
    # tangent = tank.get_tangent(y6)
    # t = tangent/np.linalg.norm(tangent)
    # print(t)
    # x4, y4 = np.array([x6, y6]) + t*1.5*t_BL
    # print(x4, y4)
    # print(tank.get_radius(y4), y4)
    p4 = add_point(x4, y4, z0, lc)
    pw4 = add_point(x4 + tw, y4, z0, lc)

    # x7, y7 = x6 - 1.5*t_BL, y_interface
    correction_factor = abs(-1 + 2 * mesh.tank.fill_level)
    d = x6 - t_BL - x10
    x7, y7 = x10 + correction_factor * d, y_interface
    p7 = add_point(x7, y7, z0, lc)

    x8, y8 = x3, y3 - t_BL
    p8 = add_point(x8, y8, z0, lc)

    x9, y9 = 0, y_outlet - t_BL
    p9 = add_point(x9, y9, z0, lc)

    x11, y11 = 0, y_bl
    p11 = add_point(x11, y11, z0, lc)

    # Compare lengths of line between

    # Add lines

    l_11_9 = add_line(p11, p9)
    l_9_2 = add_line(p9, p2)
    l_2_3 = add_line(p2, p3)
    l_8_9 = add_line(p8, p9)

    l_10_4 = add_line(p10, p4)
    l_6_7 = add_line(p6, p7)
    l_7_1 = add_line(p7, p1)
    l_1_11 = add_line(p1, p11)
    l_10_11 = add_line(p10, p11)
    l_3_3w = add_line(p3, pw3)
    l6_6w = add_line(p6, pw6)

    l_3_4 = add_ellipse(p3, origo, major_point, p4)
    l_4_6 = add_ellipse(p4, origo, major_point, p6)
    # l_4_6 = add_line(p4, p6)
    # l_10_7 = add_ellipse(p10, origo, major_point_bl, p7)
    l_10_7 = add_line(p10, p7)
    l_8_10 = add_ellipse(p8, origo, major_point_bl, p10)
    l_3w_4w = add_ellipse(pw3, origo, major_point_wall, pw4)
    l_4w_6w = add_ellipse(pw4, origo, major_point_wall, pw6)
    # l_10_6 = add_line(p10, p6)

    gmsh.model.geo.synchronize()

    cl1 = add_curve_loop([l_2_3, l_3_4, l_4_6, l_6_7, -l_10_7, -l_8_10, l_8_9, l_9_2])
    # cl1 = add_curve_loop([
    #     l_2_3, l_3_4, -l_10_4, -l_8_10, l_8_9, l_9_2])
    # cl2 = add_curve_loop([l_10_4, l_4_6, -l_10_6])
    # cl3 = add_curve_loop([-l_10_6, l_10_7, -l_6_7])

    cl2 = add_curve_loop([l_10_7, l_7_1, l_1_11, -l_10_11])
    cl3 = add_curve_loop([l_11_9, -l_8_9, l_8_10, l_10_11])
    cl_wall = add_curve_loop([l_3_3w, l_3w_4w, l_4w_6w, -l6_6w, -l_4_6, -l_3_4])

    s1 = add_surface(cl1)
    s2 = add_surface(cl2)
    s3 = add_surface(cl3)
    # s4 = add_surface(cl4)
    # s5 = add_surface(cl5)
    swall = add_surface(cl_wall)

    N_2_3 = get_N_outlet(mesh)
    print(f"{x7=}, {lc=}, {N_2_3=}")
    N_1_7 = closest_odd(x7 / lc) + 2

    y = np.linspace(y_outlet, y4)
    x = np.array([tank.get_radius(yi) for yi in y])
    L_3_4 = np.sum([np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2) for i in range(len(y) - 1)])
    N_3_4 = closest_odd(L_3_4 / lc) + 2

    gmsh.model.geo.mesh.setTransfiniteCurve(l_9_2, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_2_3, N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_8_9, N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_6_7, n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_1_11, n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_7_1, N_1_7)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_10_11, N_1_7)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_3_4, N_3_4)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_3w_4w, N_3_4)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_8_10, N_3_4)

    gmsh.model.geo.mesh.setTransfiniteCurve(l_3_3w, nw)
    gmsh.model.geo.mesh.setTransfiniteCurve(l6_6w, nw)

    # gmsh.model.geo.mesh.setTransfiniteCurve(l_6_7, n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_4_6, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_4w_6w, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_10_7, n_BL, "Progression", -r_BL)
    # gmsh.model.geo.mesh.setTransfiniteCurve(l_10_6, n_BL, "Progression", -r_BL)
    # gmsh.model.geo.mesh.setTransfiniteCurve(l_10_4, n_BL, "Progression", -r_BL)

    gmsh.model.geo.mesh.setTransfiniteSurface(s1, "Left", [p2, p6, p7, p9])
    gmsh.model.geo.mesh.setTransfiniteSurface(s2, "Left", [p1, p11, p10, p7])
    # gmsh.model.geo.mesh.setTransfiniteSurface(s3, "Left", [p11, p9, p10])
    # print(f"{s3=}")
    # gmsh.model.geo.mesh.setTransfiniteSurface(s2, "Right", [p10, p4, p6])
    # gmsh.model.geo.mesh.setTransfiniteSurface(s3, "Right", [p6, p7, p10])
    # gmsh.model.geo.mesh.setTransfiniteSurface(s4, "Left", [p1, p11, p10, p7])
    # gmsh.model.geo.mesh.setTransfiniteSurface(s5, "Left", [p11, p9, p10])
    gmsh.model.geo.mesh.setTransfiniteSurface(swall, "Left", [p3, pw3, pw6, p6])
    gmsh.model.geo.synchronize()

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [])
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [l_8_9, l_8_10, l_10_11, l_11_9])
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [])
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", mesh.wall_tan_cell_size)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", mesh.bulk_cell_size)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 2 * mesh.wall_tan_cell_size)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 3 * mesh.bulk_cell_size)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.geo.synchronize()

    for s in [s1, s2, swall]:
        gmsh.model.geo.mesh.setRecombine(2, s)
    if not mesh.tri_bulk:
        gmsh.model.geo.mesh.setRecombine(2, s3)
        gmsh.option.setNumber("Mesh.Algorithm", 8)  # 5 or 6
    else:
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # 5 or 6

    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)  # 2 or 3
    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(2)
    gmsh.model.geo.synchronize()
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.recombine()
    gmsh.model.geo.synchronize()
    angle = 2 * np.pi * revolve / 360 if revolve else wedge_angle * np.pi / 180
    n_angle = closest_odd(2 * np.pi * revolve / (360 * lc)) if revolve else 1
    _ = gmsh.model.geo.revolve(
        [(2, s1), (2, s2), (2, s3), (2, swall)],
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

    # for i in range(len(surfaces)):
    #     # Format an int string with leading zeros
    #     ind = i
    #     int_string = f"s_{ind:02d}"
    #     print(f"Adding physical surface {int_string}")
    #     add_physical_surface([ind], int_string)

    add_physical_surface([0, 1, 2, 3, 4], "cyclic_pos_gmsh")
    add_physical_surface([11, 14, 15, 20], "cyclic_neg_gmsh")
    add_physical_surface([7, 12, 19], "bottom_gmsh")
    add_physical_surface([17, 18], "walls")
    add_physical_surface([4], "outlet")
    add_physical_surface([16], "metal_outlet")

    # # if y_bl > y_cylinder:
    # #     add_physical_surface([0, 1, 2, 3, 4], "cyclic_pos_gmsh")
    # #     add_physical_surface([10, 16, 19, 20, 26], "cyclic_neg_gmsh")
    # #     add_physical_surface([22, 23, 24], "walls_gmsh")
    # #     add_physical_surface([5], "outlet")
    # #     add_physical_surface([21], "metal_outlet")
    # #     add_physical_surface([13, 17, 25], "bottom_gmsh")
    # # else:
    # #     add_physical_surface([0, 1, 2, 3, 4], "cyclic_pos_gmsh")
    # #     add_physical_surface([12, 16, 19, 20, 26], "cyclic_neg_gmsh")
    # #     add_physical_surface([22, 23, 24], "walls_gmsh")
    # #     add_physical_surface([5], "outlet")
    # #     add_physical_surface([21], "metal_outlet")
    # #     add_physical_surface([14, 17, 25], "bottom_gmsh")

    # # gmsh.model.geo.synchronize()
    # # gmsh.model.geo.synchronize()

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


def get_N_outlet(mesh: "TankMesh.TankMesh") -> int:
    if mesh.wall_tan_cell_size >= mesh.outlet_radius:
        return 2
    else:
        return closest_odd(np.ceil(mesh.outlet_radius / mesh.wall_tan_cell_size))


def gmsh_setup() -> None:
    gmsh.initialize()
    gmsh.logger.start()
    gmsh.model.add("KSite49")
    gmsh.option.setNumber("Geometry.Tolerance", 1e-9)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.TransfiniteTri", 1)


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


def get_corner_coords(mesh: "TankMesh.TankMesh") -> tuple[float, float]:
    """
    Get the corner coords of the boundary layer.
    """

    A, B = mesh.tank.cylinder_radius, mesh.tank.cap_height
    C = mesh.tank.cylinder_height
    A -= mesh.t_BL
    B -= mesh.t_BL

    print(f"{A=}, {B=}, {C=}")

    def r_ellipse(y: float) -> float:
        if y > C / 2:
            print(f"1 {y=}, {C=}")
            print(f"{1 - (y - C / 2) ** 2 / B**2=}")
            return float(A * np.sqrt(1 - (y - C / 2) ** 2 / B**2))
        elif y > -C / 2:
            print(f"2 {y=}, {C=}")
            return A
        else:
            print(f"3 {y=}, {C=}")
            return float(A * np.sqrt(1 - (y + C / 2) ** 2 / B**2))

    print(f"{mesh.tank.y_interface=}")
    y = mesh.tank.y_interface + mesh.t_BL
    print(f"{y=}")
    print(f"{r_ellipse(y)=}")
    return r_ellipse(y), y


def print_debug(mesh: "TankMesh.TankMesh", msg: str) -> None:
    if mesh.debug:
        print(msg)
