# ruff: noqa
# type: ignore
import gmsh  # type: ignore[import-untyped]
import numpy as np

from openfoam_tank_mesh import KSiteMesh


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
        y4 = y10 - normal[1] * t_BL
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
    x7, y7 = x10, y_interface
    p7 = add_point(x7, y7, z0, lc)

    x8, y8 = x3, y3 - t_BL
    p8 = add_point(x8, y8, z0, lc)

    x9, y9 = 0, y_outlet - t_BL
    p9 = add_point(x9, y9, z0, lc)

    x11, y11 = 0, y_bl
    p11 = add_point(x11, y11, z0, lc)

    # Compare lengths of line between (10, 4) and (10, 7):
    d_10_4 = np.sqrt((x10 - x4) ** 2 + (y10 - y4) ** 2)
    d_10_7 = np.sqrt((x10 - x7) ** 2 + (y10 - y7) ** 2)
    d_6_10 = np.sqrt((x6 - x10) ** 2 + (y6 - y10) ** 2)

    r_6_10 = r_BL * d_6_10 / d_10_4
    print(f"{d_10_4=}, {d_10_7=} {d_6_10=}, {r_6_10=}")

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
    l_10_6 = add_line(p10, p6)

    gmsh.model.geo.synchronize()

    # cl1 = add_curve_loop([l_2_3, l_3_4, l_4_6, l_6_7, -l_10_7, -l_8_10, l_8_9, l_9_2])
    cl1 = add_curve_loop([l_2_3, l_3_4, -l_10_4, -l_8_10, l_8_9, l_9_2])
    cl2 = add_curve_loop([l_10_4, l_4_6, -l_10_6])
    cl3 = add_curve_loop([-l_10_6, l_10_7, -l_6_7])

    cl4 = add_curve_loop([l_10_7, l_7_1, l_1_11, -l_10_11])
    cl5 = add_curve_loop([l_11_9, -l_8_9, l_8_10, l_10_11])
    cl_wall = add_curve_loop([l_3_3w, l_3w_4w, l_4w_6w, -l6_6w, -l_4_6, -l_3_4])

    s1 = add_surface(cl1)
    s2 = add_surface(cl2)
    s3 = add_surface(cl3)
    s4 = add_surface(cl4)
    s5 = add_surface(cl5)
    swall = add_surface(cl_wall)

    N_2_3 = get_N_outlet(mesh)
    N_1_7 = closest_odd(x7 / lc) + 2

    y = np.linspace(y_outlet, y4)
    x = np.array([tank.get_radius(yi) for yi in y])
    L_3_4 = np.sum([np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2) for i in range(len(y) - 1)])
    N_3_4 = closest_odd(L_3_4 / lc) + 2

    print(f"{N_2_3=}, {N_1_7=}, {N_3_4=}, {n_BL=}")

    gmsh.model.geo.mesh.setTransfiniteCurve(l_9_2, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_2_3, N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_8_9, N_2_3)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_1_11, n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_7_1, N_1_7)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_10_11, N_1_7)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_3_4, N_3_4)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_3w_4w, N_3_4)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_8_10, N_3_4)

    gmsh.model.geo.mesh.setTransfiniteCurve(l_3_3w, nw)
    gmsh.model.geo.mesh.setTransfiniteCurve(l6_6w, nw)

    # gmsh.model.geo.mesh.setTransfiniteCurve(l_6_7, n_BL, "Progression", r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_4w_6w, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_10_7, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_10_6, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_10_4, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_4_6, n_BL, "Progression", -r_BL)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_6_7, n_BL, "Progression", r_BL)

    gmsh.model.geo.mesh.setTransfiniteSurface(s1, "Left", [p2, p4, p10, p9])
    gmsh.model.geo.mesh.setTransfiniteSurface(s2, "Right", [p10, p4, p6])
    # gmsh.model.geo.mesh.setTransfiniteSurface(s2, "Left", [p6, p10, p4])
    gmsh.model.geo.mesh.setTransfiniteSurface(s3, "Right", [p6, p7, p10])
    gmsh.model.geo.mesh.setTransfiniteSurface(s4, "Left", [p1, p11, p10, p7])
    # gmsh.model.geo.mesh.setTransfiniteSurface(s5, "Left", [p11, p9, p10])
    gmsh.model.geo.mesh.setTransfiniteSurface(swall, "Left", [p3, pw3, pw6, p6])
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

    for s in [s1, s2, s3, s4, s5, swall]:
        gmsh.model.geo.mesh.setRecombine(2, s)
    gmsh.model.geo.mesh.setRecombine(2, s3)
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # 5 or 6
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)  # 2 or 3
    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(2)
    gmsh.model.geo.synchronize()

    # Get the nodes of the mesh

    gmsh.model.geo.synchronize()
    angle = 2 * np.pi * revolve / 360 if revolve else wedge_angle * np.pi / 180
    n_angle = closest_odd(2 * np.pi * revolve / (360 * lc)) if revolve else 1

    # n_angle = 5

    _ = gmsh.model.geo.revolve(
        [(2, s1), (2, s2), (2, s3), (2, s4), (2, s5), (2, swall)],
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
    gmsh.model.mesh.generate(3)
    gmsh.model.geo.synchronize()

    # volumes = gmsh.model.getEntities(dim=3)

    # # ps2 = gmsh.model.addPhysicalGroup(2, [s2])
    # # ps3 = gmsh.model.addPhysicalGroup(2, [s3])

    # # Create physical volumes, corresponding to s2 and s3 revolved
    # pv2 = gmsh.model.addPhysicalGroup(3, [volumes[1][1]])
    # pv3 = gmsh.model.addPhysicalGroup(3, [volumes[2][1]])

    # nodeTags2, nodeCoords2 = gmsh.model.mesh.getNodesForPhysicalGroup(3, pv2)
    # nodeTags3, nodeCoords3 = gmsh.model.mesh.getNodesForPhysicalGroup(3, pv3)

    # nodeCoords2 = np.array(nodeCoords2).reshape(-1, 3)
    # nodeCoords3 = np.array(nodeCoords3).reshape(-1, 3)

    def get_nodes_at_x(x: float, nc: np.ndarray, nt: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Return the tags and indices of the nodes at x-coordinate x.
        """
        # return np.where(np.isclose(nodeCoords[:, 0], x))
        indices = np.where(np.isclose(nc[:, 0], x))[0]
        tags = nt[indices]
        return tags, indices

    # Do the procedure which is commented out below,
    # but for each step in the revolution.
    # Make a for loop and create a local coordinate system for each step.

    # input(f"{n_angle=})")

    # for n in range(n_angle + 1):  # type: ignore
    if False:
        angle = 2 * np.pi * mesh.wedge_angle / 360 * (n / n_angle)  # type: ignore
        c, s = np.cos(angle), np.sin(angle)
        input(f"angle: {np.degrees(angle):.2f}")
        R = np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])

        nc2 = np.dot(nodeCoords2, R.T)
        nc3 = np.dot(nodeCoords3, R.T)

        nt2 = nodeTags2[np.where(np.isclose(nc2[:, 2], 0))]
        nc2 = nc2[np.where(np.isclose(nc2[:, 2], 0))]
        nt3 = nodeTags3[np.where(np.isclose(nc3[:, 2], 0))]
        nc3 = nc3[np.where(np.isclose(nc3[:, 2], 0))]

        tags, indices = get_nodes_at_x(x7, nc3, nt3)
        y_correct = nc3[indices, 1]
        y_correct = y_correct[np.argsort(y_correct)]

        filtered_y = nc3[np.where(np.isclose(nc3[:, 1], y6))]
        x_coords = filtered_y[:, 0]
        x_coords = x_coords[np.argsort(x_coords)]

        for x in x_coords[1:-2]:
            tags, indices = get_nodes_at_x(float(x), nc3, nt3)
            tags = tags[np.argsort(nc3[indices, 1])]
            indices = indices[np.argsort(nc3[indices, 1])]
            y = nc3[indices, 1]

            for i in range(1, len(y) - 1):
                current_coords = gmsh.model.mesh.getNode(tags[i])[0]
                xc, yc, zc = current_coords
                gmsh.model.mesh.setNode(tags[i], (xc, y_correct[i], zc), (0, 0, 0))

            # Now for s2, which is tricier. Instead of having a constant reference y=y6,
            # we need to use the line from (x6,y6) to (x10,y10) as the base line.

            a = (y6 - y10) / (x6 - x10)
            b = y6 - a * x6

            x_coords = nc2[:, 0]
            y_coords = nc2[:, 1]
            # x_coords = x_coords[np.argsort(y_coords)]
            # y_coords = y_coords[np.argsort(y_coords)]

            lower_points = nc2[np.where(np.isclose(y_coords, a * x_coords + b))]
            lower_tags = nt2[np.where(np.isclose(y_coords, a * x_coords + b))]

            # The upper coords are defined by the tank curve,
            # so we can use tank.mesh.get_radius(y) to get the upper coords.

            radius = [tank.get_radius(y) for y in y_coords]
            upper_points = nc2[np.where(np.isclose(x_coords, radius))]
            upper_tags = nt2[np.where(np.isclose(y_coords, radius))]

            import matplotlib.pyplot as plt

            plt.plot(nc2[:, 0], nc2[:, 1], "ko")
            # plt.plot(upper_points[:, 0], upper_points[:, 1], "rx")
            # plt.plot(upper_points[:, 0], upper_points[:, 1], "bx")
            # plt.plot(lower_points[:, 0], lower_points[:, 1], "o")
            # plt.show()

            # Sort the points from high to low y-values
            sort_indices = np.flip(np.argsort(upper_points[:, 1]))
            print(sort_indices)
            # upper_points = upper_points[sort_indices, :]
            # lower_points = lower_points[sort_indices, :]
            upper_points = upper_points[sort_indices]
            lower_points = lower_points[sort_indices]

            #         edited_tags = []
            #         tags = nt2.copy()
            for i in range(1, len(upper_points) - 1):
                #             # Lets go!
                #             # First, mask away points with tags in edite_tags
                #             # current_tags = nt2[np.where(~np.isin(nt2, edited_tags))]
                #             # Now we want to find the correct points in between the lower and upper line!
                #             # First, only consider points with j > lower[i, 0].

                #             current_tags = tags[np.where(nc2[:, 0] > lower_points[i, 0])]

                # Define the line from upper to lower, i.e. from upper[i] to lower[i].
                # Do not consider points with y-values lower than this line. Lets
                # make this line in the form y = A*x + B.

                # A = (upper_points[i, 1] - lower_points[i, 1])/(upper_points[i, 0] - lower_points[i, 0])
                # B = upper_points[i, 1] - A*upper_points[i, 0]

                # mask = nc2[:, 1] > A*nc2[:, 0] + B

                # Use the is_above_line function instead
                mask_above_current_line = np.array([
                    is_above_line(upper_points[i], lower_points[i], point) for point in nc2
                ])

                mask_above_previous_line = np.array([
                    is_above_line(upper_points[i - 1], lower_points[i - 1], point) for point in nc2
                ])

                mask_not_upper_or_lower = ~np.isin(nt2, list(lower_tags) + list(upper_tags))

                # Mask should be above current and NOT above previous
                mask = mask_above_current_line & ~mask_above_previous_line & mask_not_upper_or_lower

                plt.plot([lower_points[i, 0], upper_points[i, 0]], [lower_points[i, 1], upper_points[i, 1]], "k-")
                plt.plot(upper_points[:, 0], upper_points[:, 1], "rx")
                plt.plot(lower_points[:, 0], lower_points[:, 1], "bx")
                plt.plot(nc2[mask, 0], nc2[mask, 1], "gx")
                plt.show()

    #             # We know that this line should have i+1 points:
    #             # Lets only keep the i+1 one poinst with lowest y-values.

    #             current_tags = current_tags[np.argsort(nc2[np.isin(nt2, current_tags), 1])][:i+1]

    #             for t in current_tags:
    #                 if t not in edited_tags+list(lower_tags)+list(upper_tags):
    #                     current_coords = nc2[np.where(nt2 == t)][0]
    #                     xc, yc, zc = current_coords
    #                     gmsh.model.mesh.setNode(t, (xc, a*xc + b, zc), (0, 0, 0))
    #                     edited_tags.append(t)

    # print(len(current_tags))

    # print(len(mask))

    # import matplotlib.pyplot as plt
    # plt.plot(nc2[:, 0], nc2[:, 1], "o")
    # plt.plot(lower_points[:, 0], lower_points[:, 1], "o")
    # plt.plot(upper_points[:, 0], upper_points[:, 1], "o")
    # plt.show()

    # Now the

    # volumes = gmsh.model.getEntities(dim=3)
    # gas = gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes[:-1]])
    # metal = gmsh.model.addPhysicalGroup(3, [volumes[-1][1]])
    # gmsh.model.setPhysicalName(3, gas, "gas")
    # gmsh.model.setPhysicalName(3, metal, "metal")

    # surfaces: list[tuple[int, int]] = gmsh.model.getEntities(dim=2)

    # def add_physical_surface(ind: list[int], name: str) -> None:
    #     # = [surfaces[i[1]] for i in indices]
    #     m = gmsh.model.addPhysicalGroup(2, [surfaces[i][1] for i in ind])
    #     gmsh.model.setPhysicalName(2, m, name)

    # # # Recombine algorithm for quads/hexes
    # # # gmsh.option.setNumber("Mesh.RecombineAll", 1)

    # # # Generate the 3D mesh
    # gmsh.model.mesh.generate(3)
    # gmsh.model.geo.synchronize()
    # # gmsh.model.mesh.recombine()
    # gmsh.model.mesh.optimize()

    # for i in range(len(surfaces)):
    #     # Format an int string with leading zeros
    #     ind = i
    #     int_string = f"s_{ind:02d}"
    #     print(f"Adding physical surface {int_string}")
    #     add_physical_surface([ind], int_string)

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
    gmsh.option.setNumber("Mesh.TransfiniteTri", 1)
    # gmsh.option.setNumber("Mesh.Smoothing", 0)


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


def get_corner_coords(mesh: "KSiteMesh.KSiteMesh") -> tuple[float, float]:
    """
    Get the corner coords of the boundary layer.
    """

    A, B = mesh.tank.cylinder_radius, mesh.tank.cap_height
    C = mesh.tank.cylinder_height
    A -= mesh.t_BL
    B -= mesh.t_BL

    def r_ellipse(y: float) -> float:
        if y > C / 2:
            return float(A * np.sqrt(1 - (y - C / 2) ** 2 / B**2))
        elif y > -C / 2:
            return A
        else:
            return float(A * np.sqrt(1 - (y + C / 2) ** 2 / B**2))

    y = mesh.tank.y_interface + mesh.t_BL
    return r_ellipse(y), y


def is_above_line(start: np.ndarray, end: np.ndarray, point: np.ndarray) -> bool:
    """
    Check if a point is above a line defined by start -> end.
    """

    A = (end[1] - start[1]) / (end[0] - start[0])
    B = start[1] - A * start[0]

    return point[1] > A * point[0] + B


def print_debug(mesh: "KSiteMesh.KSiteMesh", msg: str) -> None:
    if mesh.debug:
        print(msg)
