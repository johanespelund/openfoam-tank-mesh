# from __future__ import annotations

import gmsh  # type: ignore[import-untyped]
import numpy as np

from openfoam_tank_mesh import TankMesh


def closest_odd(n: float) -> int:
    return max(3, int(n) // 2 * 2 + 1)


def get_N_outlet(mesh: "TankMesh.TankMesh") -> int:
    if mesh.wall_tan_cell_size >= mesh.outlet_radius:
        return 3
    else:
        temp = closest_odd(np.ceil(mesh.outlet_radius / mesh.wall_tan_cell_size) + 1)
        return temp
        # residual = abs(mesh.wall_tan_cell_size - mesh.outlet_radius / temp)
        # print(f"{residual=}")
        
        # # increase N by 2:
        # temp += 2
        # residual2 = abs(mesh.wall_tan_cell_size - mesh.outlet_radius / temp)
        # print(f"{residual2=}")

        # return temp - 2





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


def print_debug(mesh: "TankMesh.TankMesh", msg: str) -> None:
    if mesh.debug:
        print(msg)


def generate_stl2(mesh: "TankMesh.TankMesh") -> None:
    """
    Generate a stl file with named surfaces for use in cfMesh.
    """
    gmsh_setup()
    revolve = mesh.revolve
    lc = mesh.tank.outlet_radius / 2
    wedge_angle = mesh.wedge_angle
    angle = 2 * np.pi * revolve / 360 if revolve else wedge_angle * np.pi / 180
    n_angle = closest_odd(2 * np.pi * revolve / (360 * lc)) if revolve else 1

    origo = add_point(0, mesh.tank.cylinder_height / 2, 0, lc)
    major_point = add_point(mesh.tank.interface_radius, mesh.tank.cylinder_height / 2, 0, lc)

    y_mid = max(mesh.tank.y_interface, mesh.tank.cylinder_height / 2)
    x_mid = mesh.tank.get_radius(y_mid)

    points = [
        add_point(0, mesh.tank.y_outlet, 0, lc),
        add_point(mesh.tank.outlet_radius, mesh.tank.y_outlet, 0, lc),
        add_point(x_mid, y_mid, 0, lc),
        add_point(mesh.tank.interface_radius, mesh.tank.y_interface, 0, lc),
        add_point(0, mesh.tank.y_interface, 0, lc),
    ]
    gmsh.model.geo.synchronize()

    lines = [
        add_line(points[0], points[1]),
        add_ellipse(points[1], origo, major_point, points[2]),
        add_line(points[2], points[3]),
        add_line(points[3], points[4]),
        add_line(points[4], points[0]),
    ]
    gmsh.model.geo.synchronize()

    _ = gmsh.model.geo.revolve(
        [(1, line) for line in lines],
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

    cl = add_curve_loop(lines)
    cyclic_surface_pos = add_surface(cl)

    _, cyclic_surface_neg = gmsh.model.geo.copy([(2, cyclic_surface_pos)])[0]
    cyclic_surface_neg = gmsh.model.geo.rotate([(2, cyclic_surface_neg)], 0, 0, 0, 0, 1, 0, angle)

    gmsh.model.geo.rotate(gmsh.model.getEntities(dim=2), 0, 0, 0, 0, 1, 0, -angle / 2)
    gmsh.model.geo.synchronize()

    # gmsh.model.mesh.generate(2)

    s = gmsh.model.getEntities(dim=2)
    add_physical_surface([0], "outlet", s)
    add_physical_surface([1, 2], "walls", s)
    add_physical_surface([3], "bottom", s)
    add_physical_surface([4], "cyclic_pos_gmsh", s)
    add_physical_surface([5], "cyclic_neg_gmsh", s)
    # for i in range(len(s)):
    #     # Format an int string with leading zeros
    #     ind = i
    #     int_string = f"s_{ind:02d}"
    #     print(f"Adding physical surface {int_string}")
    #     add_physical_surface([ind], int_string, s)

    # Rotate the entire mesh

    # Export to STL
    gmsh.option.setNumber("Mesh.StlOneSolidPerSurface", 2)
    gmsh.model.geo.synchronize()
    gmsh.write(f"{mesh.tank.name}.stl")

    # Run the GUI
    if mesh.debug:
        gmsh.fltk.run()
        gmsh.finalize()


def add_physical_surface(ind: list[int], name: str, surfaces: list[tuple[int, int]]) -> None:
    m = gmsh.model.addPhysicalGroup(2, [surfaces[i][1] for i in ind])
    gmsh.model.setPhysicalName(2, m, name)
