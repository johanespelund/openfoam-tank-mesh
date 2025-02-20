import gmsh  # type: ignore[import-untyped]
import numpy as np

from openfoam_tank_mesh import TankMesh
from openfoam_tank_mesh.gmsh_scripts.utilities import (
    add_curve_loop,
    add_ellipse,
    add_line,
    add_physical_surface,
    add_point,
    add_surface,
    closest_odd,
    gmsh_setup,
)


def generate_3D_stl(mesh: "TankMesh.TankMesh") -> None:
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

    gmsh.model.mesh.generate(2)

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


def generate_2D_stl(mesh: "TankMesh.TankMesh") -> None:
    """
    Generate a stl file with named surfaces for use in cfMesh.
    """
    gmsh_setup()
    lc = mesh.tank.outlet_radius / 2

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
        # add_point(0, mesh.tank.y_interface + 6e-3, 0, lc),
        # add_point(mesh.tank.get_radius(mesh.tank.y_interface + 6e-3), mesh.tank.y_interface + 6e-3, 0, lc),
    ]
    gmsh.model.geo.synchronize()

    lines = [
        add_line(points[0], points[1]),
        add_ellipse(points[1], origo, major_point, points[2]),
        add_line(points[2], points[3]),
        add_line(points[3], points[4]),
        add_line(points[4], points[0]),
        # add_line(points[5], points[6]),
    ]
    gmsh.model.geo.synchronize()

    # _ = gmsh.model.geo.revolve(
    #     [(1, line) for line in lines],
    #     0,
    #     0,
    #     0,  # Point on the axis of revolution
    #     0,
    #     1,
    #     0,  # Direction of the axis of revolution
    #     angle,  # Angle of revolution
    #     numElements=[n_angle],
    #     recombine=True,
    # )

    _ = gmsh.model.geo.extrude([(1, line) for line in lines], 0, 0, 1e-3, numElements=[1])

    # cl = add_curve_loop(lines)
    # cyclic_surface_pos = add_surface(cl)

    # _, cyclic_surface_neg = gmsh.model.geo.copy([(2, cyclic_surface_pos)])[0]
    # cyclic_surface_neg = gmsh.model.geo.rotate([(2, cyclic_surface_neg)], 0, 0, 0, 0, 1, 0, angle)

    # gmsh.model.geo.rotate(gmsh.model.getEntities(dim=2), 0, 0, 0, 0, 1, 0, -angle / 2)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    s = gmsh.model.getEntities(dim=2)
    add_physical_surface([0], "outlet", s)
    add_physical_surface([1], "walls", s)
    add_physical_surface([3], "bottom", s)
    add_physical_surface([4], "axis", s)
    # add_physical_surface([5], "extra", s)
    # for i in range(len(s)):
    #     # Format an int string with leading zeros
    #     ind = i
    #     int_string = f"s_{ind:02d}"
    #     print(f"Adding physical surface {int_string}")
    #     add_physical_surface([ind], int_string, s)

    # Export to STL
    gmsh.option.setNumber("Mesh.StlOneSolidPerSurface", 2)
    gmsh.model.geo.synchronize()
    gmsh.write(f"{mesh.tank.name}.stl")

    # Run the GUI
    if mesh.debug:
        gmsh.fltk.run()
        gmsh.finalize()


def generate_2D_internal_outlet_stl(mesh: "TankMesh.TankMesh") -> None:
    """
    Generate a stl file with named surfaces for use in cfMesh.
    """
    gmsh_setup()
    lc = mesh.tank.outlet_radius / 2

    origo = add_point(0, mesh.tank.cylinder_height / 2, 0, lc)
    major_point = add_point(mesh.tank.interface_radius, mesh.tank.cylinder_height / 2, 0, lc)

    y_mid = max(mesh.tank.y_interface, mesh.tank.cylinder_height / 2)
    x_mid = mesh.tank.get_radius(y_mid)

    points = [
        add_point(0, mesh.tank.y_outlet - mesh.outlet_radius, 0, lc),
        add_point(mesh.outlet_radius, mesh.tank.y_outlet - mesh.outlet_radius, 0, lc),
        add_point(mesh.outlet_radius, mesh.tank.y_outlet, 0, lc),
        add_point(x_mid, y_mid, 0, lc),
        add_point(mesh.tank.interface_radius, mesh.tank.y_interface, 0, lc),
        add_point(0, mesh.tank.y_interface, 0, lc),
        # add_point(0, mesh.tank.y_interface + 6e-3, 0, lc),
        # add_point(mesh.tank.get_radius(mesh.tank.y_interface + 6e-3), mesh.tank.y_interface + 6e-3, 0, lc),
    ]
    gmsh.model.geo.synchronize()

    lines = [
        add_line(points[0], points[1]),
        add_line(points[1], points[2]),
        add_ellipse(points[2], origo, major_point, points[3]),
        add_line(points[3], points[4]),
        add_line(points[4], points[5]),
        add_line(points[5], points[0]),
        # add_line(points[5], points[6]),
    ]
    gmsh.model.geo.synchronize()

    # _ = gmsh.model.geo.revolve(
    #     [(1, line) for line in lines],
    #     0,
    #     0,
    #     0,  # Point on the axis of revolution
    #     0,
    #     1,
    #     0,  # Direction of the axis of revolution
    #     angle,  # Angle of revolution
    #     numElements=[n_angle],
    #     recombine=True,
    # )

    _ = gmsh.model.geo.extrude([(1, line) for line in lines], 0, 0, 1e-3, numElements=[1])

    # cl = add_curve_loop(lines)
    # cyclic_surface_pos = add_surface(cl)

    # _, cyclic_surface_neg = gmsh.model.geo.copy([(2, cyclic_surface_pos)])[0]
    # cyclic_surface_neg = gmsh.model.geo.rotate([(2, cyclic_surface_neg)], 0, 0, 0, 0, 1, 0, angle)

    # gmsh.model.geo.rotate(gmsh.model.getEntities(dim=2), 0, 0, 0, 0, 1, 0, -angle / 2)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    s = gmsh.model.getEntities(dim=2)
    add_physical_surface([0], "outlet", s)
    add_physical_surface([1], "pipe", s)
    add_physical_surface([2], "walls", s)
    add_physical_surface([4], "bottom", s)
    add_physical_surface([5], "axis", s)
    # add_physical_surface([5], "extra", s)
    # for i in range(len(s)):
    #     # Format an int string with leading zeros
    #     ind = i
    #     int_string = f"s_{ind:02d}"
    #     print(f"Adding physical surface {int_string}")
    #     add_physical_surface([ind], int_string, s)

    # Export to STL
    gmsh.option.setNumber("Mesh.StlOneSolidPerSurface", 2)
    gmsh.model.geo.synchronize()
    gmsh.write(f"{mesh.tank.name}.stl")

    # Run the GUI
    if mesh.debug:
        gmsh.fltk.run()
        gmsh.finalize()


def generate_3D_internal_outlet_stl(mesh: "TankMesh.TankMesh") -> None:
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
        add_point(0, mesh.tank.y_outlet - 2 * mesh.outlet_radius, 0, lc),
        add_point(mesh.outlet_radius, mesh.tank.y_outlet - 2 * mesh.outlet_radius, 0, lc),
        add_point(mesh.outlet_radius, mesh.tank.y_outlet, 0, lc),
        add_point(x_mid, y_mid, 0, lc),
        add_point(mesh.tank.interface_radius, mesh.tank.y_interface, 0, lc),
        add_point(0, mesh.tank.y_interface, 0, lc),
        # add_point(0, mesh.tank.y_interface + 6e-3, 0, lc),
        # add_point(mesh.tank.get_radius(mesh.tank.y_interface + 6e-3), mesh.tank.y_interface + 6e-3, 0, lc),
    ]
    gmsh.model.geo.synchronize()

    lines = [
        add_line(points[0], points[1]),
        add_line(points[1], points[2]),
        add_ellipse(points[2], origo, major_point, points[3]),
        add_line(points[3], points[4]),
        add_line(points[4], points[5]),
        add_line(points[5], points[0]),
        # add_line(points[5], points[6]),
    ]
    gmsh.model.geo.synchronize()

    # _ = gmsh.model.geo.revolve(
    #     [(1, line) for line in lines],
    #     0,
    #     0,
    #     0,  # Point on the axis of revolution
    #     0,
    #     1,
    #     0,  # Direction of the axis of revolution
    #     angle,  # Angle of revolution
    #     numElements=[n_angle],
    #     recombine=True,
    # )

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

    gmsh.model.mesh.generate(2)

    s = gmsh.model.getEntities(dim=2)
    add_physical_surface([0], "outlet", s)
    add_physical_surface([1], "pipe", s)
    add_physical_surface([2, 3], "walls", s)
    add_physical_surface([4], "bottom", s)
    add_physical_surface([5], "cyclic_pos_gmsh", s)
    add_physical_surface([6], "cyclic_neg_gmsh", s)
    # add_physical_surface([5], "extra", s)

    # for i in range(len(s)):
    #     # Format an int string with leading zeros
    #     ind = i
    #     int_string = f"s_{ind:02d}"
    #     print(f"Adding physical surface {int_string}")
    #     add_physical_surface([ind], int_string, s)

    # Export to STL
    gmsh.option.setNumber("Mesh.StlOneSolidPerSurface", 2)
    gmsh.model.geo.synchronize()
    gmsh.write(f"{mesh.tank.name}.stl")

    # Run the GUI
    if mesh.debug:
        gmsh.fltk.run()
        gmsh.finalize()
