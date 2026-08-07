"""
Example mesh script for the MHTB NASA tank.
"""

from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh

INCH = 0.0254

CYLINDER_DIAMETER = 3.05
CAP_HEIGHT = CYLINDER_DIAMETER / 4

_MESH_PARAMS = {
    "cylinder_diameter": CYLINDER_DIAMETER,
    "cylinder_height": 3.05,
    "cap_height": CAP_HEIGHT,
    "fill_level": 0.90,
    "wall_cell_size": 6.0e-3,
    "wall_tan_cell_size": 25.0e-3,
    "bulk_cell_size": 75.0e-3,
    "r_BL": 1.10,
    "tri_bulk": False,
    "outlet_radius": 0.5 * 24 * INCH,
    "internal_outlet": 0.0,
    "debug": False,
    "revolve": 0,
    "n_revolve": 0,
    "n_wall_layers": 6,
    "wall_thickness": 0.5 * INCH,
    "VoF": False,
    "smoothing": True,
    "symmetry": True,
}


if __name__ == "__main__":
    mesh = CylinderCapsMesh(input_parameters=dict(_MESH_PARAMS))
    mesh.generate()
