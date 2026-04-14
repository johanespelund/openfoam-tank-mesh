"""
Example for CylinderCapsMesh.

"""

from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh

# Common mesh parameters shared by both tests below.
_MESH_PARAMS = {
    "cylinder_radius": 1.0,
    "cylinder_height": 0.5,
    "cap_height": 0.5,
    "fill_level": 0.95,
    "wall_cell_size": 2.0e-3,
    "wall_tan_cell_size": 5.0e-3,
    "bulk_cell_size": 25e-3,
    "r_BL": 1.10,
    "tri_bulk": False,
    "outlet_radius": 0.050,
    "internal_outlet": 0,
    "debug": True,
    "revolve": 0,
    "n_revolve": 0,
    "n_wall_layers": 6,
    "wall_thickness": 1e-2,
    "VoF": False,
    "smoothing": True,
    "symmetry": True,
}


if __name__ == "__main__":
    mesh = CylinderCapsMesh(input_parameters=dict(_MESH_PARAMS))
    mesh.generate()
