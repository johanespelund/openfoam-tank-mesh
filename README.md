# openfoam-tank-mesh

This is a python package for creating meshes to use in CFD simulations with OpenFOAM.

## Example meshes

| <img width="1744" height="1144" alt="image" src="https://github.com/user-attachments/assets/9e2a3135-f6cb-4471-80b8-38bce7bfc05f" /> |
| :----------------------------------------------------------------------------------------------------------------------------------: |
|                                                         _3D cylinder mesh._                                                          |

## Installation

Before installing, [OpenFOAM 12](https://openfoam.org/version/12/) should be installed and sourced.
Install using `pip` directly

```sh
python -m venv .venv # If creating new venv
source .venv/bin/activate
pip install git+https://github.com/johanespelund/openfoam-tank-mesh.git
```

or if you are using `uv`

```sh
uv pip install git+https://github.com/johanespelund/openfoam-tank-mesh.git
```

## Usage

An axisymmetrical, vertical cylinder tank can be generated using this script (includes walls):

```python
from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh

input_parameters = dict(
    cylinder_radius=1, # Radius of vertical cylinder
    cylinder_height=1, # Height of cylinder (0 gives a spherical/ellipsoid)
    cap_height=0.5,      # Height of ellipse caps (=cylinder radius gives sphere)
    fill_level=0.35,
    wall_cell_size=15.0e-3,
    wall_tan_cell_size=15.0e-3,
    bulk_cell_size=50e-3,
    r_BL=1.05,
    tri_bulk=False,
    outlet_radius=0.0127,
    internal_outlet=0.0127 * 4,
    debug=True,
    revolve=0,
    n_revolve=0,
    n_wall_layers=6,
    VoF=False,
)

mesh = CylinderCapsMesh(input_parameters)

mesh.generate()

```

A 3D horizontal cylinder tank can be generated using this script (flat ends, no walls for now).
Yields mesh shown in figure above.

```python
from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh

input_parameters = dict(
    cylinder_radius=1, # Radius of vertical cylinder
    cylinder_height=0, # Height of cylinder (0 gives a spherical/ellipsoid)
    cap_height=1,      # Height of ellipse caps (=cylinder radius gives sphere)
    fill_level=0.25,
    wall_cell_size=15.0e-3,
    wall_tan_cell_size=15.0e-3,
    bulk_cell_size=50e-3,
    r_BL=1.05,
    tri_bulk=False,
    outlet_radius=0.0127,
    internal_outlet=0.0127 * 4,
    debug=True,
    revolve=0,
    n_revolve=0,
    n_wall_layers=6,
    VoF=False,
    empty_2d=True,
    extrude_cylinder=2.5,
    mirror=True,
)

mesh = CylinderCapsMesh(input_parameters)

mesh.generate()

```

## TODO

- [x] Add option for wall thickness (hardcoded for K-Site tank now!)
- [x] Make more robust (breaks for certain cell size and fill level combos).
- [ ] Add support for walls with extrude_cylinder.
- [x] Add Laplacian smoother to reduce non-orthogonality.
- [ ] Add support for ESI version of OpenFOAM.
