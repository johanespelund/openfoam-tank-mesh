# API reference

## High-level meshes

::: openfoam_tank_mesh.TwoPhaseMesh

## Geometry profiles

::: openfoam_tank_mesh.Profile

## Mesh base classes

::: openfoam_tank_mesh.TankMesh

::: openfoam_tank_mesh.TwoPhaseTankMesh

## Geometry base class

::: openfoam_tank_mesh.Tank

## Runtime requirements and limitations

- OpenFOAM command-line tools must be available in `PATH` for mesh generation commands.
- `gmsh` Python bindings are required for Gmsh-based mesh generation.
- Some integration tests and workflows that call OpenFOAM utilities require a fully sourced OpenFOAM environment.
