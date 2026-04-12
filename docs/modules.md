# API reference

## Canonical mesh builders

::: openfoam_tank_mesh.mesh_builders

## Geometry profiles

::: openfoam_tank_mesh.Profile

## Mesh pipeline base classes

::: openfoam_tank_mesh.mesh_pipeline

## Runtime requirements and limitations

- OpenFOAM command-line tools must be available in `PATH` for mesh generation commands.
- `gmsh` Python bindings are required for Gmsh-based mesh generation.
- Some integration tests and workflows that call OpenFOAM utilities require a fully sourced OpenFOAM environment.
