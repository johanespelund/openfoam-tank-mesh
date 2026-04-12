"""Public package API for :mod:`openfoam_tank_mesh`."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

__all__ = [
    "CylinderCapsMesh",
    "CylinderCapsTankProfile",
    "GmshMeshPipeline",
    "KSiteMesh",
    "KSiteProfile",
    "OpenFoamMeshPipeline",
    "SphereMesh",
    "SphereProfile",
    "TankProfile",
    "TwoPhaseTankMesh",
]

_EXPORT_MAP = {
    "GmshMeshPipeline": "openfoam_tank_mesh.TwoPhaseMesh",
    "CylinderCapsMesh": "openfoam_tank_mesh.mesh_builders",
    "KSiteMesh": "openfoam_tank_mesh.mesh_builders",
    "SphereMesh": "openfoam_tank_mesh.mesh_builders",
    "TankProfile": "openfoam_tank_mesh.Profile",
    "CylinderCapsTankProfile": "openfoam_tank_mesh.Profile",
    "KSiteProfile": "openfoam_tank_mesh.Profile",
    "SphereProfile": "openfoam_tank_mesh.Profile",
    "OpenFoamMeshPipeline": "openfoam_tank_mesh.mesh_pipeline",
    "TwoPhaseTankMesh": "openfoam_tank_mesh.mesh_pipeline",
}


def __getattr__(name: str) -> object:
    if name not in _EXPORT_MAP:
        msg = f"module '{__name__}' has no attribute '{name}'"
        raise AttributeError(msg)
    module = import_module(_EXPORT_MAP[name])
    return getattr(module, name)


if TYPE_CHECKING:
    from openfoam_tank_mesh.mesh_builders import CylinderCapsMesh, KSiteMesh, SphereMesh
    from openfoam_tank_mesh.mesh_pipeline import OpenFoamMeshPipeline, TwoPhaseTankMesh
    from openfoam_tank_mesh.Profile import CylinderCapsTankProfile, KSiteProfile, SphereProfile, TankProfile
    from openfoam_tank_mesh.TwoPhaseMesh import GmshMeshPipeline
