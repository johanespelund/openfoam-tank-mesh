"""Public package API for :mod:`openfoam_tank_mesh`."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

__all__ = [
    "CylinderCapsMesh",
    "CylinderCapsTankProfile",
    "KSiteMesh",
    "KSiteProfile",
    "SphereMesh",
    "SphereProfile",
    "TankMesh",
    "TankProfile",
    "TwoPhaseTankMesh",
]

_EXPORT_MAP = {
    "CylinderCapsMesh": "openfoam_tank_mesh.TwoPhaseMesh",
    "KSiteMesh": "openfoam_tank_mesh.TwoPhaseMesh",
    "SphereMesh": "openfoam_tank_mesh.TwoPhaseMesh",
    "TankProfile": "openfoam_tank_mesh.Profile",
    "CylinderCapsTankProfile": "openfoam_tank_mesh.Profile",
    "KSiteProfile": "openfoam_tank_mesh.Profile",
    "SphereProfile": "openfoam_tank_mesh.Profile",
    "TankMesh": "openfoam_tank_mesh.TankMesh",
    "TwoPhaseTankMesh": "openfoam_tank_mesh.TwoPhaseTankMesh",
}


def __getattr__(name: str) -> object:
    if name not in _EXPORT_MAP:
        msg = f"module {__name__!r} has no attribute {name!r}"
        raise AttributeError(msg)
    module = import_module(_EXPORT_MAP[name])
    return getattr(module, name)


if TYPE_CHECKING:
    from openfoam_tank_mesh.Profile import CylinderCapsTankProfile, KSiteProfile, SphereProfile, TankProfile
    from openfoam_tank_mesh.TankMesh import TankMesh
    from openfoam_tank_mesh.TwoPhaseMesh import CylinderCapsMesh, KSiteMesh, SphereMesh
    from openfoam_tank_mesh.TwoPhaseTankMesh import TwoPhaseTankMesh
