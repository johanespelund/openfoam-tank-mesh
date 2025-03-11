from __future__ import annotations

import pathlib

from openfoam_tank_mesh.gmsh_scripts.ksite83 import run as run_gmsh
from openfoam_tank_mesh.gmsh_scripts.stl import generate_3D_internal_outlet_stl, generate_3D_stl
from openfoam_tank_mesh.NASA1m3Tank import NASA1m3Tank
from openfoam_tank_mesh.TankMesh import TankMesh


class NASA1m3Mesh(TankMesh):
    def __init__(self, input_parameters: dict) -> None:
        self.tank: NASA1m3Tank = NASA1m3Tank(
            fill_level=input_parameters["fill_level"],
            outlet_radius=input_parameters["outlet_radius"],
            insulation_type=input_parameters["insulation_type"],
            cargo=input_parameters["cargo"],
        )
        super().__init__(tank=self.tank, input_parameters=input_parameters)

        return None

    def gmsh(self) -> None:
        run_gmsh(self)
        self.run_command("gmshToFoam KSite49.msh")
        self.run_command(f"transformPoints -rotate-y -{self.wedge_angle / 2}")
        if self.revolve:
            self.run_openfoam_utility(
                "createPatch -overwrite",
                "createPatchDict.gmshRevolve",
            )
            self.patch_name_pos = "cyclic_pos"
            self.patch_name_neg = "cyclic_neg"
        else:
            self.run_openfoam_utility(
                "createPatch -overwrite",
                "createPatchDict.gmsh_wedge",
            )
            self.patch_name_pos = "wedge_pos"
            self.patch_name_neg = "wedge_neg"
        self.write_mesh_parameters()

    def generate(self) -> None:
        """
        Generate the mesh.
        """
        self.gmsh()

        self.run_openfoam_utility("topoSet", "topoSetDict.createFinalFaceSets")
        self.run_command("splitMeshRegions -cellZonesOnly -overwrite")

        self.remove_wall()
        if self.internal_outlet:
            self.create_internal_outlet()

        # self.check_mesh(regions=["gas", "metal"])

        return None

    def generate_stl(self) -> None:
        """
        Generate a stl file with named surfaces for use in cfMesh.
        """
        if self.internal_outlet:
            generate_3D_internal_outlet_stl(self)
        else:
            generate_3D_stl(self)

    @property
    def dict_path(self) -> str:
        """
        The path to the OpenFOAM dict folder.
        """

        return f"{pathlib.Path(__file__).parent}/dicts/cylindrical_tanks/"

    @property
    def parameters_path(self) -> str:
        """
        The path to the mesh parameters.
        """

        return f"{pathlib.Path.cwd()}/parameters.NASA1m3Mesh"

    def q_insulation(self) -> float:
        """
        NASA 1m3 demo tank cases from [1]
        [1] doi:10.1063/1.2908497
        """
        return get_NASA1m3_heat_flux(self.tank.insulation_type, self.tank.cargo)

    def Q_parasitic(self) -> float:
        """
        Return parasitic heat loss for the tank,
        mainly support and plumbing.
        For NASA1m3, we only use average heat flux, so it is set to 0.
        """
        return 0

    def Q_outlet_pipe(self) -> float:
        """
        Return heat loss from the outlet pipe.
        For NASA1m3, we only use average heat flux, so it is set to 0.
        """
        return 0

    def Q_liquid(self) -> float:
        """
        Liquid surface heat transfer, sum of insulation and parasitic,
        i.e. assume all parasitic heat is transferred to the liquid.
        """
        return float(self.tank.area_liquid * self.q_insulation() + self.Q_parasitic())

    def Q_gas(self) -> float:
        """
        Gas surface heat transfer, sum of insulation and parasitic,
        i.e. assume all parasitic heat is transferred to the gas.
        """
        return float(self.tank.area_gas * self.q_insulation())


def get_NASA1m3_heat_flux(insulation_type: str = "bubbles", cargo: str = "LH2") -> float:
    """
    Returns the heat flux (W/m²) for the NASA 1m³ tank cases based on insulation type and cargo.

    Heat flux values are derived from [1], estimated using the script
    'cases/nasa_1m3_demo_tank.py' at commit #cce665501e3f303ebc0aa968282955f00962714f.

    Parameters:
        insulation_type (str): Type of insulation, either "bubbles" or "perlite".
        cargo (str): Type of cargo, either "LH2" or "LN2".

    Returns:
        float: Heat flux in W/m².
    """

    heat_flux_values = {
        ("bubbles", "LH2"): 2.34,
        ("bubbles", "LN2"): 2.24,
        ("perlite", "LH2"): 3.39,
        ("perlite", "LN2"): 3.2,
    }

    return heat_flux_values[(insulation_type, cargo)]
