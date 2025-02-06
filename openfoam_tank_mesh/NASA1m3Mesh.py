from __future__ import annotations

import pathlib

from openfoam_tank_mesh.gmsh_scripts.ksite83 import run as run_gmsh
from openfoam_tank_mesh.gmsh_scripts.stl import generate_stl
from openfoam_tank_mesh.NASA1m3Tank import NASA1m3Tank
from openfoam_tank_mesh.TankMesh import TankMesh


class NASA1m3Mesh(TankMesh):
    def __init__(self, input_parameters: dict) -> None:
        self.tank: NASA1m3Tank = NASA1m3Tank(
            fill_level=input_parameters["fill_level"],
            outlet_radius=input_parameters["outlet_radius"],
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
        else:
            self.run_openfoam_utility(
                "createPatch -overwrite",
                "createPatchDict.gmsh_wedge",
            )

    def generate(self) -> None:
        """
        Generate the mesh.
        """
        self.gmsh()

        self.run_openfoam_utility("topoSet", "topoSetDict.createFinalFaceSets")
        self.run_command("splitMeshRegions -cellZonesOnly -overwrite")
        # TODO: Just implement the base mesh without the wall again,
        #       and use add_wall() to add the wall.
        self.remove_wall()

        # self.check_mesh(regions=["gas", "metal"])

        return None

    def generate_stl(self) -> None:
        """
        Generate a stl file with named surfaces for use in cfMesh.
        """
        generate_stl(self)

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

    def q_MLI(self) -> float:
        """
        Return the heat flux for the MLI/insulation for T_amb = 350 K.
        From Table 2 in:
        Stochl, R. J.; Knoll, R. H. Thermal Performance of a Liquid
        Hydrogen Tank Multilayer Insulation System at Warm Boundary
        Temperatures of 630, 530, and 152 R; 1991.
        """
        return 2.97

    def Q_parasitic(self) -> float:
        """
        Return parasitic heat loss for the tank,
        mainly support and plumbing.
        """
        return 6.89

    def Q_liquid(self) -> float:
        """
        Liquid surface heat transfer, sum of insulation and parasitic,
        i.e. assume all parasitic heat is transferred to the liquid.
        """
        return float(self.tank.area_liquid * self.q_MLI() + self.Q_parasitic())

    def Q_gas(self) -> float:
        """
        Gas surface heat transfer, sum of insulation and parasitic,
        i.e. assume all parasitic heat is transferred to the gas.
        """
        return float(self.tank.area_gas * self.q_MLI())
