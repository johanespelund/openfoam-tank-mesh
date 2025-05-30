from __future__ import annotations

import pathlib

from openfoam_tank_mesh.gmsh_scripts.two_phase import run as run_gmsh
from openfoam_tank_mesh.TwoPhaseTankMesh import TwoPhaseTankMesh
from openfoam_tank_mesh.Profile import KSiteProfile


class KSiteMesh(TwoPhaseTankMesh):
    def __init__(self, input_parameters: dict) -> None:
        self.tank = KSiteProfile(
            fill_level=input_parameters["fill_level"],
            outlet_radius=input_parameters["outlet_radius"],
            bulk_cell_size=input_parameters["bulk_cell_size"],
            wall_tan_cell_size=input_parameters["wall_tan_cell_size"],
            wall_cell_size=input_parameters["wall_cell_size"],
            r_BL=input_parameters["r_BL"],
            internal_outlet=input_parameters["internal_outlet"],
        )
        super().__init__(tank=self.tank, input_parameters=input_parameters)
        self.multi_region = True

        return None

    def gmsh(self) -> None:
        """
        Generate the mesh using Gmsh.
        """

        run_gmsh(self)
        self.run_command("gmshToFoam mesh.msh")
        angle = max(self.wedge_angle, self.revolve)/2
        self.run_command(f"transformPoints -rotate-y -{angle}")  # .org
        # self.run_command(f"transformPoints \"Ry={-self.wedge_angle / 2}\"") #.com
        self.run_openfoam_utility(
            "topoSet",
            "topoSetDict.gmsh",
        )

        if self.revolve:
            self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh_symmetry")
        else:
            self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh")

    def generate(self) -> None:
        """
        Generate the mesh.
        """
        self.gmsh()

        if self.internal_outlet:
            self.create_internal_outlet()
        # self.run_openfoam_utility("topoSet", "topoSetDict.createFinalFaceSets")
        self.run_command("rm -rf 0/cellToRegion")
        # self.run_command("collapseEdges -overwrite")
        self.run_command("splitMeshRegions -cellZonesOnly -overwrite")
        self.run_command("rm -rf constant/polyMesh")
        # self.sed("metal_outlet", "outlet", "constant/metal/polyMesh/boundary")
        # self.generate_flange_boundary()
        self.check_mesh(regions=["gas", "liquid", "metal"])
        self.generate_leak_boundaries()

        self.run_command(
            "find constant/ -type f -name boundary -exec "
            + "sed -i 's/nearestPatchFace/matching/' {} \;"
        )

        return None


    def generate_leak_boundaries(self) -> None:
        """
        On the metal region, create a boundary for the flange.
        """
        self.run_openfoam_utility("topoSet -region metal", "topoSetDict.metal_patches")
        self.run_openfoam_utility(
            "createPatch -overwrite -region metal", "createPatchDict.metal_patches"
        )

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

        return f"{pathlib.Path(__file__).parent}/dicts/two_phase_tanks/"

    @property
    def parameters_path(self) -> str:
        """
        The path to the mesh parameters.
        """

        return f"{pathlib.Path.cwd()}/parameters.KSiteMesh"


    """
    The following methods use Table 1 from:
    DOI: 10.2514/6.2016-4674
    """

    def q_insulation(self) -> float:
        """
        Return the heat flux for the MLI/insulation for T_amb = 350 K.
        """
        return 2.952

    def Q_insulation(self) -> float:
        """
        Return the total heat loss from the insulation.
        """
        return 41.352

    def Q_support(self) -> float:
        """
        Return the heat loss from support struts.
        """
        return 2.813

    def Q_ducts(self) -> float:
        """
        Return  heat loss from ducting and wires.
        """
        return 3.194 + 0.879
