from __future__ import annotations

import pathlib

from openfoam_tank_mesh.gmsh_scripts.two_phase import run as run_gmsh
from openfoam_tank_mesh.TwoPhaseTankMesh import TwoPhaseTankMesh
from openfoam_tank_mesh.Profile import KSiteProfile, SphereProfile

from numpy import isclose


class KSiteMesh(TwoPhaseTankMesh):
    def __init__(self, input_parameters: dict) -> None:

        self.modify_outlet = 1
        if input_parameters["outlet_radius"] <= input_parameters["wall_tan_cell_size"]:
            input_parameters["outlet_radius"] *= 2
            self.modify_outlet = True

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
        self.n_wall_layers = input_parameters["n_wall_layers"]

        return None

    def gmsh(self) -> None:
        """
        Generate the mesh using Gmsh.
        """

        run_gmsh(self)
        if self.modify_outlet:
            self.outlet_radius /= 2
            self.write_mesh_parameters()

        self.run_command("gmshToFoam mesh.msh")
        angle = max(self.wedge_angle, self.revolve)/2
        # self.run_command(f"transformPoints -rotate-y -{angle}")  # .org
        self.run_command(f"transformPoints \"Ry={-self.wedge_angle / 2}\"") #.com
        self.run_openfoam_utility(
            "topoSet",
            "topoSetDict.gmsh",
        )

        if self.revolve and self.symmetry:
            self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh_symmetry")
        elif self.revolve and self.cyclic:
            self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh_cyclic")
        else:
            self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh")

    def generate(self) -> None:
        """
        Generate the mesh.
        """
        self.gmsh()

        if self.internal_outlet:
            self.create_internal_outlet()
        else:
            self.remove_wall_outlet()
        self.run_command("rm -rf 0/cellToRegion")
        if self.lid:
            self.regions.append("lid")
            self.write_mesh_parameters()

            y = self.tank.y_lid
            r = self.tank.r_lid
            nx, ny = self.tank.get_normal(y)/4
            r1, y1 = r - nx, y - ny
            r2, y2 = r + nx, y + ny

            topodict = self.dict("topoSetDict.splitMetalRegions")
            self.sed("radius1 .*;", f"radius1 {r1:.4f};", topodict)
            self.sed("radius2 .*;", f"radius2 {r2:.4f};", topodict)
            self.sed("point1 .*;", f"point1 (0 {y1:.4f} 0);", topodict)
            self.sed("point2 .*;", f"point2 (0 {y2:.4f} 0);", topodict)

            self.run_openfoam_utility("topoSet", "topoSetDict.splitMetalRegions")

        # self.obstacle = True
        # if self.obstacle:
        #     y = self.tank.y_lid
        #     r = self.tank.r_lid
        #     w = 0.02
        #     h = 0.01
        #     n = self.tank.get_normal(y)
        #     tr, ty = n[1], -n[0]
        #     topodict = self.dict("topoSetDict.obstacle")
        #     self.sed("n1 .*;", "n1 (-1 0 0);", topodict)
        #     self.sed("n2 .*;", f"n2 ({tr} {ty} 0);", topodict)
        #     self.sed("centre .*;", f"centre ({r} {y} 0);", topodict)
        #     self.sed("box .*;", f"box ({r-w} {y-h} -1e6) ({r} {y} 1e6);", topodict)
        #     self.run_openfoam_utility("topoSet", "topoSetDict.obstacle")

        self.run_command("splitMeshRegions -cellZonesOnly -overwrite")
        self.run_command("rm -rf constant/polyMesh")
        self.check_mesh(regions=self.regions)

        # Need to run this to convert mapped patches from .com to .org format
        # self.run_command(
        #     "find constant/ -type f -name boundary -exec "
        #     + "sed -i 's/nearestPatchFace/matching/' {} \;"
        # )

        return None


    def generate_leak_boundaries(self) -> None:
        """
        On the metal region, create a boundary for the flange.
        """
        for region in ["metal", "lid"]:
            self.run_openfoam_utility(f"topoSet -region {region}", "topoSetDict.metal_patches")
            self.run_openfoam_utility(
                f"createPatch -overwrite -region {region}", "createPatchDict.metal_patches"
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
    The following values use Table 1 from:
    DOI: 10.2514/6.2016-4674
    """

    def q_insulation(self) -> float:
        """
        Return the heat flux for the MLI/insulation for T_amb = 350 K.
        """
        return 2.952

    def Q_insulation(self) -> float:
        """
        Return the total heat loss from the insulation + support
        """
        return 41.352 + 2.813 + 3.194 + 0.879


    def Q_ducts(self) -> float:
        """
        Return  heat loss from ducting and wires.
        """
        return 0 # 3.194 + 0.879

class SphereMesh(TwoPhaseTankMesh):
    def __init__(self, input_parameters: dict) -> None:

        self.modify_outlet = 1
        if input_parameters["outlet_radius"] <= input_parameters["wall_tan_cell_size"]:
            input_parameters["outlet_radius"] *= 2
            self.modify_outlet = True

        self.tank = SphereProfile(
            radius=input_parameters["radius"],
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
        self.n_wall_layers = input_parameters["n_wall_layers"]

        return None

    def gmsh(self) -> None:
        """
        Generate the mesh using Gmsh.
        """

        run_gmsh(self)
        if self.modify_outlet:
            self.outlet_radius /= 2
            self.write_mesh_parameters()

        self.run_command("gmshToFoam mesh.msh")
        angle = max(self.wedge_angle, self.revolve)/2
        # self.run_command(f"transformPoints -rotate-y -{angle}")  # .org
        self.run_command(f"transformPoints \"Ry={-self.wedge_angle / 2}\"") #.com
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
        else:
            self.remove_wall_outlet()
        self.run_command("rm -rf 0/cellToRegion")
        self.run_command("splitMeshRegions -cellZonesOnly -overwrite")
        self.run_command("rm -rf constant/polyMesh")
        # self.generate_leak_boundaries()
        self.check_mesh(regions=["gas", "liquid", "metal"])

        # Need to run this to convert mapped patches from .com to .org format
        # self.run_command(
        #     "find constant/ -type f -name boundary -exec "
        #     + "sed -i 's/nearestPatchFace/matching/' {} \;"
        # )

        return None


    def generate_leak_boundaries(self) -> None:
        """
        On the metal region, create a boundary for the flange.
        """
        for region in ["metal", "lid"]:
            input(f"{region=}")
            self.run_openfoam_utility(f"topoSet -region {region}", "topoSetDict.metal_patches")
            self.run_openfoam_utility(
                f"createPatch -overwrite -region {region}", "createPatchDict.metal_patches"
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
    The following values use Table 1 from:
    DOI: 10.2514/6.2016-4674
    """

    def q_insulation(self) -> float:
        """
        Return the heat flux for the MLI/insulation for T_amb = 350 K.
        """
        return 2.952

    def Q_insulation(self) -> float:
        """
        Return the total heat loss from the insulation + support
        """
        return 41.352 + 2.813 + 3.194 + 0.879


    def Q_ducts(self) -> float:
        """
        Return  heat loss from ducting and wires.
        """
        return 0 # 3.194 + 0.879
# class MHTBMesh(TwoPhaseTankMesh):
#     def __init__(self, input_parameters: dict) -> None:
#         self.tank = KSiteProfile(
#             fill_level=input_parameters["fill_level"],
#             outlet_radius=input_parameters["outlet_radius"],
#             bulk_cell_size=input_parameters["bulk_cell_size"],
#             wall_tan_cell_size=input_parameters["wall_tan_cell_size"],
#             wall_cell_size=input_parameters["wall_cell_size"],
#             r_BL=input_parameters["r_BL"],
#             internal_outlet=input_parameters["internal_outlet"],
#         )
#         super().__init__(tank=self.tank, input_parameters=input_parameters)
#         self.multi_region = True

#         return None

#     def gmsh(self) -> None:
#         """
#         Generate the mesh using Gmsh.
#         """

#         run_gmsh(self)
#         self.run_command("gmshToFoam mesh.msh")
#         angle = max(self.wedge_angle, self.revolve)/2
#         self.run_command(f"transformPoints -rotate-y -{angle}")  # .org
#         # self.run_command(f"transformPoints \"Ry={-self.wedge_angle / 2}\"") #.com
#         self.run_openfoam_utility(
#             "topoSet",
#             "topoSetDict.gmsh",
#         )

#         if self.revolve:
#             self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh_symmetry")
#         else:
#             self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh")

#     def generate(self) -> None:
#         """
#         Generate the mesh.
#         """
#         self.gmsh()

#         if self.internal_outlet:
#             self.create_internal_outlet()
#         # self.run_openfoam_utility("topoSet", "topoSetDict.createFinalFaceSets")
#         self.run_command("rm -rf 0/cellToRegion")
#         # self.run_command("collapseEdges -overwrite")
#         self.run_command("splitMeshRegions -cellZonesOnly -overwrite")
#         self.run_command("rm -rf constant/polyMesh")
#         # self.sed("metal_outlet", "outlet", "constant/metal/polyMesh/boundary")
#         # self.generate_flange_boundary()
#         self.generate_leak_boundaries()
#         self.check_mesh(regions=["gas", "liquid", "metal"])
#         self.run_command(
#             "find constant/ -type f -name boundary -exec "
#             + "sed -i 's/nearestPatchFace/matching/' {} \;"
#         )

#         return None


#     def generate_leak_boundaries(self) -> None:
#         """
#         On the metal region, create a boundary for the flange.
#         """
#         self.run_openfoam_utility("topoSet -region metal", "topoSetDict.metal_patches")
#         self.run_openfoam_utility(
#             "createPatch -overwrite -region metal", "createPatchDict.metal_patches"
#         )

#     def generate_stl(self) -> None:
#         """
#         Generate a stl file with named surfaces for use in cfMesh.
#         """
#         if self.internal_outlet:
#             generate_3D_internal_outlet_stl(self)
#         else:
#             generate_3D_stl(self)

#     @property
#     def dict_path(self) -> str:
#         """
#         The path to the OpenFOAM dict folder.
#         """

#         return f"{pathlib.Path(__file__).parent}/dicts/two_phase_tanks/"

#     @property
#     def parameters_path(self) -> str:
#         """
#         The path to the mesh parameters.
#         """

#         return f"{pathlib.Path.cwd()}/parameters.KSiteMesh"


#     """
#     The following values use Table 1 from:
#     DOI: 10.2514/6.2016-4674
#     """

#     def q_insulation(self) -> float:
#         """
#         Return the heat flux for the MLI/insulation for T_amb = 350 K.
#         """
#         return 2.952

#     def Q_insulation(self) -> float:
#         """
#         Return the total heat loss from the insulation + support
#         """
#         return 41.352 + 2.813


#     def Q_ducts(self) -> float:
#         """
#         Return  heat loss from ducting and wires.
#         """
#         return 3.194 + 0.879
