from __future__ import annotations

import pathlib

from openfoam_tank_mesh.gmsh_scripts.two_phase import run as run_gmsh
from openfoam_tank_mesh.gmsh_scripts.stl import generate_3D_internal_outlet_stl, generate_3D_stl
from openfoam_tank_mesh.KSiteTank import KSiteTank
from openfoam_tank_mesh.TankMesh import TankMesh
from openfoam_tank_mesh.Profile import KSiteProfile


class KSiteMesh(TankMesh):
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
        self.run_command(f"transformPoints -rotate-y -{self.wedge_angle / 2}")
        # if self.revolve:
        #     self.run_openfoam_utility(
        #         "createPatch -overwrite",
        #         "createPatchDict.gmshRevolve",
        #     )
        self.run_openfoam_utility(
            "topoSet",
            "topoSetDict.gmsh",
        )
        self.run_openfoam_utility(
            "createPatch -overwrite",
            "createPatchDict.gmsh"
        )
        input("Press Enter to continue...")  # noqa: T201

    def generate(self) -> None:
        """
        Generate the mesh.
        """
        self.gmsh()

        if self.internal_outlet:
            self.create_internal_outlet()
        self.run_openfoam_utility("topoSet", "topoSetDict.createFinalFaceSets")
        self.run_command("rm -rf 0/cellToRegion")
        self.run_command("splitMeshRegions -cellZonesOnly -overwrite")
        self.run_command("rm -rf constant/polyMesh")
        # self.sed("metal_outlet", "outlet", "constant/metal/polyMesh/boundary")
        # self.generate_flange_boundary()
        self.check_mesh(regions=["gas", "liquid"])

        return None

    def generate_flange_boundary(self) -> None:
        """
        On the metal region, create a boundary for the flange.
        """
        self.run_openfoam_utility("topoSet -region metal", "topoSetDict.createFlange")
        self.run_openfoam_utility("createPatch -overwrite -region metal", "createPatchDict.createFlange")

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

        return f"{pathlib.Path.cwd()}/parameters.KSiteMesh"

    def q_insulation(self) -> float:
        """
        Return the heat flux for the MLI/insulation for T_amb = 350 K.
        From Table 2 in:
        Stochl, R. J.; Knoll, R. H. Thermal Performance of a Liquid
        Hydrogen Tank Multilayer Insulation System at Warm Boundary
        Temperatures of 630, 530, and 152 R; 1991.
        """
        return 2.97

    def Q_outlet_pipe(self) -> float:
        """
        Return the heat loss from the outlet pipe.
        Total heat loss from plumbing and ducts are 3.194 W,
        assume here that 50% of the heat loss is from the outlet pipe.
        """
        return 3.194 / 2

    def Q_parasitic(self) -> float:
        """
        Return parasitic heat loss for the tank,
        mainly support and plumbing. Exclude the heat loss from the outlet pipe.
        """
        return 6.89 - self.Q_outlet_pipe()

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

    def create_tank_profile(self) -> TP.TankProfile:
        A = self.cap_height
        B = self.cylinder_radius
        C = self.cylinder_height

        ellipse1 = EllipseArc(
            name="ellipse1",
            y_start=0,
            y_end=A,
            axis_major=B,
            axis_minor=A,
            length_scale=mesh.wall_tan_cell_size,
        )

        line1 = LineSegment(
            name="line1",
            y_start=A,
            y_end=A + C,
            r_start=B,
            r_end=B,
            length_scale=mesh.wall_tan_cell_size,
        )

        ellipse2 = EllipseArc(
            name="ellipse2",
            y_start=A,
            y_end=2 * A,
            axis_major=B,
            axis_minor=A,
            y_offset=C,
            length_scale=mesh.wall_tan_cell_size,
        )

        tank_profile = TankProfile(
            segments=[ellipse1, line1, ellipse2],
            fill_level=0.49,  # mesh.tank.fill_level,
            outlet_radius=mesh.tank.outlet_radius,
            internal_outlet=mesh.internal_outlet,
        )

        tank_profile.insert_interface(tol=0.02, x_wall=mesh.wall_cell_size)

        return tank_profile
