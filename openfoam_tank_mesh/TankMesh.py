from __future__ import annotations

import time
from subprocess import run
from typing import ClassVar

from rich import print as rprint
from rich.console import Console
from rich.table import Table

from openfoam_tank_mesh.exceptions import CommandFailed, MissingParameter, OpenFoamNotLoaded

# from openfoam_tank_mesh.gmsh_scripts.ksite83 import run as run_gmsh83
# from openfoam_tank_mesh.gmsh_scripts.ksite49 import run as run_gmsh49
# from openfoam_tank_mesh.gmsh_scripts.stl import generate_stl
from openfoam_tank_mesh.Tank import ABC, Tank, abstractmethod


class TankMesh(ABC):
    """
    Base class for OpenFOAM tank meshes.
    """

    REQUIRED_PARAMETERS: ClassVar[list[str]] = [
        "bulk_cell_size",
        "wall_cell_size",
        "outlet_radius",
        "fill_level",
        "debug",
    ]

    def __init__(self, tank: Tank, input_parameters: dict) -> None:
        self.tank = tank

        self.name = tank.name
        self.fill_level = tank.fill_level
        self.y_interface = tank.y_interface
        self.y_outlet = tank.y_outlet
        self.outlet_radius = tank.outlet_radius
        self.interface_radius = tank.interface_radius
        self.area_liquid = tank.area_liquid
        self.area_gas = tank.area_gas

        # Will be defined in set_parameters()
        self.bulk_cell_size: float = 0  # Bulk cell size
        self.wall_cell_size: float = 0  # Wall cell size
        self.wall_tan_cell_size: float = 0  # Wall tangential cell size
        self.debug: bool = False
        self.tri_bulk: bool = False  # Use triangle cells for bulk mesh

        # Default parameters (can be overwritten by input_parameters)
        self.r_BL: float = 1.1  # Boundary layer growth rate
        self.revolve: float = 0  # Revolution angle (0 means 2D)
        self.wedge_angle: float = 1  # Revolution angle if 2D
        self.surface_file = f"{self.name}.stl"
        self.non_coupled_cyclic = False
        self.patch_name_pos = "wedge_pos"
        self.patch_name_neg = "wedge_neg"

        self.check_openfoam_loaded(version="com")
        self.validate_parameters(input_parameters)
        self.set_parameters(input_parameters)
        self.n_BL, self.t_BL, self.e_BL = self.calculate_boundary_layer()
        self.wedge_angle = self.revolve if self.revolve else self.wedge_angle
        super().__init__()
        self.write_mesh_parameters()

    @property
    @abstractmethod
    def dict_path(self) -> str:
        """
        The path to the OpenFOAM dict folder.
        """

    @property
    @abstractmethod
    def parameters_path(self) -> str:
        """
        The path to the mesh parameters.
        """

    @abstractmethod
    def generate(self) -> None:
        """
        Use gmsh to create a mesh in .msh format,
        and then use gmshToFoam to convert it to OpenFOAM.
        """

    @abstractmethod
    def generate_stl(self) -> None:
        """
        Generate an STL file of the tank geometry for use with cfMesh.
        """

    def validate_parameters(self, input_parameters: dict) -> None:
        for param in self.REQUIRED_PARAMETERS:
            if param not in input_parameters:
                raise MissingParameter(param)

    def set_parameters(self, input_parameters: dict) -> None:
        """
        Set the mesh parameters.
        """
        for key, value in input_parameters.items():
            setattr(self, key, value)
        if "wall_tan_cell_size" not in input_parameters:
            self.wall_tan_cell_size = self.bulk_cell_size

    def write_mesh_parameters(self) -> None:
        """
        Write the mesh parameters to a file.
        """
        with open(self.parameters_path, "w") as f:
            for key, value in self.__dict__.items():
                if type(value) in [int, float, str]:
                    f.write(f"{key} {value};\n")
        self.run_command(f"cp {self.parameters_path} system/meshdata")

    def run_command(self, command: str, no_output: bool = False, return_exception: bool = False) -> None | Exception:
        """
        Method to run shell commands. The result should always be captured,
        and an error should be raised if the command fails (even if no_output is True).
        """
        t = time.time()

        command_words = command.split()
        for i, word in enumerate(command_words):
            command_words[i] = word.replace(self.dict_path, f"<{self.name}.dict_path>")
        rprint(" ".join(command_words), end="")

        capture_output = True
        result = run(command, shell=True, capture_output=capture_output)  # noqa: S602
        dt = time.time() - t
        if not no_output:
            color = "green" if result.returncode == 0 else "red"
            rprint(f"[{color}] ({dt:.6f} s)[/{color}]")
        if result.returncode != 0:
            exception = CommandFailed(command)
            if return_exception:
                return exception
            rprint(exception)
            if self.debug:
                rprint(result.stdout.decode())
            raise CommandFailed(command)

        if self.debug:
            rprint(result.stdout.decode())

        return None

    def check_openfoam_loaded(self, version: str = "com") -> None:
        """
        Check if OpenFOAM is loaded.
        version: str ("org" or "com")
        """
        command = "simpleFoam -help"
        result = run(command, shell=True, capture_output=True)  # noqa: S602
        if f"openfoam.{version}" not in result.stdout.decode():
            raise OpenFoamNotLoaded

    def run_openfoam_utility(
        self, utility: str, foam_dict: str = "", return_exception: bool = False
    ) -> None | Exception:
        """
        Run an OpenFOAM utility with optional input dictionary.
        """
        command = f"{utility} -dict {self.dict_path}/{foam_dict}"
        self.sed("include .*", f'include "{self.parameters_path}"', f"{self.dict_path}/{foam_dict}")
        return self.run_command(command, return_exception=return_exception)

    def calculate_boundary_layer(self) -> tuple[int, float, float]:
        """
        Calculate boundary layer parameters.
        returns:
            n_bl: int = number of boundary layer cells
            t_bl: float = thickness of boundary layer
            e_bl: float = boundary layer expansion
        """
        if self.r_BL == 1 or self.wall_cell_size == self.bulk_cell_size:
            self.r_BL = 1.0
            return 4, self.wall_cell_size * 4, 1
        n = 1
        x = self.wall_cell_size
        t = self.wall_cell_size
        outer_cell_size = self.wall_tan_cell_size * (1 - 0.2 * (self.tri_bulk))
        while x <= outer_cell_size / self.r_BL:
            n += 1
            x *= self.r_BL
            t += x

        return n, t, x / self.wall_cell_size

    def check_mesh(self, regions: None | list = None) -> None:
        console = Console()
        regions = [""] if not regions else [f"-region {region}" for region in regions]
        for region in regions:
            table = Table(title=f"Mesh Summary ({region})", show_header=False)
            output = run(["checkMesh", *region.split()], capture_output=True, text=True)  # noqa: S603 S607
            if "FAILED" in output.stdout:
                rprint(output.stdout)
                raise CommandFailed("checkMesh")
            for line in output.stdout.split("\n"):
                if "Mesh non-orthogonality Max" in line:
                    angle = float(line.split()[3])
                    table.add_row("Max. non-orthogonality", f"{angle:.2f}°")
                if "Max skewness" in line:
                    skewness = float(line.split()[3])
                    table.add_row("Max. skewness", f"{skewness:.2f}")
                if "  cells: " in line:
                    cells = int(line.split()[1])
                    table.add_row("Cells", f"{cells:,}")
            console.print(table)

    def remove_wall(self) -> None:
        self.run_command("rm -r constant/polyMesh")
        self.run_command("rm -r constant/metal/polyMesh")
        self.run_command("mv constant/gas/polyMesh constant/polyMesh")
        self.sed("gas_to_metal", "walls", "constant/polyMesh/boundary")
        self.sed("mappedWall", "wall", "constant/polyMesh/boundary")
        self.run_command("rm -f 0/cellToRegion")

    def add_wall(self, wall_thickness: float, n_layers: int) -> None:
        """
        Add a wall region to the mesh.
        """
        extrude_mesh_dict = "extrudeMeshDict.addWall"
        topo_set_dict = "topoSetDict.addWall"
        self.sed("thickness.*;", f"thickness {wall_thickness};", self.dict(extrude_mesh_dict))
        self.sed("nLayers.*;", f"nLayers {n_layers};", self.dict(extrude_mesh_dict))
        self.run_openfoam_utility("extrudeMesh", extrude_mesh_dict)
        self.run_openfoam_utility("topoSet", topo_set_dict)
        self.run_command("splitMeshRegions -cellZonesOnly -overwrite")

    def cfMesh(self, nLayers: int = 0) -> None:
        """
        Use cfMesh to create a 3D mesh.
        """

        self.generate_stl()
        self.run_command(f"cp {self.dict_path}/meshDict system/meshDict")
        self.sed("nLayers.*;", f"nLayers {nLayers};", "system/meshDict")
        self.run_command("cartesianMesh")
        result = self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.cfMesh", return_exception=True)
        if result:
            self.run_openfoam_utility(
                "createPatch -overwrite",
                "createPatchDict.cfMeshNonConformal",
            )
            self.non_coupled_cyclic = True
            self.patch_name_pos = "couple1"
            self.patch_name_neg = "couple2"
            self.run_command(f"cp {self.dict_path}/createNonConformalCouplesDict system/")
            self.non_coupled_info()
        else:
            self.patch_name_pos = "cyclic_pos"
            self.patch_name_neg = "cyclic_neg"

        self.write_mesh_parameters()

        self.run_command(f"transformPoints -rotate-y -{self.wedge_angle / 2}")

        # self.check_mesh()

    def sed(self, orig: str, new: str, path: str) -> None | Exception:
        """
        Run sed, e.g
        self.run_command("sed -i 's/gas_to_metal/walls/g' constant/polyMesh/boundary")
        """
        return self.run_command(f"sed -i 's|{orig}|{new}|g' {path}")

    def dict(self, name: str) -> str:
        """
        Return the path to a dictionary file.
        """
        return f"{self.dict_path}/{name}"

    def non_coupled_info(self) -> None:
        rprint("\n[bold yellow]Non-coupled cyclic boundary conditions detected[/bold yellow]")
        rprint("    [yellow]Adding createNonConformalCouplesDict to system/[/yellow]")
        rprint("    [yellow]Load OpenFOAM 11 and run the following command:[/yellow]")
        rprint("    [bold yellow]createNonConformalCouples -overwrite[/yellow bold]\n")
