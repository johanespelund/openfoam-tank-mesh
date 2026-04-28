from __future__ import annotations

import logging
import shlex
import shutil
import time
from abc import ABC, abstractmethod
from subprocess import CompletedProcess, run
from typing import ClassVar

import numpy as np
from rich.console import Console
from rich.progress import Progress, TaskID
from rich.table import Table

from openfoam_tank_mesh.exceptions import (
    CommandFailed,
    ExtrudeCylinderRequiresEmpty2D,
    MirrorRequiresEmpty2D,
    MissingParameter,
    OpenFoamNotLoaded,
    WallMeshOutletRequiresNoInternalOutlet,
)
from openfoam_tank_mesh.Profile import TankProfile

console = Console()
logger = logging.getLogger(__name__)


def _setup_file_logging(log_path: str = "mesh_generation.log") -> None:
    """Configure file-based logging for the *openfoam_tank_mesh* package.

    Log records accumulate across runs within the same process (append mode).
    This is intentional so that a session with multiple mesh objects produces a
    single, continuous log file.  Delete or rotate ``mesh_generation.log``
    between invocations if a clean log per run is preferred.
    """
    pkg_logger = logging.getLogger("openfoam_tank_mesh")
    if any(isinstance(h, logging.FileHandler) for h in pkg_logger.handlers):
        return  # already configured
    pkg_logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
            datefmt="%H:%M:%S",
        )
    )
    pkg_logger.addHandler(handler)
    pkg_logger.propagate = False


class OpenFoamMeshPipeline(ABC):
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
    _ALLOWED_COMMANDS: ClassVar[set[str]] = {
        "cartesianMesh",
        "checkMesh",
        "cp",
        "createPatch",
        "extrudeMesh",
        "gmshToFoam",
        "laplacianMeshSmoother",
        "mirrorMesh",
        "mv",
        "rm",
        "sed",
        "simpleFoam",
        "splitMeshRegions",
        "subsetMesh",
        "topoSet",
        "transformPoints",
    }

    def __init__(self, tank: TankProfile, input_parameters: dict) -> None:
        _setup_file_logging()

        # Rich Progress context set by generate() so that run_command() can
        # update the active task description instead of creating a nested Live.
        self._progress: Progress | None = None
        self._progress_task: TaskID | None = None

        self.tank: TankProfile = tank

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
        self.internal_outlet: float = 0  # Create internal outlet
        self.debug: bool = False
        self.tri_bulk: bool = False  # Use triangle cells for bulk mesh
        self.multi_region: bool = False  # Use multiple regions
        self.n_wall_layers: int = 0  # Number of wall layers
        self.wall_thickness: float = 2.08e-3  # Default tank wall thickness

        # Default parameters (can be overwritten by input_parameters)
        self.r_BL: float = 1.1  # Boundary layer growth rate
        self.revolve: float = 0  # Revolution angle (0 means 2D)
        self.n_revolve: float = 0  # Revolution angle (0 means 2D)
        self.symmetry: bool = True  # Revolution angle (0 means 2D)
        self.wedge_angle: float = 1  # Revolution angle if 2D
        self.wedge_pos_normal: list = [0, 0, 0]
        self.wedge_neg_normal: list = [0, 0, 0]
        self.surface_file = f"{self.name}.stl"
        self.non_coupled_cyclic = False
        self.patch_name_pos = "wedgePos"
        self.patch_name_neg = "wedgeNeg"
        self.ymax = tank.ymax()
        self.lid = False
        self.obstacle = False
        self.regions = ["gas", "liquid", "metal"]
        self.VoF = False
        self.empty_2d: bool = False  # Use 2D planar (empty BC) instead of wedge
        self.empty_2d_thickness: float = 0.0  # z-extrusion thickness; defaults to bulk_cell_size
        self.symmetry_normal: list = [0, 0, 0]  # Normal for symmetry plane (set when empty_2d=True)
        self.mirror: bool = False  # Mirror the mesh about the symmetry axis (requires empty_2d=True)
        self.extrude_cylinder: float = 0  # Extrude in z to create 3D cylinder (requires empty_2d=True)
        self.smoothing: bool = False  # Run laplacianMeshSmoother after splitMeshRegions
        self.wall_mesh_outlet: bool = True  # Wall mesh covers the outlet boundary

        if self.VoF:
            self.regions.remove("liquid")

        # self.check_openfoam_loaded(version="org")
        self.validate_parameters(input_parameters)
        self.set_parameters(input_parameters)
        if self.internal_outlet > 0:
            if "wall_mesh_outlet" in input_parameters and self.wall_mesh_outlet:
                raise WallMeshOutletRequiresNoInternalOutlet()
            self.wall_mesh_outlet = False
        if self.mirror and not self.empty_2d:
            raise MirrorRequiresEmpty2D()
        if self.extrude_cylinder and not self.empty_2d:
            raise ExtrudeCylinderRequiresEmpty2D()
        self.n_BL, self.t_BL, self.e_BL = self.calculate_boundary_layer()
        self.wedge_angle = self.revolve if self.revolve else self.wedge_angle
        if self.empty_2d and not self.empty_2d_thickness:
            self.empty_2d_thickness = self.bulk_cell_size
        self.calc_wedge_normal()
        super().__init__()
        self.cyclic = not self.symmetry
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
                elif type(value) in [list]:
                    content = [f"{v}" if type(v) is str else f"{v: .6g}" for v in value]
                    f.write(f"{key} ({' '.join(content)});\n")

            for key, value in self.tank.__dict__.items():
                if key not in self.__dict__ and type(value) in [int, float, str]:
                    f.write(f"{key} {value};\n")

    def run_command(self, command: str, no_output: bool = False, return_exception: bool = False) -> None | Exception:
        """Run a shell command with uv-style console feedback and file logging.

        While the command executes the terminal shows a transient spinner (or
        the active :class:`~rich.progress.Progress` task description is
        updated).  On success the spinner disappears cleanly; on failure a red
        error block is printed and a :exc:`~openfoam_tank_mesh.exceptions.CommandFailed`
        exception is raised (or returned when *return_exception* is ``True``).

        All commands, timings and output are written to the log file regardless
        of *no_output* or *debug*.
        """
        display_cmd = " ".join(word.replace(self.dict_path, f"<{self.name}.dict_path>") for word in command.split())
        logger.info("$ %s", display_cmd)

        t = time.time()
        if self._progress is not None and self._progress_task is not None:
            # Inside a Progress context - update description in place, no
            # nested Live display.
            self._progress.update(self._progress_task, description=f"[dim cyan]{display_cmd}[/dim cyan]")
            result = self._run_subprocess(command, capture_output=True)
        else:
            # Standalone call - transient spinner that disappears on success.
            with console.status(f"[cyan]{display_cmd}[/cyan]"):
                result = self._run_subprocess(command, capture_output=True)
        dt = time.time() - t

        if result.returncode == 0:
            logger.info("  ✓ %.3f s", dt)
            if self.debug:
                stdout = result.stdout.decode()
                if stdout.strip():
                    logger.debug(stdout)
        else:
            stderr = result.stderr.decode()
            stdout = result.stdout.decode()
            logger.error("  ✗ FAILED (%.3f s)  cmd=%s", dt, display_cmd)
            if stderr.strip():
                logger.error("stderr: %s", stderr)
            if stdout.strip():
                logger.error("stdout: %s", stdout)

            exception = CommandFailed(command, stderr)
            if return_exception:
                return exception

            out = self._progress.console if self._progress is not None else console
            out.print(f"[bold red]✗ {display_cmd}[/bold red] [red]failed ({dt:.3f} s)[/red]")
            if not no_output and stderr.strip():
                out.print(f"[red]{stderr}[/red]")
            raise exception

        return None

    def check_openfoam_loaded(self, version: str = "com") -> None:
        """
        Check if OpenFOAM is loaded.
        version: str ("org" or "com")
        """
        result = self._run_subprocess("simpleFoam -help", capture_output=True)
        if f"openfoam.{version}" not in result.stdout.decode():
            raise OpenFoamNotLoaded

    def _run_subprocess(self, command: str, capture_output: bool = True) -> CompletedProcess[bytes]:
        command_parts = shlex.split(command)
        if not command_parts:
            raise CommandFailed(command, "Empty command")
        if command_parts[0] not in self._ALLOWED_COMMANDS:
            raise CommandFailed(command, f"Disallowed command: {command_parts[0]}")
        return run(command_parts, capture_output=capture_output)  # noqa: S603

    def run_openfoam_utility(
        self, utility: str, foam_dict: str = "", return_exception: bool = False
    ) -> None | Exception:
        """
        Run an OpenFOAM utility with optional input dictionary.
        """
        command = f"{utility} -dict {self.dict_path}/{foam_dict}"
        self.sed(
            "include .*",
            f'include "{self.parameters_path}"',
            f"{self.dict_path}/{foam_dict}",
        )
        return self.run_command(command, return_exception=return_exception)

    def check_smoothing_utility(self) -> bool:
        """
        Check that ``laplacianMeshSmoother`` is available in ``PATH``.
        """
        if shutil.which("laplacianMeshSmoother"):
            return True

        msg = "laplacianMeshSmoother not found in PATH. Disabling smoothing."
        logger.error(msg)
        console.print(f"[bold red]ERROR:[/bold red] {msg}")
        self.smoothing = False
        return False

    def smooth_mesh(self) -> None:
        """
        Run Laplacian mesh smoothing when enabled and available.
        """
        if not self.smoothing:
            return
        if not self.check_smoothing_utility():
            return
        for region in ("gas", "liquid"):
            if region in self.regions:
                self.run_openfoam_utility(
                    f"laplacianMeshSmoother -overwrite -region {region}",
                    "laplacianSmoothDict",
                )

    def calculate_boundary_layer(self) -> tuple[int, float, float]:
        """
        Calculate boundary layer parameters.
        returns:
            n_bl: int = number of boundary layer cells
            t_bl: float = thickness of boundary layer
            e_bl: float = boundary layer expansion
        """
        if self.r_BL == 1 or self.wall_cell_size == self.wall_tan_cell_size:
            self.r_BL = 1.0
            return 4, self.wall_cell_size * 4, 1
        n = 1
        x = self.wall_cell_size
        t = self.wall_cell_size
        outer_cell_size = self.wall_tan_cell_size * (1 - 0.25 * (self.tri_bulk))
        while x <= outer_cell_size / self.r_BL:
            n += 1
            x *= self.r_BL
            t += x

        return n, t, x / self.wall_cell_size

    def check_mesh(self, regions: None | list = None) -> None:
        regions = [""] if not regions else [f"-region {region}" for region in regions]
        for region in regions:
            table = Table(title=f"Mesh Summary ({region})", show_header=False)
            output = run(["checkMesh", *region.split()], capture_output=True, text=True)  # noqa: S603, S607
            logger.info("checkMesh %s → returncode=%d", region, output.returncode)
            logger.debug(output.stdout)
            if "FAILED" in output.stdout:
                console.print(output.stdout)
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
        self.run_command("rm -r constant/metal/polyMesh")
        self.regions.remove("metal")
        for region in ("gas", "liquid"):
            # self.sed(f"{region}_to_metal", "walls", f"constant/{region}/polyMesh/boundary")
            # self.sed("mappedWall", "wall", f"constant/{region}/polyMesh/boundary")
            # Do the same sed command, but only for the LAST match!
            self.run_command(f"sed -i 'H;$!d;g;s/\\(.*\\)mappedWall/\\1wall/' constant/{region}/polyMesh/boundary")

    def add_wall(self, wall_thickness: float, n_layers: int) -> None:
        """
        Add a wall region to the mesh.
        """
        self.regions.append("metal")
        extrude_mesh_dict = "extrudeMeshDict.addWall"
        topo_set_dict = "topoSetDict.addWall"
        self.sed("thickness.*;", f"thickness {wall_thickness};", self.dict(extrude_mesh_dict))
        self.sed("nLayers.*;", f"nLayers {n_layers};", self.dict(extrude_mesh_dict))
        self.run_openfoam_utility("extrudeMesh", extrude_mesh_dict)
        self.run_openfoam_utility("topoSet", topo_set_dict)
        self.run_command("splitMeshRegions -cellZonesOnly -overwrite")
        self.sed("pipe", "outlet", "constant/metal/polyMesh/boundary")
        self.multi_region = True

    def create_internal_outlet(self) -> None:
        self.run_command("rm -rf 0/cellToRegion")
        self.run_openfoam_utility("topoSet", "topoSetDict.subsetMesh")
        self.run_command("cp -r 0 0.temp")
        self.run_command("subsetMesh cellsToKeep -overwrite -patch outlet")
        self.run_command("rm -r 0")
        self.run_command("mv 0.temp 0")
        self.run_openfoam_utility("topoSet", "topoSetDict.pipe2outlet")
        self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.pipe2outlet")
        # self.y_outlet = self.tank.y2 - self.internal_outlet
        self.write_mesh_parameters()

    def remove_wall_outlet(self) -> None:
        self.run_command("rm -rf 0/cellToRegion")
        self.run_openfoam_utility("topoSet", "topoSetDict.removeWallOutlet")
        self.run_command("cp -r 0 0.temp")
        self.run_command("subsetMesh cellsToKeep -overwrite -patch outlet")
        self.run_command("rm -r 0")
        self.run_command("mv 0.temp 0")
        self.run_openfoam_utility("topoSet", "topoSetDict.pipe2outlet")
        self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.pipe2outlet")
        # self.y_outlet = self.tank.y2 - self.internal_outlet
        self.write_mesh_parameters()

    def remove_outlet(self) -> None:
        """Remove fluid cells at the outlet while keeping wall (metal) cells.

        Used when ``wall_mesh_outlet=True``: the metal region is preserved at
        the outlet so that it covers the outlet boundary after
        ``splitMeshRegions``.
        """
        self.run_command("rm -rf 0/cellToRegion")
        self.run_openfoam_utility("topoSet", "topoSetDict.removeOutletFluidOnly")
        self.run_command("cp -r 0 0.temp")
        self.run_command("subsetMesh cellsToKeep -overwrite -patch outlet")
        self.run_command("rm -r 0")
        self.run_command("mv 0.temp 0")
        self.run_openfoam_utility("topoSet", "topoSetDict.pipe2outlet")
        self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.pipe2outlet")
        self.write_mesh_parameters()

    def add_wall_thickness(
        self, region: str, patchName: str, ranges: list[tuple[float, float]], thicknesses: list[float]
    ) -> None:
        """
        Extrude the wall outwards in region. Extrudes between y_start and y_end,
        specified in the list of ranges, each with their own thickness.
        """

        topo_dict_path = self.dict("topoSetDict.add_wall_thickness")
        create_dict_path = self.dict("createPatchDict.add_wall_thickness")
        extrude_dict_path = self.dict("extrudeMeshDict.add_wall_thickness")

        for p in [topo_dict_path, create_dict_path, extrude_dict_path]:
            self.sed("^patchName .*;", f"patchName {patchName};", p)
        for r, t in zip(ranges, thicknesses, strict=True):
            logger.info("Adding wall thickness %s in range %s", t, r)
            ys, ye = r
            n = max(3, int(self.n_wall_layers * (t / self.wall_thickness)))
            self.sed("(-1e6 .* 1e6)", f"(-1e6 {ys} -1e6) (1e6 {ye} 1e6)", topo_dict_path)
            self.run_openfoam_utility(f"topoSet -region {region}", "topoSetDict.add_wall_thickness")
            self.sed("^extrude .*;", "extrude true;", create_dict_path)
            self.run_openfoam_utility(f"createPatch -overwrite -region {region}", "createPatchDict.add_wall_thickness")
            self.sed("^thickness .*;", f"thickness {t};", extrude_dict_path)
            self.sed("^nLayers .*;", f"nLayers {n};", extrude_dict_path)
            self.run_openfoam_utility(f"extrudeMesh  -region {region}", "extrudeMeshDict.add_wall_thickness")
            self.sed("^extrude .*;", "extrude false;", create_dict_path)
            self.run_openfoam_utility(f"createPatch -overwrite -region {region}", "createPatchDict.add_wall_thickness")

    def extrude_outlet(self, length: float) -> None:
        n_layers = int(length / self.wall_tan_cell_size)
        dict_path = self.dict("extrudeMeshDict.outlet")
        self.sed("nLayers.*;", f"nLayers {n_layers};", dict_path)
        self.sed("thickness.*;", f"thickness {length};", dict_path)
        region = "gas" if self.multi_region else "region0"
        self.run_openfoam_utility(f"extrudeMesh -region {region}", "extrudeMeshDict.outlet")
        self.y_outlet += length
        self.write_mesh_parameters()

    def mirror_mesh(self) -> None:
        """Mirror the mesh about the symmetry plane (x = 0).

        Requires ``empty_2d=True``.  Runs ``mirrorMesh`` for each mesh region
        so that the full (un-halved) geometry is obtained.  The symmetry
        boundary at x = 0 is merged into the interior and disappears from the
        boundary list after this step.
        """
        mirror_dict = self.dict("mirrorMeshDict")
        for region in self.regions:
            self.run_command(f"mirrorMesh -dict {mirror_dict} -region {region} -overwrite")

    def do_extrude_cylinder(self) -> None:
        """Extrude the 2D planar mesh in the z-direction to create a 3D cylinder mesh.

        Requires ``empty_2d=True`` and ``extrude_cylinder > 0``.  For each
        mesh region the ``empty_pos`` patch (front face, outward normal +z) is
        extruded linearly by ``extrude_cylinder`` metres.  The number of layers
        is chosen so that each cell matches ``bulk_cell_size``.

        After extrusion:
        * The original ``empty_pos`` boundary (front face, z = dz) is renamed
          to ``front`` and its type is changed from ``empty`` to ``wall``.
        * The original ``empty_neg`` boundary (back face, z = 0) is renamed to
          ``back`` and its type is changed from ``empty`` to ``wall``.

        Call :meth:`remove_wall` before this method to delete the metal region.
        """
        n_layers = max(1, round(self.extrude_cylinder / self.bulk_cell_size))
        dict_path = self.dict("extrudeMeshDict.cylinder")
        self.sed("nLayers.*;", f"nLayers {n_layers};", dict_path)
        self.sed("thickness.*;", f"thickness {self.extrude_cylinder};", dict_path)
        for region in self.regions:
            self.run_openfoam_utility(f"extrudeMesh -region {region}", "extrudeMeshDict.cylinder")
            boundary_file = f"constant/{region}/polyMesh/boundary"
            # Rename empty_neg → back, empty_pos → front
            self.sed("empty_neg", "back", boundary_file)
            self.sed("empty_pos", "front", boundary_file)
            # Change both (formerly empty) patch types to wall
            self.sed("type.*empty;", "type wall;", boundary_file)

    def cfMesh(self, nLayers: int = 0) -> None:
        """
        Use cfMesh to create a 3D mesh.
        """

        self.generate_stl()
        self.run_command("rm -rf constant/polyMesh constant/metal/polyMesh constant/gas/polyMesh")
        self.run_command(f"cp {self.dict_path}/meshDict system/meshDict")
        self.sed("nLayers.*;", f"nLayers {nLayers};", "system/meshDict")
        self.sed("include .*", f'include "{self.parameters_path}"', "system/meshDict")
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

        # self.run_command(f"transformPoints -rotate-y -{self.wedge_angle / 2}")
        self.run_command(f'transformPoints "Ry=-{self.wedge_angle / 2}"')

        self.check_mesh()

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
        logger.warning("Non-coupled cyclic boundary conditions detected")
        console.print("\n[bold yellow]Non-coupled cyclic boundary conditions detected[/bold yellow]")
        console.print("    [yellow]Adding createNonConformalCouplesDict to system/[/yellow]")
        console.print("    [yellow]Load OpenFOAM 11 and run the following command:[/yellow]")
        console.print("    [bold yellow]createNonConformalCouples -overwrite[/bold yellow]\n")

    def calc_wedge_normal(self) -> None:
        """
        Calculate normals for the lateral boundary faces.

        For the axisymmetric wedge case the normals are derived from the wedge
        half-angle (rotation around y-axis).

        For the 2D planar case (``empty_2d=True``) the mesh is extruded one
        cell in the z-direction instead of being revolved.  The front/back
        faces therefore have normals along ±z, and the axis boundary at x=0
        becomes a symmetry plane with normal along -x.  The patch names are
        updated accordingly.
        """
        if self.empty_2d:
            self.wedge_pos_normal = [0, 0, 1]
            self.wedge_neg_normal = [0, 0, -1]
            self.symmetry_normal = [-1, 0, 0]
            self.patch_name_pos = "empty_pos"
            self.patch_name_neg = "empty_neg"
        else:
            angle = max(self.revolve, self.wedge_angle) / 2
            alpha = np.radians(angle)
            dx = np.cos(alpha)
            dz = np.sin(alpha)
            self.wedge_pos_normal = [-dz, 0, dx]
            self.wedge_neg_normal = [-dz, 0, -dx]


# Backward-compatible alias
TwoPhaseTankMesh = OpenFoamMeshPipeline
