from __future__ import annotations

import logging
import pathlib
import shutil
from typing import ClassVar

import numpy as np
from rich.panel import Panel
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from openfoam_tank_mesh.gmsh_scripts.stl import generate_3D_internal_outlet_stl, generate_3D_stl
from openfoam_tank_mesh.gmsh_scripts.two_phase import run as run_gmsh
from openfoam_tank_mesh.Profile import CylinderCapsTankProfile, KSiteProfile, SphereProfile, TankProfile
from openfoam_tank_mesh.TwoPhaseTankMesh import OpenFoamMeshPipeline, console


def _make_progress() -> Progress:
    """Return a :class:`~rich.progress.Progress` bar styled for mesh generation."""
    return Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
    )


class GmshMeshPipeline(OpenFoamMeshPipeline):
    """
    Base class for two-phase tank meshes that use the ``two_phase.py`` Gmsh script.

    Provides the standard mesh-generation pipeline
    (gmsh → gmshToFoam → topoSet → createPatch → splitMeshRegions)
    and eliminates boilerplate from the individual tank-mesh subclasses.

    Subclasses **must** implement
    --------------------------------
    ``_create_profile(input_parameters)``
        Create and return the :class:`~openfoam_tank_mesh.Profile.TankProfile`
        that describes the tank geometry.

    Subclasses **may** override
    ----------------------------
    ``_pre_init_setup()``
        Called at the very start of ``__init__``, before the profile is
        created.  Use for setup that must precede
        :class:`TwoPhaseTankMesh.__init__` (e.g. copying dict templates to a
        writable directory on disk).
    ``_pre_split_setup()``
        Called inside :meth:`generate` immediately before
        ``splitMeshRegions``.  Use to inject tank-specific OpenFOAM
        operations (e.g. lid / obstacle extrusion for the K-Site tank).
    ``dict_path``
        Override when the dict templates live outside the default package
        location.
    ``parameters_path``
        Override when a custom parameters-file name is required.
    """

    REQUIRED_PARAMETERS: ClassVar[list[str]] = [
        *OpenFoamMeshPipeline.REQUIRED_PARAMETERS,
        "r_BL",
        "internal_outlet",
        "n_wall_layers",
    ]

    def __init__(self, input_parameters: dict) -> None:
        # Allow subclasses to do setup before profile creation (e.g. writable dicts).
        self._pre_init_setup()

        # Ensure the outlet opening is large enough to be meshed.
        self.modify_outlet = False
        if input_parameters["outlet_radius"] <= input_parameters["wall_tan_cell_size"]:
            input_parameters["outlet_radius"] *= 2
            self.modify_outlet = True

        # Create the tank profile (subclass-specific geometry).
        self.tank = self._create_profile(input_parameters)

        super().__init__(tank=self.tank, input_parameters=input_parameters)
        self.multi_region = True
        self.n_wall_layers = input_parameters["n_wall_layers"]

    # ------------------------------------------------------------------
    # Hooks for subclasses
    # ------------------------------------------------------------------

    def _pre_init_setup(self) -> None:
        """Copy package dict templates to a writable working-directory location.

        This ensures the dict files can be modified in-place by
        :meth:`run_openfoam_utility` even when the package is installed
        read-only (e.g. in a site-packages directory).
        """
        self._work_dict_path = pathlib.Path.cwd() / "system" / "openfoam-tank-mesh" / "dicts"
        self._work_dict_path.mkdir(parents=True, exist_ok=True)
        pkg_dicts = pathlib.Path(__file__).parent / "dicts" / "two_phase_tanks"
        shutil.copytree(pkg_dicts, self._work_dict_path, dirs_exist_ok=True)

    def _create_profile(self, input_parameters: dict) -> TankProfile:
        """Create and return the tank profile.

        Subclasses must implement this method.
        """
        raise NotImplementedError(f"{type(self).__name__} must implement _create_profile()")

    def _pre_split_setup(self) -> None:
        """Run just before ``splitMeshRegions`` inside :meth:`generate`.

        Default implementation does nothing.
        """

    # ------------------------------------------------------------------
    # Mesh-generation pipeline
    # ------------------------------------------------------------------

    def gmsh(self) -> None:
        """Generate a Gmsh mesh and convert it to OpenFOAM format."""
        run_gmsh(self)
        if self.modify_outlet:
            self.outlet_radius /= 2
            self.write_mesh_parameters()

        self.run_command("gmshToFoam mesh.msh")

        if self.empty_2d:
            # 2D planar mesh: no rotation needed; use dedicated topoSet / createPatch dicts.
            self.run_openfoam_utility("topoSet", "topoSetDict.gmsh_empty")
            self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh_empty")
            self.patch_name_pos = "empty_pos"
            self.patch_name_neg = "empty_neg"
        elif self.revolve and self.symmetry:
            self.run_command(f'transformPoints "Ry={-self.wedge_angle / 2}"')
            self.run_openfoam_utility("topoSet", "topoSetDict.gmsh")
            self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh_symmetry")
            self.patch_name_pos = "symmetryPos"
            self.patch_name_neg = "symmetryNeg"
        elif self.revolve and self.cyclic:
            self.run_command(f'transformPoints "Ry={-self.wedge_angle / 2}"')
            self.run_openfoam_utility("topoSet", "topoSetDict.gmsh")
            self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh_cyclic")
            self.patch_name_pos = "cyclic_pos"
            self.patch_name_neg = "cyclic_neg"
        else:
            self.run_command(f'transformPoints "Ry={-self.wedge_angle / 2}"')
            self.run_openfoam_utility("topoSet", "topoSetDict.gmsh")
            self.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh")
            self.patch_name_pos = "wedgePos"
            self.patch_name_neg = "wedgeNeg"
        self.write_mesh_parameters()

    def generate(self) -> None:
        """Standard two-phase mesh-generation pipeline.

        Calls :meth:`gmsh`, handles the outlet, invokes the
        :meth:`_pre_split_setup` hook, then splits and checks the mesh.

        When ``mirror=True`` an extra stage mirrors each mesh region about the
        x = 0 symmetry plane using :meth:`mirror_mesh`.

        When ``extrude_cylinder > 0`` an extra stage removes the wall region
        and extrudes each remaining region in the z-direction to produce a 3D
        cylinder mesh (requires ``empty_2d=True``).

        A :class:`~rich.progress.Progress` bar tracks the high-level pipeline
        stages.  While each stage runs, the active sub-command is shown as the
        task description so the terminal never looks idle.
        """
        if self.smoothing:
            self.check_smoothing_utility()
        _STAGES = 5 + (1 if self.smoothing else 0) + (1 if self.mirror else 0) + (1 if self.extrude_cylinder else 0)
        _n = str(_STAGES)
        with _make_progress() as progress:
            task = progress.add_task("Mesh generation pipeline", total=_STAGES)
            self._progress = progress
            self._progress_task = task

            # ── Stage 1 ──────────────────────────────────────────────
            progress.update(task, description=f"[bold]Stage 1/{_n}[/bold] Gmsh mesh generation")
            self.gmsh()
            progress.advance(task)

            # ── Stage 2 ──────────────────────────────────────────────
            stage = 2
            progress.update(task, description=f"[bold]Stage {stage}/{_n}[/bold] Outlet configuration")
            if self.internal_outlet:
                self.create_internal_outlet()
            else:
                self.remove_wall_outlet()
            self.run_command("rm -rf 0/cellToRegion")
            progress.advance(task)

            stage += 1
            # ── Stage 3/4 ────────────────────────────────────────────
            progress.update(task, description=f"[bold]Stage {stage}/{_n}[/bold] Pre-split setup")
            self._pre_split_setup()
            progress.advance(task)

            stage += 1
            # ── Stage 4/5 ────────────────────────────────────────────
            progress.update(task, description=f"[bold]Stage {stage}/{_n}[/bold] Splitting mesh regions")
            self.run_command("splitMeshRegions -cellZonesOnly -overwrite")
            self.run_command("rm -rf constant/polyMesh")
            progress.advance(task)

            stage += 1
            # ── Optional smoothing stage ─────────────────────────────
            if self.smoothing:
                progress.update(task, description=f"[bold]Stage {stage}/{_n}[/bold] Laplacian smoothing")
                self.smooth_mesh()
                progress.advance(task)
                stage += 1

            # ── Optional extrusion stage ─────────────────────────────
            if self.extrude_cylinder:
                progress.update(task, description=f"[bold]Stage {stage}/{_n}[/bold] Cylinder extrusion")
                self.remove_wall()
                self.do_extrude_cylinder()
                progress.advance(task)
                stage += 1

            # ── Optional mirror stage ────────────────────────────────
            if self.mirror:
                progress.update(task, description=f"[bold]Stage {stage}/{_n}[/bold] Mirroring mesh")
                self.mirror_mesh()
                progress.advance(task)
                stage += 1

            # ── Final stage ───────────────────────────────────────────
            progress.update(task, description=f"[bold]Stage {_n}/{_n}[/bold] Checking mesh quality")
            progress.advance(task)

            self._progress = None
            self._progress_task = None

        # Render the quality table after the progress bar has finished.
        self.check_mesh(regions=self.regions)

        log_path = pathlib.Path.cwd() / "mesh_generation.log"
        logging.getLogger(__name__).info("Mesh generation complete")
        console.print(
            Panel(
                f"[bold green]Mesh generation complete[/bold green]\n"
                f"[dim]Full log written to [underline]{log_path}[/underline][/dim]",
                expand=False,
            )
        )

    def generate_stl(self) -> None:
        """Generate an STL file for use with cfMesh."""
        if self.internal_outlet:
            generate_3D_internal_outlet_stl(self)
        else:
            generate_3D_stl(self)

    # ------------------------------------------------------------------
    # Paths (may be overridden by subclasses)
    # ------------------------------------------------------------------

    @property
    def dict_path(self) -> str:
        """Path to the writable copy of the OpenFOAM dict templates."""
        return str(self._work_dict_path)

    @property
    def parameters_path(self) -> str:
        """Path to the mesh-parameters file (named after the tank profile)."""
        return f"{pathlib.Path.cwd()}/parameters.{self.name}"


# ---------------------------------------------------------------------------
# Concrete mesh classes
# ---------------------------------------------------------------------------


class KSiteMesh(GmshMeshPipeline):
    """Two-phase mesh for the NASA K-Site cryogenic tank (Stochl & Knoll, 1991).

    The tank geometry is taken from
    :class:`~openfoam_tank_mesh.Profile.KSiteProfile` (hard-coded K-Site
    dimensions).  Dict templates are copied to a writable sub-directory of
    the OpenFOAM case before any meshing is performed.

    .. seealso:: :class:`CylinderCapsMesh` for the equivalent mesh with
        user-specified geometry.
    """

    # -- hooks ----------------------------------------------------------

    def _create_profile(self, input_parameters: dict) -> KSiteProfile:
        return KSiteProfile(
            fill_level=input_parameters["fill_level"],
            outlet_radius=input_parameters["outlet_radius"],
            bulk_cell_size=input_parameters["bulk_cell_size"],
            wall_tan_cell_size=input_parameters["wall_tan_cell_size"],
            wall_cell_size=input_parameters["wall_cell_size"],
            r_BL=input_parameters["r_BL"],
            internal_outlet=input_parameters["internal_outlet"],
            wall_thickness=input_parameters.get("wall_thickness", 2.08e-3),
        )

    def _pre_split_setup(self) -> None:
        """K-Site-specific obstacle and lid extrusion before ``splitMeshRegions``."""
        if self.obstacle:
            self.add_wall_thickness("region0", "walls", [(1.849, 1e6)], [self.wall_thickness])
            self.add_wall_thickness("region0", "walls", [(1.859, 1e6)], [10e-3 - self.wall_thickness])
            self.add_wall_thickness("region0", "walls", [(0.9136, 0.9605)], [self.wall_thickness])

        if self.lid:
            self.regions.append("lid")
            self.write_mesh_parameters()

            y = self.tank.y_lid  # type: ignore[attr-defined]
            r = self.tank.r_lid  # type: ignore[attr-defined]
            nx, ny = self.tank.get_normal(y)
            r1, y1 = r - nx / 4, y - ny / 4
            r2, y2 = r + nx, y + ny
            y_duct = y - ny * self.wall_thickness

            topodict = self.dict("topoSetDict.splitMetalRegions")
            self.sed("radius1 .*;", f"radius1 {r1:.4f};", topodict)
            self.sed("radius2 .*;", f"radius2 {r2:.4f};", topodict)
            self.sed("point1 .*;", f"point1 (0 {y1:.4f} 0);", topodict)
            self.sed("point2 .*;", f"point2 (0 {y2:.4f} 0);", topodict)

            topodict = self.dict("topoSetDict.metal_patches")
            self.sed("p2 .*;", f"p2 (0 {y_duct:.4f} 0);", topodict)

            self.run_openfoam_utility("topoSet", "topoSetDict.splitMetalRegions")

        # Build the ring of y/r/w/h values used by the obstacle extrusion.
        ys = [self.tank.y_lid]  # type: ignore[attr-defined]
        rs = [self.tank.r_lid]  # type: ignore[attr-defined]
        ws, hs = [0.05], [0.02]
        r0 = 0.16
        n_subdivisions = 10
        for _i in range(n_subdivisions):
            _w = (0.16 - self.outlet_radius) / n_subdivisions
            _r = r0 - _i * _w
            if _i > 1:
                _r += self.wall_tan_cell_size
                _w += self.wall_tan_cell_size
            elif _i < n_subdivisions - 1:
                _w += self.wall_tan_cell_size
            _y = self.tank.get_y(_r, 1.5, self.y_outlet)
            ys.append(_y)
            rs.append(_r)
            ws.append(_w)
            hs.append(0.008)

        if self.obstacle:
            for y, r, w, h in zip(ys, rs, ws, hs, strict=True):
                y_average = (y + self.tank.get_y(r - w, 1.5, self.y_outlet)) / 2
                n = self.tank.get_normal(y_average)
                n = np.array([n[0], n[1], 0])
                tr, ty = n[1], -n[0]
                topodict = self.dict("topoSetDict.obstacle")
                region = "lid" if self.lid else "metal"
                self.sed("^obstacleRegion .*;", f"obstacleRegion {region};", topodict)

                origin = np.array([r, y, -1e6]) - n * 5 * h
                iHat = [tr * w, ty * w, 0]
                jHat = np.array([n[0] * h, n[1] * h, 0]) * 6
                kHat = [0, 0, 2e6]

                self.sed("origin .*;", f"origin ({origin[0]:.6f} {origin[1]:.6f} {origin[2]:.6f});", topodict)
                self.sed("i .*(.*;", f"i ({iHat[0]:.6f} {iHat[1]:.6f} {iHat[2]:.6f});", topodict)
                self.sed("j .*(.*;", f"j ({jHat[0]:.6f} {jHat[1]:.6f} {jHat[2]:.6f});", topodict)
                self.sed("k .*(.*);", f"k ({kHat[0]:.6f} {kHat[1]:.6f} {kHat[2]:.6f});", topodict)

                self.run_openfoam_utility("topoSet", "topoSetDict.obstacle")

    # -- overrides ------------------------------------------------------

    def generate(self) -> None:
        """K-Site generation: verify OpenFOAM is loaded, then run the standard pipeline."""
        self.check_openfoam_loaded(version="org")
        super().generate()

    def generate_leak_boundaries(self, region: str) -> None:
        """Create flange boundaries on the given metal region."""
        self.run_openfoam_utility(f"topoSet -region {region}", "topoSetDict.metal_patches")
        self.run_openfoam_utility(f"createPatch -overwrite -region {region}", "createPatchDict.metal_patches")

    @property
    def parameters_path(self) -> str:
        return f"{pathlib.Path.cwd()}/parameters.KSiteMesh"

    # -- K-Site heat-transfer helpers -----------------------------------
    # Values from Table 1 in DOI: 10.2514/6.2016-4674

    def q_insulation(self) -> float:
        """Heat flux for the MLI/insulation at T_amb = 350 K."""
        return 2.952

    def Q_insulation(self) -> float:
        """Total heat loss from insulation + support."""
        return 41.352 + 2.813 + 3.194 + 0.879

    def Q_ducts(self) -> float:
        """Heat loss from ducting and wires."""
        return 0


class SphereMesh(GmshMeshPipeline):
    """Two-phase mesh for a spherical tank."""

    def _create_profile(self, input_parameters: dict) -> SphereProfile:
        return SphereProfile(
            radius=input_parameters["radius"],
            fill_level=input_parameters["fill_level"],
            outlet_radius=input_parameters["outlet_radius"],
            bulk_cell_size=input_parameters["bulk_cell_size"],
            wall_tan_cell_size=input_parameters["wall_tan_cell_size"],
            wall_cell_size=input_parameters["wall_cell_size"],
            r_BL=input_parameters["r_BL"],
            internal_outlet=input_parameters["internal_outlet"],
            wall_thickness=input_parameters.get("wall_thickness", 2.08e-3),
        )

    def generate_leak_boundaries(self) -> None:
        """Create leak boundaries on metal and lid regions."""
        for region in ["metal", "lid"]:
            self.run_openfoam_utility(f"topoSet -region {region}", "topoSetDict.metal_patches")
            self.run_openfoam_utility(f"createPatch -overwrite -region {region}", "createPatchDict.metal_patches")

    @property
    def parameters_path(self) -> str:
        return f"{pathlib.Path.cwd()}/parameters.SphereMesh"


class CylinderCapsMesh(GmshMeshPipeline):
    """Two-phase mesh for a general tank with a cylindrical midsection and ellipsoidal caps.

    This is the general-purpose mesh class.  Supplying the K-Site dimensions
    (see :class:`~openfoam_tank_mesh.Profile.KSiteProfile`) produces a mesh
    that is geometrically equivalent to :class:`KSiteMesh`.

    Required keys in ``input_parameters``
    --------------------------------------
    In addition to the base-class requirements:

    * ``cylinder_radius`` *or* ``cylinder_diameter`` — radius / diameter of
      the cylindrical section.
    * ``cylinder_height`` — axial height of the cylindrical section.
    * ``cap_height`` — axial height of each ellipsoidal cap.
    """

    def _create_profile(self, input_parameters: dict) -> CylinderCapsTankProfile:
        if "cylinder_radius" not in input_parameters and "cylinder_diameter" not in input_parameters:
            raise ValueError(  # noqa: TRY003
                "CylinderCapsMesh requires either 'cylinder_radius' or 'cylinder_diameter' in input_parameters."
            )
        return CylinderCapsTankProfile(
            cylinder_radius=input_parameters.get("cylinder_radius", 0),
            cylinder_height=input_parameters["cylinder_height"],
            cap_height=input_parameters["cap_height"],
            fill_level=input_parameters["fill_level"],
            outlet_radius=input_parameters["outlet_radius"],
            bulk_cell_size=input_parameters["bulk_cell_size"],
            wall_tan_cell_size=input_parameters["wall_tan_cell_size"],
            wall_cell_size=input_parameters["wall_cell_size"],
            r_BL=input_parameters["r_BL"],
            internal_outlet=input_parameters["internal_outlet"],
            cylinder_diameter=input_parameters.get("cylinder_diameter"),
            wall_thickness=input_parameters.get("wall_thickness", 2.08e-3),
        )

    @property
    def parameters_path(self) -> str:
        return f"{pathlib.Path.cwd()}/parameters.CylinderCapsMesh"


# Backward-compatible alias
TwoPhaseGmshMesh = GmshMeshPipeline
