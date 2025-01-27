import time
from subprocess import run

from rich import print as rprint

from openfoam_tank_mesh.Tank import Tank


class TankMesh:
    """
    Base class for OpenFOAM tank meshes.
    """

    def __init__(self, tank: Tank, debug: bool = False) -> None:
        self.tank = tank
        self.debug = debug
        self.mesh_parameters: dict = {}

        self.name = tank.name
        self.fill_level = tank.fill_level
        self.y_interface = tank.y_interface
        self.y_outlet = tank.y_outlet
        self.area_liquid = tank.area_liquid
        self.area_gas = tank.area_gas

        self.check_openfoam_loaded()

    def write_mesh_parameters(self, filename: str) -> None:
        """
        Write the mesh parameters to a file.
        """
        with open(filename, "w") as f:
            for key, value in self.mesh_parameters.items():
                f.write(f"{key} {value};\n")

    def run_command(self, command: str, no_output: bool = False) -> None:
        """
        Method to run shell commands. The result should always be captured,
        and an error should be raised if the command fails (even if no_output is True).
        """
        t = time.time()
        if self.debug:
            rprint(f"Running command: {command}")

        result = run(command, shell=True, capture_output=no_output)  # noqa: S602
        dt = time.time() - t
        if result.returncode != 0:
            rprint(result.stderr.decode())
            rprint(f"Command failed: {command}")
            exit(1)

        if not no_output:
            rprint(f" ({dt:.6f} s)")

            if self.debug:
                rprint(result.stdout.decode())

    def check_openfoam_loaded(self, version: str = "org") -> bool:
        """
        Check if OpenFOAM is loaded.
        version: str ("org" or "com")
        """
        command = "simpleFoam -help"
        result = run(command, shell=True, capture_output=True)  # noqa: S602
        return f"openfoam.{version}" in result.stdout.decode()


if __name__ == "__main__":
    # Test the run_command method and check_openfoam_loaded function
    from openfoam_tank_mesh.SphericalTank import SphericalTank

    tank = SphericalTank("test", 0.5, 0.1, 1.0)
    mesh = TankMesh(tank, debug=False)
    mesh.run_command("simpleFoam -help", no_output=True)
