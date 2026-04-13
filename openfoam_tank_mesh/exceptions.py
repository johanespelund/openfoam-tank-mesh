import logging

from rich.console import Console

console = Console()
logger = logging.getLogger(__name__)


class OutOfRange(Exception):
    def __init__(self, value: float) -> None:
        logger.error("y = %s is out of range", value)
        console.print(f"[bold red]y = {value} is out of range.[/bold red]")


class SegmentNotInitialized(RuntimeError):
    def __init__(self) -> None:
        super().__init__("Segment not fully initialized.")


class BoundaryLayerTooThick(RuntimeError):
    def __init__(self) -> None:
        super().__init__("Boundary layer too thick!")


class SegmentsNotConnected(ValueError):
    def __init__(self, name1: str, name2: str) -> None:
        super().__init__(f"Segments {name1} and {name2} are not connected.")


class MissingParameter(Exception):
    def __init__(self, missing_parameter: str) -> None:
        logger.error("Missing parameter: %s", missing_parameter)
        console.print(f"[bold red]Missing parameter: {missing_parameter}[/bold red]")


class OpenFoamNotLoaded(Exception):
    def __init__(self) -> None:
        logger.error("OpenFOAM has not been loaded")
        console.print("[bold red]OpenFOAM has not been loaded.[/bold red]")


class MirrorRequiresEmpty2D(ValueError):
    def __init__(self) -> None:
        msg = "mirror=True requires empty_2d=True"
        logger.error(msg)
        console.print(f"[bold red]{msg}[/bold red]")
        super().__init__(msg)


class ExtrudeCylinderRequiresEmpty2D(ValueError):
    def __init__(self) -> None:
        msg = "extrude_cylinder > 0 requires empty_2d=True"
        logger.error(msg)
        console.print(f"[bold red]{msg}[/bold red]")
        super().__init__(msg)


class CommandFailed(Exception):
    def __init__(self, command: str, output: str = "") -> None:
        self.command = command
        logger.error("Command failed: %s", command)
        if output.strip():
            logger.error("Output: %s", output)
        super().__init__()

    def __rich__(self) -> str:
        return f"[bold red]Command {self.command} failed.[/bold red]"
