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


class CommandFailed(Exception):
    def __init__(self, command: str, output: str = "") -> None:
        self.command = command
        logger.error("Command failed: %s", command)
        if output.strip():
            logger.error("Output: %s", output)
        super().__init__()

    def __rich__(self) -> str:
        return f"[bold red]Command {self.command} failed.[/bold red]"
