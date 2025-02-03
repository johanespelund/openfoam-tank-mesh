from rich import print as rprint


class OutOfRange(Exception):
    def __init__(self, y: float) -> None:
        # User rprint() instead of print() to get red and bold text.
        rprint(f"[bold red]y = {y} is out of range.[/bold red]")

        # super().__init__(f"y = {y} is out of range.")


class MissingParameter(Exception):
    def __init__(self, missing_parameter: str) -> None:
        # super().__init__(f"Missing parameter: {missing_parameter}")
        rprint(f"[bold red]Missing parameter: {missing_parameter}[/bold red]")


class OpenFoamNotLoaded(Exception):
    def __init__(self) -> None:
        # super().__init__("OpenFOAM has not been loaded.")
        rprint("[bold red]OpenFOAM has not been loaded.[/bold red]")


class CommandFailed(Exception):
    def __init__(self, command: str) -> None:
        # super().__init__(f"Command '{command}' failed.")
        rprint(f"[bold red]Command '{command}' failed.[/bold red]")
