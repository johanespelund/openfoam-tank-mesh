from rich import print as rprint


class OutOfRange(Exception):
    def __init__(self, y: float) -> None:
        rprint(f"[bold red]y = {y} is out of range.[/bold red]")


class MissingParameter(Exception):
    def __init__(self, missing_parameter: str) -> None:
        rprint(f"[bold red]Missing parameter: {missing_parameter}[/bold red]")


class OpenFoamNotLoaded(Exception):
    def __init__(self) -> None:
        rprint("[bold red]OpenFOAM has not been loaded.[/bold red]")


class CommandFailed(Exception):
    def __init__(self, command: str) -> None:
        self.command = command
        super().__init__()
        # Q: How can I achieve the f

    def __rich__(self) -> str:
        return f"[bold red]Command {self.command} failed.[/bold red]"
