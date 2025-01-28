class OutOfRange(Exception):
    def __init__(self, y: float) -> None:
        super().__init__(f"y = {y} is out of range.")


class MissingParameter(Exception):
    def __init__(self, missing_parameter: str) -> None:
        super().__init__(f"Missing parameter: {missing_parameter}")


class OpenFoamNotLoaded(Exception):
    def __init__(self) -> None:
        super().__init__("OpenFOAM has not been loaded.")


class CommandFailed(Exception):
    def __init__(self, command: str) -> None:
        super().__init__(f"Command '{command}' failed.")
