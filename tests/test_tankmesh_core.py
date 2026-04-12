from __future__ import annotations

from types import SimpleNamespace

import pytest

from openfoam_tank_mesh.exceptions import CommandFailed, MissingParameter
from openfoam_tank_mesh.mesh_pipeline import OpenFoamMeshPipeline


class _DummyPipeline(OpenFoamMeshPipeline):
    def __init__(self, input_parameters: dict, parameters_path: str, dict_path: str):
        self._parameters_path = parameters_path
        self._dict_path = dict_path
        tank = SimpleNamespace(
            name="dummy",
            fill_level=0.5,
            y_interface=0.0,
            y_outlet=0.1,
            outlet_radius=0.01,
            interface_radius=0.2,
            area_liquid=1.0,
            area_gas=1.0,
            ymax=lambda: 1.0,
        )
        super().__init__(tank=tank, input_parameters=input_parameters)

    @property
    def dict_path(self) -> str:
        return self._dict_path

    @property
    def parameters_path(self) -> str:
        return self._parameters_path

    def generate(self) -> None:
        return None

    def generate_stl(self) -> None:
        return None


def _valid_parameters() -> dict:
    return {
        "bulk_cell_size": 0.1,
        "wall_cell_size": 0.01,
        "outlet_radius": 0.01,
        "fill_level": 0.5,
        "debug": False,
    }


def test_validate_parameters_missing_required_key(tmp_path):
    params = _valid_parameters()
    del params["fill_level"]
    with pytest.raises(MissingParameter):
        _DummyPipeline(params, str(tmp_path / "parameters"), str(tmp_path / "dicts"))


def test_run_openfoam_utility_builds_command_and_updates_include(monkeypatch: pytest.MonkeyPatch, tmp_path):
    mesh = _DummyPipeline(_valid_parameters(), str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    sed_calls: list[tuple[str, str, str]] = []
    commands: list[str] = []

    monkeypatch.setattr(mesh, "sed", lambda orig, new, path: sed_calls.append((orig, new, path)))
    monkeypatch.setattr(mesh, "run_command", lambda command, **kwargs: commands.append(command))

    mesh.run_openfoam_utility("createPatch -overwrite", "createPatchDict.gmsh")

    assert sed_calls == [
        (
            "include .*",
            f'include "{mesh.parameters_path}"',
            f"{mesh.dict_path}/createPatchDict.gmsh",
        )
    ]
    assert commands == [f"createPatch -overwrite -dict {mesh.dict_path}/createPatchDict.gmsh"]


def test_run_command_raises_command_failed(monkeypatch: pytest.MonkeyPatch, tmp_path):
    mesh = _DummyPipeline(_valid_parameters(), str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    fake_result = SimpleNamespace(returncode=1, stderr=b"boom", stdout=b"")

    monkeypatch.setattr(mesh, "_run_subprocess", lambda *args, **kwargs: fake_result)

    with pytest.raises(CommandFailed):
        mesh.run_command("bad-command")
