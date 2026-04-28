from __future__ import annotations

from types import SimpleNamespace

import pytest

from openfoam_tank_mesh.exceptions import CommandFailed, MissingParameter, WallMeshOutletRequiresNoInternalOutlet
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


def test_wall_mesh_outlet_defaults_to_true(tmp_path):
    """wall_mesh_outlet should be True by default."""
    mesh = _DummyPipeline(_valid_parameters(), str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    assert mesh.wall_mesh_outlet is True


def test_wall_mesh_outlet_can_be_set_false(tmp_path):
    """wall_mesh_outlet=False should be accepted."""
    params = {**_valid_parameters(), "wall_mesh_outlet": False}
    mesh = _DummyPipeline(params, str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    assert mesh.wall_mesh_outlet is False


def test_wall_mesh_outlet_true_with_internal_outlet_raises(tmp_path):
    """wall_mesh_outlet=True is not possible when internal_outlet > 0."""
    params = {**_valid_parameters(), "internal_outlet": 0.05, "wall_mesh_outlet": True}
    with pytest.raises(WallMeshOutletRequiresNoInternalOutlet):
        _DummyPipeline(params, str(tmp_path / "parameters"), str(tmp_path / "dicts"))


def test_wall_mesh_outlet_auto_false_when_internal_outlet_set(tmp_path):
    """When internal_outlet > 0 and wall_mesh_outlet is not explicitly set, it should be forced False."""
    params = {**_valid_parameters(), "internal_outlet": 0.05}
    mesh = _DummyPipeline(params, str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    assert mesh.wall_mesh_outlet is False


def test_wall_mesh_outlet_false_explicit_with_internal_outlet(tmp_path):
    """wall_mesh_outlet=False is allowed with internal_outlet > 0."""
    params = {**_valid_parameters(), "internal_outlet": 0.05, "wall_mesh_outlet": False}
    mesh = _DummyPipeline(params, str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    assert mesh.wall_mesh_outlet is False


def test_lid_defaults_to_zero(tmp_path):
    """lid should default to 0.0 (no lid) and has_lid should be False."""
    mesh = _DummyPipeline(_valid_parameters(), str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    assert mesh.lid == 0.0
    assert mesh.has_lid is False


def test_lid_positive_float(tmp_path):
    """lid > 0 should be accepted as the y-position of the lid boundary."""
    params = {**_valid_parameters(), "lid": 1.5}
    mesh = _DummyPipeline(params, str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    assert mesh.lid == 1.5
    assert mesh.has_lid is True


def test_lid_zero_means_no_lid(tmp_path):
    """lid = 0.0 means no lid region is created."""
    params = {**_valid_parameters(), "lid": 0.0}
    mesh = _DummyPipeline(params, str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    assert mesh.lid == 0.0
    assert mesh.has_lid is False
    assert "lid" not in mesh.regions


def test_lid_negative_means_no_lid(tmp_path):
    """lid <= 0 means no lid region is created."""
    params = {**_valid_parameters(), "lid": -1.0}
    mesh = _DummyPipeline(params, str(tmp_path / "parameters"), str(tmp_path / "dicts"))
    assert mesh.lid == -1.0
    assert mesh.has_lid is False
    assert "lid" not in mesh.regions


def test_lid_written_to_parameters_file_on_init_default(tmp_path):
    """lid 0.0 must appear in the parameters file after initialisation (default case).

    The topoSetDict.splitMetalAtYLid dict resolves ``$lid`` via ``#include``.
    If lid is a bool (the old default ``False``) write_mesh_parameters skips it
    and OpenFOAM cannot resolve the variable.  This test guards against that.
    """
    params_path = str(tmp_path / "parameters")
    mesh = _DummyPipeline(_valid_parameters(), params_path, str(tmp_path / "dicts"))
    content = (tmp_path / "parameters").read_text()
    assert "lid 0.0;" in content


def test_lid_written_to_parameters_file_on_init_positive(tmp_path):
    """A positive lid y-position must appear in the parameters file after init."""
    params_path = str(tmp_path / "parameters")
    params = {**_valid_parameters(), "lid": 1.5}
    _DummyPipeline(params, params_path, str(tmp_path / "dicts"))
    content = (tmp_path / "parameters").read_text()
    assert "lid 1.5;" in content

