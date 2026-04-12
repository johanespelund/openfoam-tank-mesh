import pytest

from openfoam_tank_mesh.TwoPhaseTankMesh import TwoPhaseTankMesh


class _DummyMesh(TwoPhaseTankMesh):
    @property
    def dict_path(self) -> str:
        return ""

    @property
    def parameters_path(self) -> str:
        return ""

    def generate(self) -> None:
        return None

    def generate_stl(self) -> None:
        return None


def test_smooth_mesh_runs_only_for_gas_and_liquid(monkeypatch: pytest.MonkeyPatch):
    mesh = object.__new__(_DummyMesh)
    mesh.smoothing = True
    mesh.regions = ["gas", "liquid", "metal"]
    calls: list[str] = []

    monkeypatch.setattr(mesh, "check_smoothing_utility", lambda: True)
    monkeypatch.setattr(mesh, "run_openfoam_utility", lambda utility, foam_dict="": calls.append(utility))

    _DummyMesh.smooth_mesh(mesh)

    assert calls == [
        "laplacianMeshSmoother -overwrite -region gas",
        "laplacianMeshSmoother -overwrite -region liquid",
    ]
