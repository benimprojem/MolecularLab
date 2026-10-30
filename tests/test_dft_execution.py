import subprocess

import pytest

from molecular_lab.calculations.dft_settings import (
    DFTCalculationSettings, DFTOutput, DFTTask, EnergyCutoff, EnergyUnit,
    KPointMesh, KPointScheme, PseudopotentialReference,
)
from molecular_lab.calculations import quantum_espresso as qe_module
from molecular_lab.calculations.execution import (
    DFTExecutionStatus, execute_quantum_espresso_scf,
)
from molecular_lab.structures.crystal import CrystalSite, CrystalStructure


@pytest.fixture(autouse=True)
def fake_ase_atoms(monkeypatch):
    class Cell:
        array = ((5.6, 0.0, 0.0), (0.0, 5.6, 0.0), (0.0, 0.0, 5.6))

    class Atoms:
        cell = Cell()

        def get_chemical_symbols(self):
            return ["Na", "Cl"]

        def get_masses(self):
            return [22.989769, 35.45]

        def get_scaled_positions(self, wrap=False):
            assert wrap is False
            return ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5))

    monkeypatch.setattr(qe_module, "to_ase_atoms", lambda structure: Atoms())


def make_structure():
    return CrystalStructure(
        structure_id="NaCl",
        lattice_vectors_angstrom=((5.6, 0, 0), (0, 5.6, 0), (0, 0, 5.6)),
        sites=(
            CrystalSite("Na1", "Na", (0, 0, 0)),
            CrystalSite("Cl1", "Cl", (0.5, 0.5, 0.5)),
        ),
    )


def make_settings(**overrides):
    pseudo = {
        element: PseudopotentialReference(
            element, f"{element}.UPF", "PBE", "test",
        )
        for element in ("Na", "Cl")
    }
    values = dict(
        task=DFTTask.SCF,
        exchange_correlation="PBE",
        energy_cutoff=EnergyCutoff(60, EnergyUnit.RYDBERG),
        kpoint_mesh=KPointMesh((2, 2, 2), KPointScheme.MONKHORST_PACK),
        pseudopotentials=pseudo,
        requested_outputs=(DFTOutput.TOTAL_ENERGY,),
    )
    values.update(overrides)
    return DFTCalculationSettings(**values)


def test_qe_missing_executable_is_reported(tmp_path, monkeypatch):
    monkeypatch.setattr("molecular_lab.calculations.execution.shutil.which", lambda _: None)
    result = execute_quantum_espresso_scf(
        make_structure(), make_settings(), working_directory=tmp_path,
        executable="not-a-real-pw.x",
    )
    assert result.status is DFTExecutionStatus.BACKEND_UNAVAILABLE
    assert result.total_energy_ev is None


def test_qe_missing_pseudopotential_is_reported(tmp_path, monkeypatch):
    monkeypatch.setattr("molecular_lab.calculations.execution.shutil.which", lambda _: "/usr/bin/pw.x")
    result = execute_quantum_espresso_scf(
        make_structure(), make_settings(), working_directory=tmp_path,
    )
    assert result.status is DFTExecutionStatus.BACKEND_UNAVAILABLE
    assert "missing" in result.message


def test_qe_converged_output_is_parsed_and_energy_converted(tmp_path, monkeypatch):
    (tmp_path / "pseudo").mkdir()
    for name in ("Na.UPF", "Cl.UPF"):
        (tmp_path / "pseudo" / name).write_text("test pseudo")
    monkeypatch.setattr("molecular_lab.calculations.execution.shutil.which", lambda _: "/usr/bin/pw.x")

    def fake_run(*args, **kwargs):
        assert kwargs["cwd"] == tmp_path.resolve()
        assert "calculation = 'scf'" in kwargs["input"]
        return subprocess.CompletedProcess(
            args[0], 0,
            "Program PWSCF v.7.3 starts\n"
            "!    total energy              =   -10.00000000 Ry\n"
            "convergence has been achieved\n"
            "JOB DONE.\n",
            "",
        )

    monkeypatch.setattr("molecular_lab.calculations.execution.subprocess.run", fake_run)
    result = execute_quantum_espresso_scf(
        make_structure(), make_settings(), working_directory=tmp_path,
    )
    assert result.status is DFTExecutionStatus.SUCCEEDED
    assert result.converged is True
    assert result.backend_version == "7.3"
    assert result.total_energy_ev == pytest.approx(-136.05693122994)


def test_qe_nonconverged_output_never_reports_energy_as_success(tmp_path, monkeypatch):
    (tmp_path / "pseudo").mkdir()
    for name in ("Na.UPF", "Cl.UPF"):
        (tmp_path / "pseudo" / name).write_text("test pseudo")
    monkeypatch.setattr("molecular_lab.calculations.execution.shutil.which", lambda _: "/usr/bin/pw.x")
    monkeypatch.setattr(
        "molecular_lab.calculations.execution.subprocess.run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0, "! total energy = -10 Ry\nJOB DONE.\n", "",
        ),
    )
    result = execute_quantum_espresso_scf(
        make_structure(), make_settings(), working_directory=tmp_path,
    )
    assert result.status is DFTExecutionStatus.NOT_CONVERGED
    assert result.total_energy_ev is None


def test_qe_timeout_is_explicit(tmp_path, monkeypatch):
    (tmp_path / "pseudo").mkdir()
    for name in ("Na.UPF", "Cl.UPF"):
        (tmp_path / "pseudo" / name).write_text("test pseudo")
    monkeypatch.setattr("molecular_lab.calculations.execution.shutil.which", lambda _: "/usr/bin/pw.x")

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], kwargs["timeout"], output="partial")

    monkeypatch.setattr("molecular_lab.calculations.execution.subprocess.run", timeout)
    result = execute_quantum_espresso_scf(
        make_structure(), make_settings(), working_directory=tmp_path,
        timeout_seconds=5,
    )
    assert result.status is DFTExecutionStatus.TIMED_OUT
    assert result.total_energy_ev is None


def test_only_scf_total_energy_is_supported(tmp_path):
    with pytest.raises(ValueError, match="SCF"):
        execute_quantum_espresso_scf(
            make_structure(), make_settings(task=DFTTask.RELAX),
            working_directory=tmp_path,
        )
    with pytest.raises(ValueError, match="total energy only"):
        execute_quantum_espresso_scf(
            make_structure(),
            make_settings(requested_outputs=(DFTOutput.FORCES,)),
            working_directory=tmp_path,
        )
