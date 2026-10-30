import subprocess

import pytest

from molecular_lab.calculations.dft_settings import (
    DFTCalculationSettings, DFTOutput, DFTTask, EnergyCutoff, EnergyUnit,
    KPointMesh, KPointScheme, PseudopotentialReference,
)
from molecular_lab.calculations import quantum_espresso as qe_module
from molecular_lab.calculations.execution import (
    DFTExecutionStatus, execute_quantum_espresso_scf, execute_quantum_espresso, execute_dft_calculation,
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



def test_qe_relax_returns_final_structure_forces_and_energy(tmp_path, monkeypatch):
    (tmp_path / "pseudo").mkdir()
    for name in ("Na.UPF", "Cl.UPF"):
        (tmp_path / "pseudo" / name).write_text("test pseudo")
    monkeypatch.setattr("molecular_lab.calculations.execution.shutil.which", lambda _: "/usr/bin/pw.x")
    output = (
        "Program PWSCF v.7.3 starts\\n"
        "! total energy = -12.00000000 Ry\\n"
        "Forces acting on atoms (cartesian axes, Ry/au):\\n"
        " atom 1 type 1 force = 0.001 0.000 0.000\\n"
        " atom 2 type 2 force = -0.001 0.000 0.000\\n"
        "bfgs converged in 4 scf cycles and 3 bfgs steps\\n"
        "Begin final coordinates\\n"
        "CELL_PARAMETERS (angstrom)\\n"
        "5.6 0 0\\n0 5.6 0\\n0 0 5.6\\n"
        "ATOMIC_POSITIONS (crystal)\\n"
        "Na 0.01 0.00 0.00\\n"
        "Cl 0.49 0.50 0.50\\n"
        "End final coordinates\\n"
        "JOB DONE.\\n"
    )
    monkeypatch.setattr(
        "molecular_lab.calculations.execution.subprocess.run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 0, output, ""),
    )
    value = make_settings(
        task=DFTTask.RELAX,
        convergence=__import__(
            "molecular_lab.calculations.dft_settings",
            fromlist=["ConvergenceSettings"],
        ).ConvergenceSettings(force_tolerance_ev_per_angstrom=0.03),
        requested_outputs=(DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES),
    )
    result = execute_quantum_espresso(
        make_structure(), value, working_directory=tmp_path,
    )
    assert result.status is DFTExecutionStatus.SUCCEEDED
    assert result.total_energy_ev == pytest.approx(-163.268317475928)
    assert result.optimized_structure is not None
    assert result.optimized_structure.sites[0].fractional_coordinates == (0.01, 0.0, 0.0)
    assert len(result.forces_ev_per_angstrom) == 2
    assert result.forces_ev_per_angstrom[0][0] == pytest.approx(
        0.001 * 13.605693122994 / 0.529177210903
    )


def test_common_dispatch_rejects_unknown_backend():
    result = execute_dft_calculation("missing-engine", make_structure(), make_settings())
    assert result.status is DFTExecutionStatus.BACKEND_UNAVAILABLE
    assert "Unknown" in result.message


def test_relax_requires_explicit_force_tolerance(tmp_path):
    with pytest.raises(ValueError, match="force_tolerance"):
        execute_quantum_espresso(
            make_structure(),
            make_settings(task=DFTTask.RELAX),
            working_directory=tmp_path,
        )


def test_scf_compatibility_wrapper_rejects_relax(tmp_path):
    with pytest.raises(ValueError, match="SCF tasks only"):
        execute_quantum_espresso_scf(
            make_structure(),
            make_settings(task=DFTTask.RELAX),
            working_directory=tmp_path,
        )
