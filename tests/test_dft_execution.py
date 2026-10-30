import subprocess
import sys
import types

import pytest

from molecular_lab.calculations.dft_settings import (
    ConvergenceSettings, DFTCalculationSettings, DFTOutput, DFTTask, EnergyCutoff, EnergyUnit,
    KPointMesh, KPointScheme, PseudopotentialReference,
)
from molecular_lab.calculations import quantum_espresso as qe_module
from molecular_lab.calculations import gpaw as gpaw_module
from molecular_lab.calculations.execution import (
    DFTExecutionStatus, execute_quantum_espresso_scf, execute_quantum_espresso, execute_dft_calculation, execute_gpaw,
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
        "Program PWSCF v.7.3 starts\n"
        "! total energy = -12.00000000 Ry\n"
        "Forces acting on atoms (cartesian axes, Ry/au):\n"
        " atom 1 type 1 force = 0.001 0.000 0.000\n"
        " atom 2 type 2 force = -0.001 0.000 0.000\n"
        "bfgs converged in 4 scf cycles and 3 bfgs steps\n"
        "Begin final coordinates\n"
        "CELL_PARAMETERS (angstrom)\n"
        "5.6 0 0\n0 5.6 0\n0 0 5.6\n"
        "ATOMIC_POSITIONS (crystal)\n"
        "Na 0.01 0.00 0.00\n"
        "Cl 0.49 0.50 0.50\n"
        "End final coordinates\n"
        "JOB DONE.\n"
    )
    monkeypatch.setattr(
        "molecular_lab.calculations.execution.subprocess.run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 0, output, ""),
    )
    value = make_settings(
        task=DFTTask.RELAX,
        convergence=ConvergenceSettings(force_tolerance_ev_per_angstrom=0.03),
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



def test_gpaw_relax_returns_optimized_structure_and_forces(monkeypatch):
    class Cell:
        array = ((5.6, 0.0, 0.0), (0.0, 5.6, 0.0), (0.0, 0.0, 5.6))

    class FakeAtoms:
        cell = Cell()
        calc = None

        def get_potential_energy(self):
            return -4.5

        def get_forces(self):
            return ((0.001, 0.0, 0.0), (-0.001, 0.0, 0.0))

        def get_scaled_positions(self, wrap=False):
            assert wrap is False
            return ((0.02, 0.0, 0.0), (0.48, 0.5, 0.5))

    class FakeCalculator:
        scf = types.SimpleNamespace(converged=True)

    class FakeBFGS:
        def __init__(self, atoms, logfile=None):
            self.atoms = atoms
            self.logfile = logfile

        def run(self, fmax):
            assert fmax == pytest.approx(0.025)
            return True

    fake_gpaw = types.ModuleType("gpaw")
    fake_gpaw.__version__ = "test-gpaw"
    fake_gpaw.GPAW = lambda **kwargs: FakeCalculator()
    fake_gpaw.PW = lambda cutoff: ("PW", cutoff)
    fake_ase_optimize = types.ModuleType("ase.optimize")
    fake_ase_optimize.BFGS = FakeBFGS
    monkeypatch.setitem(sys.modules, "gpaw", fake_gpaw)
    monkeypatch.setitem(sys.modules, "ase.optimize", fake_ase_optimize)
    monkeypatch.setattr(gpaw_module, "to_ase_atoms", lambda structure: FakeAtoms())

    value = make_settings(
        task=DFTTask.RELAX,
        convergence=ConvergenceSettings(force_tolerance_ev_per_angstrom=0.025),
        requested_outputs=(DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES),
    )
    result = execute_gpaw(make_structure(), value)
    assert result.status is DFTExecutionStatus.SUCCEEDED
    assert result.backend_version == "test-gpaw"
    assert result.total_energy_ev == pytest.approx(-4.5)
    assert result.optimized_structure is not None
    assert result.optimized_structure.sites[0].fractional_coordinates == (0.02, 0.0, 0.0)
    assert result.forces_ev_per_angstrom[1] == (-0.001, 0.0, 0.0)



def test_qe_vc_relax_normalizes_stress_and_requires_pressure_tolerance(tmp_path, monkeypatch):
    (tmp_path / "pseudo").mkdir()
    for name in ("Na.UPF", "Cl.UPF"):
        (tmp_path / "pseudo" / name).write_text("test pseudo")
    monkeypatch.setattr("molecular_lab.calculations.execution.shutil.which", lambda _: "/usr/bin/pw.x")
    output = (
        "Program PWSCF v.7.5 starts\\n"
        "! total energy = -13.00000000 Ry\\n"
        "Forces acting on atoms (cartesian axes, Ry/au):\\n"
        " atom 1 type 1 force = 0.001 0.000 0.000\\n"
        " atom 2 type 2 force = -0.001 0.000 0.000\\n"
        "total   stress  (Ry/bohr**3)                   (kbar)     P=  10.00\\n"
        " 0.00000000 0.00000000 0.00000000 10.00 0.00 0.00\\n"
        " 0.00000000 0.00000000 0.00000000 0.00 20.00 0.00\\n"
        " 0.00000000 0.00000000 0.00000000 0.00 0.00 30.00\\n"
        "bfgs converged in 5 scf cycles and 4 bfgs steps\\n"
        "Begin final coordinates\\n"
        "CELL_PARAMETERS (angstrom)\\n"
        "5.7 0 0\\n0 5.7 0\\n0 0 5.7\\n"
        "ATOMIC_POSITIONS (crystal)\\n"
        "Na 0.01 0.00 0.00\\n"
        "Cl 0.49 0.50 0.50\\n"
        "End final coordinates\\n"
        "JOB DONE.\\n"
    )

    def fake_run(*args, **kwargs):
        assert "calculation = 'vc-relax'" in kwargs["input"]
        assert "press_conv_thr = 2.5" in kwargs["input"]
        return subprocess.CompletedProcess(args[0], 0, output, "")

    monkeypatch.setattr("molecular_lab.calculations.execution.subprocess.run", fake_run)
    value = make_settings(
        task=DFTTask.VC_RELAX,
        pressure_gpa=1.0,
        convergence=ConvergenceSettings(
            force_tolerance_ev_per_angstrom=0.03,
            stress_tolerance_gpa=0.25,
        ),
        requested_outputs=(DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES, DFTOutput.STRESS),
    )
    result = execute_quantum_espresso(
        make_structure(), value, working_directory=tmp_path,
    )
    assert result.status is DFTExecutionStatus.SUCCEEDED
    assert result.optimized_structure.lattice_vectors_angstrom[0][0] == pytest.approx(5.7)
    assert result.stress_gpa == (
        (-1.0, 0.0, 0.0),
        (0.0, -2.0, 0.0),
        (0.0, 0.0, -3.0),
    )


def test_vc_relax_requires_explicit_stress_tolerance():
    value = make_settings(
        task=DFTTask.VC_RELAX,
        convergence=ConvergenceSettings(force_tolerance_ev_per_angstrom=0.03),
        requested_outputs=(DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES, DFTOutput.STRESS),
    )
    with pytest.raises(ValueError, match="stress_tolerance_gpa"):
        execute_dft_calculation("qe", make_structure(), value, working_directory=".")


def test_gpaw_vc_relax_uses_cell_filter_and_checks_stress(monkeypatch):
    from molecular_lab.calculations.execution import execute_gpaw

    class Cell:
        array = ((5.7, 0.0, 0.0), (0.0, 5.7, 0.0), (0.0, 0.0, 5.7))

    class FakeAtoms:
        cell = Cell()

        def get_forces(self):
            return ((0.001, 0.0, 0.0), (-0.001, 0.0, 0.0))

        def get_stress(self, voigt=False):
            assert voigt is False
            value = 0.01 / 160.21766208
            return ((value, 0.0, 0.0), (0.0, value, 0.0), (0.0, 0.0, value))

        def get_scaled_positions(self, wrap=False):
            assert wrap is False
            return ((0.01, 0.0, 0.0), (0.49, 0.5, 0.5))

        def get_potential_energy(self):
            return -5.25

    atoms = FakeAtoms()
    prepared = types.SimpleNamespace(atoms=atoms)
    fake_calculator = types.SimpleNamespace()
    fake_gpaw = types.ModuleType("gpaw")
    fake_gpaw.__version__ = "test-gpaw-vc"
    fake_ase_optimize = types.ModuleType("ase.optimize")

    class FakeBFGS:
        def __init__(self, target, logfile=None):
            self.target = target

        def run(self, fmax):
            assert fmax == pytest.approx(0.02)
            return True

    fake_ase_optimize.BFGS = FakeBFGS
    monkeypatch.setitem(sys.modules, "gpaw", fake_gpaw)
    monkeypatch.setitem(sys.modules, "ase.optimize", fake_ase_optimize)
    monkeypatch.setattr(gpaw_module, "prepare_gpaw_input", lambda structure, settings: prepared)
    monkeypatch.setattr(gpaw_module, "create_gpaw_calculator", lambda prepared, txt=None: fake_calculator)
    monkeypatch.setattr(gpaw_module, "build_optimization_target", lambda prepared: prepared.atoms)

    value = make_settings(
        task=DFTTask.VC_RELAX,
        convergence=ConvergenceSettings(
            force_tolerance_ev_per_angstrom=0.02,
            stress_tolerance_gpa=0.05,
        ),
        requested_outputs=(DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES, DFTOutput.STRESS),
    )
    result = execute_gpaw(make_structure(), value)
    assert result.status is DFTExecutionStatus.SUCCEEDED
    assert result.optimized_structure.lattice_vectors_angstrom[0][0] == pytest.approx(5.7)
    assert result.stress_gpa[0][0] == pytest.approx(0.01)
