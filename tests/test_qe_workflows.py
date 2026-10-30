from types import SimpleNamespace

import pytest

from molecular_lab.calculations import execution
from molecular_lab.calculations.dft_settings import (
    BandPath, BandPathPoint, DOSSettings, DFTCalculationSettings,
    DFTOutput, DFTTask, EnergyCutoff, EnergyUnit, KPointMesh,
    KPointScheme, PseudopotentialReference,
)
from molecular_lab.calculations.execution import (
    DFTExecutionResult, DFTExecutionStatus, execute_qe_bands, execute_qe_dos,
)
from molecular_lab.structures.crystal import CrystalSite, CrystalStructure


def crystal():
    return CrystalStructure(
        structure_id="Si",
        lattice_vectors_angstrom=((5.4, 0, 0), (0, 5.4, 0), (0, 0, 5.4)),
        sites=(CrystalSite("Si1", "Si", (0, 0, 0)),),
    )


def pseudo():
    return PseudopotentialReference(
        element="Si", file_name="Si.UPF", family="PBE",
        source="test", valence_electrons=4,
    )


def bands_settings():
    return DFTCalculationSettings(
        task=DFTTask.BANDS, exchange_correlation="PBE",
        energy_cutoff=EnergyCutoff(50, EnergyUnit.RYDBERG),
        kpoint_mesh=KPointMesh((6, 6, 6), KPointScheme.MONKHORST_PACK),
        pseudopotentials={"Si": pseudo()},
        band_path=BandPath((
            BandPathPoint((0, 0, 0), "GAMMA", 10),
            BandPathPoint((0.5, 0, 0), "X", 0),
        )),
        requested_outputs=(DFTOutput.BAND_STRUCTURE,),
    )


def dos_settings():
    return DFTCalculationSettings(
        task=DFTTask.DOS, exchange_correlation="PBE",
        energy_cutoff=EnergyCutoff(50, EnergyUnit.RYDBERG),
        kpoint_mesh=KPointMesh((8, 8, 8), KPointScheme.MONKHORST_PACK),
        pseudopotentials={"Si": pseudo()},
        dos=DOSSettings(-5, 5, points=101),
        requested_outputs=(DFTOutput.DENSITY_OF_STATES,),
    )


def scf_success():
    return DFTExecutionResult(
        "quantum_espresso", DFTExecutionStatus.SUCCEEDED,
        "SCF done", total_energy_ev=-10.0, converged=True,
        backend_version="7.5", stdout="the Fermi energy is 0.25 ev",
    )


def patch_scf(monkeypatch, result=None):
    monkeypatch.setattr(
        execution, "execute_quantum_espresso",
        lambda *args, **kwargs: result or scf_success(),
    )
    monkeypatch.setattr(
        execution, "prepare_pw_input",
        lambda *args, **kwargs: SimpleNamespace(input_text="pw stage"),
    )
    monkeypatch.setattr(
        execution, "prepare_qe_bands_postprocess_input",
        lambda *args, **kwargs: SimpleNamespace(input_text="bands post"),
    )
    monkeypatch.setattr(
        execution, "prepare_qe_dos_postprocess_input",
        lambda *args, **kwargs: SimpleNamespace(input_text="dos post"),
    )


def test_qe_bands_workflow_runs_prerequisite_and_parses_output(tmp_path, monkeypatch):
    patch_scf(monkeypatch)

    def run(executable, input_text, *, workdir, **kwargs):
        if executable == "bands.x":
            (workdir / "molecular_lab.bands.gnu").write_text(
                "0 -1\n0.5 -0.5\n\n0 1\n0.5 2\n",
                encoding="utf-8",
            )
        return SimpleNamespace(
            returncode=0, stdout="JOB DONE.", stderr="",
        ), None

    monkeypatch.setattr(execution, "_run_qe_input", run)
    result = execute_qe_bands(crystal(), bands_settings(), working_directory=tmp_path)
    assert result.status is DFTExecutionStatus.SUCCEEDED
    assert result.total_energy_ev == pytest.approx(-10)
    assert result.band_structure.channels[0].band_count == 2
    assert result.band_structure.fermi_level_ev == pytest.approx(0.25)


def test_qe_dos_workflow_runs_scf_nscf_and_parses_output(tmp_path, monkeypatch):
    patch_scf(monkeypatch)

    def run(executable, input_text, *, workdir, **kwargs):
        if executable == "dos.x":
            (workdir / "molecular_lab.dos").write_text(
                "-1 0.2 0.0\n0 0.5 0.3\n1 0.1 0.6\n",
                encoding="utf-8",
            )
        return SimpleNamespace(
            returncode=0, stdout="JOB DONE.", stderr="",
        )

    monkeypatch.setattr(execution, "_run_qe_input", run)
    result = execute_qe_dos(crystal(), dos_settings(), working_directory=tmp_path)
    assert result.status is DFTExecutionStatus.SUCCEEDED
    assert result.total_energy_ev == pytest.approx(-10)
    assert result.density_of_states.energies_ev == (-1, 0, 1)
    assert result.density_of_states.fermi_level_ev == pytest.approx(0.25)


def test_qe_workflow_stops_when_scf_fails(tmp_path, monkeypatch):
    failed = DFTExecutionResult(
        "quantum_espresso", DFTExecutionStatus.NOT_CONVERGED,
        "SCF did not converge",
    )
    patch_scf(monkeypatch, failed)
    called = []
    monkeypatch.setattr(
        execution, "_run_qe_input",
        lambda *args, **kwargs: called.append(True),
    )
    result = execute_qe_bands(crystal(), bands_settings(), working_directory=tmp_path)
    assert result.status is DFTExecutionStatus.NOT_CONVERGED
    assert "prerequisite SCF" in result.message
    assert called == []


def test_spin_polarized_bands_is_explicitly_rejected(tmp_path, monkeypatch):
    patch_scf(monkeypatch)
    value = bands_settings()
    from dataclasses import replace
    value = replace(value, spin_polarized=True)
    result = execute_qe_bands(crystal(), value, working_directory=tmp_path)
    assert result.status is DFTExecutionStatus.BACKEND_UNAVAILABLE
