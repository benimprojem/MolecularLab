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
            BandPathPoint((0, 0, 0), "GAMMA", 1),
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
    # Existing workflow tests isolate the stages after prerequisite SCF.
    monkeypatch.setattr(execution, "_preflight_qe_bands", lambda *args, **kwargs: None)
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
    monkeypatch.setattr(
        execution, "prepare_qe_projected_dos_postprocess_input",
        lambda *args, **kwargs: SimpleNamespace(input_text="pdos post"),
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
        ), None

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



def test_qe_projected_dos_workflow_runs_projwfc_and_returns_atom_orbitals(tmp_path, monkeypatch):
    from dataclasses import replace
    patch_scf(monkeypatch)
    settings = replace(
        dos_settings(),
        requested_outputs=(DFTOutput.PROJECTED_DOS,),
    )

    def run(executable, input_text, *, workdir, **kwargs):
        if executable == "projwfc.x":
            (workdir / "molecular_lab.pdos.pdos_atm#1(Si)_wfc#1(p)").write_text(
                "-1 1.0 0.2 0.3 0.5\n0 1.2 0.4 0.5 0.6\n",
                encoding="utf-8",
            )
        return SimpleNamespace(returncode=0, stdout="JOB DONE.", stderr=""), None

    monkeypatch.setattr(execution, "_run_qe_input", run)
    result = execute_qe_dos(
        crystal(), settings, working_directory=tmp_path,
        projwfc_executable="projwfc.x",
    )
    assert result.status is DFTExecutionStatus.SUCCEEDED, result.message
    assert result.density_of_states is None
    assert result.projected_density_of_states.records[0].element == "Si"
    assert result.projected_density_of_states.records[0].angular_momentum == 1


def test_qe_fatband_workflow_runs_projwfc_and_joins_original_state_weights(tmp_path, monkeypatch):
    from dataclasses import replace

    patch_scf(monkeypatch)
    settings = replace(
        bands_settings(),
        requested_outputs=(DFTOutput.BAND_STRUCTURE, DFTOutput.FAT_BAND_STRUCTURE),
    )
    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="2" NUMBER_OF_K-POINTS="2"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES>
        <PROJS>
          <ATOMIC_WFC index="1">0.5 0.0 0.25 0.0</ATOMIC_WFC>
          <ATOMIC_WFC index="2">0.0 0.5 0.0 0.25</ATOMIC_WFC>
          <ATOMIC_WFC index="3">0.1 0.2 0.3 0.4</ATOMIC_WFC>
        </PROJS>
        <PROJS>
          <ATOMIC_WFC index="1">0.4 0.0 0.2 0.0</ATOMIC_WFC>
          <ATOMIC_WFC index="2">0.0 0.4 0.0 0.2</ATOMIC_WFC>
          <ATOMIC_WFC index="3">0.2 0.1 0.4 0.3</ATOMIC_WFC>
        </PROJS>
      </EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    projection_stdout = """
    Atomic states used for projection
    state # 1: atom 1 (Si), wfc 1 (l=0 m=0)
    state # 2: atom 1 (Si), wfc 2 (l=2 m=0)
    state # 3: atom 2 (O), wfc 1 (l=1 m=1)
    k = 0.0000000000 0.0000000000 0.0000000000
    ==== e(   1) = -1.00000 eV ====
    ==== e(   2) =  1.00000 eV ====
    k = 0.5000000000 0.0000000000 0.0000000000
    ==== e(   1) = -0.80000 eV ====
    ==== e(   2) =  1.20000 eV ====
    JOB DONE.
    """

    def run(executable, input_text, *, workdir, **kwargs):
        if executable == "bands.x":
            (workdir / "molecular_lab.bands.gnu").write_text(
                "0 -1\n0.5 -0.8\n\n0 1\n0.5 1.2\n",
                encoding="utf-8",
            )
            stdout = "JOB DONE."
        elif executable == "projwfc.x":
            save_dir = workdir / "tmp" / "molecular_lab.save"
            save_dir.mkdir(parents=True, exist_ok=True)
            (save_dir / "atomic_proj.xml").write_text(xml, encoding="utf-8")
            stdout = projection_stdout
        else:
            stdout = "JOB DONE."
        return SimpleNamespace(returncode=0, stdout=stdout, stderr=""), None

    monkeypatch.setattr(execution, "_run_qe_input", run)
    result = execute_qe_bands(
        crystal(), settings, working_directory=tmp_path,
        projwfc_executable="projwfc.x",
    )
    assert result.status is DFTExecutionStatus.SUCCEEDED, result.message
    assert result.band_structure is not None
    assert result.fat_band_structure is not None
    assert len(result.fat_band_structure.points) == 4
    first = next(
        point for point in result.fat_band_structure.points
        if point.kpoint_index == 0 and point.band_index == 1
    )
    assert first.energy_ev == pytest.approx(-1.0)
    si_s = next(
        item for item in first.contributions
        if item.element == "Si" and item.angular_momentum == 0
    )
    assert si_s.raw_weight == pytest.approx(0.25)


def test_qe_fatband_workflow_fails_if_atomic_projection_xml_is_missing(tmp_path, monkeypatch):
    from dataclasses import replace

    patch_scf(monkeypatch)
    settings = replace(
        bands_settings(),
        requested_outputs=(DFTOutput.BAND_STRUCTURE, DFTOutput.FAT_BAND_STRUCTURE),
    )
    stale_dir = tmp_path / "tmp" / "molecular_lab.save"
    stale_dir.mkdir(parents=True)
    (stale_dir / "atomic_proj.xml").write_text("<stale/>", encoding="utf-8")

    def run(executable, input_text, *, workdir, **kwargs):
        if executable == "bands.x":
            (workdir / "molecular_lab.bands.gnu").write_text(
                "0 -1\n0.5 -0.8\n\n0 1\n0.5 1.2\n",
                encoding="utf-8",
            )
        return SimpleNamespace(returncode=0, stdout="JOB DONE.", stderr=""), None

    monkeypatch.setattr(execution, "_run_qe_input", run)
    result = execute_qe_bands(crystal(), settings, working_directory=tmp_path)
    assert result.status is DFTExecutionStatus.FAILED
    assert "atomic_proj.xml" in result.message



def test_qe_fatband_preflight_detects_missing_projwfc_before_scf(tmp_path, monkeypatch):
    from dataclasses import replace

    settings = replace(
        bands_settings(),
        requested_outputs=(DFTOutput.BAND_STRUCTURE, DFTOutput.FAT_BAND_STRUCTURE),
    )
    (tmp_path / "pseudo").mkdir()
    (tmp_path / "pseudo" / "Si.UPF").write_text("test pseudo", encoding="utf-8")
    monkeypatch.setattr(execution, "ase_availability", lambda: (True, ""))
    monkeypatch.setattr(
        execution, "_resolve_qe_executable",
        lambda name: None if name == "projwfc-missing" else f"/fake/{name}",
    )
    scf_called = []
    monkeypatch.setattr(
        execution, "execute_quantum_espresso",
        lambda *args, **kwargs: scf_called.append(True) or scf_success(),
    )

    result = execute_qe_bands(
        crystal(), settings, working_directory=tmp_path,
        projwfc_executable="projwfc-missing",
    )

    assert result.status is DFTExecutionStatus.BACKEND_UNAVAILABLE
    assert "projwfc.x executable was not found" in result.message
    assert "no calculation was started" in result.message
    assert scf_called == []


def test_qe_bands_preflight_reports_missing_pseudopotential(tmp_path, monkeypatch):
    monkeypatch.setattr(execution, "ase_availability", lambda: (True, ""))
    monkeypatch.setattr(
        execution, "_resolve_qe_executable",
        lambda name: f"/fake/{name}",
    )
    (tmp_path / "pseudo").mkdir()

    result = execution._preflight_qe_bands(
        crystal(), bands_settings(),
        working_directory=tmp_path,
        executable="pw.x",
        bands_executable="bands.x",
        projwfc_executable="projwfc.x",
        wants_fatband=False,
    )

    assert result is not None
    assert result.status is DFTExecutionStatus.BACKEND_UNAVAILABLE
    assert "Missing Si pseudopotential" in result.message


def test_qe_bands_preflight_does_not_require_projwfc_for_plain_bands(tmp_path, monkeypatch):
    monkeypatch.setattr(execution, "ase_availability", lambda: (True, ""))
    monkeypatch.setattr(
        execution, "_resolve_qe_executable",
        lambda name: f"/fake/{name}",
    )
    (tmp_path / "pseudo").mkdir()
    (tmp_path / "pseudo" / "Si.UPF").write_text("test pseudo", encoding="utf-8")

    result = execution._preflight_qe_bands(
        crystal(), bands_settings(),
        working_directory=tmp_path,
        executable="pw.x",
        bands_executable="bands.x",
        projwfc_executable="projwfc-missing",
        wants_fatband=False,
    )

    assert result is None
