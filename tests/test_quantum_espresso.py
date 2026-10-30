import pytest

from molecular_lab.calculations import quantum_espresso as qe
from molecular_lab.calculations.dft_settings import (
    DFTCalculationSettings, DFTOutput, DFTTask, EnergyCutoff, EnergyUnit,
    KPointMesh, KPointScheme, PseudopotentialReference, SmearingScheme,
    SmearingSettings, ConvergenceSettings,
)
from molecular_lab.calculations.quantum_espresso import (
    QuantumEspressoInputError, prepare_pw_input,
)
from molecular_lab.structures.crystal import CrystalSite, CrystalStructure


class FakeCell:
    array = ((5.6, 0, 0), (0, 5.6, 0), (0, 0, 5.6))


class FakeAtoms:
    cell = FakeCell()

    def get_chemical_symbols(self):
        return ["Na", "Cl"]

    def get_masses(self):
        return [22.989769, 35.45]

    def get_scaled_positions(self, wrap=False):
        assert wrap is False
        return ((0, 0, 0), (0.5, 0.5, 0.5))


@pytest.fixture(autouse=True)
def fake_ase(monkeypatch):
    monkeypatch.setattr(qe, "to_ase_atoms", lambda structure: FakeAtoms())


def crystal(pbc=(True, True, True), moments=(None, None)):
    return CrystalStructure(
        structure_id="NaCl",
        lattice_vectors_angstrom=((5.6, 0, 0), (0, 5.6, 0), (0, 0, 5.6)),
        sites=(
            CrystalSite("Na1", "Na", (0, 0, 0), magnetic_moment_mu_b=moments[0]),
            CrystalSite("Cl1", "Cl", (0.5, 0.5, 0.5), magnetic_moment_mu_b=moments[1]),
        ),
        periodic_boundary_conditions=pbc,
    )


def pseudo(element):
    return PseudopotentialReference(
        element=element, file_name=f"{element}.UPF", family="PBE",
        source="test source", valence_electrons=1,
    )


def settings(**overrides):
    values = dict(
        task=DFTTask.SCF,
        exchange_correlation="PBE",
        energy_cutoff=EnergyCutoff(60, EnergyUnit.RYDBERG),
        kpoint_mesh=KPointMesh((6, 6, 6), KPointScheme.MONKHORST_PACK),
        pseudopotentials={"Na": pseudo("Na"), "Cl": pseudo("Cl")},
    )
    values.update(overrides)
    return DFTCalculationSettings(**values)


def test_scf_deck_contains_species_cell_positions_and_kmesh():
    result = prepare_pw_input(crystal(), settings())
    deck = result.input_text
    assert "calculation = 'scf'" in deck
    assert "ecutwfc = 60" in deck
    assert "ATOMIC_SPECIES" in deck
    assert "Na" in deck and "Na.UPF" in deck
    assert "ATOMIC_POSITIONS crystal" in deck
    assert "K_POINTS automatic\n  6 6 6 0 0 0" in deck
    assert result.required_pseudopotentials == ("Na.UPF", "Cl.UPF")
    assert any("not downloaded" in warning for warning in result.warnings)


def test_cutoff_and_smearing_are_converted_to_rydberg():
    value = settings(
        energy_cutoff=EnergyCutoff(27.211386245988, EnergyUnit.EV),
        smearing=SmearingSettings(SmearingScheme.COLD, EnergyCutoff(0.02, EnergyUnit.EV)),
    )
    deck = prepare_pw_input(crystal(), value).input_text
    assert "ecutwfc = 2" in deck
    assert "smearing = 'marzari-vanderbilt'" in deck
    assert "degauss = 0.001470997" in deck


def test_relax_maps_ionic_dynamics_and_force_threshold():
    value = settings(
        task=DFTTask.RELAX,
        convergence=ConvergenceSettings(force_tolerance_ev_per_angstrom=0.05),
    )
    deck = prepare_pw_input(crystal(), value).input_text
    assert "calculation = 'relax'" in deck
    assert "ion_dynamics = 'bfgs'" in deck
    assert "forc_conv_thr =" in deck
    assert "tprnfor = .true." in deck


def test_vc_relax_maps_pressure_gpa_to_kbar():
    value = settings(task=DFTTask.VC_RELAX, pressure_gpa=1.5)
    deck = prepare_pw_input(crystal(), value).input_text
    assert "calculation = 'vc-relax'" in deck
    assert "cell_dynamics = 'bfgs'" in deck
    assert "press = 15" in deck
    assert "tstress = .true." in deck


@pytest.mark.parametrize("task", [DFTTask.BANDS, DFTTask.DOS])
def test_bands_and_dos_are_rejected(task):
    value = settings(task=task, band_path_ref="path:test") if task is DFTTask.BANDS else settings(task=task)
    with pytest.raises(QuantumEspressoInputError, match="supports SCF"):
        prepare_pw_input(crystal(), value)


def test_soc_is_rejected_until_relativistic_pseudo_metadata_exists():
    with pytest.raises(QuantumEspressoInputError, match="fully relativistic"):
        prepare_pw_input(crystal(), settings(spin_orbit_coupling=True))


def test_site_magnetic_moments_are_not_silently_dropped():
    with pytest.raises(QuantumEspressoInputError, match="site-resolved"):
        prepare_pw_input(crystal(moments=(1.0, -1.0)), settings(spin_polarized=True))


def test_partial_pbc_and_gamma_mesh_are_rejected():
    with pytest.raises(QuantumEspressoInputError, match="Monkhorst-Pack"):
        prepare_pw_input(crystal(), settings(kpoint_mesh=KPointMesh((6, 6, 6))))
    with pytest.raises(QuantumEspressoInputError, match="three-dimensional"):
        prepare_pw_input(crystal((True, True, False)), settings(
            task=DFTTask.VC_RELAX,
            kpoint_mesh=KPointMesh((6, 6, 1), KPointScheme.MONKHORST_PACK),
        ))


def test_invalid_qe_tokens_and_paths_are_rejected():
    with pytest.raises(QuantumEspressoInputError, match="prefix"):
        prepare_pw_input(crystal(), settings(), prefix="bad name")
    with pytest.raises(QuantumEspressoInputError, match="outdir"):
        prepare_pw_input(crystal(), settings(), outdir="./tmp\n&CONTROL")


def test_unsupported_outputs_are_rejected():
    with pytest.raises(QuantumEspressoInputError, match="not implemented"):
        prepare_pw_input(crystal(), settings(
            requested_outputs=(DFTOutput.TOTAL_ENERGY, DFTOutput.BAND_STRUCTURE)
        ))


def test_spin_and_temperature_warnings_are_explicit():
    result = prepare_pw_input(crystal(), settings(
        spin_polarized=True, temperature_kelvin=300
    ))
    assert "nspin = 2" in result.input_text
    assert any("starting magnetization" in warning for warning in result.warnings)
    assert any("finite-temperature" in warning for warning in result.warnings)
