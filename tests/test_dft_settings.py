import pytest

from molecular_lab.calculations.dft_settings import (
    BandPath, BandPathPoint, DOSSettings, ConvergenceSettings,
    DFTCalculationSettings, DFTOutput, DFTTask,
    EnergyCutoff, EnergyUnit, KPointMesh, KPointScheme, PseudopotentialReference,
    SmearingScheme, SmearingSettings,
)
from molecular_lab.structures.crystal import CrystalSite, CrystalStructure


def crystal(pbc=(True, True, True)):
    return CrystalStructure(
        structure_id="NaCl",
        lattice_vectors_angstrom=((5.6, 0, 0), (0, 5.6, 0), (0, 0, 5.6)),
        sites=(
            CrystalSite("Na1", "Na", (0, 0, 0)),
            CrystalSite("Cl1", "Cl", (0.5, 0.5, 0.5)),
        ),
        periodic_boundary_conditions=pbc,
    )


def pseudo(element):
    return PseudopotentialReference(
        element=element, file_name=f"{element}.UPF", family="test-family",
        source="test fixture", valence_electrons=1,
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


def test_energy_cutoff_unit_conversion_round_trip():
    cutoff = EnergyCutoff(1, EnergyUnit.HARTREE)
    assert cutoff.value_ev == pytest.approx(27.211386245988)
    assert cutoff.to(EnergyUnit.RYDBERG).value == pytest.approx(2.0)
    assert cutoff.to(EnergyUnit.EV).to(EnergyUnit.HARTREE).value == pytest.approx(1.0)


@pytest.mark.parametrize("value,unit", [(0, EnergyUnit.EV), (-1, EnergyUnit.EV), (float("inf"), EnergyUnit.EV)])
def test_invalid_energy_cutoff(value, unit):
    with pytest.raises(ValueError):
        EnergyCutoff(value, unit)


def test_kpoint_mesh_requires_three_positive_integers():
    with pytest.raises(ValueError):
        KPointMesh((2, 0, 2))
    with pytest.raises(ValueError):
        KPointMesh((2, 2, 2, 2))


def test_nonperiodic_axis_requires_single_kpoint():
    slab = crystal((True, True, False))
    settings(kpoint_mesh=KPointMesh((6, 6, 1))).validate_for_structure(slab)
    with pytest.raises(ValueError, match="non-periodic axis 2"):
        settings().validate_for_structure(slab)


def test_structure_validation_requires_exact_pseudopotential_coverage():
    with pytest.raises(ValueError, match="missing pseudopotential"):
        settings(pseudopotentials={"Na": pseudo("Na")}).validate_for_structure(crystal())
    with pytest.raises(ValueError, match="absent from structure"):
        settings(pseudopotentials={
            "Na": pseudo("Na"), "Cl": pseudo("Cl"), "K": pseudo("K")
        }).validate_for_structure(crystal())


def test_settings_mapping_is_immutable():
    value = settings()
    with pytest.raises(TypeError):
        value.pseudopotentials["K"] = pseudo("K")


def test_smearing_width_and_scheme_are_explicit():
    value = settings(smearing=SmearingSettings(
        SmearingScheme.GAUSSIAN, EnergyCutoff(0.02, EnergyUnit.EV)
    ))
    assert value.smearing.width.value_ev == pytest.approx(0.02)
    assert value.smearing.scheme is SmearingScheme.GAUSSIAN


def test_convergence_settings_validate_positive_values():
    assert ConvergenceSettings(max_scf_iterations=250).max_scf_iterations == 250
    assert ConvergenceSettings(stress_tolerance_gpa=0.2).stress_tolerance_gpa == 0.2
    with pytest.raises(ValueError):
        ConvergenceSettings(energy_tolerance_ev=0)
    with pytest.raises(ValueError):
        ConvergenceSettings(stress_tolerance_gpa=0)
    with pytest.raises(ValueError):
        ConvergenceSettings(max_scf_iterations=0)


def test_band_task_requires_explicit_path_reference():
    with pytest.raises(ValueError, match="band_path_ref"):
        settings(task=DFTTask.BANDS)
    value = settings(task=DFTTask.BANDS, band_path_ref="path:NaCl-standard")
    assert value.band_path_ref == "path:NaCl-standard"


def test_band_path_reference_is_rejected_for_non_band_task():
    with pytest.raises(ValueError, match="only valid"):
        settings(band_path_ref="path:unexpected")


def test_outputs_are_unique_and_normalized():
    value = settings(requested_outputs=("total_energy", "forces"))
    assert value.requested_outputs == (DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES)
    with pytest.raises(ValueError):
        settings(requested_outputs=(DFTOutput.TOTAL_ENERGY, DFTOutput.TOTAL_ENERGY))


def test_spin_orbit_is_not_forced_to_mean_magnetic_spin_polarization():
    value = settings(spin_orbit_coupling=True, spin_polarized=False)
    assert value.spin_orbit_coupling is True



def test_structured_band_path_and_dos_settings_validation():
    path = BandPath((
        BandPathPoint((0.0, 0.0, 0.0), "GAMMA", 20),
        BandPathPoint((0.5, 0.0, 0.0), "X", 0),
    ))
    assert path.sample_count == 21
    with pytest.raises(ValueError):
        BandPath((BandPathPoint((0, 0, 0), "GAMMA"),))
    with pytest.raises(ValueError):
        DOSSettings(5.0, -5.0)
    with pytest.raises(ValueError):
        DOSSettings(-5.0, 5.0, points=1)
    with pytest.raises(ValueError):
        DOSSettings(-5.0, 5.0, broadening_ev=0)
