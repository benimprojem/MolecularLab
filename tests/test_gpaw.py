import pytest

from molecular_lab.calculations import gpaw as gpaw_module
from molecular_lab.calculations.dft_settings import (
    ConvergenceSettings, DFTCalculationSettings, DFTOutput, DFTTask,
    EnergyCutoff, EnergyUnit, KPointMesh, KPointScheme,
    PseudopotentialReference, SmearingScheme, SmearingSettings,
)
from molecular_lab.calculations.gpaw import (
    GPAWInputError, GPAWUnavailableError, PreparedGPAWInput,
    build_optimization_target, create_gpaw_calculator, prepare_gpaw_input,
)
from molecular_lab.structures.crystal import CrystalSite, CrystalStructure


class FakeAtoms:
    def __init__(self):
        self.calc = None


@pytest.fixture(autouse=True)
def fake_ase_bridge(monkeypatch):
    monkeypatch.setattr(gpaw_module, "to_ase_atoms", lambda structure: FakeAtoms())


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


def settings(**overrides):
    values = dict(
        task=DFTTask.SCF,
        exchange_correlation="PBE",
        energy_cutoff=EnergyCutoff(500, EnergyUnit.EV),
        kpoint_mesh=KPointMesh((6, 6, 6), KPointScheme.GAMMA_CENTERED),
        # GPAW uses its local PAW setup library; these QE-style references are ignored.
        pseudopotentials={},
    )
    values.update(overrides)
    return DFTCalculationSettings(**values)


def test_scf_prepares_gpaw_parameters_without_external_pseudos():
    result = prepare_gpaw_input(crystal(), settings())
    assert result.task is DFTTask.SCF
    assert result.plane_wave_cutoff_ev == 500
    assert result.calculator_parameters["xc"] == "PBE"
    assert result.calculator_parameters["kpts"] == {
        "size": (6, 6, 6), "gamma": True,
    }
    assert result.calculator_parameters["maxiter"] == 100
    assert any("not used" in warning for warning in result.warnings)


def test_monkhorst_pack_and_nonzero_charge_are_mapped():
    result = prepare_gpaw_input(crystal(), settings(
        kpoint_mesh=KPointMesh((4, 4, 4), KPointScheme.MONKHORST_PACK),
    ))
    assert result.calculator_parameters["kpts"]["gamma"] is False


def test_energy_cutoff_converts_to_ev():
    result = prepare_gpaw_input(crystal(), settings(
        energy_cutoff=EnergyCutoff(2, EnergyUnit.HARTREE),
    ))
    assert result.plane_wave_cutoff_ev == pytest.approx(54.422772491976)


def test_smearing_and_convergence_are_mapped():
    result = prepare_gpaw_input(crystal(), settings(
        smearing=SmearingSettings(
            SmearingScheme.FERMI_DIRAC,
            EnergyCutoff(0.1, EnergyUnit.EV),
        ),
        convergence=ConvergenceSettings(
            energy_tolerance_ev=1e-5,
            force_tolerance_ev_per_angstrom=0.03,
            max_scf_iterations=250,
            density_tolerance=1e-7,
        ),
    ))
    assert result.calculator_parameters["occupations"] == {
        "name": "fermi-dirac", "width": 0.1,
    }
    assert result.calculator_parameters["convergence"] == {
        "energy": 1e-5, "density": 1e-7, "forces": 0.03,
    }
    assert result.optimizer_fmax_ev_per_angstrom == 0.03


def test_relaxation_target_is_atoms_and_fmax_is_separate():
    result = prepare_gpaw_input(crystal(), settings(
        task=DFTTask.RELAX,
        convergence=ConvergenceSettings(force_tolerance_ev_per_angstrom=0.04),
    ))
    assert result.optimizer_fmax_ev_per_angstrom == 0.04
    assert build_optimization_target(result) is result.atoms


def test_vc_relax_pressure_converts_to_ev_per_angstrom3(monkeypatch):
    class Filter:
        def __init__(self, atoms, **kwargs):
            self.atoms = atoms
            self.kwargs = kwargs

    import sys
    import types
    fake_filters = types.ModuleType("ase.filters")
    fake_filters.FrechetCellFilter = Filter
    fake_ase = types.ModuleType("ase")
    fake_ase.filters = fake_filters
    monkeypatch.setitem(sys.modules, "ase", fake_ase)
    monkeypatch.setitem(sys.modules, "ase.filters", fake_filters)

    result = prepare_gpaw_input(crystal(), settings(
        task=DFTTask.VC_RELAX, pressure_gpa=1.0,
    ))
    assert result.cell_filter_parameters["scalar_pressure"] == pytest.approx(
        0.006241509074
    )
    target = build_optimization_target(result)
    assert target.kwargs["scalar_pressure"] == pytest.approx(0.006241509074)


@pytest.mark.parametrize("task", [DFTTask.BANDS, DFTTask.DOS])
def test_bands_and_dos_are_rejected(task):
    value = settings(task=task, band_path_ref="path:test") if task is DFTTask.BANDS else settings(task=task)
    with pytest.raises(GPAWInputError, match="supports SCF"):
        prepare_gpaw_input(crystal(), value)


def test_soc_partial_pbc_and_site_moments_are_rejected():
    with pytest.raises(GPAWInputError, match="SOC"):
        prepare_gpaw_input(crystal(), settings(spin_orbit_coupling=True))
    with pytest.raises(GPAWInputError, match="3D PBC"):
        prepare_gpaw_input(crystal((True, True, False)), settings())
    with pytest.raises(GPAWInputError, match="site-resolved"):
        prepare_gpaw_input(crystal(moments=(1.0, -1.0)), settings(spin_polarized=True))


def test_unsupported_smearing_shift_and_outputs_are_rejected():
    with pytest.raises(GPAWInputError, match="Fermi-Dirac"):
        prepare_gpaw_input(crystal(), settings(
            smearing=SmearingSettings(
                SmearingScheme.GAUSSIAN, EnergyCutoff(0.1, EnergyUnit.EV)
            )
        ))
    with pytest.raises(GPAWInputError, match="zero explicit"):
        prepare_gpaw_input(crystal(), settings(
            kpoint_mesh=KPointMesh(
                (4, 4, 4), KPointScheme.GAMMA_CENTERED, (0.5, 0, 0)
            )
        ))
    with pytest.raises(GPAWInputError, match="not implemented"):
        prepare_gpaw_input(crystal(), settings(
            requested_outputs=(DFTOutput.TOTAL_ENERGY, DFTOutput.BAND_STRUCTURE)
        ))


def test_calculator_factory_reports_missing_gpaw(monkeypatch):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "gpaw":
            raise ImportError("not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    prepared = PreparedGPAWInput(
        atoms=FakeAtoms(), calculator_parameters={}, plane_wave_cutoff_ev=300,
        task=DFTTask.SCF, optimizer_fmax_ev_per_angstrom=None,
        cell_filter_parameters={},
    )
    with pytest.raises(GPAWUnavailableError, match="gpaw"):
        create_gpaw_calculator(prepared)


def test_scf_cannot_build_optimization_target():
    prepared = prepare_gpaw_input(crystal(), settings())
    with pytest.raises(GPAWInputError, match="SCF"):
        build_optimization_target(prepared)
