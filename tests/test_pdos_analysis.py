import pytest

from molecular_lab.calculations.pdos_analysis import (
    aggregate_projected_dos, integrate_pdos_window,
)
from molecular_lab.materials.electronic_structure import (
    ProjectedDOSComponent, ProjectedDOSRecord, ProjectedDensityOfStatesData,
)


def _record(atom, element, l, spin, local, components, wfc=1):
    return ProjectedDOSRecord(
        atom_index=atom, element=element, wavefunction_index=wfc,
        angular_momentum=l, spin_channel=spin,
        local_dos_states_per_ev=tuple(local),
        components=tuple(
            ProjectedDOSComponent(i, f"m{i}", tuple(values))
            for i, values in enumerate(components)
        ),
    )


def test_element_aggregation_sums_atoms_and_wavefunctions_but_keeps_spin():
    data = ProjectedDensityOfStatesData(
        energies_ev=(-1.0, 0.0),
        records=(
            _record(1, "Fe", 0, "spin_up", (1, 2), ((1, 2),)),
            _record(2, "Fe", 0, "spin_up", (3, 4), ((3, 4),)),
            _record(1, "Fe", 1, "spin_down", (5, 6),
                    ((1, 1), (2, 2), (2, 3))),
            _record(1, "O", 0, "spin_up", (7, 8), ((7, 8),)),
        ),
    )
    result = aggregate_projected_dos(data, group_by="element")
    fe_s_up = next(x for x in result if x.element == "Fe" and x.angular_momentum == 0)
    fe_p_down = next(x for x in result if x.element == "Fe" and x.angular_momentum == 1)
    assert fe_s_up.local_dos_states_per_ev == (4, 6)
    assert fe_s_up.component_dos_states_per_ev == (4, 6)
    assert fe_p_down.component_dos_states_per_ev == (5, 6)
    assert len(result) == 3


def test_atom_aggregation_keeps_atoms_separate():
    data = ProjectedDensityOfStatesData(
        energies_ev=(-1.0, 0.0),
        records=(
            _record(1, "Fe", 0, "total", (1, 1), ((1, 1),)),
            _record(2, "Fe", 0, "total", (2, 2), ((2, 2),)),
        ),
    )
    result = aggregate_projected_dos(data, group_by="atom")
    assert [x.atom_index for x in result] == [1, 2]
    assert [x.local_dos_states_per_ev for x in result] == [(1, 1), (2, 2)]


def test_aggregation_rejects_unknown_grouping():
    data = ProjectedDensityOfStatesData(
        energies_ev=(-1.0, 0.0),
        records=(_record(1, "Fe", 0, "total", (1, 1), ((1, 1),)),),
    )
    with pytest.raises(ValueError, match="group_by"):
        aggregate_projected_dos(data, group_by="orbital")


def test_integrates_exact_window_with_interpolated_boundaries():
    # y=x+2 over [-1, 0, 1], integrated from -0.5 to 0.5 gives 2.0.
    value = integrate_pdos_window((-1, 0, 1), (1, 2, 3), -0.5, 0.5)
    assert value == pytest.approx(2.0)


def test_integral_rejects_window_outside_sampled_grid():
    with pytest.raises(ValueError, match="within"):
        integrate_pdos_window((0, 1), (1, 1), -0.1, 0.5)



def test_window_summary_integrates_element_orbital_and_spin_contributions():
    from molecular_lab.calculations.pdos_analysis import summarize_pdos_window

    data = ProjectedDensityOfStatesData(
        energies_ev=(-1.0, 0.0, 1.0),
        records=(
            _record(1, "Fe", 0, "spin_up", (1, 2, 3), ((1, 2, 3),)),
            _record(2, "Fe", 0, "spin_up", (2, 3, 4), ((2, 3, 4),)),
            _record(1, "Fe", 1, "spin_down", (1, 1, 1),
                    ((0.2, 0.2, 0.2), (0.3, 0.3, 0.3), (0.5, 0.5, 0.5))),
        ),
    )
    result = summarize_pdos_window(data, -0.5, 0.5)
    fe_s_up = next(x for x in result if x.element == "Fe" and x.angular_momentum == 0)
    fe_p_down = next(x for x in result if x.element == "Fe" and x.angular_momentum == 1)
    assert fe_s_up.integrated_local_states == pytest.approx(5.0)
    assert fe_s_up.integrated_component_states == pytest.approx(5.0)
    assert fe_p_down.integrated_component_states == pytest.approx(1.0)
    assert fe_s_up.spin_channel == "spin_up"
    assert fe_p_down.spin_channel == "spin_down"


def test_window_summary_rejects_energy_window_outside_grid():
    from molecular_lab.calculations.pdos_analysis import summarize_pdos_window

    data = ProjectedDensityOfStatesData(
        energies_ev=(-1.0, 1.0),
        records=(_record(1, "Fe", 0, "total", (1, 1), ((1, 1),)),),
    )
    with pytest.raises(ValueError, match="within"):
        summarize_pdos_window(data, -2.0, 0.0)


def test_band_edge_analysis_summarizes_orbitals_near_vbm_and_cbm():
    from molecular_lab.calculations.pdos_analysis import analyze_band_edge_pdos
    from molecular_lab.materials.electronic_structure import BandChannelData, BandStructureData, KPoint

    bands = BandStructureData(
        (KPoint((0, 0, 0), 1.0),),
        (BandChannelData('total', ((-1.0,), (1.0,)), ((2.0,), (0.0,))),),
    )
    energies = (-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
    curve = tuple(float(i) for i in range(len(energies)))
    zero = (0.0,) * len(energies)
    pdos = ProjectedDensityOfStatesData(
        energies,
        (_record(1, 'Fe', 2, 'total', curve, (curve, zero, zero, zero, zero)),),
    )
    result = analyze_band_edge_pdos(bands, pdos, window_ev=0.5)
    assert result.band_gap.indirect_gap_ev == pytest.approx(2.0)
    assert result.vbm_window[0].lower_ev == pytest.approx(-1.5)
    assert result.vbm_window[0].upper_ev == pytest.approx(-1.0)
    assert result.cbm_window[0].lower_ev == pytest.approx(1.0)
    assert result.cbm_window[0].upper_ev == pytest.approx(1.5)


def test_band_edge_analysis_applies_explicit_pdos_energy_shift():
    from molecular_lab.calculations.pdos_analysis import analyze_band_edge_pdos
    from molecular_lab.materials.electronic_structure import BandChannelData, BandStructureData, KPoint

    bands = BandStructureData(
        (KPoint((0, 0, 0), 1.0),),
        (BandChannelData('total', ((-1.0,), (1.0,)), ((2.0,), (0.0,))),),
    )
    energies = (-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0)
    zero = (0.0,) * len(energies)
    pdos = ProjectedDensityOfStatesData(
        energies, (_record(1, 'O', 1, 'total', (1.0,)*7, ((1.0,)*7, zero, zero)),),
    )
    result = analyze_band_edge_pdos(bands, pdos, window_ev=0.5, pdos_energy_shift_ev=1.0)
    assert result.pdos_energy_shift_ev == 1.0
    assert result.vbm_window[0].lower_ev == pytest.approx(-1.5)


def test_band_edge_analysis_rejects_metal_and_out_of_grid_window():
    from molecular_lab.calculations.pdos_analysis import analyze_band_edge_pdos
    from molecular_lab.materials.electronic_structure import BandChannelData, BandStructureData, KPoint

    metal = BandStructureData(
        (KPoint((0, 0, 0), 1.0),),
        (BandChannelData('total', ((0.0,),), ((1.0,),)),),
        fermi_level_ev=0.0,
    )
    pdos = ProjectedDensityOfStatesData(
        (-0.5, 0.0, 0.5), (_record(1, 'Fe', 0, 'total', (1, 1, 1), ((1, 1, 1),)),),
    )
    with pytest.raises(ValueError, match='resolved gap'):
        analyze_band_edge_pdos(metal, pdos)
