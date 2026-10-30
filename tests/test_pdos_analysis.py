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
