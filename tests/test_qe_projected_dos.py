import pytest

from molecular_lab.calculations.qe_output import (
    QEOutputParseError, parse_qe_projected_dos,
)


def test_parse_unpolarized_atomic_pdos_keeps_atom_wavefunction_and_m_components():
    files = {
        "sample.pdos_atm#2(Si)_wfc#1(s)": "-2 0.4 0.4\n-1 0.8 0.8\n",
        "sample.pdos_atm#2(Si)_wfc#2(p)": (
            "-2 1.0 0.2 0.3 0.5\n"
            "-1 1.5 0.4 0.5 0.6\n"
        ),
    }
    result = parse_qe_projected_dos(files, fermi_level_ev=0.1, backend_version="7.5")
    assert result.energies_ev == (-2.0, -1.0)
    assert result.fermi_level_ev == pytest.approx(0.1)
    assert len(result.records) == 2
    p_record = result.records[1]
    assert (p_record.atom_index, p_record.element, p_record.wavefunction_index) == (2, "Si", 2)
    assert p_record.angular_momentum == 1
    assert [c.component_label for c in p_record.components] == [
        "p_z", "p_x", "p_y"
    ]
    assert p_record.components[2].density_states_per_ev == (0.5, 0.6)


def test_parse_spin_polarized_pdos_keeps_up_and_down_separate():
    files = {
        "sample.pdos_atm#1(Fe)_wfc#3(d)": (
            "-1 1.0 0.8 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0\n"
            "0 1.1 0.9 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1\n"
        ),
    }
    result = parse_qe_projected_dos(files, spin_polarized=True)
    assert len(result.records) == 2
    assert result.records[0].spin_channel == "spin_up"
    assert result.records[1].spin_channel == "spin_down"
    assert result.records[0].components[0].density_states_per_ev == (0.1, 0.2)
    assert result.records[1].components[0].density_states_per_ev == (0.2, 0.3)


def test_pdos_parser_rejects_inconsistent_energy_grids():
    files = {
        "sample.pdos_atm#1(Fe)_wfc#1(s)": "-1 0.1 0.1\n0 0.2 0.2\n",
        "sample.pdos_atm#2(Fe)_wfc#1(s)": "-1 0.1 0.1\n1 0.2 0.2\n",
    }
    with pytest.raises(QEOutputParseError, match="common energy grid"):
        parse_qe_projected_dos(files)


def test_pdos_parser_rejects_spin_column_mismatch():
    files = {"sample.pdos_atm#1(Fe)_wfc#1(s)": "-1 0.1 0.1\n0 0.2 0.2\n"}
    with pytest.raises(QEOutputParseError, match="expected at least 4"):
        parse_qe_projected_dos(files, spin_polarized=True)


def test_pdos_parser_rejects_j_resolved_filename_instead_of_mislabeling():
    files = {"sample.pdos_atm#1(Pt)_wfc#1(p_j)": "-1 0.1 0.1\n0 0.2 0.2\n"}
    with pytest.raises(QEOutputParseError, match="Unsupported QE PDOS filename"):
        parse_qe_projected_dos(files)
