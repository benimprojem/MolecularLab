import pytest

from molecular_lab.calculations.qe_output import (
    QEOutputParseError,
    parse_qe_bands_gnu,
    parse_qe_fermi_level,
    parse_qe_total_dos,
)


def test_parse_qe_bands_gnu_converts_to_common_band_model():
    source = """0.0 -2.0
0.5 -1.0
1.0 -2.0

0.0 1.0
0.5 2.0
1.0 1.0
"""
    result = parse_qe_bands_gnu(
        source, fermi_level_ev=0.0, backend_version="7.5",
        exchange_correlation="PBE",
    )
    assert result.backend == "quantum_espresso"
    assert result.backend_version == "7.5"
    assert result.fermi_level_ev == pytest.approx(0.0)
    assert result.kpoints[1].coordinates == (0.5, 0.0, 0.0)
    assert result.kpoints[1].coordinate_system == "path_distance"
    assert result.channels[0].energies_ev == ((-2.0, -1.0, -2.0), (1.0, 2.0, 1.0))


def test_bands_parser_rejects_mismatched_path_grids():
    source = "0.0 -1\n0.5 -0.5\n\n0.1 1\n0.5 2\n"
    with pytest.raises(QEOutputParseError, match="same k-path grid"):
        parse_qe_bands_gnu(source)


def test_parse_qe_total_dos_unpolarized_and_ignore_integrated_column():
    source = "-2.0 0.1 0.0\n-1.0 0.4 0.2\n0.0 0.8 0.7\n"
    result = parse_qe_total_dos(source)
    assert result.energies_ev == (-2.0, -1.0, 0.0)
    assert result.channels[0].spin_channel == "total"
    assert result.channels[0].density_states_per_ev == (0.1, 0.4, 0.8)


def test_parse_qe_total_dos_spin_channels():
    source = "-1 0.2 0.3\n0 0.4 0.5\n1 0.1 0.2\n"
    result = parse_qe_total_dos(source, spin_polarized=True)
    assert [channel.spin_channel for channel in result.channels] == ["spin_up", "spin_down"]
    assert result.channels[1].density_states_per_ev == (0.3, 0.5, 0.2)


def test_dos_parser_rejects_nonmonotonic_energy_grid():
    with pytest.raises(QEOutputParseError, match="strictly increasing"):
        parse_qe_total_dos("0 1\n-1 2\n")


def test_fermi_energy_parser_only_uses_explicit_value():
    assert parse_qe_fermi_level("the Fermi energy is 5.2500 ev") == pytest.approx(5.25)
    assert parse_qe_fermi_level("highest occupied, lowest unoccupied level: 2.0 3.0 ev") is None
