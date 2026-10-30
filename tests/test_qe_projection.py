import pytest

from molecular_lab.calculations.fatband import build_fatband_dataset
from molecular_lab.calculations.qe_projection import (
    QEProjectionParseError,
    parse_qe_atomic_proj_xml,
)
from molecular_lab.materials.electronic_structure import (
    BandChannelData, BandStructureData, KPoint,
)


STATE_OUTPUT = """
state # 1: atom 1 (Fe), wfc 1 (l=0 m=0)
state # 2: atom 1 (Fe), wfc 2 (l=2 m=0)
state # 3: atom 2 (O), wfc 1 (l=1 m=1)
"""


def test_parse_new_qe_eigenstates_projection_layout():
    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="2" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES><PROJS>
        <ATOMIC_WFC index="1">0.5 0.0 0.25 0.0</ATOMIC_WFC>
        <ATOMIC_WFC index="2">0.0 0.5 0.0 0.25</ATOMIC_WFC>
        <ATOMIC_WFC index="3">0.1 0.2 0.3 0.4</ATOMIC_WFC>
      </PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    data = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    assert data.format_generation == "eigenstates_projs"
    assert data.number_of_bands == 2
    assert data.number_of_kpoints == 1
    assert data.spin_channels == ("total",)
    item = next(x for x in data.coefficients
                if x.atomic_wfc_index == 3 and x.band_index == 2)
    assert item.real == pytest.approx(0.3)
    assert item.imaginary == pytest.approx(0.4)
    assert item.weight == pytest.approx(0.25)


def test_parse_legacy_qe_projection_layout_with_spin_channels():
    xml = """<ATOMIC_PROJECTIONS>
      <HEADER>
        <NUMBER_OF_BANDS>1</NUMBER_OF_BANDS>
        <NUMBER_OF_K-POINTS>1</NUMBER_OF_K-POINTS>
        <NUMBER_OF_SPIN_COMPONENTS>2</NUMBER_OF_SPIN_COMPONENTS>
        <NUMBER_OF_ATOMIC_WFC>3</NUMBER_OF_ATOMIC_WFC>
      </HEADER>
      <PROJECTIONS><K-POINT.1>
        <SPIN.1>
          <ATMWFC.1>0.5, 0.0</ATMWFC.1>
          <ATMWFC.2>0.0, 0.5</ATMWFC.2>
          <ATMWFC.3>0.1, 0.2</ATMWFC.3>
        </SPIN.1>
        <SPIN.2>
          <ATMWFC.1>0.25, 0.0</ATMWFC.1>
          <ATMWFC.2>0.0, 0.25</ATMWFC.2>
          <ATMWFC.3>0.2, 0.1</ATMWFC.3>
        </SPIN.2>
      </K-POINT.1></PROJECTIONS>
    </ATOMIC_PROJECTIONS>"""
    data = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    assert data.format_generation == "legacy_projections"
    assert data.spin_channels == ("spin_up", "spin_down")
    assert len(data.coefficients) == 6
    down = next(x for x in data.coefficients if x.spin_channel == "spin_down")
    assert down.weight == pytest.approx(0.25 ** 2)


def test_parser_rejects_spinor_data_until_supported():
    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="4" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES/>
    </ATOMIC_PROJECTIONS>"""
    with pytest.raises(QEProjectionParseError, match="spinor"):
        parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)


def test_parser_rejects_missing_state_metadata():
    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="1"/>
      <EIGENSTATES><PROJS><ATOMIC_WFC index="1">1 0</ATOMIC_WFC></PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    with pytest.raises(QEProjectionParseError, match="No state"):
        parse_qe_atomic_proj_xml(xml, "projwfc output without state records")


def test_fatband_builder_joins_weights_to_band_energy_and_orbitals():
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
    projections = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    bands = BandStructureData(
        kpoints=(
            KPoint((0.0, 0.0, 0.0), 1.0),
            KPoint((0.5, 0.0, 0.0), 1.0),
        ),
        channels=(BandChannelData("total", ((-1.0, -0.8), (1.0, 1.2))),),
        fermi_level_ev=0.0,
    )
    result = build_fatband_dataset(bands, projections)
    assert len(result.points) == 4
    point = next(x for x in result.points if x.kpoint_index == 0 and x.band_index == 1)
    fe_s = next(x for x in point.contributions
                if x.element == "Fe" and x.angular_momentum == 0)
    assert fe_s.raw_weight == pytest.approx(0.25)
    assert point.energy_ev == pytest.approx(-1.0)
    assert point.total_projected_weight >= fe_s.raw_weight
    assert "not normalized" in " ".join(result.notes)


def test_fatband_builder_rejects_kpoint_and_spin_mismatch():
    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES><PROJS>
        <ATOMIC_WFC index="1">1 0</ATOMIC_WFC>
        <ATOMIC_WFC index="2">0 0</ATOMIC_WFC>
        <ATOMIC_WFC index="3">0 0</ATOMIC_WFC>
      </PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    projections = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    bands = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0),),
        channels=(BandChannelData("total", ((0.0,),)),),
    )
    assert len(build_fatband_dataset(bands, projections).points) == 1
    wrong = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0),),
        channels=(BandChannelData("spin_up", ((0.0,),)),),
    )
    with pytest.raises(ValueError, match="spin channel"):
        build_fatband_dataset(wrong, projections)


def test_parse_projwfc_stdout_preserves_original_band_order_and_fractional_kpoints():
    from molecular_lab.calculations.qe_projection import parse_qe_projwfc_band_output

    output = """
    Atomic states used for projection
    state # 1: atom 1 (Si), wfc 1 (l=0 m=0)
    k = 0.0000000000 0.0000000000 0.0000000000
    ==== e(   1) = -1.25000 eV ====
    ==== e(   2) =  0.75000 eV ====
    k = 0.5000000000 0.0000000000 0.0000000000
    ==== e(   1) = -1.10000 eV ====
    ==== e(   2) =  0.90000 eV ====
    """
    from molecular_lab.calculations.dft_settings import BandPath, BandPathPoint
    path = BandPath((
        BandPathPoint((0.0, 0.0, 0.0), "GAMMA", 1),
        BandPathPoint((0.5, 0.0, 0.0), "X", 0),
    ))
    data = parse_qe_projwfc_band_output(output, fermi_level_ev=0.0, band_path=path)
    assert data.channels[0].spin_channel == "total"
    assert data.channels[0].energies_ev == (
        (-1.25, -1.1),
        (0.75, 0.9),
    )
    assert data.kpoints[1].coordinates == (0.5, 0.0, 0.0)
    assert data.kpoints[1].coordinate_system == "fractional"


def test_parse_projwfc_stdout_rejects_incomplete_band_indices():
    from molecular_lab.calculations.qe_projection import (
        QEProjectionParseError, parse_qe_projwfc_band_output,
    )

    output = """
    k = 0 0 0
    ==== e( 2) = 0.5 eV ====
    """
    with pytest.raises(QEProjectionParseError, match="contiguous"):
        parse_qe_projwfc_band_output(output)
