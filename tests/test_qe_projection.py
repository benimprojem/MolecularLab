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


def test_projection_data_rejects_inconsistent_weight_and_nonfinite_values():
    from dataclasses import replace

    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES><PROJS>
        <ATOMIC_WFC index="1">0.5 0.0</ATOMIC_WFC>
        <ATOMIC_WFC index="2">0.0 0.5</ATOMIC_WFC>
        <ATOMIC_WFC index="3">0.1 0.2</ATOMIC_WFC>
      </PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    data = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    item = data.coefficients[0]
    with pytest.raises(ValueError, match="real²"):
        replace(item, weight=item.weight + 0.1)
    with pytest.raises(ValueError, match="finite"):
        replace(item, real=float("nan"))


def test_projection_data_rejects_missing_duplicate_and_out_of_range_keys():
    from dataclasses import replace

    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES><PROJS>
        <ATOMIC_WFC index="1">0.5 0.0</ATOMIC_WFC>
        <ATOMIC_WFC index="2">0.0 0.5</ATOMIC_WFC>
        <ATOMIC_WFC index="3">0.1 0.2</ATOMIC_WFC>
      </PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    data = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    with pytest.raises(ValueError, match="grid"):
        replace(data, coefficients=data.coefficients[:-1])
    with pytest.raises(ValueError, match="duplicate"):
        replace(data, coefficients=data.coefficients + (data.coefficients[0],))
    with pytest.raises(ValueError, match="grid"):
        replace(data, coefficients=(replace(data.coefficients[0], band_index=2),) + data.coefficients[1:])


def test_fatband_dataset_rejects_duplicate_state_and_orbital_records():
    from dataclasses import replace

    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES><PROJS>
        <ATOMIC_WFC index="1">0.5 0.0</ATOMIC_WFC>
        <ATOMIC_WFC index="2">0.0 0.5</ATOMIC_WFC>
        <ATOMIC_WFC index="3">0.1 0.2</ATOMIC_WFC>
      </PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    projections = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    bands = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0),),
        channels=(BandChannelData("total", ((0.0,),)),),
    )
    dataset = build_fatband_dataset(bands, projections)
    with pytest.raises(ValueError, match="duplicate state"):
        replace(dataset, points=dataset.points + (dataset.points[0],))
    point = dataset.points[0]
    duplicate = replace(point, contributions=point.contributions + (point.contributions[0],))
    with pytest.raises(ValueError, match="duplicate atom/orbital"):
        replace(dataset, points=(duplicate,))


def test_fatband_builder_checks_projection_atom_identity_against_crystal():
    from molecular_lab.structures.crystal import CrystalSite, CrystalStructure

    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES><PROJS>
        <ATOMIC_WFC index="1">0.5 0.0</ATOMIC_WFC>
        <ATOMIC_WFC index="2">0.0 0.5</ATOMIC_WFC>
        <ATOMIC_WFC index="3">0.1 0.2</ATOMIC_WFC>
      </PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    projections = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    bands = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0),),
        channels=(BandChannelData("total", ((0.0,),)),),
    )
    structure = CrystalStructure(
        structure_id="FeO-test",
        lattice_vectors_angstrom=((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
        sites=(
            CrystalSite("Fe1", "Fe", (0.0, 0.0, 0.0)),
            CrystalSite("O1", "O", (0.5, 0.5, 0.5)),
        ),
    )
    result = build_fatband_dataset(bands, projections, structure=structure)
    assert len(result.points) == 1

    wrong_order = CrystalStructure(
        structure_id="OFe-test",
        lattice_vectors_angstrom=structure.lattice_vectors_angstrom,
        sites=(
            CrystalSite("O1", "O", (0.0, 0.0, 0.0)),
            CrystalSite("Fe1", "Fe", (0.5, 0.5, 0.5)),
        ),
    )
    with pytest.raises(ValueError, match="input structure site 1"):
        build_fatband_dataset(bands, projections, structure=wrong_order)


def test_band_energy_spectrum_validation_ignores_band_reordering():
    from molecular_lab.calculations.fatband import validate_band_energy_spectra

    original_order = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0), KPoint((0.5, 0.0, 0.0), 1.0)),
        channels=(BandChannelData("total", ((-1.0, -0.8), (1.0, 1.2))),),
    )
    reordered = BandStructureData(
        kpoints=original_order.kpoints,
        channels=(BandChannelData("non_spin_polarized", ((1.0, 1.2), (-1.0, -0.8))),),
    )
    validate_band_energy_spectra(original_order, reordered)


def test_band_energy_spectrum_validation_detects_wrong_kpoint_spectrum():
    from molecular_lab.calculations.fatband import validate_band_energy_spectra

    original = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0), KPoint((0.5, 0.0, 0.0), 1.0)),
        channels=(BandChannelData("total", ((-1.0, -0.8), (1.0, 1.2))),),
    )
    mismatched = BandStructureData(
        kpoints=original.kpoints,
        channels=(BandChannelData("non_spin_polarized", ((-1.0, -0.8), (1.0, 1.3))),),
    )
    with pytest.raises(ValueError, match="spectra disagree at k-point 1"):
        validate_band_energy_spectra(original, mismatched)


def test_band_energy_spectrum_validation_checks_tolerance_and_dimensions():
    from molecular_lab.calculations.fatband import validate_band_energy_spectra

    reference = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0),),
        channels=(BandChannelData("total", ((-1.0,), (1.0,))),),
    )
    slightly_shifted = BandStructureData(
        kpoints=reference.kpoints,
        channels=(BandChannelData("total", ((-0.99995,), (1.00005,))),),
    )
    validate_band_energy_spectra(reference, slightly_shifted)
    with pytest.raises(ValueError, match="tolerance"):
        validate_band_energy_spectra(reference, slightly_shifted, tolerance_ev=1e-6)
    wrong_count = BandStructureData(
        kpoints=reference.kpoints,
        channels=(BandChannelData("total", ((-1.0,),)),),
    )
    with pytest.raises(ValueError, match="Band count"):
        validate_band_energy_spectra(reference, wrong_count)


def test_fatband_projection_diagnostics_report_zero_and_high_weights_without_normalizing():
    from dataclasses import replace
    from molecular_lab.calculations.fatband import diagnose_fatband_projection_weights

    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES><PROJS>
        <ATOMIC_WFC index="1">0.5 0.0</ATOMIC_WFC>
        <ATOMIC_WFC index="2">0.0 0.5</ATOMIC_WFC>
        <ATOMIC_WFC index="3">0.1 0.2</ATOMIC_WFC>
      </PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    projections = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    bands = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0),),
        channels=(BandChannelData("total", ((0.0,),)),),
    )
    dataset = build_fatband_dataset(bands, projections)
    point = dataset.points[0]
    high = replace(
        point,
        contributions=tuple(
            replace(item, raw_weight=1.25 if item.angular_momentum == 0 else item.raw_weight)
            for item in point.contributions
        ),
    )
    empty = replace(point, contributions=())
    changed = replace(dataset, points=(high,))
    warnings = diagnose_fatband_projection_weights(changed)
    assert any(item.code == "high_total_weight" for item in warnings)
    assert any(item.code == "high_orbital_weight" for item in warnings)
    assert high.total_projected_weight > 1.0
    assert not diagnose_fatband_projection_weights(replace(dataset, points=(empty,)))[0].code == "high_total_weight"
    assert diagnose_fatband_projection_weights(replace(dataset, points=(empty,)))[0].code == "zero_total_weight"


def test_fatband_projection_diagnostics_validate_thresholds():
    from molecular_lab.calculations.fatband import diagnose_fatband_projection_weights

    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES><PROJS>
        <ATOMIC_WFC index="1">0.5 0.0</ATOMIC_WFC>
        <ATOMIC_WFC index="2">0.0 0.5</ATOMIC_WFC>
        <ATOMIC_WFC index="3">0.1 0.2</ATOMIC_WFC>
      </PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    projections = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    bands = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0),),
        channels=(BandChannelData("total", ((0.0,),)),),
    )
    dataset = build_fatband_dataset(bands, projections)
    with pytest.raises(ValueError, match="thresholds"):
        diagnose_fatband_projection_weights(dataset, high_total_weight_threshold=float("nan"))
    with pytest.raises(ValueError, match="thresholds"):
        diagnose_fatband_projection_weights(dataset, zero_weight_tolerance=True)


def test_fatband_builder_attaches_projection_warnings_automatically():
    from molecular_lab.calculations.fatband import diagnose_fatband_projection_weights

    xml = """<ATOMIC_PROJECTIONS>
      <HEADER NUMBER_OF_BANDS="1" NUMBER_OF_K-POINTS="1"
              NUMBER_OF_SPIN_COMPONENTS="1" NUMBER_OF_ATOMIC_WFC="3"/>
      <EIGENSTATES><PROJS>
        <ATOMIC_WFC index="1">0.0 0.0</ATOMIC_WFC>
        <ATOMIC_WFC index="2">0.0 0.0</ATOMIC_WFC>
        <ATOMIC_WFC index="3">0.0 0.0</ATOMIC_WFC>
      </PROJS></EIGENSTATES>
    </ATOMIC_PROJECTIONS>"""
    projections = parse_qe_atomic_proj_xml(xml, STATE_OUTPUT)
    bands = BandStructureData(
        kpoints=(KPoint((0.0, 0.0, 0.0), 1.0),),
        channels=(BandChannelData("total", ((0.0,),)),),
    )
    dataset = build_fatband_dataset(bands, projections)
    assert len(dataset.projection_warnings) == 1
    assert dataset.projection_warnings[0].code == "zero_total_weight"
    assert dataset.projection_warnings == diagnose_fatband_projection_weights(dataset)
