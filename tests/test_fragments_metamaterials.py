"""Tests for fragment and metamaterial domain records."""

import pytest

from molecular_lab.core.models import (
    FragmentPort,
    MolecularAtom,
    MolecularBond,
    MolecularFragment,
    MolecularGraph,
)
from molecular_lab.materials.metamaterials import (
    MetamaterialDesign,
    MetamaterialMechanism,
    MetamaterialUnitCell,
    UnitCellParameter,
)


def _methyl_like_fragment():
    graph = MolecularGraph(
        atoms=(MolecularAtom(0, "C"), MolecularAtom(1, "H"),
               MolecularAtom(2, "H"), MolecularAtom(3, "H")),
        bonds=(MolecularBond(0, 1), MolecularBond(0, 2), MolecularBond(0, 3)),
        method="unit_test_fixture",
    )
    return MolecularFragment(
        fragment_id="test.methyl",
        name="Methyl attachment fragment",
        version="1",
        graph=graph,
        ports=(FragmentPort("carbon-attachment", 0, (1,)),),
        source="unit test fixture",
    )


def test_fragment_accepts_valid_attachment_port():
    fragment = _methyl_like_fragment()
    assert fragment.ports[0].atom_id == 0
    assert fragment.ports[0].allowed_bond_orders == (1,)


def test_fragment_rejects_unknown_port_atom():
    graph = MolecularGraph(
        atoms=(MolecularAtom(0, "C"),),
        bonds=(),
        method="unit_test_fixture",
    )
    with pytest.raises(ValueError, match="unknown atom"):
        MolecularFragment(
            fragment_id="bad",
            name="Bad fragment",
            version="1",
            graph=graph,
            ports=(FragmentPort("missing", 9),),
            source="unit test",
        )


def test_fragment_rejects_duplicate_port_ids():
    graph = MolecularGraph(
        atoms=(MolecularAtom(0, "C"),),
        bonds=(),
        method="unit_test_fixture",
    )
    with pytest.raises(ValueError, match="unique"):
        MolecularFragment(
            fragment_id="bad",
            name="Bad fragment",
            version="1",
            graph=graph,
            ports=(FragmentPort("same", 0), FragmentPort("same", 0)),
            source="unit test",
        )


def test_fragment_port_rejects_invalid_bond_order():
    with pytest.raises(ValueError):
        FragmentPort("bad", 0, (4,))


def test_metamaterial_unit_cell_accepts_2d_geometry():
    cell = MetamaterialUnitCell(
        name="Split-ring example",
        geometry_family="split_ring_resonator",
        dimension=2,
        lattice_period_m=(300e-9, 300e-9),
        parameters=(UnitCellParameter("ring_width", 20e-9, "m"),),
        constituent_materials=("Au", "SiO2"),
    )
    assert cell.dimension == 2
    assert cell.lattice_period_m == (300e-9, 300e-9)


def test_metamaterial_unit_cell_checks_dimension_and_periods():
    with pytest.raises(ValueError, match="dimension"):
        MetamaterialUnitCell("bad", "generic", 1, (1.0,))
    with pytest.raises(ValueError, match="positive"):
        MetamaterialUnitCell("bad", "generic", 2, (1.0, 0.0))


def test_metamaterial_design_records_band_and_targets():
    cell = MetamaterialUnitCell(
        name="Photonic lattice",
        geometry_family="dielectric_rods",
        dimension=2,
        lattice_period_m=(500e-9, 500e-9),
    )
    design = MetamaterialDesign(
        design_id="meta-001",
        mechanism=MetamaterialMechanism.PHOTONIC_CRYSTAL,
        unit_cell=cell,
        wavelength_band_m=(1.2e-6, 1.8e-6),
        target_properties=("photonic_band_gap", "low_loss_transmission"),
        operating_conditions={"temperature_K": 300.0},
    )
    assert design.mechanism is MetamaterialMechanism.PHOTONIC_CRYSTAL
    assert design.target_properties[0] == "photonic_band_gap"
    assert design.operating_conditions["temperature_K"] == 300.0


def test_metamaterial_design_rejects_invalid_wavelength_band():
    cell = MetamaterialUnitCell("cell", "generic", 2, (1e-6, 1e-6))
    with pytest.raises(ValueError, match="smaller"):
        MetamaterialDesign(
            design_id="bad",
            mechanism=MetamaterialMechanism.HYBRID,
            unit_cell=cell,
            wavelength_band_m=(2e-6, 1e-6),
            target_properties=("target",),
        )
