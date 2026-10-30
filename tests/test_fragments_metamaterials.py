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



def test_multiphysics_property_targets_cover_optical_magnetic_and_electrical():
    from molecular_lab.materials.metamaterials import (
        MetamaterialPropertyDomain,
        MetamaterialPropertyTarget,
        TargetRelation,
    )

    cell = MetamaterialUnitCell("Coupled cell", "hybrid_resonator", 3, (1e-6, 1e-6, 1e-6))
    design = MetamaterialDesign(
        design_id="meta-multiphysics-001",
        mechanism=MetamaterialMechanism.HYBRID,
        unit_cell=cell,
        wavelength_band_m=(1.0e-6, 2.0e-6),
        property_targets=(
            MetamaterialPropertyTarget(
                MetamaterialPropertyDomain.OPTICAL,
                "effective_refractive_index",
                TargetRelation.APPROACH,
                "1",
                target_value=-1.0,
                spectral_band_m=(1.2e-6, 1.6e-6),
            ),
            MetamaterialPropertyTarget(
                MetamaterialPropertyDomain.MAGNETIC,
                "relative_permeability",
                TargetRelation.MAXIMIZE,
                "1",
                minimum=1.5,
            ),
            MetamaterialPropertyTarget(
                MetamaterialPropertyDomain.ELECTRICAL,
                "electrical_conductivity",
                TargetRelation.WITHIN_RANGE,
                "S/m",
                minimum=1e3,
                maximum=1e5,
            ),
            MetamaterialPropertyTarget(
                MetamaterialPropertyDomain.ELECTRONIC_STRUCTURE,
                "band_gap",
                TargetRelation.WITHIN_RANGE,
                "eV",
                minimum=0.5,
                maximum=2.0,
            ),
        ),
    )
    assert len(design.property_targets) == 4
    assert design.property_targets[0].target_value == -1.0
    assert design.property_targets[2].unit == "S/m"


def test_property_target_requires_valid_target_range():
    from molecular_lab.materials.metamaterials import (
        MetamaterialPropertyDomain,
        MetamaterialPropertyTarget,
        TargetRelation,
    )

    with pytest.raises(ValueError, match="within_range"):
        MetamaterialPropertyTarget(
            MetamaterialPropertyDomain.ELECTRICAL,
            "electrical_conductivity",
            TargetRelation.WITHIN_RANGE,
            "S/m",
            minimum=10.0,
        )
    with pytest.raises(ValueError, match="target_value"):
        MetamaterialPropertyTarget(
            MetamaterialPropertyDomain.OPTICAL,
            "refractive_index",
            TargetRelation.APPROACH,
            "1",
        )
