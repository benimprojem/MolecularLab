"""
Tests for the restricted saturated acyclic molecular graph prototype.

Version: 0.1
Rev.: 01
"""

from collections import Counter

import pytest

from molecular_lab.chemistry.formula import parse_formula
from molecular_lab.core.models import MolecularAtom, MolecularBond, MolecularGraph
from molecular_lab.engines.molecular_graphs import iter_saturated_acyclic_graphs


def atom_counts(graph):
    return Counter(atom.element for atom in graph.atoms)


def test_water_graph_is_proposed():
    graphs = list(iter_saturated_acyclic_graphs(parse_formula("H2O")))
    assert len(graphs) == 1
    assert atom_counts(graphs[0]) == {"H": 2, "O": 1}
    assert len(graphs[0].bonds) == 2


def test_hydrogen_peroxide_has_oxygen_oxygen_connectivity():
    graphs = list(iter_saturated_acyclic_graphs(parse_formula("H2O2")))
    assert len(graphs) == 1
    graph = graphs[0]
    oxygen_ids = {atom.atom_id for atom in graph.atoms if atom.element == "O"}
    assert any({bond.atom_a, bond.atom_b} == oxygen_ids for bond in graph.bonds)
    assert atom_counts(graph) == {"H": 2, "O": 2}


def test_methane_graph_is_proposed():
    graphs = list(iter_saturated_acyclic_graphs(parse_formula("CH4")))
    assert len(graphs) == 1
    assert atom_counts(graphs[0]) == {"C": 1, "H": 4}


def test_c2h6o_connectivity_isomers_are_proposed():
    graphs = list(iter_saturated_acyclic_graphs(parse_formula("C2H6O")))
    assert len(graphs) == 2
    oxygen_heavy_degrees = []
    for graph in graphs:
        oxygen_id = next(atom.atom_id for atom in graph.atoms if atom.element == "O")
        heavy_ids = {atom.atom_id for atom in graph.atoms if atom.element != "H"}
        oxygen_heavy_degrees.append(sum(
            1 for bond in graph.bonds
            if oxygen_id in (bond.atom_a, bond.atom_b)
            and bond.atom_a in heavy_ids and bond.atom_b in heavy_ids
        ))
    # Ethanol has one heavy-atom neighbour; dimethyl ether has two.
    assert sorted(oxygen_heavy_degrees) == [1, 2]


def test_unsupported_formula_is_not_misreported_as_impossible():
    # A formula requiring a double bond is outside this saturated-tree model.
    assert list(iter_saturated_acyclic_graphs(parse_formula("C2H4"))) == []


def test_charged_species_are_not_generated():
    charged = parse_formula("H2O")
    from molecular_lab.core.models import Composition
    charged = Composition(counts=charged.counts, net_charge=1)
    assert list(iter_saturated_acyclic_graphs(charged)) == []


def test_unknown_element_raises_from_identity_validation():
    from molecular_lab.core.models import Composition
    with pytest.raises(ValueError, match="Unknown element"):
        list(iter_saturated_acyclic_graphs(Composition(counts={"Xx": 1})))


def test_graph_rejects_duplicate_atom_ids():
    with pytest.raises(ValueError, match="unique"):
        MolecularGraph(
            atoms=(MolecularAtom(0, "H"), MolecularAtom(0, "H")),
            bonds=(),
            method="test",
        )


def test_graph_rejects_unknown_bond_endpoint():
    with pytest.raises(ValueError, match="unknown atom"):
        MolecularGraph(
            atoms=(MolecularAtom(0, "H"),),
            bonds=(MolecularBond(0, 1),),
            method="test",
        )


def test_graph_rejects_duplicate_bonds():
    with pytest.raises(ValueError, match="duplicate bond"):
        MolecularGraph(
            atoms=(MolecularAtom(0, "H"), MolecularAtom(1, "H")),
            bonds=(MolecularBond(0, 1), MolecularBond(1, 0)),
            method="test",
        )
