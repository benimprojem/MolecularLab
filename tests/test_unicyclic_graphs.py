"""
Tests for bounded unicyclic molecular graph proposals.

Version: 0.1
Rev.: 00
"""

from collections import Counter

import pytest

from molecular_lab.chemistry.formula import parse_formula
from molecular_lab.engines.molecular_graphs import iter_unicyclic_graphs


def heavy_graph_cycle_rank(graph):
    heavy_ids = {atom.atom_id for atom in graph.atoms if atom.element != "H"}
    heavy_edges = [
        bond for bond in graph.bonds
        if bond.atom_a in heavy_ids and bond.atom_b in heavy_ids
    ]
    return len(heavy_edges) - len(heavy_ids) + 1


def test_cyclopropane_ring_is_generated():
    graphs = list(iter_unicyclic_graphs(parse_formula("C3H6")))
    assert graphs
    assert len(graphs) == 1
    assert all(heavy_graph_cycle_rank(graph) == 1 for graph in graphs)
    assert any(
        Counter(atom.element for atom in graph.atoms) == Counter({"C": 3, "H": 6})
        for graph in graphs
    )


def test_cyclobutane_and_methylcyclopropane_topologies_are_in_scope():
    graphs = list(iter_unicyclic_graphs(parse_formula("C4H8")))
    assert graphs
    assert all(heavy_graph_cycle_rank(graph) == 1 for graph in graphs)
    assert all(
        Counter(atom.element for atom in graph.atoms) == Counter({"C": 4, "H": 8})
        for graph in graphs
    )
    heavy_degrees = set()
    for graph in graphs:
        heavy_ids = {a.atom_id for a in graph.atoms if a.element != "H"}
        degrees = tuple(sorted(
            sum(
                1 for bond in graph.bonds
                if atom_id in (bond.atom_a, bond.atom_b)
                and bond.atom_a in heavy_ids and bond.atom_b in heavy_ids
            )
            for atom_id in heavy_ids
        ))
        heavy_degrees.add(degrees)
    assert len(heavy_degrees) >= 2


def test_ring_bond_orders_obey_fixed_carbon_valence():
    graphs = list(iter_unicyclic_graphs(parse_formula("C3H6")))
    for graph in graphs:
        valence_use = {a.atom_id: 0 for a in graph.atoms}
        for bond in graph.bonds:
            valence_use[bond.atom_a] += bond.order
            valence_use[bond.atom_b] += bond.order
        assert all(
            valence_use[a.atom_id] == (1 if a.element == "H" else 4)
            for a in graph.atoms
        )


def test_too_few_heavy_atoms_cannot_form_a_simple_cycle():
    assert list(iter_unicyclic_graphs(parse_formula("CH4"))) == []
    assert list(iter_unicyclic_graphs(parse_formula("C2H4"))) == []


def test_assignment_and_candidate_limits_are_validated():
    with pytest.raises(ValueError, match="max_assignments"):
        list(iter_unicyclic_graphs(parse_formula("C3H6"), max_assignments=0))
    with pytest.raises(ValueError, match="max_candidates"):
        list(iter_unicyclic_graphs(parse_formula("C3H6"), max_candidates=0))
