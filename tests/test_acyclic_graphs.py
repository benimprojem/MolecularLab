"""
Tests for bond-order-aware acyclic molecular graph proposals.

Version: 0.1
Rev.: 00
"""

from collections import Counter

import pytest

from molecular_lab.chemistry.formula import parse_formula
from molecular_lab.engines.molecular_graphs import iter_acyclic_graphs


def graph_bond_orders(graph):
    return sorted(bond.order for bond in graph.bonds)


@pytest.mark.parametrize(
    ("formula", "expected_order"),
    [("C2H4", 2), ("C2H2", 3)],
)
def test_carbon_carbon_unsaturation_is_generated(formula, expected_order):
    graphs = list(iter_acyclic_graphs(parse_formula(formula)))
    assert graphs
    assert any(expected_order in graph_bond_orders(graph) for graph in graphs)
    assert all(Counter(a.element for a in graph.atoms) == Counter(
        parse_formula(formula).counts
    ) for graph in graphs)


def test_formaldehyde_double_bond_is_generated():
    graphs = list(iter_acyclic_graphs(parse_formula("CH2O")))
    assert len(graphs) == 1
    graph = graphs[0]
    carbon = next(a.atom_id for a in graph.atoms if a.element == "C")
    oxygen = next(a.atom_id for a in graph.atoms if a.element == "O")
    assert any(
        {bond.atom_a, bond.atom_b} == {carbon, oxygen} and bond.order == 2
        for bond in graph.bonds
    )


def test_carbon_dioxide_has_two_double_bonds():
    graphs = list(iter_acyclic_graphs(parse_formula("CO2")))
    assert len(graphs) == 1
    assert graph_bond_orders(graphs[0]).count(2) == 2


def test_single_bond_saturated_examples_remain_available():
    graphs = list(iter_acyclic_graphs(parse_formula("C2H6O")))
    assert len(graphs) == 2
    assert all(set(graph_bond_orders(graph)) == {1} for graph in graphs)


def test_assignment_limit_is_validated():
    with pytest.raises(ValueError, match="max_assignments"):
        list(iter_acyclic_graphs(parse_formula("C2H4"), max_assignments=0))


def test_candidate_limit_is_validated():
    with pytest.raises(ValueError, match="max_candidates"):
        list(iter_acyclic_graphs(parse_formula("C2H4"), max_candidates=0))
