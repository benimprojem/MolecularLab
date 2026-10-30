"""
Tests for bounded linearly fused acene graph proposals.

Version: 0.1
Rev.: 00
"""

from collections import Counter

import pytest

from molecular_lab.chemistry.formula import parse_formula
from molecular_lab.engines.molecular_graphs import iter_linear_fused_acene_graphs
from molecular_lab.engines.structure_identification import identify_graph_against_catalogue


@pytest.mark.parametrize(
    ("formula", "expected_cycle_rank", "expected_aromatic_carbons"),
    [
        ("C10H8", 2, 10),  # naphthalene
        ("C14H10", 3, 14),  # anthracene
    ],
)
def test_linear_fused_acene_topology_and_valence(
    formula, expected_cycle_rank, expected_aromatic_carbons
):
    graphs = list(iter_linear_fused_acene_graphs(parse_formula(formula)))
    assert graphs

    for graph in graphs:
        heavy_ids = {atom.atom_id for atom in graph.atoms if atom.element != "H"}
        heavy_edges = [
            bond for bond in graph.bonds
            if bond.atom_a in heavy_ids and bond.atom_b in heavy_ids
        ]
        assert len(heavy_edges) - len(heavy_ids) + 1 == expected_cycle_rank
        assert Counter(atom.element for atom in graph.atoms) == Counter(
            {"C": expected_aromatic_carbons, "H": int(formula.split("H")[1])}
        )

        valence_use = {atom.atom_id: 0 for atom in graph.atoms}
        for bond in graph.bonds:
            valence_use[bond.atom_a] += bond.order
            valence_use[bond.atom_b] += bond.order
        assert all(
            valence_use[atom.atom_id] == (1 if atom.element == "H" else 4)
            for atom in graph.atoms
        )

        assert sum(bond.order == 2 for bond in heavy_edges) == expected_aromatic_carbons // 2
        assert all(bond.order in (1, 2) for bond in heavy_edges)


def test_fused_acene_generator_respects_candidate_and_atom_bounds():
    composition = parse_formula("C10H8")
    assert len(list(iter_linear_fused_acene_graphs(composition, max_candidates=1))) == 1
    assert list(iter_linear_fused_acene_graphs(
        composition, max_heavy_atoms=8
    )) == []


def test_fused_acene_generator_rejects_out_of_scope_formulas():
    for formula in ("C6H6", "C10H10", "C10H8O", "C14H10"):
        if formula == "C14H10":
            continue
        assert list(iter_linear_fused_acene_graphs(parse_formula(formula))) == []


def test_naphthalene_and_anthracene_are_recognized_as_aromatic_when_rdkit_exists():
    pytest.importorskip("rdkit")
    for formula in ("C10H8", "C14H10"):
        graphs = list(iter_linear_fused_acene_graphs(parse_formula(formula)))
        assert graphs
        identified = [identify_graph_against_catalogue(graph) for graph in graphs]
        assert any(item.aromaticity_status == "aromatic_detected" for item in identified)


def test_fused_acene_search_limit_is_validated():
    with pytest.raises(ValueError, match="max_assignments"):
        list(iter_linear_fused_acene_graphs(
            parse_formula("C10H8"), max_assignments=0
        ))
