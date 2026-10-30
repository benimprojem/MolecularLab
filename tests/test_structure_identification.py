"""
Tests for molecular-graph reference identification.

Version: 0.1
Rev.: 00
"""

import pytest

from molecular_lab.chemistry.formula import parse_formula
from molecular_lab.engines.molecular_graphs import iter_saturated_acyclic_graphs
from molecular_lab.engines import structure_identification


def _ethanol_and_ether_graphs():
    return list(iter_saturated_acyclic_graphs(parse_formula("C2H6O")))


def test_graph_identification_falls_back_to_formula_when_rdkit_missing(monkeypatch):
    monkeypatch.setattr(structure_identification, "_load_rdkit", lambda: None)
    results = [
        structure_identification.identify_graph_against_catalogue(graph)
        for graph in _ethanol_and_ether_graphs()
    ]

    assert all(result.status == "rdkit_unavailable" for result in results)
    assert all(result.match_level.value == "formula_only" for result in results)
    assert all(result.canonical_smiles is None for result in results)
    assert {
        match.preferred_name
        for result in results
        for match in result.matches
    } == {"Ethanol", "Dimethyl ether"}


def test_graph_identification_distinguishes_c2h6o_isomers_when_rdkit_available():
    pytest.importorskip("rdkit")
    results = [
        structure_identification.identify_graph_against_catalogue(graph)
        for graph in _ethanol_and_ether_graphs()
    ]

    assert all(result.status == "identified" for result in results)
    assert all(result.inchi_key for result in results)
    assert all(result.canonical_smiles for result in results)
    assert all(result.match_level.value == "exact_structure" for result in results)
    assert {
        result.matches[0].preferred_name
        for result in results
    } == {"Ethanol", "Dimethyl ether"}


def test_graph_identification_handles_known_water():
    pytest.importorskip("rdkit")
    graph = next(iter_saturated_acyclic_graphs(parse_formula("H2O")))
    result = structure_identification.identify_graph_against_catalogue(graph)

    assert result.formula == "H2O"
    assert result.match_level.value == "exact_structure"
    assert result.matches[0].preferred_name == "Water"


def test_missing_formula_match_is_not_claimed_to_be_new(monkeypatch):
    monkeypatch.setattr(structure_identification, "_load_rdkit", lambda: None)
    graph = next(iter_saturated_acyclic_graphs(parse_formula("C2H6")))
    result = structure_identification.identify_graph_against_catalogue(graph)

    assert result.match_level.value == "not_found"
    assert "does not prove" in result.explanation


def test_benzene_is_reported_as_aromatic_when_rdkit_is_available():
    pytest.importorskip("rdkit")
    from molecular_lab.engines.molecular_graphs import iter_unicyclic_graphs

    candidates = iter_unicyclic_graphs(
        parse_formula("C6H6"), max_candidates=1000, max_assignments=100000
    )
    graph = next(
        candidate for candidate in candidates
        if all(
            sum(
                1 for bond in candidate.bonds
                if bond.order == 2
                and atom.atom_id in (bond.atom_a, bond.atom_b)
            ) == 1
            for atom in candidate.atoms if atom.element == "C"
        )
        and all(bond.order in (1, 2) for bond in candidate.bonds
                if bond.atom_a < 6 and bond.atom_b < 6)
    )
    result = structure_identification.identify_graph_against_catalogue(graph)
    assert result.status == "identified"
    assert result.aromaticity_status == "aromatic_detected"
    assert "c1ccccc1" in result.canonical_smiles
