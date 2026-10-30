"""
Tests for the integrated forward molecular discovery workflow.

Version: 0.1
Rev.: 00
"""

import pytest

from molecular_lab.core.models import SearchConstraints, SearchDomain
from molecular_lab.engines.forward_discovery import (
    iter_forward_discovery,
    run_forward_discovery,
)


def test_hydrogen_oxygen_search_runs_composition_to_graph_pipeline(monkeypatch):
    from molecular_lab.engines import structure_identification

    monkeypatch.setattr(
        structure_identification,
        "_load_rdkit",
        lambda: None,
    )
    constraints = SearchConstraints(
        allowed_elements=frozenset({"H", "O"}),
        min_total_atoms=1,
        max_total_atoms=4,
        domain=SearchDomain.MOLECULE,
    )
    results = run_forward_discovery(
        constraints,
        required_elements={"H", "O"},
        max_compositions=20,
    )

    water = next(
        result for result in results
        if dict(result.composition.counts) == {"H": 2, "O": 1}
    )
    peroxide = next(
        result for result in results
        if dict(result.composition.counts) == {"H": 2, "O": 2}
    )

    assert water.graph_generation_status == "graphs_proposed"
    assert len(water.graphs) == 1
    assert water.reference_lookup.match_level.value == "formula_only"
    assert peroxide.graph_generation_status == "graphs_proposed"
    assert len(peroxide.graphs) == 1
    assert peroxide.reference_lookup.match_level.value == "formula_only"
    assert len(peroxide.graph_identifications) == 1


def test_formula_matches_are_reported_alongside_graph_proposals(monkeypatch):
    from molecular_lab.engines import structure_identification

    monkeypatch.setattr(structure_identification, "_load_rdkit", lambda: None)
    constraints = SearchConstraints(
        allowed_elements=frozenset({"C", "H", "O"}),
        min_total_atoms=9,
        max_total_atoms=9,
        domain=SearchDomain.MOLECULE,
        per_element_max={"C": 2, "H": 6, "O": 1},
    )
    results = run_forward_discovery(
        constraints,
        required_elements={"C", "H", "O"},
        max_compositions=10,
    )
    assert results
    result = results[0]
    assert result.reference_lookup.match_level.value == "formula_only"
    assert result.graph_generation_status == "graphs_proposed"


def test_periodic_formula_unit_is_not_sent_to_molecular_graph_generator():
    constraints = SearchConstraints(
        allowed_elements=frozenset({"H", "O"}),
        min_total_atoms=3,
        max_total_atoms=3,
        domain=SearchDomain.MIXED,
    )
    results = run_forward_discovery(
        constraints,
        required_elements={"H", "O"},
        max_compositions=2,
    )

    periodic = next(
        result for result in results
        if result.composition.kind.value == "formula_unit"
    )
    assert periodic.graphs == ()
    assert periodic.graph_generation_status == "not_attempted_non_molecular"
    assert periodic.reference_lookup.match_level.value == "not_applicable"


def test_charged_compositions_are_visible_but_not_graph_generated():
    constraints = SearchConstraints(
        allowed_elements=frozenset({"H", "O"}),
        min_total_atoms=3,
        max_total_atoms=3,
        domain=SearchDomain.MOLECULE,
        allowed_charges=(1,),
    )
    result = next(iter_forward_discovery(
        constraints,
        required_elements={"H", "O"},
        max_compositions=1,
    ))
    assert result.graphs == ()
    assert result.graph_generation_status == "not_attempted_charged"


def test_composition_limit_is_enforced():
    constraints = SearchConstraints(
        allowed_elements=frozenset({"H", "O"}),
        min_total_atoms=1,
        max_total_atoms=5,
        domain=SearchDomain.MOLECULE,
    )
    assert len(run_forward_discovery(constraints, max_compositions=3)) == 3


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_compositions": 0},
        {"max_graphs_per_composition": 0},
        {"max_heavy_atoms": 0},
    ],
)
def test_invalid_workflow_limits_are_rejected(kwargs):
    constraints = SearchConstraints(
        allowed_elements=frozenset({"H", "O"}),
        min_total_atoms=1,
        max_total_atoms=2,
    )
    with pytest.raises(ValueError, match="positive integer"):
        run_forward_discovery(constraints, **kwargs)


def test_periodic_only_domain_is_rejected_explicitly():
    constraints = SearchConstraints(
        allowed_elements=frozenset({"Na", "Cl"}),
        min_total_atoms=2,
        max_total_atoms=2,
        domain=SearchDomain.PERIODIC,
    )
    with pytest.raises(ValueError, match="Periodic-only"):
        run_forward_discovery(constraints)


def test_forward_pipeline_includes_unicyclic_ring_proposals(monkeypatch):
    from molecular_lab.engines import structure_identification

    monkeypatch.setattr(structure_identification, "_load_rdkit", lambda: None)
    constraints = SearchConstraints(
        allowed_elements=frozenset({"C", "H"}),
        min_total_atoms=9,
        max_total_atoms=9,
        domain=SearchDomain.MOLECULE,
        per_element_max={"C": 3, "H": 6},
    )
    results = run_forward_discovery(
        constraints,
        required_elements={"C", "H"},
        max_compositions=10,
    )
    cyclopropane_formula = next(
        result for result in results
        if dict(result.composition.counts) == {"C": 3, "H": 6}
    )
    assert any(
        graph.method == "bounded_unicyclic_bond_order_graph_v1"
        for graph in cyclopropane_formula.graphs
    )
