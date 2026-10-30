"""
Tests for foundational Molecular Lab domain models.

Version: 0.1
Rev.: 00
"""

import pytest

from molecular_lab.core.models import (
    Candidate,
    Composition,
    CompositionKind,
    EvidenceCategory,
    EvidenceRecord,
    SearchConstraints,
    SearchDomain,
)


def test_composition_preserves_molecular_atom_counts():
    water_peroxide = Composition(
        counts={"H": 2, "O": 2},
        kind=CompositionKind.MOLECULAR,
        input_formula="H2O2",
    )
    assert water_peroxide.total_atoms == 4
    assert dict(water_peroxide.counts) == {"H": 2, "O": 2}
    assert water_peroxide.input_formula == "H2O2"


def test_composition_rejects_zero_or_negative_counts():
    with pytest.raises(ValueError):
        Composition(counts={"H": 2, "O": 0})


def test_composition_rejects_boolean_count():
    with pytest.raises(ValueError):
        Composition(counts={"H": True})


def test_search_constraints_check_bounds_only():
    constraints = SearchConstraints(
        allowed_elements=frozenset({"H", "O"}),
        min_total_atoms=2,
        max_total_atoms=6,
        domain=SearchDomain.MOLECULE,
        per_element_max={"H": 4, "O": 3},
        allowed_charges=(0,),
    )
    candidate = Composition(counts={"H": 2, "O": 2}, net_charge=0)
    outside_bounds = Composition(counts={"H": 5, "O": 2}, net_charge=0)
    disallowed_element = Composition(counts={"H": 2, "N": 1}, net_charge=0)

    assert constraints.accepts_composition(candidate)
    assert not constraints.accepts_composition(outside_bounds)
    assert not constraints.accepts_composition(disallowed_element)


def test_search_constraints_reject_invalid_limits():
    with pytest.raises(ValueError):
        SearchConstraints(allowed_elements=frozenset({"H"}), max_total_atoms=0)


def test_candidate_keeps_evidence_separate():
    composition = Composition(counts={"H": 2, "O": 1})
    evidence = EvidenceRecord(
        category=EvidenceCategory.EXPERIMENTAL,
        claim="Illustrative evidence record",
        source="Example source",
    )
    candidate = Candidate(candidate_id="candidate-001", composition=composition, evidence=(evidence,))

    assert candidate.composition.total_atoms == 3
    assert candidate.evidence[0].category is EvidenceCategory.EXPERIMENTAL
