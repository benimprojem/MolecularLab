"""
Tests for bounded composition enumeration.

Version: 0.1
Rev.: 00
"""

import pytest

from molecular_lab.core.models import CompositionKind, SearchConstraints, SearchDomain
from molecular_lab.engines.composition_enumerator import iter_compositions


def _constraints(**overrides):
    values = {
        "allowed_elements": frozenset({"H", "O"}),
        "min_total_atoms": 1,
        "max_total_atoms": 4,
        "domain": SearchDomain.MOLECULE,
    }
    values.update(overrides)
    return SearchConstraints(**values)


def test_enumerator_is_bounded_and_deterministic():
    first = list(iter_compositions(_constraints()))
    second = list(iter_compositions(_constraints()))

    assert first == second
    assert all(1 <= item.total_atoms <= 4 for item in first)
    assert all(set(item.counts) <= {"H", "O"} for item in first)


def test_hydrogen_oxygen_compositions_include_water_and_peroxide():
    compositions = {
        tuple(sorted(item.counts.items()))
        for item in iter_compositions(_constraints())
    }

    assert (("H", 2), ("O", 1)) in compositions
    assert (("H", 2), ("O", 2)) in compositions


def test_per_element_maximum_is_respected():
    constraints = _constraints(per_element_max={"H": 2, "O": 1})
    results = list(iter_compositions(constraints))

    assert all(item.counts.get("H", 0) <= 2 for item in results)
    assert all(item.counts.get("O", 0) <= 1 for item in results)


def test_required_elements_must_be_present():
    results = list(iter_compositions(_constraints(), required_elements={"O"}))

    assert results
    assert all("O" in item.counts for item in results)


def test_required_element_must_be_in_allowed_set():
    with pytest.raises(ValueError, match="not allowed"):
        list(iter_compositions(_constraints(), required_elements={"C"}))


def test_charge_options_are_enumerated():
    constraints = _constraints(allowed_charges=(-1, 0, 1))
    results = list(iter_compositions(constraints))

    assert {item.net_charge for item in results} == {-1, 0, 1}


def test_unspecified_charge_defaults_to_neutral():
    results = list(iter_compositions(_constraints()))
    assert {item.net_charge for item in results} == {0}


def test_periodic_domain_uses_formula_unit_kind():
    constraints = _constraints(domain=SearchDomain.PERIODIC)
    results = list(iter_compositions(constraints))

    assert results
    assert {item.kind for item in results} == {CompositionKind.FORMULA_UNIT}


def test_mixed_domain_emits_both_interpretations():
    constraints = _constraints(domain=SearchDomain.MIXED)
    results = list(iter_compositions(constraints))

    assert {item.kind for item in results} == {
        CompositionKind.MOLECULAR,
        CompositionKind.FORMULA_UNIT,
    }


def test_required_all_selected_elements():
    results = list(
        iter_compositions(_constraints(), required_elements={"H", "O"})
    )

    assert results
    assert all(set(item.counts) == {"H", "O"} for item in results)
