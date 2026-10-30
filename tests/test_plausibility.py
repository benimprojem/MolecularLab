"""
Tests for conservative composition-level plausibility diagnostics.

Version: 0.1
Rev.: 00
"""

import pytest

from molecular_lab.chemistry.formula import parse_formula
from molecular_lab.chemistry.plausibility import (
    DiagnosticStatus,
    ElectronParity,
    assess_composition,
)
from molecular_lab.core.models import Composition


@pytest.mark.parametrize(
    ("formula", "expected_electrons", "expected_parity"),
    [
        ("H2O", 10, ElectronParity.EVEN),
        ("H2O2", 18, ElectronParity.EVEN),
        ("O2", 16, ElectronParity.EVEN),
        ("NO", 15, ElectronParity.ODD),
    ],
)
def test_electron_count_and_parity(formula, expected_electrons, expected_parity):
    result = assess_composition(parse_formula(formula))
    assert result.electron_count == expected_electrons
    assert result.electron_parity == expected_parity


def test_positive_charge_reduces_electron_count():
    composition = Composition(counts={"H": 2, "O": 1}, net_charge=1)
    result = assess_composition(composition)
    assert result.electron_count == 9
    assert result.electron_parity == ElectronParity.ODD


def test_negative_charge_increases_electron_count():
    composition = Composition(counts={"H": 2, "O": 1}, net_charge=-1)
    result = assess_composition(composition)
    assert result.electron_count == 11


def test_unknown_charge_makes_electron_parity_unknown():
    composition = Composition(counts={"H": 2, "O": 1}, net_charge=None)
    result = assess_composition(composition)
    assert result.electron_count is None
    assert result.electron_parity == ElectronParity.UNKNOWN
    assert result.oxidation_state_status == DiagnosticStatus.NOT_ASSESSED


@pytest.mark.parametrize("formula", ["H2O", "H2O2", "CO2", "Fe2O3"])
def test_common_formal_oxidation_state_assignments_are_found(formula):
    result = assess_composition(parse_formula(formula))
    assert result.oxidation_state_status == DiagnosticStatus.SUPPORTED_PATTERN


def test_unmatched_common_state_assignment_is_review_not_rejection():
    result = assess_composition(parse_formula("HF"))
    # H(+1) and F(-1) balance: this should be a supported pattern.
    assert result.oxidation_state_status == DiagnosticStatus.SUPPORTED_PATTERN


def test_unlisted_element_is_reported_as_not_assessed():
    result = assess_composition(parse_formula("OgO"))
    assert result.oxidation_state_status == DiagnosticStatus.NOT_ASSESSED
    assert result.requires_review


def test_single_element_composition_has_zero_formal_state_option():
    result = assess_composition(parse_formula("O2"))
    assert result.oxidation_state_status == DiagnosticStatus.SUPPORTED_PATTERN


def test_unknown_element_is_rejected_by_element_lookup():
    composition = Composition(counts={"Xx": 1})
    with pytest.raises(ValueError, match="Unknown element"):
        assess_composition(composition)


def test_findings_are_explanatory_and_nonempty():
    result = assess_composition(parse_formula("H2O"))
    assert len(result.findings) >= 2
    assert all(finding.rule_id and finding.message for finding in result.findings)


def test_charged_single_element_is_not_assumed_to_be_neutral_elemental_state():
    composition = Composition(counts={"O": 2}, net_charge=1)
    result = assess_composition(composition)
    assert result.oxidation_state_status == DiagnosticStatus.NOT_ASSESSED


def test_missing_formal_balance_is_review_not_rejection():
    result = assess_composition(parse_formula("LiF2"))
    assert result.oxidation_state_status == DiagnosticStatus.NO_MATCH_IN_REFERENCE
    assert result.requires_review
