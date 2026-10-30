"""
Tests for periodic-table lookup and formula parsing.

Version: 0.1
Rev.: 00
"""

import pytest

from molecular_lab.chemistry.elements import (
    all_elements,
    find_element_by_name,
    get_element,
    get_element_by_atomic_number,
)
from molecular_lab.chemistry.formula import FormulaParseError, parse_formula


def test_periodic_table_contains_118_elements():
    elements = all_elements()
    assert len(elements) == 118
    assert elements[0].symbol == "H"
    assert elements[-1].symbol == "Og"


def test_element_lookup_by_symbol_number_and_name():
    assert get_element("Fe").atomic_number == 26
    assert get_element_by_atomic_number(8).symbol == "O"
    assert find_element_by_name("oxygen").symbol == "O"


@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        ("H2O", {"H": 2, "O": 1}),
        ("H2O2", {"H": 2, "O": 2}),
        ("CO", {"C": 1, "O": 1}),
        ("Co", {"Co": 1}),
        ("Fe2O3", {"Fe": 2, "O": 3}),
        ("Fe2(SO4)3", {"Fe": 2, "S": 3, "O": 12}),
        ("K4[Fe(CN)6]", {"K": 4, "Fe": 1, "C": 6, "N": 6}),
        ("Al2(SO4)3", {"Al": 2, "S": 3, "O": 12}),
    ],
)
def test_formula_parser(formula, expected):
    composition = parse_formula(formula)
    assert dict(composition.counts) == expected
    assert composition.input_formula == formula
    assert composition.total_atoms == sum(expected.values())


@pytest.mark.parametrize(
    "formula",
    [
        "",
        "H 2O",
        "h2o",
        "Xx2",
        "H0",
        "H-2",
        "2H2O",
        "H(OH",
        "H)O",
        "Mg(OH]2",
        "H()",
        "CuSO4·5H2O",
        "NH4+",
    ],
)
def test_formula_parser_rejects_unsupported_or_invalid_syntax(formula):
    with pytest.raises(FormulaParseError):
        parse_formula(formula)
