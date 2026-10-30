"""
Tests for the offline reference-molecule catalogue.

Version: 0.1
Rev.: 00
"""

from molecular_lab.database.reference_molecules import (
    MatchLevel,
    all_reference_molecules,
    identify_candidate,
    lookup_by_formula,
    lookup_by_name,
)


def test_seed_catalogue_loads_with_unique_records():
    records = all_reference_molecules()
    assert len(records) == 5
    assert len({record.record_id for record in records}) == len(records)


def test_formula_lookup_returns_all_c2h6o_isomers():
    matches = lookup_by_formula("C2H6O")
    assert {record.preferred_name for record in matches} == {
        "Ethanol",
        "Dimethyl ether",
    }


def test_formula_lookup_finds_water_and_peroxide():
    assert lookup_by_formula("H2O")[0].preferred_name == "Water"
    assert lookup_by_formula("H2O2")[0].preferred_name == "Hydrogen peroxide"


def test_name_lookup_accepts_synonyms_case_insensitively():
    assert lookup_by_name("METHOXYMETHANE")[0].preferred_name == "Dimethyl ether"
    assert lookup_by_name("ethyl alcohol")[0].preferred_name == "Ethanol"


def test_inchikey_establishes_exact_catalogue_identity():
    result = identify_candidate(
        "C2H6O",
        inchi_key="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    )
    assert result.match_level is MatchLevel.EXACT_STRUCTURE
    assert [record.preferred_name for record in result.matches] == ["Ethanol"]


def test_formula_only_match_does_not_confuse_isomers():
    result = identify_candidate("C2H6O")
    assert result.match_level is MatchLevel.FORMULA_ONLY
    assert len(result.matches) == 2


def test_unknown_formula_is_not_claimed_to_be_a_new_molecule():
    result = identify_candidate("C2H4O")
    assert result.match_level is MatchLevel.NOT_FOUND
    assert "does not prove" in result.explanation


def test_incorrect_inchikey_falls_back_to_formula_only():
    result = identify_candidate("C2H6O", inchi_key="NOT-A-REAL-KEY")
    assert result.match_level is MatchLevel.FORMULA_ONLY
    assert len(result.matches) == 2
