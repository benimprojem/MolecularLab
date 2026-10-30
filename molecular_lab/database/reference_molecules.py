"""
Local reference-molecule lookup and identity matching.

Version: 0.1
Rev.: 00
Status: Initial offline reference catalogue; no universal structure canonicalizer.

Formula matches identify possible known records, not exact molecular identity.
Exact identity matching currently requires an externally supplied InChIKey.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
from importlib.resources import files

from molecular_lab.chemistry.formula import parse_formula
from molecular_lab.core.models import Composition


class MatchLevel(str, Enum):
    EXACT_STRUCTURE = "exact_structure"
    FORMULA_ONLY = "formula_only"
    NOT_FOUND = "not_found"
    NOT_APPLICABLE = "not_applicable"


@dataclass(frozen=True, slots=True)
class ReferenceMolecule:
    record_id: str
    preferred_name: str
    formula: str
    synonyms: tuple[str, ...]
    pubchem_cid: int
    canonical_smiles: str
    inchi_key: str
    cas_rn: str | None
    source: str
    source_url: str
    retrieved_on: str


@dataclass(frozen=True, slots=True)
class MoleculeLookupResult:
    match_level: MatchLevel
    matches: tuple[ReferenceMolecule, ...]
    explanation: str


@lru_cache(maxsize=1)
def all_reference_molecules() -> tuple[ReferenceMolecule, ...]:
    """Load and validate the bundled offline reference catalogue."""
    path = files("molecular_lab.database").joinpath("reference_molecules.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != "0.1":
        raise RuntimeError("Unsupported reference-molecule schema version")

    records = tuple(
        ReferenceMolecule(
            record_id=item["record_id"],
            preferred_name=item["preferred_name"],
            formula=item["formula"],
            synonyms=tuple(item.get("synonyms", ())),
            pubchem_cid=int(item["pubchem_cid"]),
            canonical_smiles=item["canonical_smiles"],
            inchi_key=item["inchi_key"],
            cas_rn=item.get("cas_rn"),
            source=item["source"],
            source_url=item["source_url"],
            retrieved_on=item["retrieved_on"],
        )
        for item in payload["molecules"]
    )

    record_ids = [record.record_id for record in records]
    cids = [record.pubchem_cid for record in records]
    if len(record_ids) != len(set(record_ids)):
        raise RuntimeError("Reference-molecule record IDs must be unique")
    if len(cids) != len(set(cids)):
        raise RuntimeError("PubChem CIDs must be unique in this seed catalogue")
    for record in records:
        parsed = parse_formula(record.formula)
        if not record.preferred_name.strip() or not record.inchi_key.strip():
            raise RuntimeError(f"Incomplete reference record: {record.record_id}")
        if not record.source_url.startswith("https://"):
            raise RuntimeError(f"Invalid source URL: {record.record_id}")
        if not parsed.counts:
            raise RuntimeError(f"Invalid formula in reference record: {record.record_id}")
    return records


def lookup_by_formula(
    formula_or_composition: str | Composition,
) -> tuple[ReferenceMolecule, ...]:
    """Return all records with the same element counts, preserving isomers."""
    composition = (
        parse_formula(formula_or_composition)
        if isinstance(formula_or_composition, str)
        else formula_or_composition
    )
    if composition.net_charge != 0:
        return ()
    target = dict(composition.counts)
    return tuple(
        record
        for record in all_reference_molecules()
        if dict(parse_formula(record.formula).counts) == target
    )


def lookup_by_name(name: str) -> tuple[ReferenceMolecule, ...]:
    """Case-insensitive lookup by preferred name or recorded synonym."""
    normalized = name.strip().casefold()
    if not normalized:
        return ()
    return tuple(
        record
        for record in all_reference_molecules()
        if normalized == record.preferred_name.casefold()
        or normalized in {synonym.casefold() for synonym in record.synonyms}
    )


def identify_candidate(
    formula_or_composition: str | Composition,
    *,
    inchi_key: str | None = None,
) -> MoleculeLookupResult:
    """
    Identify a candidate against the local catalogue.

    An exact structure match is reported only when InChIKeys match. Formula-only
    matches are deliberately returned as a set because multiple isomers can
    share the same elemental composition.
    """
    formula_matches = lookup_by_formula(formula_or_composition)
    if inchi_key:
        exact = tuple(
            record for record in formula_matches
            if record.inchi_key.casefold() == inchi_key.strip().casefold()
        )
        if exact:
            return MoleculeLookupResult(
                MatchLevel.EXACT_STRUCTURE,
                exact,
                "The supplied InChIKey matches a local reference record.",
            )
    if formula_matches:
        return MoleculeLookupResult(
            MatchLevel.FORMULA_ONLY,
            formula_matches,
            "The formula exists in the reference catalogue, but exact structure identity is not established.",
        )
    return MoleculeLookupResult(
        MatchLevel.NOT_FOUND,
        (),
        "No matching formula was found in the local reference catalogue. This does not prove the molecule is new or unknown.",
    )
