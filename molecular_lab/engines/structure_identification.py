"""
Identify restricted molecular graph proposals against the local reference catalogue.

Version: 0.1
Rev.: 01
Status: Optional RDKit adapter with aromaticity reporting; formula-only fallback remains available.

RDKit is imported lazily. If it is unavailable, graph composition is still
checked against the catalogue, but exact structure identity is not claimed.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Any

from molecular_lab.core.models import Composition, MolecularGraph
from molecular_lab.database.reference_molecules import (
    MatchLevel,
    ReferenceMolecule,
    identify_candidate,
)


@dataclass(frozen=True, slots=True)
class GraphIdentificationResult:
    """Reference-catalogue result plus the identifiers actually computed."""

    formula: str
    canonical_smiles: str | None
    inchi_key: str | None
    backend: str | None
    status: str
    match_level: MatchLevel
    matches: tuple[ReferenceMolecule, ...]
    explanation: str
    aromaticity_status: str = "not_assessed"


def _load_rdkit() -> tuple[Any, Any] | None:
    """Load RDKit only when this adapter is used."""
    try:
        from rdkit import Chem
        from rdkit.Chem import inchi
    except ImportError:
        return None
    return Chem, inchi


def _format_formula(counts: Counter[str]) -> str:
    """Format a simple neutral molecular formula using Hill ordering."""
    symbols = list(counts)
    if "C" in counts:
        ordered = ["C"]
        if "H" in counts:
            ordered.append("H")
        ordered.extend(sorted(symbol for symbol in symbols if symbol not in {"C", "H"}))
    else:
        ordered = sorted(symbols)
    return "".join(
        symbol + (str(counts[symbol]) if counts[symbol] != 1 else "")
        for symbol in ordered
    )


def _graph_to_identifiers(
    graph: MolecularGraph,
    chem: Any,
    inchi: Any,
) -> tuple[str, str]:
    """Convert an explicit graph to canonical SMILES and an InChIKey."""
    editable = chem.RWMol()
    atom_indices: dict[int, int] = {}
    periodic_table = chem.GetPeriodicTable()

    for atom in graph.atoms:
        atomic_number = periodic_table.GetAtomicNumber(atom.element)
        if atomic_number <= 0:
            raise ValueError(f"RDKit does not recognize element {atom.element!r}")
        atom_indices[atom.atom_id] = editable.AddAtom(chem.Atom(atomic_number))

    bond_types = {
        1: chem.BondType.SINGLE,
        2: chem.BondType.DOUBLE,
        3: chem.BondType.TRIPLE,
    }
    for bond in graph.bonds:
        editable.AddBond(
            atom_indices[bond.atom_a],
            atom_indices[bond.atom_b],
            bond_types[bond.order],
        )

    molecule = editable.GetMol()
    chem.SanitizeMol(molecule)
    smiles_molecule = chem.RemoveHs(molecule)
    smiles = chem.MolToSmiles(smiles_molecule, canonical=True)
    inchi_key = inchi.MolToInchiKey(molecule)
    aromatic_bond_count = sum(1 for bond in molecule.GetBonds() if bond.GetIsAromatic())
    aromatic_atom_count = sum(1 for atom in molecule.GetAtoms() if atom.GetIsAromatic())
    if not smiles or not inchi_key:
        raise ValueError("RDKit did not produce a canonical identifier")
    aromaticity = "aromatic_detected" if aromatic_atom_count and aromatic_bond_count else "no_aromaticity_detected"
    return smiles, inchi_key, aromaticity


def identify_graph_against_catalogue(
    graph: MolecularGraph,
) -> GraphIdentificationResult:
    """
    Derive the graph formula and identify it against the offline catalogue.

    Without RDKit, the result is limited to formula matching. Conversion or
    sanitization errors also fall back to formula-level evidence; they never
    become an exact match.
    """
    counts = Counter(atom.element for atom in graph.atoms)
    composition = Composition(counts=counts, net_charge=0)
    formula = _format_formula(counts)
    rdkit = _load_rdkit()

    if rdkit is None:
        lookup = identify_candidate(composition)
        return GraphIdentificationResult(
            formula=formula,
            canonical_smiles=None,
            inchi_key=None,
            backend=None,
            status="rdkit_unavailable",
            match_level=lookup.match_level,
            matches=lookup.matches,
            explanation=(
                "RDKit is not installed. Only formula-level catalogue matching "
                "was performed; exact structure identity is undetermined. "
                + lookup.explanation
            ),
        )

    chem, inchi = rdkit
    try:
        smiles, inchi_key, aromaticity_status = _graph_to_identifiers(graph, chem, inchi)
    except Exception as exc:
        lookup = identify_candidate(composition)
        return GraphIdentificationResult(
            formula=formula,
            canonical_smiles=None,
            inchi_key=None,
            backend="RDKit",
            status="structure_conversion_failed",
            match_level=lookup.match_level,
            matches=lookup.matches,
            explanation=(
                f"RDKit could not validate or identify this graph ({type(exc).__name__}: {exc}). "
                "No exact structure match is claimed. "
                + lookup.explanation
            ),
        )

    lookup = identify_candidate(composition, inchi_key=inchi_key)
    return GraphIdentificationResult(
        formula=formula,
        canonical_smiles=smiles,
        inchi_key=inchi_key,
        backend="RDKit",
        status="identified",
        match_level=lookup.match_level,
        matches=lookup.matches,
        explanation=lookup.explanation,
        aromaticity_status=aromaticity_status,
    )
