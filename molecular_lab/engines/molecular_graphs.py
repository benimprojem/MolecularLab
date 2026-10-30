"""
Restricted saturated acyclic molecular graph proposal engine.

Version: 0.1
Rev.: 00
Status: Prototype; explicit scope and limitations are returned with each graph.

This module proposes single-bond, connected, acyclic graphs for a small set of
neutral main-group elements. It is not a universal molecule generator.
"""

from __future__ import annotations

from itertools import permutations, product
from math import factorial
from typing import Iterator

from molecular_lab.chemistry.elements import get_element
from molecular_lab.core.models import (
    Composition,
    MolecularAtom,
    MolecularBond,
    MolecularGraph,
)

# This is a deliberately narrow, fixed-valence model for the first prototype.
# It does not cover radicals, hypervalence, aromaticity, multiple bonds,
# coordination chemistry, charged species, or variable-valence cases.
_SUPPORTED_VALENCES = {"H": 1, "C": 4, "N": 3, "O": 2, "F": 1, "Cl": 1}
_METHOD = "bounded_saturated_acyclic_tree_v1"
_LIMITATIONS = (
    "Neutral species only; charged and radical states are not generated.",
    "Single bonds only; double bonds, triple bonds, aromaticity, and resonance are outside scope.",
    "Acyclic connected graphs only; rings and cages are outside scope.",
    "Fixed common valences are used for H, C, N, O, F, and Cl only.",
    "A graph proposal is not a 3D geometry, stability result, or synthesis prediction.",
)


def _prufer_trees(node_count: int) -> Iterator[tuple[tuple[int, int], ...]]:
    """Yield labeled trees using Prüfer sequences; node_count must be >= 2."""
    if node_count == 2:
        yield ((0, 1),)
        return
    for sequence in product(range(node_count), repeat=node_count - 2):
        degrees = [1] * node_count
        for node in sequence:
            degrees[node] += 1
        edges: list[tuple[int, int]] = []
        for node in sequence:
            leaf = next(i for i, degree in enumerate(degrees) if degree == 1)
            edges.append((min(leaf, node), max(leaf, node)))
            degrees[leaf] -= 1
            degrees[node] -= 1
        left = [i for i, degree in enumerate(degrees) if degree == 1]
        edges.append((min(left), max(left)))
        yield tuple(sorted(edges))


def _canonical_signature(elements: tuple[str, ...], edges: tuple[tuple[int, int], ...],
                          hydrogen_counts: tuple[int, ...]) -> tuple:
    """Canonicalize a labeled heavy-atom tree under same-element permutations."""
    groups = [tuple(i for i, symbol in enumerate(elements) if symbol == element)
              for element in sorted(set(elements))]
    permutation_count = 1
    for group in groups:
        permutation_count *= factorial(len(group))
    if permutation_count > 50000:
        # The caller's candidate cap still bounds work; retain deterministic
        # representation rather than attempting an expensive canonical search.
        return (elements, edges, hydrogen_counts)

    best = None
    group_permutations = [tuple(permutations(group)) for group in groups]
    for choices in product(*group_permutations):
        old_to_new: dict[int, int] = {}
        offset = 0
        for group, ordered in zip(groups, choices):
            for local_index, old_index in enumerate(ordered):
                old_to_new[old_index] = offset + local_index
            offset += len(group)
        # Rebuild element and attached-H arrays in canonical index order.
        new_elements = [""] * len(elements)
        new_hydrogens = [0] * len(elements)
        for old_index, new_index in old_to_new.items():
            new_elements[new_index] = elements[old_index]
            new_hydrogens[new_index] = hydrogen_counts[old_index]
        new_edges = tuple(sorted(
            (min(old_to_new[a], old_to_new[b]), max(old_to_new[a], old_to_new[b]))
            for a, b in edges
        ))
        signature = (tuple(new_elements), new_edges, tuple(new_hydrogens))
        if best is None or signature < best:
            best = signature
    return best


def iter_saturated_acyclic_graphs(
    composition: Composition,
    *,
    max_candidates: int = 1000,
    max_heavy_atoms: int = 8,
) -> Iterator[MolecularGraph]:
    """
    Yield bounded neutral saturated acyclic graph proposals for a composition.

    Hydrogen atoms are explicit in the returned graph. The heavy-atom
    connectivity is a tree; hydrogens fill the remaining fixed valences.
    """
    if isinstance(max_candidates, bool) or max_candidates < 1:
        raise ValueError("max_candidates must be a positive integer")
    if isinstance(max_heavy_atoms, bool) or max_heavy_atoms < 1:
        raise ValueError("max_heavy_atoms must be a positive integer")
    if composition.net_charge != 0:
        return

    for symbol in composition.counts:
        get_element(symbol)  # Validate against the periodic-table identity data.
        if symbol not in _SUPPORTED_VALENCES:
            return

    heavy_elements = tuple(
        symbol
        for symbol in sorted(composition.counts)
        if symbol != "H"
        for _ in range(composition.counts[symbol])
    )
    if not heavy_elements:
        # The first prototype recognizes only the simplest diatomic hydrogen.
        if composition.counts.get("H") == 2:
            yield MolecularGraph(
                atoms=(MolecularAtom(0, "H"), MolecularAtom(1, "H")),
                bonds=(MolecularBond(0, 1),),
                method=_METHOD,
                limitations=_LIMITATIONS,
            )
        return

    heavy_count = len(heavy_elements)
    if heavy_count > max_heavy_atoms:
        return

    hydrogen_target = composition.counts.get("H", 0)
    emitted = 0
    seen: set[tuple] = set()

    for heavy_edges in _prufer_trees(heavy_count) if heavy_count > 1 else ((),):
        degrees = [0] * heavy_count
        for a, b in heavy_edges:
            degrees[a] += 1
            degrees[b] += 1
        if any(degrees[i] > _SUPPORTED_VALENCES[element]
               for i, element in enumerate(heavy_elements)):
            continue

        attached_h = tuple(
            _SUPPORTED_VALENCES[element] - degrees[i]
            for i, element in enumerate(heavy_elements)
        )
        if sum(attached_h) != hydrogen_target:
            continue

        signature = _canonical_signature(heavy_elements, heavy_edges, attached_h)
        if signature in seen:
            continue
        seen.add(signature)

        atoms = [MolecularAtom(i, element) for i, element in enumerate(heavy_elements)]
        bonds = [MolecularBond(a, b, 1) for a, b in heavy_edges]
        next_id = heavy_count
        for heavy_index, hydrogen_count in enumerate(attached_h):
            for _ in range(hydrogen_count):
                atoms.append(MolecularAtom(next_id, "H"))
                bonds.append(MolecularBond(heavy_index, next_id, 1))
                next_id += 1

        yield MolecularGraph(
            atoms=tuple(atoms),
            bonds=tuple(bonds),
            method=_METHOD,
            limitations=_LIMITATIONS,
        )
        emitted += 1
        if emitted >= max_candidates:
            return
