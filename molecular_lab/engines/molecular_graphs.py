"""
Restricted acyclic molecular graph proposal engines.

Version: 0.1
Rev.: 01
Status: Acyclic, unicyclic, and bounded linear fused-acene graph prototypes.

The generators propose connectivity only. They do not establish geometry,
energetics, stability, experimental existence, or synthesis feasibility.
"""

from __future__ import annotations

from itertools import permutations, product
from math import factorial
from typing import Iterator

from molecular_lab.chemistry.elements import get_element
from molecular_lab.core.models import Composition, MolecularAtom, MolecularBond, MolecularGraph

_SUPPORTED_VALENCES = {"H": 1, "C": 4, "N": 3, "O": 2, "F": 1, "Cl": 1}
_LIMITATIONS = (
    "Neutral species only; charged and radical states are not generated.",
    "Fixed common valences are used for H, C, N, O, F, and Cl only.",
    "Acyclic, unicyclic, and linear fused-carbon acene skeletons are supported by separate generators; other fused, bridged, spiro, and general multi-ring systems remain outside scope.",
    "Aromaticity is not an explicit graph bond type; optional RDKit perception is limited to supported Kekule graphs. General resonance enumeration, hypervalence, coordination chemistry, and variable valence are outside scope.",
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


def _canonical_signature(
    elements: tuple[str, ...],
    edges: tuple[tuple[int, int, int], ...],
    hydrogen_counts: tuple[int, ...],
) -> tuple:
    """Canonicalize a labeled heavy-atom tree under same-element permutations."""
    groups = [tuple(i for i, symbol in enumerate(elements) if symbol == element)
              for element in sorted(set(elements))]
    permutation_count = 1
    for group in groups:
        permutation_count *= factorial(len(group))
    if permutation_count > 50000:
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
        new_elements = [""] * len(elements)
        new_hydrogens = [0] * len(elements)
        for old_index, new_index in old_to_new.items():
            new_elements[new_index] = elements[old_index]
            new_hydrogens[new_index] = hydrogen_counts[old_index]
        new_edges = tuple(sorted(
            (min(old_to_new[a], old_to_new[b]),
             max(old_to_new[a], old_to_new[b]), order)
            for a, b, order in edges
        ))
        signature = (tuple(new_elements), new_edges, tuple(new_hydrogens))
        if best is None or signature < best:
            best = signature
    return best


def _validate_limits(max_candidates: int, max_heavy_atoms: int) -> None:
    for name, value in (("max_candidates", max_candidates),
                        ("max_heavy_atoms", max_heavy_atoms)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")


def _build_graph(
    heavy_elements: tuple[str, ...],
    heavy_edges: tuple[tuple[int, int, int], ...],
    attached_h: tuple[int, ...],
    method: str,
) -> MolecularGraph:
    atoms = [MolecularAtom(i, element) for i, element in enumerate(heavy_elements)]
    bonds = [MolecularBond(a, b, order) for a, b, order in heavy_edges]
    next_id = len(heavy_elements)
    for heavy_index, hydrogen_count in enumerate(attached_h):
        for _ in range(hydrogen_count):
            atoms.append(MolecularAtom(next_id, "H"))
            bonds.append(MolecularBond(heavy_index, next_id, 1))
            next_id += 1
    return MolecularGraph(
        atoms=tuple(atoms),
        bonds=tuple(bonds),
        method=method,
        limitations=_LIMITATIONS,
    )


def _heavy_elements(composition: Composition, max_heavy_atoms: int) -> tuple[str, ...] | None:
    if composition.net_charge != 0:
        return None
    for symbol in composition.counts:
        get_element(symbol)
        if symbol not in _SUPPORTED_VALENCES:
            return None
    elements = tuple(
        symbol for symbol in sorted(composition.counts) if symbol != "H"
        for _ in range(composition.counts[symbol])
    )
    if len(elements) > max_heavy_atoms:
        return None
    return elements


def iter_saturated_acyclic_graphs(
    composition: Composition,
    *,
    max_candidates: int = 1000,
    max_heavy_atoms: int = 8,
) -> Iterator[MolecularGraph]:
    """Yield neutral connected acyclic proposals with single heavy-atom bonds."""
    _validate_limits(max_candidates, max_heavy_atoms)
    elements = _heavy_elements(composition, max_heavy_atoms)
    if elements is None:
        return
    if not elements:
        if composition.counts.get("H") == 2:
            yield MolecularGraph(
                atoms=(MolecularAtom(0, "H"), MolecularAtom(1, "H")),
                bonds=(MolecularBond(0, 1),),
                method="bounded_saturated_acyclic_tree_v1",
                limitations=_LIMITATIONS,
            )
        return

    target_h = composition.counts.get("H", 0)
    emitted = 0
    seen: set[tuple] = set()
    topologies = _prufer_trees(len(elements)) if len(elements) > 1 else ((),)
    for topology in topologies:
        degrees = [0] * len(elements)
        for a, b in topology:
            degrees[a] += 1
            degrees[b] += 1
        if any(degrees[i] > _SUPPORTED_VALENCES[e] for i, e in enumerate(elements)):
            continue
        attached_h = tuple(_SUPPORTED_VALENCES[e] - degrees[i]
                           for i, e in enumerate(elements))
        if sum(attached_h) != target_h:
            continue
        signature = _canonical_signature(
            elements, tuple((a, b, 1) for a, b in topology), attached_h
        )
        if signature in seen:
            continue
        seen.add(signature)
        yield _build_graph(
            elements, tuple((a, b, 1) for a, b in topology),
            attached_h, "bounded_saturated_acyclic_tree_v1"
        )
        emitted += 1
        if emitted >= max_candidates:
            return


def iter_acyclic_graphs(
    composition: Composition,
    *,
    max_candidates: int = 1000,
    max_heavy_atoms: int = 8,
    max_assignments: int = 100000,
) -> Iterator[MolecularGraph]:
    """
    Yield neutral acyclic connectivity proposals with single, double, or triple
    heavy-atom bonds. Each heavy-atom skeleton is a tree; explicit hydrogens
    fill unused valence. max_assignments bounds bond-order combinations examined.
    """
    _validate_limits(max_candidates, max_heavy_atoms)
    if isinstance(max_assignments, bool) or not isinstance(max_assignments, int) or max_assignments < 1:
        raise ValueError("max_assignments must be a positive integer")
    elements = _heavy_elements(composition, max_heavy_atoms)
    if elements is None:
        return
    if not elements:
        if composition.counts.get("H") == 2:
            yield MolecularGraph(
                atoms=(MolecularAtom(0, "H"), MolecularAtom(1, "H")),
                bonds=(MolecularBond(0, 1),),
                method="bounded_acyclic_bond_order_tree_v1",
                limitations=_LIMITATIONS,
            )
        return

    target_h = composition.counts.get("H", 0)
    emitted = 0
    examined = 0
    seen: set[tuple] = set()
    topologies = _prufer_trees(len(elements)) if len(elements) > 1 else ((),)
    for topology in topologies:
        order_options = product((1, 2, 3), repeat=len(topology))
        for orders in order_options:
            examined += 1
            if examined > max_assignments:
                return
            used = [0] * len(elements)
            weighted_edges = tuple(
                (a, b, order) for (a, b), order in zip(topology, orders)
            )
            for a, b, order in weighted_edges:
                used[a] += order
                used[b] += order
            if any(used[i] > _SUPPORTED_VALENCES[element]
                   for i, element in enumerate(elements)):
                continue
            attached_h = tuple(_SUPPORTED_VALENCES[element] - used[i]
                               for i, element in enumerate(elements))
            if sum(attached_h) != target_h:
                continue
            signature = _canonical_signature(elements, weighted_edges, attached_h)
            if signature in seen:
                continue
            seen.add(signature)
            yield _build_graph(
                elements, weighted_edges, attached_h,
                "bounded_acyclic_bond_order_tree_v1"
            )
            emitted += 1
            if emitted >= max_candidates:
                return


def iter_unicyclic_graphs(
    composition: Composition,
    *,
    max_candidates: int = 1000,
    max_heavy_atoms: int = 8,
    max_assignments: int = 100000,
) -> Iterator[MolecularGraph]:
    """
    Yield neutral connected heavy-atom graphs containing exactly one independent
    cycle. The heavy-atom skeleton is constructed by adding one chord to a tree;
    single, double, and triple bond orders are considered under fixed valences.

    This covers simple rings and ring-plus-acyclic-substituent topologies, but
    not fused, bridged, or multi-ring systems.
    """
    _validate_limits(max_candidates, max_heavy_atoms)
    if isinstance(max_assignments, bool) or not isinstance(max_assignments, int) or max_assignments < 1:
        raise ValueError("max_assignments must be a positive integer")
    elements = _heavy_elements(composition, max_heavy_atoms)
    if elements is None or len(elements) < 3:
        return

    target_h = composition.counts.get("H", 0)
    emitted = 0
    examined = 0
    seen: set[tuple] = set()

    # Enumerate simple cycles first. This ensures compact rings (including
    # Kekule representations of aromatic systems) are not hidden behind the
    # much larger spanning-tree-plus-chord search space.
    for ring_elements in sorted(set(permutations(elements))):
        ring_topology = tuple(
            (i, (i + 1) % len(ring_elements))
            for i in range(len(ring_elements))
        )
        for orders in product((1, 2, 3), repeat=len(ring_topology)):
            examined += 1
            if examined > max_assignments:
                return
            used = [0] * len(ring_elements)
            weighted_edges = tuple(
                (a, b, order)
                for (a, b), order in zip(ring_topology, orders)
            )
            for a, b, order in weighted_edges:
                used[a] += order
                used[b] += order
            if any(
                used[i] > _SUPPORTED_VALENCES[element]
                for i, element in enumerate(ring_elements)
            ):
                continue
            attached_h = tuple(
                _SUPPORTED_VALENCES[element] - used[i]
                for i, element in enumerate(ring_elements)
            )
            if sum(attached_h) != target_h:
                continue
            signature = _canonical_signature(
                ring_elements, weighted_edges, attached_h
            )
            if signature in seen:
                continue
            seen.add(signature)
            yield _build_graph(
                ring_elements,
                weighted_edges,
                attached_h,
                "bounded_unicyclic_bond_order_graph_v1",
            )
            emitted += 1
            if emitted >= max_candidates:
                return

    # A connected unicyclic graph is a spanning tree plus one non-tree edge.
    # Repeated constructions are removed by the element-aware canonical key.
    for tree in _prufer_trees(len(elements)):
        tree_pairs = {tuple(sorted(edge)) for edge in tree}
        missing_edges = tuple(
            (a, b)
            for a in range(len(elements))
            for b in range(a + 1, len(elements))
            if (a, b) not in tree_pairs
        )
        for chord in missing_edges:
            topology = tuple(sorted((*tree, chord)))
            order_options = product((1, 2, 3), repeat=len(topology))
            for orders in order_options:
                examined += 1
                if examined > max_assignments:
                    return

                used = [0] * len(elements)
                weighted_edges = tuple(
                    (a, b, order) for (a, b), order in zip(topology, orders)
                )
                for a, b, order in weighted_edges:
                    used[a] += order
                    used[b] += order
                if any(
                    used[i] > _SUPPORTED_VALENCES[element]
                    for i, element in enumerate(elements)
                ):
                    continue

                attached_h = tuple(
                    _SUPPORTED_VALENCES[element] - used[i]
                    for i, element in enumerate(elements)
                )
                if sum(attached_h) != target_h:
                    continue

                signature = _canonical_signature(
                    elements, weighted_edges, attached_h
                )
                if signature in seen:
                    continue
                seen.add(signature)

                yield _build_graph(
                    elements,
                    weighted_edges,
                    attached_h,
                    "bounded_unicyclic_bond_order_graph_v1",
                )
                emitted += 1
                if emitted >= max_candidates:
                    return



def _linear_acene_topology(ring_count: int) -> tuple[tuple[int, int], ...]:
    """Build the carbon skeleton of a linearly fused sequence of six-membered rings."""
    if ring_count < 2:
        raise ValueError("ring_count must be at least 2")

    edges: set[tuple[int, int]] = {
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5)
    }
    shared_edge = (2, 3)
    next_atom = 6

    for _ in range(1, ring_count):
        left, right = shared_edge
        new_path = (left, next_atom, next_atom + 1, next_atom + 2,
                    next_atom + 3, right)
        for a, b in zip(new_path, new_path[1:]):
            edges.add((min(a, b), max(a, b)))
        # The next ring shares the edge opposite the previous shared edge.
        shared_edge = (next_atom + 1, next_atom + 2)
        next_atom += 4

    return tuple(sorted(edges))


def _perfect_matchings(
    node_count: int,
    edges: tuple[tuple[int, int], ...],
    *,
    max_assignments: int,
) -> Iterator[tuple[tuple[int, int], ...]]:
    """Yield perfect matchings, bounding recursive search nodes examined."""
    adjacency = {node: set() for node in range(node_count)}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)

    examined = 0

    def search(
        remaining: frozenset[int],
        chosen: tuple[tuple[int, int], ...],
    ) -> Iterator[tuple[tuple[int, int], ...]]:
        nonlocal examined
        examined += 1
        if examined > max_assignments:
            return
        if not remaining:
            yield tuple(sorted(chosen))
            return

        # Minimum-remaining-values ordering reduces branching and is stable.
        node = min(
            remaining,
            key=lambda item: (len(adjacency[item] & remaining), item),
        )
        for neighbor in sorted(adjacency[node] & remaining):
            yield from search(
                remaining - {node, neighbor},
                chosen + ((min(node, neighbor), max(node, neighbor)),),
            )

    yield from search(frozenset(range(node_count)), ())


def iter_linear_fused_acene_graphs(
    composition: Composition,
    *,
    max_candidates: int = 100,
    max_heavy_atoms: int = 16,
    max_assignments: int = 100000,
) -> Iterator[MolecularGraph]:
    """
    Yield Kekule graph representations for neutral linearly fused benzene-ring
    systems (acenes), beginning with naphthalene and anthracene.

    This is a deliberately narrow multi-ring prototype: carbon-only heavy-atom
    skeletons, consecutive six-membered rings sharing one edge, and one double
    bond incident to each carbon in every Kekule assignment. Different yielded
    Kekule assignments are resonance representations, not distinct constitutional
    isomers or independent evidence of stability.
    """
    _validate_limits(max_candidates, max_heavy_atoms)
    if isinstance(max_assignments, bool) or not isinstance(max_assignments, int) or max_assignments < 1:
        raise ValueError("max_assignments must be a positive integer")
    if composition.net_charge != 0:
        return
    if set(composition.counts) - {"C", "H"}:
        return

    carbon_count = composition.counts.get("C", 0)
    hydrogen_count = composition.counts.get("H", 0)
    if carbon_count < 10 or (carbon_count - 2) % 4:
        return
    ring_count = (carbon_count - 2) // 4
    if ring_count < 2 or carbon_count > max_heavy_atoms:
        return
    if hydrogen_count != 2 * ring_count + 4:
        return

    edges = _linear_acene_topology(ring_count)
    if len(edges) - carbon_count + 1 != ring_count:
        return

    degree = [0] * carbon_count
    for a, b in edges:
        degree[a] += 1
        degree[b] += 1
    if any(value not in (2, 3) for value in degree):
        return
    attached_h = tuple(3 - value for value in degree)
    if sum(attached_h) != hydrogen_count:
        return

    emitted = 0
    for matching in _perfect_matchings(
        carbon_count, edges, max_assignments=max_assignments
    ):
        matched = set(matching)
        weighted_edges = tuple(
            (a, b, 2 if (a, b) in matched else 1) for a, b in edges
        )
        yield _build_graph(
            tuple("C" for _ in range(carbon_count)),
            weighted_edges,
            attached_h,
            "bounded_linear_fused_acene_kekule_v1",
        )
        emitted += 1
        if emitted >= max_candidates:
            return
