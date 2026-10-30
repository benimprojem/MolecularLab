"""
Bounded composition enumeration.

Version: 0.1
Rev.: 01
Status: First deterministic enumerator; no chemical plausibility scoring.

This module enumerates element-count combinations inside explicit user bounds.
It does not generate molecular bonds or crystal structures and does not assess
existence, stability, synthesizability, or chemical plausibility.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from molecular_lab.chemistry.elements import get_element
from molecular_lab.core.models import (
    Composition,
    CompositionKind,
    SearchConstraints,
    SearchDomain,
)


def iter_compositions(
    constraints: SearchConstraints,
    *,
    required_elements: Iterable[str] = (),
) -> Iterator[Composition]:
    """
    Yield bounded element-count compositions in deterministic order.

    Selected elements are optional by default. Pass required_elements when
    every listed element must occur.

    Counts are enumerated by increasing total atom count, then by sorted
    element symbols and ascending counts. If no charge list is supplied,
    neutral compositions are emitted. For a mixed domain, both molecular and
    formula-unit interpretations are emitted; neither is a structural claim.
    """
    required = frozenset(required_elements)
    unknown_required = required - constraints.allowed_elements
    if unknown_required:
        names = ", ".join(sorted(unknown_required))
        raise ValueError(f"required elements are not allowed: {names}")

    for symbol in constraints.allowed_elements:
        try:
            get_element(symbol)
        except ValueError as exc:
            raise ValueError(f"unknown allowed element: {symbol}") from exc

    elements = tuple(sorted(constraints.allowed_elements))
    limits = tuple(
        min(
            constraints.per_element_max.get(symbol, constraints.max_total_atoms),
            constraints.max_total_atoms,
        )
        for symbol in elements
    )
    charges = constraints.allowed_charges
    if charges is None:
        charges = (0,)
    else:
        charges = tuple(sorted(set(charges)))

    if constraints.domain == SearchDomain.MOLECULE:
        kinds = (CompositionKind.MOLECULAR,)
    elif constraints.domain == SearchDomain.PERIODIC:
        kinds = (CompositionKind.FORMULA_UNIT,)
    else:
        kinds = (CompositionKind.MOLECULAR, CompositionKind.FORMULA_UNIT)

    required_indices = frozenset(
        index for index, symbol in enumerate(elements) if symbol in required
    )

    for total_atoms in range(
        constraints.min_total_atoms, constraints.max_total_atoms + 1
    ):
        counts = [0] * len(elements)

        def visit(index: int, remaining: int) -> Iterator[tuple[int, ...]]:
            if index == len(elements):
                if remaining == 0:
                    yield tuple(counts)
                return

            minimum = 1 if index in required_indices else 0
            maximum = min(limits[index], remaining)
            for count in range(minimum, maximum + 1):
                counts[index] = count
                yield from visit(index + 1, remaining - count)
            counts[index] = 0

        for count_vector in visit(0, total_atoms):
            element_counts = {
                symbol: count
                for symbol, count in zip(elements, count_vector)
                if count > 0
            }
            for charge in charges:
                for kind in kinds:
                    yield Composition(
                        counts=element_counts,
                        net_charge=charge,
                        kind=kind,
                    )
