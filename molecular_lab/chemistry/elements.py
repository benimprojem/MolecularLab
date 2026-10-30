"""
Periodic-table reference data access.

Version: 0.1
Rev.: 00
Status: Initial element lookup; no chemical properties or valence rules inferred.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files


@dataclass(frozen=True, slots=True)
class Element:
    atomic_number: int
    symbol: str
    name: str


@lru_cache(maxsize=1)
def all_elements() -> tuple[Element, ...]:
    """Load the bundled 118-element reference table."""
    data_path = files("molecular_lab.chemistry").joinpath("data/elements.json")
    payload = json.loads(data_path.read_text(encoding="utf-8"))
    elements = tuple(Element(**item) for item in payload["elements"])

    if len(elements) != 118:
        raise RuntimeError("Periodic-table data must contain exactly 118 elements")
    if [element.atomic_number for element in elements] != list(range(1, 119)):
        raise RuntimeError("Element atomic numbers must be unique and sequential")
    if len({element.symbol for element in elements}) != len(elements):
        raise RuntimeError("Element symbols must be unique")

    return elements


@lru_cache(maxsize=1)
def _indexes() -> tuple[dict[str, Element], dict[int, Element], dict[str, Element]]:
    elements = all_elements()
    return (
        {element.symbol: element for element in elements},
        {element.atomic_number: element for element in elements},
        {element.name.casefold(): element for element in elements},
    )


def get_element(symbol: str) -> Element:
    """Return an element by its case-sensitive chemical symbol."""
    try:
        return _indexes()[0][symbol]
    except KeyError as exc:
        raise ValueError(f"Unknown element symbol: {symbol!r}") from exc


def get_element_by_atomic_number(atomic_number: int) -> Element:
    try:
        return _indexes()[1][atomic_number]
    except KeyError as exc:
        raise ValueError(f"Unknown atomic number: {atomic_number!r}") from exc


def find_element_by_name(name: str) -> Element:
    try:
        return _indexes()[2][name.casefold()]
    except KeyError as exc:
        raise ValueError(f"Unknown element name: {name!r}") from exc
