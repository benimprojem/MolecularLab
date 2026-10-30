"""
Chemical formula parsing for ordinary element-count formulae.

Version: 0.1
Rev.: 00
Status: Initial parser; does not infer bonding, charge, phase, or stability.

Supported:
- Element symbols with case-sensitive capitalization (H, He, Fe).
- Positive integer subscripts (H2, Fe2).
- Nested round and square groups (Ca(OH)2, K4[Fe(CN)6]).
- Adjacent groups and elements.

Not yet supported:
- Ionic charge suffixes, isotope notation, hydrate separators, leading
  coefficients, dot-disconnected mixtures, and fractional occupancies.
"""

from __future__ import annotations

from collections import defaultdict

from molecular_lab.chemistry.elements import get_element
from molecular_lab.core.models import Composition, CompositionKind


_OPEN_TO_CLOSE = {"(": ")", "[": "]"}
_CLOSE_TO_OPEN = {close: open_ for open_, close in _OPEN_TO_CLOSE.items()}


class FormulaParseError(ValueError):
    """Raised when a formula is malformed or contains an unknown element."""


def parse_formula(formula: str) -> Composition:
    """
    Parse a conventional formula into an element-count Composition.

    The parser validates syntax and element identity only. It does not determine
    whether the resulting composition is a realizable or stable compound.
    """
    if not isinstance(formula, str) or not formula:
        raise FormulaParseError("Formula must be a non-empty string")
    if any(character.isspace() for character in formula):
        raise FormulaParseError("Whitespace is not allowed inside a formula")

    stack: list[tuple[str | None, defaultdict[str, int]]] = [
        (None, defaultdict(int))
    ]
    index = 0
    token_count = 0

    while index < len(formula):
        character = formula[index]

        if character in _OPEN_TO_CLOSE:
            stack.append((character, defaultdict(int)))
            index += 1
            continue

        if character in _CLOSE_TO_OPEN:
            if len(stack) == 1:
                raise FormulaParseError(
                    f"Unexpected closing bracket {character!r} at position {index}"
                )
            opener, group_counts = stack.pop()
            if _OPEN_TO_CLOSE[opener] != character:
                raise FormulaParseError(
                    f"Mismatched bracket {character!r} at position {index}"
                )
            index += 1
            multiplier, index = _read_number(formula, index)
            for symbol, count in group_counts.items():
                stack[-1][1][symbol] += count * multiplier
            token_count += 1
            continue

        if character.isupper():
            start = index
            index += 1
            if index < len(formula) and formula[index].islower():
                index += 1
            symbol = formula[start:index]
            try:
                get_element(symbol)
            except ValueError as exc:
                raise FormulaParseError(str(exc)) from exc
            count, index = _read_number(formula, index)
            stack[-1][1][symbol] += count
            token_count += 1
            continue

        if character.isdigit():
            raise FormulaParseError(
                f"Unexpected number at position {index}; counts must follow an element or group"
            )

        raise FormulaParseError(
            f"Unsupported character {character!r} at position {index}"
        )

    if len(stack) != 1:
        opener = stack[-1][0]
        raise FormulaParseError(f"Unclosed group starting with {opener!r}")
    if token_count == 0 or not stack[0][1]:
        raise FormulaParseError("Formula contains no elements")

    return Composition(
        counts=dict(stack[0][1]),
        net_charge=0,
        kind=CompositionKind.MOLECULAR,
        input_formula=formula,
    )


def _read_number(formula: str, index: int) -> tuple[int, int]:
    """Read an optional positive integer multiplier."""
    start = index
    while index < len(formula) and formula[index].isdigit():
        index += 1
    if start == index:
        return 1, index
    value = int(formula[start:index])
    if value < 1:
        raise FormulaParseError("Subscripts must be positive integers")
    return value, index
