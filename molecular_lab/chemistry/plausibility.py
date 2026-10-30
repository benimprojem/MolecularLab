"""
Conservative composition-level plausibility diagnostics.

Version: 0.1
Rev.: 00
Status: Initial heuristic diagnostics; no structure or stability calculations.

The checks in this module are advisory. They must never be interpreted as a
probability of existence, a stability verdict, or a universal chemical rule.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
from itertools import product

from molecular_lab.chemistry.elements import get_element
from molecular_lab.core.models import Composition


class DiagnosticStatus(str, Enum):
    SUPPORTED_PATTERN = "supported_pattern"
    NO_MATCH_IN_REFERENCE = "no_match_in_limited_reference"
    NOT_ASSESSED = "not_assessed"


class ElectronParity(str, Enum):
    EVEN = "even"
    ODD = "odd"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class PlausibilityFinding:
    rule_id: str
    status: str
    message: str
    details: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class PlausibilityAssessment:
    composition: Composition
    electron_count: int | None
    electron_parity: ElectronParity
    oxidation_state_status: DiagnosticStatus
    findings: tuple[PlausibilityFinding, ...]

    @property
    def requires_review(self) -> bool:
        """True when a diagnostic found a gap that deserves later review."""
        return any(
            finding.status in {"review", "not_assessed"}
            for finding in self.findings
        )


# Selected common formal oxidation states, not exhaustive allowed-state lists.
# Zero is handled separately for single-element elemental compositions.
_COMMON_OXIDATION_STATES: dict[str, tuple[int, ...]] = {
    "H": (-1, 1),
    "Li": (1,), "Na": (1,), "K": (1,), "Rb": (1,), "Cs": (1,),
    "Be": (2,), "Mg": (2,), "Ca": (2,), "Sr": (2,), "Ba": (2,),
    "B": (-3, 3),
    "C": (-4, -3, -2, -1, 1, 2, 3, 4),
    "N": (-3, -2, -1, 1, 2, 3, 4, 5),
    "O": (-2, -1, 1, 2),
    "F": (-1,),
    "Si": (-4, 2, 4),
    "P": (-3, 1, 3, 5),
    "S": (-2, 2, 4, 6),
    "Cl": (-1, 1, 3, 5, 7),
    "Br": (-1, 1, 3, 5, 7),
    "I": (-1, 1, 3, 5, 7),
    "Al": (3,), "Ga": (1, 3), "In": (1, 3),
    "Zn": (2,), "Cd": (2,), "Hg": (1, 2),
    "Ag": (1,), "Cu": (1, 2), "Au": (1, 3),
    "Fe": (2, 3), "Co": (2, 3), "Ni": (2, 3),
    "Mn": (2, 3, 4, 6, 7), "Cr": (2, 3, 6),
    "Ti": (2, 3, 4), "V": (2, 3, 4, 5),
    "Zr": (4,), "Hf": (4,),
    "Sn": (2, 4), "Pb": (2, 4),
    "As": (-3, 3, 5), "Se": (-2, 2, 4, 6),
    "Te": (-2, 2, 4, 6),
    "Xe": (2, 4, 6, 8),
}

_MAX_ASSIGNMENT_STATES = 200_000


def assess_composition(composition: Composition) -> PlausibilityAssessment:
    """
    Run conservative composition-level diagnostics.

    The electron count is exact for the supplied neutral-atom composition and
    net charge. Oxidation-state matching uses only the small curated table
    above and is a formal bookkeeping check, not a bonding model.
    """
    atomic_numbers: dict[str, int] = {}
    for symbol in composition.counts:
        atomic_numbers[symbol] = get_element(symbol).atomic_number

    if composition.net_charge is None:
        electron_count = None
        parity = ElectronParity.UNKNOWN
    else:
        electron_count = sum(
            atomic_numbers[symbol] * count
            for symbol, count in composition.counts.items()
        ) - composition.net_charge
        parity = (
            ElectronParity.EVEN if electron_count % 2 == 0
            else ElectronParity.ODD
        )

    findings: list[PlausibilityFinding] = [
        PlausibilityFinding(
            rule_id="electron_count.parity",
            status="informational",
            message=(
                "Net electron count is unspecified because net charge is unknown."
                if electron_count is None
                else f"Total electron count is {electron_count} ({parity.value})."
            ),
            details=(
                "Electron parity can inform later spin-state enumeration.",
                "Even/odd parity alone does not establish stability or a valid structure.",
            ),
        )
    ]

    oxidation_status, oxidation_detail = _check_oxidation_state_balance(
        composition
    )
    if oxidation_status == DiagnosticStatus.SUPPORTED_PATTERN:
        findings.append(
            PlausibilityFinding(
                rule_id="oxidation_state.formal_balance",
                status="supporting",
                message="At least one formal oxidation-state assignment matches the net charge.",
                details=(oxidation_detail,),
            )
        )
    elif oxidation_status == DiagnosticStatus.NO_MATCH_IN_REFERENCE:
        findings.append(
            PlausibilityFinding(
                rule_id="oxidation_state.formal_balance",
                status="review",
                message=(
                    "No matching assignment was found in the limited common-state "
                    "reference; this is not proof of chemical impossibility."
                ),
                details=(oxidation_detail,),
            )
        )
    else:
        findings.append(
            PlausibilityFinding(
                rule_id="oxidation_state.formal_balance",
                status="not_assessed",
                message=(
                    "Oxidation-state matching was not assessed because at least "
                    "one element is outside the curated reference."
                ),
                details=(oxidation_detail,),
            )
        )

    return PlausibilityAssessment(
        composition=composition,
        electron_count=electron_count,
        electron_parity=parity,
        oxidation_state_status=oxidation_status,
        findings=tuple(findings),
    )


def _check_oxidation_state_balance(
    composition: Composition,
) -> tuple[DiagnosticStatus, str]:
    if composition.net_charge is None:
        return (
            DiagnosticStatus.NOT_ASSESSED,
            "Net charge is unspecified.",
        )

    symbols = tuple(sorted(composition.counts))
    if len(symbols) == 1:
        return (
            DiagnosticStatus.SUPPORTED_PATTERN,
            "Single-element composition: formal oxidation state zero balances the composition.",
        )

    missing = [symbol for symbol in symbols if symbol not in _COMMON_OXIDATION_STATES]
    if missing:
        return (
            DiagnosticStatus.NOT_ASSESSED,
            "No curated common-state list for: " + ", ".join(missing),
        )

    target = composition.net_charge
    states = tuple(_COMMON_OXIDATION_STATES[symbol] for symbol in symbols)
    counts = tuple(composition.counts[symbol] for symbol in symbols)

    possible, explored = _has_charge_balanced_assignment(
        states, counts, target
    )
    if possible:
        return (
            DiagnosticStatus.SUPPORTED_PATTERN,
            "A formal charge-balanced assignment exists using the curated state lists.",
        )
    if not explored:
        return (
            DiagnosticStatus.NOT_ASSESSED,
            "Assignment search exceeded its configured internal safety limit.",
        )
    return (
        DiagnosticStatus.NO_MATCH_IN_REFERENCE,
        "Curated oxidation-state combinations did not sum to the supplied net charge.",
    )


@lru_cache(maxsize=4096)
def _has_charge_balanced_assignment(
    states: tuple[tuple[int, ...], ...],
    counts: tuple[int, ...],
    target: int,
) -> tuple[bool, bool]:
    combinations = 1
    for options in states:
        combinations *= len(options)
        if combinations > _MAX_ASSIGNMENT_STATES:
            return False, False

    for assignment in product(*states):
        if sum(count * state for count, state in zip(counts, assignment)) == target:
            return True, True
    return False, True
