"""Explicit-criterion screening labels for electrical behavior.

This module does not impose universal conductor/semiconductor/insulator cutoffs.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from math import isfinite

from molecular_lab.core.models import EvidenceRecord


class ElectricalBehavior(str, Enum):
    CONDUCTOR = "conductor"
    SEMICONDUCTOR = "semiconductor"
    INSULATOR = "insulator"
    INDETERMINATE = "indeterminate"


class ElectricalMetric(str, Enum):
    BAND_GAP = "band_gap"
    CONDUCTIVITY = "conductivity"


@dataclass(frozen=True, slots=True)
class ElectricalClassificationCriteria:
    """Two explicit cutoffs defining lower/upper screening regions.

    For band gap (eV): <= lower is conductor-like; >= upper is insulator-like.
    For conductivity (S/m): <= lower is insulator-like; >= upper is conductor-like.
    Values between the cutoffs receive a semiconductor-like screening label.
    """
    metric: ElectricalMetric
    lower_threshold: float
    upper_threshold: float
    source: str
    notes: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.metric, ElectricalMetric):
            object.__setattr__(self, "metric", ElectricalMetric(self.metric))
        if not isfinite(self.lower_threshold) or not isfinite(self.upper_threshold):
            raise ValueError("classification thresholds must be finite")
        if self.lower_threshold < 0 or self.upper_threshold <= self.lower_threshold:
            raise ValueError("thresholds must satisfy 0 <= lower < upper")
        if not self.source.strip():
            raise ValueError("classification criterion source cannot be empty")


@dataclass(frozen=True, slots=True)
class ElectricalBehaviorAssessment:
    """A condition-bound classification with explicit criterion provenance."""
    behavior: ElectricalBehavior
    metric: ElectricalMetric
    measured_or_calculated_value: float
    unit: str
    temperature_K: float
    frequency_hz: float
    criteria: ElectricalClassificationCriteria
    rationale: str
    evidence: tuple[EvidenceRecord, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.behavior, ElectricalBehavior):
            object.__setattr__(self, "behavior", ElectricalBehavior(self.behavior))
        if not isinstance(self.metric, ElectricalMetric):
            object.__setattr__(self, "metric", ElectricalMetric(self.metric))
        if self.metric is not self.criteria.metric:
            raise ValueError("assessment metric must match classification criteria")
        if not isfinite(self.measured_or_calculated_value) or self.measured_or_calculated_value < 0:
            raise ValueError("measured_or_calculated_value must be finite and non-negative")
        expected_unit = "eV" if self.metric is ElectricalMetric.BAND_GAP else "S/m"
        if self.unit != expected_unit:
            raise ValueError(f"{self.metric.value} classification requires unit {expected_unit}")
        if not isfinite(self.temperature_K) or self.temperature_K <= 0:
            raise ValueError("temperature_K must be finite and positive")
        if not isfinite(self.frequency_hz) or self.frequency_hz < 0:
            raise ValueError("frequency_hz must be finite and non-negative (0 means DC)")
        if not self.rationale.strip():
            raise ValueError("classification rationale cannot be empty")
        object.__setattr__(self, "evidence", tuple(self.evidence))


def assess_electrical_behavior(
    value: float,
    *,
    temperature_K: float,
    criteria: ElectricalClassificationCriteria,
    frequency_hz: float = 0.0,
    evidence: tuple[EvidenceRecord, ...] = (),
) -> ElectricalBehaviorAssessment:
    """Apply caller-supplied thresholds; this is not a universal material label."""
    if not isfinite(value) or value < 0:
        raise ValueError("classification value must be finite and non-negative")

    low = criteria.lower_threshold
    high = criteria.upper_threshold
    if criteria.metric is ElectricalMetric.BAND_GAP:
        if value <= low:
            behavior = ElectricalBehavior.CONDUCTOR
            rationale = f"Band gap {value:g} eV is at or below the configured conductor-like cutoff {low:g} eV."
        elif value >= high:
            behavior = ElectricalBehavior.INSULATOR
            rationale = f"Band gap {value:g} eV is at or above the configured insulator-like cutoff {high:g} eV."
        else:
            behavior = ElectricalBehavior.SEMICONDUCTOR
            rationale = f"Band gap {value:g} eV lies between the configured cutoffs {low:g} and {high:g} eV."
        unit = "eV"
    else:
        if value <= low:
            behavior = ElectricalBehavior.INSULATOR
            rationale = f"Conductivity {value:g} S/m is at or below the configured insulator-like cutoff {low:g} S/m."
        elif value >= high:
            behavior = ElectricalBehavior.CONDUCTOR
            rationale = f"Conductivity {value:g} S/m is at or above the configured conductor-like cutoff {high:g} S/m."
        else:
            behavior = ElectricalBehavior.SEMICONDUCTOR
            rationale = f"Conductivity {value:g} S/m lies between the configured cutoffs {low:g} and {high:g} S/m."
        unit = "S/m"

    return ElectricalBehaviorAssessment(
        behavior=behavior,
        metric=criteria.metric,
        measured_or_calculated_value=value,
        unit=unit,
        temperature_K=temperature_K,
        frequency_hz=frequency_hz,
        criteria=criteria,
        rationale=rationale,
        evidence=evidence,
    )
