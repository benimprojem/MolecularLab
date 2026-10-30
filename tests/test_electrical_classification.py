"""Tests for criterion-bound conductor/semiconductor/insulator screening."""

import pytest

from molecular_lab.materials.electrical_classification import (
    ElectricalBehavior, ElectricalClassificationCriteria, ElectricalMetric,
    assess_electrical_behavior,
)


@pytest.mark.parametrize(
    "value,expected",
    [(0.0, ElectricalBehavior.CONDUCTOR), (0.5, ElectricalBehavior.SEMICONDUCTOR),
     (2.0, ElectricalBehavior.INSULATOR)],
)
def test_band_gap_screening_uses_explicit_cutoffs(value, expected):
    criteria = ElectricalClassificationCriteria(ElectricalMetric.BAND_GAP, 0.1, 1.5, "test protocol")
    result = assess_electrical_behavior(value, temperature_K=300.0, criteria=criteria)
    assert result.behavior is expected
    assert result.unit == "eV"
    assert result.frequency_hz == 0.0


@pytest.mark.parametrize(
    "value,expected",
    [(1e-8, ElectricalBehavior.INSULATOR), (1e2, ElectricalBehavior.SEMICONDUCTOR),
     (1e6, ElectricalBehavior.CONDUCTOR)],
)
def test_conductivity_screening_uses_explicit_cutoffs(value, expected):
    criteria = ElectricalClassificationCriteria(ElectricalMetric.CONDUCTIVITY, 1e-4, 1e4, "test protocol")
    result = assess_electrical_behavior(
        value, temperature_K=77.0, frequency_hz=1e9, criteria=criteria,
    )
    assert result.behavior is expected
    assert result.unit == "S/m"
    assert result.temperature_K == 77.0
    assert result.frequency_hz == 1e9


def test_classification_requires_explicit_criterion_source():
    with pytest.raises(ValueError, match="source"):
        ElectricalClassificationCriteria(ElectricalMetric.BAND_GAP, 0.1, 1.0, " ")


def test_classification_rejects_invalid_threshold_order():
    with pytest.raises(ValueError, match="lower < upper"):
        ElectricalClassificationCriteria(ElectricalMetric.CONDUCTIVITY, 2.0, 1.0, "test")


def test_assessment_rejects_wrong_unit_and_temperature():
    criteria = ElectricalClassificationCriteria(ElectricalMetric.BAND_GAP, 0.1, 1.0, "test")
    with pytest.raises(ValueError, match="requires unit"):
        from molecular_lab.materials.electrical_classification import ElectricalBehaviorAssessment
        ElectricalBehaviorAssessment(
            ElectricalBehavior.SEMICONDUCTOR, ElectricalMetric.BAND_GAP, 0.5, "S/m",
            300.0, 0.0, criteria, "test rationale",
        )
    with pytest.raises(ValueError, match="temperature_K"):
        assess_electrical_behavior(0.5, temperature_K=0.0, criteria=criteria)
