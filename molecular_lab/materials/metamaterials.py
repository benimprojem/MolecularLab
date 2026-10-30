"""
Metamaterial design records.

Version: 0.1
Rev.: 00
Status: Foundational geometry and target-property records only.

These records do not calculate electromagnetic response or validate a design.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from math import isfinite
from types import MappingProxyType
from typing import Mapping

from molecular_lab.core.models import EvidenceRecord


class MetamaterialPropertyDomain(str, Enum):
    """Physical response family requested or evaluated for a design."""
    OPTICAL = "optical"
    MAGNETIC = "magnetic"
    ELECTRICAL = "electrical"
    ELECTRONIC_STRUCTURE = "electronic_structure"
    MULTIPHYSICS = "multiphysics"


class TargetRelation(str, Enum):
    """How a requested value or interval should be interpreted."""
    MINIMIZE = "minimize"
    MAXIMIZE = "maximize"
    APPROACH = "approach"
    WITHIN_RANGE = "within_range"


@dataclass(frozen=True, slots=True)
class MetamaterialPropertyTarget:
    """A measurable property objective with units and optional operating band.

    Values are goals/constraints, not computed results. The metric name and unit
    are explicit so the same schema can represent conductivity, band gap,
    permeability, refractive index, optical loss, and related quantities.
    """
    domain: MetamaterialPropertyDomain
    metric: str
    relation: TargetRelation
    unit: str
    target_value: float | None = None
    minimum: float | None = None
    maximum: float | None = None
    spectral_band_m: tuple[float, float] | None = None
    notes: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.domain, MetamaterialPropertyDomain):
            object.__setattr__(self, "domain", MetamaterialPropertyDomain(self.domain))
        if not isinstance(self.relation, TargetRelation):
            object.__setattr__(self, "relation", TargetRelation(self.relation))
        if not self.metric.strip() or not self.unit.strip():
            raise ValueError("property metric and unit cannot be empty")
        values = (self.target_value, self.minimum, self.maximum)
        if any(value is not None and not isfinite(value) for value in values):
            raise ValueError("property target values must be finite")
        if self.relation is TargetRelation.APPROACH and self.target_value is None:
            raise ValueError("approach relation requires target_value")
        if self.relation is TargetRelation.WITHIN_RANGE:
            if self.minimum is None or self.maximum is None or self.minimum > self.maximum:
                raise ValueError("within_range requires minimum <= maximum")
        if self.minimum is not None and self.maximum is not None and self.minimum > self.maximum:
            raise ValueError("minimum cannot exceed maximum")
        if self.spectral_band_m is not None:
            band = tuple(self.spectral_band_m)
            if len(band) != 2 or any(not isfinite(v) or v <= 0 for v in band) or band[0] >= band[1]:
                raise ValueError("spectral_band_m must be an increasing positive metre interval")
            object.__setattr__(self, "spectral_band_m", band)


class MetamaterialMechanism(str, Enum):
    RESONANT = "resonant"
    GEOMETRIC_PHASE = "geometric_phase"
    PHOTONIC_CRYSTAL = "photonic_crystal"
    EFFECTIVE_MEDIUM = "effective_medium"
    HYBRID = "hybrid"


@dataclass(frozen=True, slots=True)
class UnitCellParameter:
    """One geometric/material parameter with an explicit unit."""
    name: str
    value: float
    unit: str

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.unit.strip():
            raise ValueError("parameter name and unit cannot be empty")
        if not isfinite(self.value):
            raise ValueError("parameter value must be finite")
        if self.value <= 0:
            raise ValueError("geometric parameter values must be positive")


@dataclass(frozen=True, slots=True)
class MetamaterialUnitCell:
    """Parameterized 2D or 3D repeating cell; no physical response is implied."""
    name: str
    geometry_family: str
    dimension: int
    lattice_period_m: tuple[float, ...]
    parameters: tuple[UnitCellParameter, ...] = ()
    constituent_materials: tuple[str, ...] = ()
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.geometry_family.strip():
            raise ValueError("unit-cell name and geometry_family cannot be empty")
        if self.dimension not in (2, 3):
            raise ValueError("metamaterial unit-cell dimension must be 2 or 3")
        periods = tuple(self.lattice_period_m)
        if len(periods) != self.dimension:
            raise ValueError("lattice_period_m must match the unit-cell dimension")
        if any(not isfinite(value) or value <= 0 for value in periods):
            raise ValueError("lattice periods must be finite positive values in metres")
        parameters = tuple(self.parameters)
        if len({item.name for item in parameters}) != len(parameters):
            raise ValueError("unit-cell parameter names must be unique")
        materials = tuple(self.constituent_materials)
        if any(not item.strip() for item in materials):
            raise ValueError("constituent material names cannot be empty")
        object.__setattr__(self, "lattice_period_m", periods)
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "constituent_materials", materials)


@dataclass(frozen=True, slots=True)
class MetamaterialDesign:
    """A proposed periodic design and its requested spectral/functional targets."""
    design_id: str
    mechanism: MetamaterialMechanism
    unit_cell: MetamaterialUnitCell
    wavelength_band_m: tuple[float, float]
    target_properties: tuple[str, ...] = ()
    property_targets: tuple[MetamaterialPropertyTarget, ...] = ()
    operating_conditions: Mapping[str, str | float] | None = None
    evidence: tuple[EvidenceRecord, ...] = ()
    limitations: tuple[str, ...] = (
        "Design record only; no electromagnetic simulation or experimental validation is implied.",
    )

    def __post_init__(self) -> None:
        if not self.design_id.strip():
            raise ValueError("design_id cannot be empty")
        if not isinstance(self.mechanism, MetamaterialMechanism):
            object.__setattr__(self, "mechanism", MetamaterialMechanism(self.mechanism))
        band = tuple(self.wavelength_band_m)
        if len(band) != 2 or any(not isfinite(value) or value <= 0 for value in band):
            raise ValueError("wavelength_band_m must contain two finite positive metre values")
        if band[0] >= band[1]:
            raise ValueError("wavelength band minimum must be smaller than maximum")
        targets = tuple(self.target_properties)
        property_targets = tuple(self.property_targets)
        if any(not target.strip() for target in targets):
            raise ValueError("target property labels cannot be empty")
        if not targets and not property_targets:
            raise ValueError("at least one target property or typed property target is required")
        conditions = dict(self.operating_conditions or {})
        if any(not key.strip() for key in conditions):
            raise ValueError("operating condition names cannot be empty")
        object.__setattr__(self, "wavelength_band_m", band)
        object.__setattr__(self, "target_properties", targets)
        object.__setattr__(self, "property_targets", property_targets)
        object.__setattr__(self, "operating_conditions", MappingProxyType(conditions))
        object.__setattr__(self, "evidence", tuple(self.evidence))
        object.__setattr__(self, "limitations", tuple(self.limitations))
