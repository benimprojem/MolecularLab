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
    target_properties: tuple[str, ...]
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
        if not targets or any(not target.strip() for target in targets):
            raise ValueError("at least one non-empty target property is required")
        conditions = dict(self.operating_conditions or {})
        if any(not key.strip() for key in conditions):
            raise ValueError("operating condition names cannot be empty")
        object.__setattr__(self, "wavelength_band_m", band)
        object.__setattr__(self, "target_properties", targets)
        object.__setattr__(self, "operating_conditions", MappingProxyType(conditions))
        object.__setattr__(self, "evidence", tuple(self.evidence))
        object.__setattr__(self, "limitations", tuple(self.limitations))
