"""Frequency-dependent material-response result records.

Version: 0.1 | Rev.: 00
Status: Data schema only; no solver is invoked by these records.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from math import isfinite, sqrt
from types import MappingProxyType
from typing import Mapping

from molecular_lab.core.models import EvidenceRecord
from molecular_lab.materials.metamaterials import MetamaterialPropertyDomain


class ConvergenceStatus(str, Enum):
    NOT_CHECKED = "not_checked"
    CONVERGED = "converged"
    NOT_CONVERGED = "not_converged"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class ComplexValue:
    """Real and imaginary components of one complex-valued response entry."""
    real: float
    imaginary: float = 0.0

    def __post_init__(self) -> None:
        if not isfinite(self.real) or not isfinite(self.imaginary):
            raise ValueError("complex response components must be finite")


@dataclass(frozen=True, slots=True)
class ConvergenceInfo:
    """Numerical convergence metadata; absent checks remain explicitly unknown."""
    status: ConvergenceStatus = ConvergenceStatus.NOT_CHECKED
    relative_error: float | None = None
    iterations: int | None = None
    discretization: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.status, ConvergenceStatus):
            object.__setattr__(self, "status", ConvergenceStatus(self.status))
        if self.relative_error is not None and (
            not isfinite(self.relative_error) or self.relative_error < 0
        ):
            raise ValueError("relative_error must be finite and non-negative")
        if self.iterations is not None and (
            isinstance(self.iterations, bool) or not isinstance(self.iterations, int)
            or self.iterations < 0
        ):
            raise ValueError("iterations must be a non-negative integer")


@dataclass(frozen=True, slots=True)
class SolverProvenance:
    """Identifies the method and exact software context behind a result."""
    method: str
    backend: str
    backend_version: str = ""
    settings: Mapping[str, str | float | int | bool] | None = None
    source: str = ""

    def __post_init__(self) -> None:
        if not self.method.strip() or not self.backend.strip():
            raise ValueError("solver method and backend cannot be empty")
        settings = dict(self.settings or {})
        if any(not key.strip() for key in settings):
            raise ValueError("solver setting names cannot be empty")
        for value in settings.values():
            if isinstance(value, float) and not isfinite(value):
                raise ValueError("solver settings must not contain non-finite floats")
        object.__setattr__(self, "settings", MappingProxyType(settings))


@dataclass(frozen=True, slots=True)
class ResponseSample:
    """One scalar or tensor response at a specified frequency and conditions.

    A 1x1 tensor represents an isotropic/scalar response. 2x2 and 3x3 tensors
    represent in-plane or full 3D tensor responses. Entries retain real and
    imaginary parts separately; the sign convention belongs in provenance.
    """
    frequency_hz: float
    tensor: tuple[tuple[ComplexValue, ...], ...]
    wavelength_m: float | None = None
    temperature_K: float | None = None
    propagation_direction: tuple[float, ...] | None = None
    polarization: str = ""
    uncertainty: tuple[tuple[float, ...], ...] | None = None
    convergence: ConvergenceInfo = ConvergenceInfo()

    def __post_init__(self) -> None:
        if not isfinite(self.frequency_hz) or self.frequency_hz <= 0:
            raise ValueError("frequency_hz must be finite and positive")
        if self.wavelength_m is not None and (
            not isfinite(self.wavelength_m) or self.wavelength_m <= 0
        ):
            raise ValueError("wavelength_m must be finite and positive")
        if self.temperature_K is not None and (
            not isfinite(self.temperature_K) or self.temperature_K <= 0
        ):
            raise ValueError("temperature_K must be finite and positive")

        tensor = tuple(tuple(row) for row in self.tensor)
        size = len(tensor)
        if size not in (1, 2, 3) or any(len(row) != size for row in tensor):
            raise ValueError("response tensor must be square with dimension 1, 2, or 3")
        if any(not isinstance(value, ComplexValue) for row in tensor for value in row):
            raise TypeError("tensor entries must be ComplexValue instances")
        object.__setattr__(self, "tensor", tensor)

        if self.propagation_direction is not None:
            direction = tuple(self.propagation_direction)
            if len(direction) != 3 or any(not isfinite(v) for v in direction):
                raise ValueError("propagation_direction must contain three finite components")
            norm = sqrt(sum(v * v for v in direction))
            if abs(norm - 1.0) > 1e-6:
                raise ValueError("propagation_direction must be a unit vector")
            object.__setattr__(self, "propagation_direction", direction)

        if self.uncertainty is not None:
            uncertainty = tuple(tuple(row) for row in self.uncertainty)
            if len(uncertainty) != size or any(len(row) != size for row in uncertainty):
                raise ValueError("uncertainty tensor must match response tensor dimensions")
            if any(not isfinite(v) or v < 0 for row in uncertainty for v in row):
                raise ValueError("uncertainties must be finite and non-negative")
            object.__setattr__(self, "uncertainty", uncertainty)


@dataclass(frozen=True, slots=True)
class MaterialResponseDataset:
    """Frequency sweep for one material property and explicit provenance."""
    domain: MetamaterialPropertyDomain
    metric: str
    unit: str
    samples: tuple[ResponseSample, ...]
    provenance: SolverProvenance
    evidence: tuple[EvidenceRecord, ...] = ()
    limitations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.domain, MetamaterialPropertyDomain):
            object.__setattr__(self, "domain", MetamaterialPropertyDomain(self.domain))
        if not self.metric.strip() or not self.unit.strip():
            raise ValueError("response metric and unit cannot be empty")
        samples = tuple(self.samples)
        if not samples:
            raise ValueError("response dataset must contain at least one sample")
        frequencies = [sample.frequency_hz for sample in samples]
        if len(frequencies) != len(set(frequencies)):
            raise ValueError("response dataset cannot contain duplicate frequencies")
        object.__setattr__(self, "samples", tuple(sorted(samples, key=lambda x: x.frequency_hz)))
        object.__setattr__(self, "evidence", tuple(self.evidence))
        object.__setattr__(self, "limitations", tuple(self.limitations))
