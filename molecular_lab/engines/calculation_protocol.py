"""Common calculation request, solver capability, and execution protocol.

This module defines orchestration contracts only. It does not perform physics
calculations and does not assume that any external solver is installed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from math import isfinite
from types import MappingProxyType
from typing import Mapping, Protocol, runtime_checkable

from molecular_lab.materials.electronic_structure import BandStructureData, DensityOfStatesData
from molecular_lab.materials.response import MaterialResponseDataset


class CalculationDomain(str, Enum):
    ELECTRONIC_STRUCTURE = "electronic_structure"
    ELECTRICAL_TRANSPORT = "electrical_transport"
    MAGNETIC = "magnetic"
    ELECTROMAGNETIC = "electromagnetic"


class CalculationStatus(str, Enum):
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    UNSUPPORTED = "unsupported"
    BACKEND_UNAVAILABLE = "backend_unavailable"


@dataclass(frozen=True, slots=True)
class CalculationRequest:
    """A reproducible request passed to a domain-specific solver adapter."""
    request_id: str
    domain: CalculationDomain
    method: str
    parameters: Mapping[str, str | int | float | bool]
    conditions: Mapping[str, str | int | float | bool] = field(default_factory=dict)
    requested_outputs: tuple[str, ...] = ()
    structure_ref: str = ""

    def __post_init__(self) -> None:
        if not self.request_id.strip() or not self.method.strip():
            raise ValueError("request_id and method cannot be empty")
        if not isinstance(self.domain, CalculationDomain):
            object.__setattr__(self, "domain", CalculationDomain(self.domain))
        parameters = dict(self.parameters)
        conditions = dict(self.conditions)
        for mapping_name, mapping in (("parameters", parameters), ("conditions", conditions)):
            if any(not isinstance(key, str) or not key.strip() for key in mapping):
                raise ValueError(f"{mapping_name} keys must be non-empty strings")
            for value in mapping.values():
                if isinstance(value, float) and not isfinite(value):
                    raise ValueError(f"{mapping_name} cannot contain non-finite floats")
                if not isinstance(value, (str, int, float, bool)):
                    raise TypeError(f"{mapping_name} values must be scalar JSON-compatible values")
        outputs = tuple(self.requested_outputs)
        if any(not item.strip() for item in outputs):
            raise ValueError("requested output names cannot be empty")
        if len(outputs) != len(set(outputs)):
            raise ValueError("requested output names must be unique")
        object.__setattr__(self, "parameters", MappingProxyType(parameters))
        object.__setattr__(self, "conditions", MappingProxyType(conditions))
        object.__setattr__(self, "requested_outputs", outputs)


@dataclass(frozen=True, slots=True)
class SolverCapability:
    """Declared solver support; declaration is not proof of installation."""
    solver_id: str
    display_name: str
    domains: frozenset[CalculationDomain]
    methods: tuple[str, ...]
    optional_dependencies: tuple[str, ...] = ()
    description: str = ""

    def __post_init__(self) -> None:
        if not self.solver_id.strip() or not self.display_name.strip():
            raise ValueError("solver identity fields cannot be empty")
        domains = frozenset(
            item if isinstance(item, CalculationDomain) else CalculationDomain(item)
            for item in self.domains
        )
        if not domains:
            raise ValueError("solver must declare at least one calculation domain")
        methods = tuple(item.strip() for item in self.methods)
        if not methods or any(not item for item in methods):
            raise ValueError("solver must declare non-empty methods")
        object.__setattr__(self, "domains", domains)
        object.__setattr__(self, "methods", methods)
        object.__setattr__(self, "optional_dependencies", tuple(self.optional_dependencies))


@dataclass(frozen=True, slots=True)
class SolverAvailability:
    capability: SolverCapability
    available: bool
    version: str = ""
    reason: str = ""


@dataclass(frozen=True, slots=True)
class CalculationResult:
    """Normalized orchestration result; datasets retain their own provenance."""
    request_id: str
    status: CalculationStatus
    solver_id: str
    message: str = ""
    response_datasets: tuple[MaterialResponseDataset, ...] = ()
    band_structures: tuple[BandStructureData, ...] = ()
    density_of_states: tuple[DensityOfStatesData, ...] = ()
    scalar_outputs: Mapping[str, str | int | float | bool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.request_id.strip():
            raise ValueError("result request_id cannot be empty")
        if not isinstance(self.status, CalculationStatus):
            object.__setattr__(self, "status", CalculationStatus(self.status))
        outputs = dict(self.scalar_outputs)
        if any(not key.strip() for key in outputs):
            raise ValueError("scalar output names cannot be empty")
        for value in outputs.values():
            if isinstance(value, float) and not isfinite(value):
                raise ValueError("scalar outputs cannot contain non-finite floats")
            if not isinstance(value, (str, int, float, bool)):
                raise TypeError("scalar outputs must be scalar JSON-compatible values")
        object.__setattr__(self, "scalar_outputs", MappingProxyType(outputs))
        object.__setattr__(self, "response_datasets", tuple(self.response_datasets))
        object.__setattr__(self, "band_structures", tuple(self.band_structures))
        object.__setattr__(self, "density_of_states", tuple(self.density_of_states))


@runtime_checkable
class SolverAdapter(Protocol):
    """Contract implemented by each optional physics backend adapter."""

    @property
    def capability(self) -> SolverCapability: ...

    def availability(self) -> SolverAvailability: ...

    def run(self, request: CalculationRequest) -> CalculationResult: ...


class SolverRegistry:
    """Small explicit registry; no automatic package imports or fake fallbacks."""

    def __init__(self) -> None:
        self._adapters: dict[str, SolverAdapter] = {}

    def register(self, adapter: SolverAdapter) -> None:
        solver_id = adapter.capability.solver_id
        if solver_id in self._adapters:
            raise ValueError(f"solver already registered: {solver_id}")
        self._adapters[solver_id] = adapter

    def availability(self) -> tuple[SolverAvailability, ...]:
        return tuple(
            self._adapters[key].availability()
            for key in sorted(self._adapters)
        )

    def matching(self, request: CalculationRequest) -> tuple[SolverAdapter, ...]:
        return tuple(
            adapter
            for _, adapter in sorted(self._adapters.items())
            if request.domain in adapter.capability.domains
            and request.method in adapter.capability.methods
        )

    def execute(self, request: CalculationRequest, solver_id: str | None = None) -> CalculationResult:
        if solver_id is not None:
            adapter = self._adapters.get(solver_id)
            if adapter is None:
                return CalculationResult(
                    request.request_id, CalculationStatus.UNSUPPORTED, solver_id,
                    "Requested solver is not registered.",
                )
            candidates = (adapter,)
        else:
            candidates = self.matching(request)

        adapter = next(
            (
                candidate for candidate in candidates
                if request.domain in candidate.capability.domains
                and request.method in candidate.capability.methods
            ),
            None,
        )
        if adapter is None:
            return CalculationResult(
                request.request_id, CalculationStatus.UNSUPPORTED, solver_id or "",
                "No registered solver declares support for this domain and method.",
            )

        state = adapter.availability()
        if not state.available:
            return CalculationResult(
                request.request_id, CalculationStatus.BACKEND_UNAVAILABLE,
                adapter.capability.solver_id, state.reason or "Solver backend is unavailable.",
            )

        try:
            result = adapter.run(request)
        except Exception as exc:  # adapter boundary: convert backend errors to explicit failure
            return CalculationResult(
                request.request_id, CalculationStatus.FAILED,
                adapter.capability.solver_id, f"{type(exc).__name__}: {exc}",
            )

        if result.request_id != request.request_id:
            return CalculationResult(
                request.request_id, CalculationStatus.FAILED,
                adapter.capability.solver_id,
                "Solver returned a result for a different request ID.",
            )
        if result.solver_id != adapter.capability.solver_id:
            return CalculationResult(
                request.request_id, CalculationStatus.FAILED,
                adapter.capability.solver_id,
                "Solver result identity does not match the registered adapter.",
            )
        return result
