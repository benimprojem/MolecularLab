"""Tests for the common physics solver orchestration protocol."""

import pytest

from molecular_lab.materials.electronic_structure import (\n    BandChannelData, BandStructureData, DensityOfStatesChannel, DensityOfStatesData, KPoint,\n)\nfrom molecular_lab.engines.calculation_protocol import (
    CalculationDomain, CalculationRequest, CalculationResult, CalculationStatus,
    SolverAvailability, SolverCapability, SolverRegistry,
)


class FakeAdapter:
    def __init__(self, available=True, fail=False, wrong_id=False):
        self._available = available
        self._fail = fail
        self._wrong_id = wrong_id
        self._capability = SolverCapability(
            solver_id="fake-fem",
            display_name="Fake FEM test adapter",
            domains=frozenset({CalculationDomain.ELECTROMAGNETIC}),
            methods=("frequency_domain_fem",),
            optional_dependencies=("fake-solver-package",),
        )

    @property
    def capability(self):
        return self._capability

    def availability(self):
        return SolverAvailability(
            self.capability, self._available, "0.test",
            "" if self._available else "test backend intentionally unavailable",
        )

    def run(self, request):
        if self._fail:
            raise RuntimeError("deliberate test failure")
        return CalculationResult(
            request.request_id if not self._wrong_id else "wrong-id",
            CalculationStatus.SUCCEEDED,
            self.capability.solver_id,
            scalar_outputs={"test_scalar": 1.25},
        )


def _request():
    return CalculationRequest(
        request_id="req-001",
        domain=CalculationDomain.ELECTROMAGNETIC,
        method="frequency_domain_fem",
        parameters={"mesh_order": 2, "periodic": True},
        conditions={"temperature_K": 300.0},
        requested_outputs=("reflection", "transmission"),
        structure_ref="unit-cell-001",
    )


def test_request_freezes_parameters_and_conditions():
    request = _request()
    with pytest.raises(TypeError):
        request.parameters["mesh_order"] = 3
    with pytest.raises(TypeError):
        request.conditions["temperature_K"] = 77.0


def test_request_rejects_nonfinite_parameter():
    with pytest.raises(ValueError, match="non-finite"):
        CalculationRequest("r", CalculationDomain.MAGNETIC, "method", {"x": float("nan")})


def test_registry_reports_declared_capabilities_and_matches_request():
    registry = SolverRegistry()
    adapter = FakeAdapter()
    registry.register(adapter)
    assert registry.matching(_request()) == (adapter,)
    assert registry.availability()[0].available


def test_registry_returns_unsupported_when_no_solver_matches():
    result = SolverRegistry().execute(_request())
    assert result.status is CalculationStatus.UNSUPPORTED


def test_registry_does_not_claim_unavailable_backend_ran():
    registry = SolverRegistry()
    registry.register(FakeAdapter(available=False))
    result = registry.execute(_request())
    assert result.status is CalculationStatus.BACKEND_UNAVAILABLE
    assert "unavailable" in result.message


def test_registry_executes_matching_adapter_and_returns_result():
    registry = SolverRegistry()
    registry.register(FakeAdapter())
    result = registry.execute(_request())
    assert result.status is CalculationStatus.SUCCEEDED
    assert result.scalar_outputs["test_scalar"] == 1.25


def test_registry_converts_backend_exception_to_failed_result():
    registry = SolverRegistry()
    registry.register(FakeAdapter(fail=True))
    result = registry.execute(_request())
    assert result.status is CalculationStatus.FAILED
    assert "deliberate test failure" in result.message


def test_registry_rejects_result_for_another_request():
    registry = SolverRegistry()
    registry.register(FakeAdapter(wrong_id=True))
    result = registry.execute(_request())
    assert result.status is CalculationStatus.FAILED
    assert "different request ID" in result.message


def test_registry_rejects_duplicate_solver_ids():
    registry = SolverRegistry()
    registry.register(FakeAdapter())
    with pytest.raises(ValueError, match="already registered"):
        registry.register(FakeAdapter())


def test_explicit_unknown_solver_is_unsupported():
    result = SolverRegistry().execute(_request(), solver_id="missing")
    assert result.status is CalculationStatus.UNSUPPORTED


def test_calculation_result_can_carry_band_structure_and_dos():
    band = BandStructureData(
        (KPoint((0.0, 0.0, 0.0), 1.0),),
        (BandChannelData("total", ((-1.0,), (1.0,)), ((2.0,), (0.0,))),),
    )
    dos = DensityOfStatesData(
        (-1.0, 0.0, 1.0),
        (DensityOfStatesChannel("total", (0.1, 0.2, 0.1)),),
    )
    result = CalculationResult(
        "req-electronic", CalculationStatus.SUCCEEDED, "fake-dft",
        band_structures=(band,), density_of_states=(dos,),
    )
    assert result.band_structures == (band,)
    assert result.density_of_states == (dos,)
