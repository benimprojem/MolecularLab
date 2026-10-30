"""Tests for frequency-dependent material response records."""

import pytest

from molecular_lab.core.models import EvidenceCategory, EvidenceRecord
from molecular_lab.materials.metamaterials import MetamaterialPropertyDomain
from molecular_lab.materials.response import (
    ComplexValue, ConvergenceInfo, ConvergenceStatus, MaterialResponseDataset,
    ResponseSample, SolverProvenance,
)


def _sample(frequency=1.0e12):
    return ResponseSample(
        frequency_hz=frequency, wavelength_m=3.0e-4, temperature_K=300.0,
        tensor=((ComplexValue(-2.0, 0.15), ComplexValue(0.1, -0.02)),
                (ComplexValue(0.1, -0.02), ComplexValue(-2.0, 0.15))),
        propagation_direction=(0.0, 0.0, 1.0), polarization="TE",
        uncertainty=((0.02, 0.01), (0.01, 0.02)),
        convergence=ConvergenceInfo(ConvergenceStatus.CONVERGED, 1e-5, 120, "FEM mesh"),
    )


def test_response_sample_preserves_complex_tensor_and_conditions():
    sample = _sample()
    assert sample.tensor[0][0] == ComplexValue(-2.0, 0.15)
    assert sample.temperature_K == 300.0
    assert sample.propagation_direction == (0.0, 0.0, 1.0)
    assert sample.convergence.status is ConvergenceStatus.CONVERGED


def test_response_sample_accepts_scalar_as_one_by_one_tensor():
    sample = ResponseSample(1e9, ((ComplexValue(4.0, -0.2),),))
    assert sample.tensor[0][0].imaginary == -0.2


@pytest.mark.parametrize("frequency", [0.0, -1.0, float("inf")])
def test_response_sample_rejects_invalid_frequency(frequency):
    with pytest.raises(ValueError, match="frequency_hz"):
        ResponseSample(frequency, ((ComplexValue(1.0),),))


def test_response_sample_rejects_non_square_tensor():
    with pytest.raises(ValueError, match="square"):
        ResponseSample(1e9, ((ComplexValue(1.0), ComplexValue(0.0)),))


def test_response_sample_rejects_non_unit_direction():
    with pytest.raises(ValueError, match="unit vector"):
        ResponseSample(1e9, ((ComplexValue(1.0),),), propagation_direction=(0.0, 0.0, 2.0))


def test_response_sample_rejects_mismatched_uncertainty():
    with pytest.raises(ValueError, match="match response tensor"):
        ResponseSample(
            1e9,
            ((ComplexValue(1.0), ComplexValue(0.0)), (ComplexValue(0.0), ComplexValue(1.0))),
            uncertainty=((0.1,),),
        )


def test_dataset_sorts_samples_and_retains_provenance_and_evidence():
    dataset = MaterialResponseDataset(
        domain=MetamaterialPropertyDomain.OPTICAL, metric="relative_permittivity", unit="1",
        samples=(_sample(2e12), _sample(1e12)),
        provenance=SolverProvenance(
            "frequency-domain FEM", "example-solver", "1.2",
            {"mesh_order": 2, "periodic_boundary": True}, "simulation input archive",
        ),
        evidence=(EvidenceRecord(EvidenceCategory.THEORETICAL, "Illustrative test result", "unit test fixture"),),
    )
    assert [s.frequency_hz for s in dataset.samples] == [1e12, 2e12]
    assert dataset.provenance.backend_version == "1.2"
    assert dataset.evidence[0].category is EvidenceCategory.THEORETICAL
    with pytest.raises(TypeError):
        dataset.provenance.settings["mesh_order"] = 3


def test_dataset_rejects_duplicate_frequencies():
    with pytest.raises(ValueError, match="duplicate frequencies"):
        MaterialResponseDataset(
            MetamaterialPropertyDomain.MAGNETIC, "relative_permeability", "1",
            (_sample(), _sample()), SolverProvenance("test", "test-backend"),
        )


def test_convergence_metadata_rejects_negative_error():
    with pytest.raises(ValueError, match="relative_error"):
        ConvergenceInfo(ConvergenceStatus.CONVERGED, relative_error=-0.1)
