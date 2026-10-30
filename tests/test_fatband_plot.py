import math

import pytest

from molecular_lab.calculations.dft_settings import BandPath, BandPathPoint
from molecular_lab.calculations.fatband import (
    FatBandContribution, FatBandDataset, FatBandPoint,
)
from molecular_lab.calculations.fatband_plot import prepare_fatband_plot_data
from molecular_lab.materials.electronic_structure import KPoint


def sample_dataset():
    return FatBandDataset(
        kpoints=(
            KPoint((0.0, 0.0, 0.0), 1.0, "fractional"),
            KPoint((0.5, 0.0, 0.0), 1.0, "fractional"),
        ),
        spin_channels=("total",),
        number_of_bands=1,
        points=(
            FatBandPoint(
                0, 1, "total", -1.0,
                (
                    FatBandContribution(1, "Fe", 0, 0, None, "s", (1,), 0.25),
                    FatBandContribution(2, "O", 1, 0, None, "p", (2,), 0.10),
                ),
            ),
            FatBandPoint(
                1, 1, "total", -0.8,
                (
                    FatBandContribution(1, "Fe", 0, 0, None, "s", (1,), 0.50),
                    FatBandContribution(2, "O", 1, 0, None, "p", (2,), 0.20),
                ),
            ),
        ),
    )


def sample_path():
    return BandPath((
        BandPathPoint((0.0, 0.0, 0.0), "GAMMA", 1),
        BandPathPoint((0.5, 0.0, 0.0), "X", 0),
    ))


def test_plot_data_uses_reciprocal_lattice_distance_and_symmetry_ticks():
    result = prepare_fatband_plot_data(
        sample_dataset(),
        band_path=sample_path(),
        lattice_vectors_angstrom=((2.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.0)),
        group_by="element_orbital",
        normalization="global_max",
    )
    assert result.kpoint_distances_inverse_angstrom == pytest.approx((0.0, math.pi / 2))
    assert [(tick.label, tick.kpoint_index) for tick in result.symmetry_ticks] == [
        ("GAMMA", 0), ("X", 1),
    ]
    fe = next(item for item in result.series if item.label == "Fe · s")
    assert [point.raw_weight for point in fe.points] == pytest.approx((0.25, 0.50))
    assert [point.display_weight for point in fe.points] == pytest.approx((0.5, 1.0))
    assert fe.points[0].energy_ev == pytest.approx(-1.0)


def test_plot_filters_and_minimum_weight_keep_raw_values_separate():
    result = prepare_fatband_plot_data(
        sample_dataset(),
        band_path=sample_path(),
        lattice_vectors_angstrom=((2.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.0)),
        group_by="element",
        normalization="raw",
        minimum_raw_weight=0.15,
        elements=("O",),
    )
    assert len(result.series) == 1
    assert result.series[0].label == "O"
    assert result.series[0].points[0].raw_weight == pytest.approx(0.10)
    assert result.series[0].points[0].display_weight == 0.0
    assert result.series[0].points[1].display_weight == pytest.approx(0.20)


def test_plot_rejects_path_order_or_lattice_mismatch():
    with pytest.raises(ValueError, match="does not match"):
        prepare_fatband_plot_data(
            sample_dataset(),
            band_path=BandPath((
                BandPathPoint((0.0, 0.0, 0.0), "GAMMA", 1),
                BandPathPoint((0.0, 0.5, 0.0), "Y", 0),
            )),
            lattice_vectors_angstrom=((2.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.0)),
        )
    with pytest.raises(ValueError, match="non-zero"):
        prepare_fatband_plot_data(
            sample_dataset(),
            band_path=sample_path(),
            lattice_vectors_angstrom=((1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        )


def test_plot_normalization_does_not_mutate_raw_fatband_dataset():
    dataset = sample_dataset()
    before = tuple(
        contribution.raw_weight
        for point in dataset.points
        for contribution in point.contributions
    )
    prepare_fatband_plot_data(
        dataset,
        band_path=sample_path(),
        lattice_vectors_angstrom=((2.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.0)),
        normalization="per_band_max",
    )
    after = tuple(
        contribution.raw_weight
        for point in dataset.points
        for contribution in point.contributions
    )
    assert after == before
