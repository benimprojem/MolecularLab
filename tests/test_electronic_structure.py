"""Tests for electronic-structure records and band-gap screening."""
import pytest
from molecular_lab.materials.electronic_structure import (
    BandChannelData, BandGapStatus, BandStructureData, DensityOfStatesChannel,
    DensityOfStatesData, KPoint, assess_band_gap,
)

def _kpoints():
    return (KPoint((0.0, 0.0, 0.0), 0.5), KPoint((0.5, 0.0, 0.0), 0.5))

def test_indirect_gap_is_distinguished_from_direct_gap():
    channel = BandChannelData("total", ((-1.0, -0.4), (1.2, 0.7)), ((2.0, 2.0), (0.0, 0.0)))
    result = assess_band_gap(BandStructureData(_kpoints(), (channel,)))
    assert result.status is BandGapStatus.GAP
    assert result.indirect_gap_ev == pytest.approx(1.1)
    assert result.direct_gap_ev == pytest.approx(1.1)
    assert result.vbm_kpoint_index == 1
    assert result.cbm_kpoint_index == 1

def test_partial_occupation_is_classified_as_metallic():
    channel = BandChannelData("spin-up", ((-0.2, 0.1), (1.0, 1.2)), ((0.5, 0.5), (0.0, 0.0)), maximum_occupation=1.0)
    result = assess_band_gap(BandStructureData(_kpoints(), (channel,)))
    assert result.status is BandGapStatus.METALLIC
    assert "Partial" in result.reason

def test_overlapping_band_edges_are_metallic():
    channel = BandChannelData("total", ((-0.5, 0.2), (0.1, 0.8)), ((2.0, 2.0), (0.0, 0.0)))
    result = assess_band_gap(BandStructureData(_kpoints(), (channel,)))
    assert result.status is BandGapStatus.METALLIC
    assert result.indirect_gap_ev < 0

def test_missing_occupation_and_fermi_information_is_undetermined():
    channel = BandChannelData("total", ((-1.0, -0.8), (1.0, 1.2)))
    assert assess_band_gap(BandStructureData(_kpoints(), (channel,))).status is BandGapStatus.UNDETERMINED

def test_fermi_level_can_support_gap_estimate_when_occupations_absent():
    channel = BandChannelData("total", ((-1.0, -0.8), (1.0, 1.2)))
    result = assess_band_gap(BandStructureData(_kpoints(), (channel,), fermi_level_ev=0.0))
    assert result.status is BandGapStatus.GAP
    assert result.indirect_gap_ev == pytest.approx(1.8)

def test_band_data_validates_dimensions_and_occupancy():
    with pytest.raises(ValueError, match="same number"):
        BandChannelData("total", ((-1.0, -0.5), (1.0,)))
    with pytest.raises(ValueError, match="maximum_occupation"):
        BandChannelData("total", ((-1.0,),), maximum_occupation=0.0)
    with pytest.raises(ValueError, match="between zero"):
        BandChannelData("total", ((-1.0,),), ((3.0,),))

def test_band_structure_validates_kpoint_count_and_unique_channels():
    with pytest.raises(ValueError, match="k-point count"):
        BandStructureData(_kpoints(), (BandChannelData("total", ((-1.0,), (1.0,))),))
    channel = BandChannelData("same", ((-1.0, -0.8), (1.0, 1.2)))
    with pytest.raises(ValueError, match="unique"):
        BandStructureData(_kpoints(), (channel, channel))

def test_dos_validates_grid_and_channel_shapes():
    dos = DensityOfStatesData((-2.0, -1.0, 0.0, 1.0), (DensityOfStatesChannel("total", (0.1, 0.5, 0.2, 0.0)),), fermi_level_ev=0.0)
    assert dos.energies_ev[0] == -2.0
    with pytest.raises(ValueError, match="strictly increasing"):
        DensityOfStatesData((0.0, 1.0, 1.0), (DensityOfStatesChannel("total", (0.0, 1.0, 0.0)),))
    with pytest.raises(ValueError, match="energy grid length"):
        DensityOfStatesData((0.0, 1.0), (DensityOfStatesChannel("total", (0.0,)),))

def test_dos_rejects_negative_density():
    with pytest.raises(ValueError, match="non-negative"):
        DensityOfStatesChannel("total", (0.0, -0.1))
