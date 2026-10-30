"""Electronic-structure result records and conservative band-gap analysis.

Schema and analysis only: this module does not run DFT or infer experimental truth.
Version: 0.1 | Rev.: 01
"""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from math import isfinite

@dataclass(frozen=True, slots=True)
class KPoint:
    """Fractional reciprocal-space coordinate and integration weight."""
    coordinates: tuple[float, float, float]
    weight: float
    def __post_init__(self):
        coordinates = tuple(self.coordinates)
        if len(coordinates) != 3 or any(not isfinite(v) for v in coordinates):
            raise ValueError("k-point coordinates must contain three finite values")
        if not isfinite(self.weight) or self.weight <= 0:
            raise ValueError("k-point weight must be finite and positive")
        object.__setattr__(self, "coordinates", coordinates)

@dataclass(frozen=True, slots=True)
class BandChannelData:
    """Band energies in eV, arranged as bands × k-points for one spin channel."""
    spin_channel: str
    energies_ev: tuple[tuple[float, ...], ...]
    occupations: tuple[tuple[float, ...], ...] | None = None
    maximum_occupation: float = 2.0
    def __post_init__(self):
        if not self.spin_channel.strip():
            raise ValueError("spin_channel cannot be empty")
        energies = tuple(tuple(row) for row in self.energies_ev)
        if not energies or not energies[0]:
            raise ValueError("at least one band and one k-point are required")
        point_count = len(energies[0])
        if any(len(row) != point_count for row in energies):
            raise ValueError("all bands must contain the same number of k-point energies")
        if any(not isfinite(v) for row in energies for v in row):
            raise ValueError("band energies must be finite")
        object.__setattr__(self, "energies_ev", energies)
        if not isfinite(self.maximum_occupation) or self.maximum_occupation <= 0:
            raise ValueError("maximum_occupation must be finite and positive")
        if self.occupations is not None:
            occupations = tuple(tuple(row) for row in self.occupations)
            if len(occupations) != len(energies) or any(len(row) != point_count for row in occupations):
                raise ValueError("occupation array must match the band-energy dimensions")
            if any(not isfinite(v) or v < 0 or v > self.maximum_occupation
                   for row in occupations for v in row):
                raise ValueError("occupations must lie between zero and maximum_occupation")
            object.__setattr__(self, "occupations", occupations)
    @property
    def band_count(self):
        return len(self.energies_ev)
    @property
    def kpoint_count(self):
        return len(self.energies_ev[0])

@dataclass(frozen=True, slots=True)
class BandStructureData:
    """Spin-resolved band energies with k-point path/grid and calculation context."""
    kpoints: tuple[KPoint, ...]
    channels: tuple[BandChannelData, ...]
    fermi_level_ev: float | None = None
    occupation_tolerance: float = 1e-3
    method: str = ""
    backend: str = ""
    backend_version: str = ""
    exchange_correlation: str = ""
    pseudopotential_set: str = ""
    notes: tuple[str, ...] = ()
    def __post_init__(self):
        kpoints, channels = tuple(self.kpoints), tuple(self.channels)
        if not kpoints or not channels:
            raise ValueError("band structure requires k-points and at least one spin channel")
        if len({c.spin_channel for c in channels}) != len(channels):
            raise ValueError("spin channel names must be unique")
        if any(c.kpoint_count != len(kpoints) for c in channels):
            raise ValueError("band data k-point count must match the k-point list")
        if self.fermi_level_ev is not None and not isfinite(self.fermi_level_ev):
            raise ValueError("Fermi level must be finite")
        if not isfinite(self.occupation_tolerance) or not 0 <= self.occupation_tolerance < 0.5:
            raise ValueError("occupation_tolerance must be in [0, 0.5)")
        object.__setattr__(self, "kpoints", kpoints)
        object.__setattr__(self, "channels", channels)
        object.__setattr__(self, "notes", tuple(self.notes))

class BandGapStatus(str, Enum):
    GAP = "gap"
    METALLIC = "metallic"
    UNDETERMINED = "undetermined"

@dataclass(frozen=True, slots=True)
class BandGapAssessment:
    status: BandGapStatus
    indirect_gap_ev: float | None = None
    direct_gap_ev: float | None = None
    valence_band_maximum_ev: float | None = None
    conduction_band_minimum_ev: float | None = None
    vbm_kpoint_index: int | None = None
    cbm_kpoint_index: int | None = None
    reason: str = ""

def assess_band_gap(data: BandStructureData, *, gap_tolerance_ev: float = 1e-5) -> BandGapAssessment:
    """Estimate direct/indirect gaps; missing evidence yields undetermined."""
    if not isfinite(gap_tolerance_ev) or gap_tolerance_ev < 0:
        raise ValueError("gap_tolerance_ev must be finite and non-negative")
    global_vbm = None
    global_cbm = None
    direct_candidates = []
    for channel in data.channels:
        if channel.occupations is None and data.fermi_level_ev is None:
            return BandGapAssessment(BandGapStatus.UNDETERMINED,
                reason="Neither occupation data nor a Fermi level is available.")
        for k_index in range(len(data.kpoints)):
            occupied, empty = [], []
            for band_index, band in enumerate(channel.energies_ev):
                energy = band[k_index]
                if channel.occupations is not None:
                    occupation = channel.occupations[band_index][k_index]
                    full = channel.maximum_occupation
                    if data.occupation_tolerance < occupation < full - data.occupation_tolerance:
                        if data.fermi_level_ev is not None and abs(energy - data.fermi_level_ev) <= gap_tolerance_ev:
                            return BandGapAssessment(BandGapStatus.METALLIC,
                                reason=f"Partially occupied state lies at the Fermi level in {channel.spin_channel} at k-point {k_index}.")
                        return BandGapAssessment(BandGapStatus.UNDETERMINED,
                            reason="Partial occupation detected; finite-temperature or smearing effects cannot be distinguished from metallic occupation without additional context.")
                    if occupation >= full - data.occupation_tolerance:
                        occupied.append(energy)
                    elif occupation <= data.occupation_tolerance:
                        empty.append(energy)
                else:
                    delta = energy - data.fermi_level_ev
                    if abs(delta) <= gap_tolerance_ev:
                        return BandGapAssessment(BandGapStatus.METALLIC,
                            reason=f"A band lies at the Fermi level in {channel.spin_channel} at k-point {k_index}.")
                    (occupied if delta < 0 else empty).append(energy)
            if not occupied or not empty:
                continue
            local_vbm, local_cbm = max(occupied), min(empty)
            direct_candidates.append(local_cbm - local_vbm)
            if global_vbm is None or local_vbm > global_vbm[0]:
                global_vbm = (local_vbm, k_index)
            if global_cbm is None or local_cbm < global_cbm[0]:
                global_cbm = (local_cbm, k_index)
    if global_vbm is None or global_cbm is None or not direct_candidates:
        return BandGapAssessment(BandGapStatus.UNDETERMINED,
            reason="Could not identify both occupied and unoccupied band edges across the supplied channels.")
    indirect_gap = global_cbm[0] - global_vbm[0]
    direct_gap = min(direct_candidates)
    fields = dict(indirect_gap_ev=indirect_gap, direct_gap_ev=direct_gap,
        valence_band_maximum_ev=global_vbm[0], conduction_band_minimum_ev=global_cbm[0],
        vbm_kpoint_index=global_vbm[1], cbm_kpoint_index=global_cbm[1])
    if indirect_gap <= gap_tolerance_ev:
        return BandGapAssessment(BandGapStatus.METALLIC, **fields,
            reason="Valence and conduction edges overlap or are unresolved within tolerance.")
    return BandGapAssessment(BandGapStatus.GAP, **fields,
        reason="Gap estimated from supplied occupations/Fermi level and sampled k-points.")

@dataclass(frozen=True, slots=True)
class DensityOfStatesChannel:
    spin_channel: str
    density_states_per_ev: tuple[float, ...]
    def __post_init__(self):
        values = tuple(self.density_states_per_ev)
        if not self.spin_channel.strip() or not values:
            raise ValueError("DOS channel requires a name and at least one value")
        if any(not isfinite(v) or v < 0 for v in values):
            raise ValueError("DOS values must be finite and non-negative")
        object.__setattr__(self, "density_states_per_ev", values)

@dataclass(frozen=True, slots=True)
class DensityOfStatesData:
    energies_ev: tuple[float, ...]
    channels: tuple[DensityOfStatesChannel, ...]
    fermi_level_ev: float | None = None
    projection: str = "total"
    normalization: str = "per_cell"
    method: str = ""
    backend: str = ""
    backend_version: str = ""
    notes: tuple[str, ...] = ()
    def __post_init__(self):
        energies, channels = tuple(self.energies_ev), tuple(self.channels)
        if len(energies) < 2 or any(not isfinite(v) for v in energies):
            raise ValueError("DOS energy grid requires at least two finite values")
        if any(b <= a for a, b in zip(energies, energies[1:])):
            raise ValueError("DOS energy grid must be strictly increasing")
        if not channels:
            raise ValueError("DOS requires at least one channel")
        if len({c.spin_channel for c in channels}) != len(channels):
            raise ValueError("DOS spin channel names must be unique")
        if any(len(c.density_states_per_ev) != len(energies) for c in channels):
            raise ValueError("each DOS channel must match the energy grid length")
        if self.fermi_level_ev is not None and not isfinite(self.fermi_level_ev):
            raise ValueError("Fermi level must be finite")
        if not self.projection.strip() or not self.normalization.strip():
            raise ValueError("projection and normalization descriptions cannot be empty")
        object.__setattr__(self, "energies_ev", energies)
        object.__setattr__(self, "channels", channels)
        object.__setattr__(self, "notes", tuple(self.notes))
