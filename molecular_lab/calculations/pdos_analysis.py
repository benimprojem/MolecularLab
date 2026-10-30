"""Aggregation helpers for atom/orbital projected density of states.

These helpers sum supplied projector values; they do not make the projection
unique or guarantee completeness of the atomic projector basis.
Version: 0.1 | Rev.: 01
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from math import isfinite

from molecular_lab.materials.electronic_structure import (
    BandGapAssessment, BandGapStatus, BandStructureData,
    ProjectedDensityOfStatesData, assess_band_gap,
)


@dataclass(frozen=True, slots=True)
class ProjectedDOSAggregate:
    """Summed PDOS for one selected atom/element, angular momentum and spin."""
    element: str | None
    atom_index: int | None
    angular_momentum: int
    spin_channel: str
    local_dos_states_per_ev: tuple[float, ...]
    component_dos_states_per_ev: tuple[float, ...]

    def __post_init__(self):
        local = tuple(self.local_dos_states_per_ev)
        components = tuple(self.component_dos_states_per_ev)
        if self.angular_momentum < 0 or not self.spin_channel.strip():
            raise ValueError("angular momentum and spin channel must be valid")
        if (self.element is None) == (self.atom_index is None):
            raise ValueError("aggregate must identify exactly one element or atom")
        if not local or len(local) != len(components):
            raise ValueError("aggregate arrays must be non-empty and have equal lengths")
        if any(not isfinite(v) or v < 0 for v in local + components):
            raise ValueError("aggregate DOS values must be finite and non-negative")
        object.__setattr__(self, "local_dos_states_per_ev", local)
        object.__setattr__(self, "component_dos_states_per_ev", components)


def aggregate_projected_dos(
    data: ProjectedDensityOfStatesData,
    *,
    group_by: str = "element",
) -> tuple[ProjectedDOSAggregate, ...]:
    """Sum PDOS over atoms or wavefunctions, retaining l and spin separately.

    group_by='element' sums same-element atoms; group_by='atom' keeps each
    atom separate. For each l channel, angular components are summed across
    wavefunctions and (when requested) atoms. Spin channels are never merged.
    """
    if group_by not in {"element", "atom"}:
        raise ValueError("group_by must be 'element' or 'atom'")
    buckets = {}
    for record in data.records:
        identity = record.element if group_by == "element" else record.atom_index
        key = (identity, record.angular_momentum, record.spin_channel)
        if key not in buckets:
            buckets[key] = ([0.0] * len(data.energies_ev),
                            [0.0] * len(data.energies_ev))
        local_sum, component_sum = buckets[key]
        for i, value in enumerate(record.local_dos_states_per_ev):
            local_sum[i] += value
        for component in record.components:
            for i, value in enumerate(component.density_states_per_ev):
                component_sum[i] += value

    result = []
    for (identity, angular_momentum, spin), (local, components) in sorted(
        buckets.items(), key=lambda item: (str(item[0][0]), item[0][1], item[0][2])
    ):
        result.append(ProjectedDOSAggregate(
            element=identity if group_by == "element" else None,
            atom_index=identity if group_by == "atom" else None,
            angular_momentum=angular_momentum,
            spin_channel=spin,
            local_dos_states_per_ev=tuple(local),
            component_dos_states_per_ev=tuple(components),
        ))
    return tuple(result)


def integrate_pdos_window(
    energies_ev: tuple[float, ...],
    density_states_per_ev: tuple[float, ...],
    lower_ev: float,
    upper_ev: float,
) -> float:
    """Trapezoid-integrate a DOS curve over an energy window in eV.

    Window boundaries are linearly interpolated when they fall between
    sampled points. The returned value is an integrated projected-state
    count under the supplied DOS normalization.
    """
    energies = tuple(energies_ev)
    density = tuple(density_states_per_ev)
    if len(energies) < 2 or len(energies) != len(density):
        raise ValueError("energy and DOS arrays must have equal length >= 2")
    if any(not isfinite(v) for v in energies) or any(
        not isfinite(v) or v < 0 for v in density
    ):
        raise ValueError("energies and DOS values must be finite; DOS must be non-negative")
    if any(b <= a for a, b in zip(energies, energies[1:])):
        raise ValueError("energy grid must be strictly increasing")
    if not isfinite(lower_ev) or not isfinite(upper_ev) or lower_ev >= upper_ev:
        raise ValueError("energy window must have finite lower < upper bounds")
    if lower_ev < energies[0] or upper_ev > energies[-1]:
        raise ValueError("energy window must lie within the sampled energy grid")

    def interpolate(x):
        for i in range(len(energies) - 1):
            left, right = energies[i], energies[i + 1]
            if left <= x <= right:
                fraction = (x - left) / (right - left)
                return density[i] + fraction * (density[i + 1] - density[i])
        return density[-1]

    points = [(lower_ev, interpolate(lower_ev))]
    points.extend((e, d) for e, d in zip(energies, density) if lower_ev < e < upper_ev)
    points.append((upper_ev, interpolate(upper_ev)))
    return sum(
        (x1 - x0) * (y0 + y1) / 2
        for (x0, y0), (x1, y1) in zip(points, points[1:])
    )


@dataclass(frozen=True, slots=True)
class PDOSWindowSummary:
    """Integrated PDOS weights for one atom/element, l channel and spin."""
    element: str | None
    atom_index: int | None
    angular_momentum: int
    spin_channel: str
    lower_ev: float
    upper_ev: float
    integrated_local_states: float
    integrated_component_states: float

    def __post_init__(self):
        if (self.element is None) == (self.atom_index is None):
            raise ValueError("summary must identify exactly one element or atom")
        if self.angular_momentum < 0 or not self.spin_channel.strip():
            raise ValueError("angular momentum and spin channel must be valid")
        if (not isfinite(self.lower_ev) or not isfinite(self.upper_ev)
                or self.lower_ev >= self.upper_ev):
            raise ValueError("summary energy bounds must be finite and increasing")
        if (not isfinite(self.integrated_local_states)
                or not isfinite(self.integrated_component_states)
                or self.integrated_local_states < 0
                or self.integrated_component_states < 0):
            raise ValueError("integrated PDOS values must be finite and non-negative")


def summarize_pdos_window(
    data: ProjectedDensityOfStatesData,
    lower_ev: float,
    upper_ev: float,
    *,
    group_by: str = "element",
) -> tuple[PDOSWindowSummary, ...]:
    """Integrate each aggregated l/spin contribution over an energy interval.

    Returns both the sum of QE local-DOS columns and the sum of angular
    components. The two are kept separate as a consistency check; neither is
    assumed to equal the total DOS or to be a complete projector basis.
    """
    aggregates = aggregate_projected_dos(data, group_by=group_by)
    summaries = []
    for item in aggregates:
        summaries.append(PDOSWindowSummary(
            element=item.element,
            atom_index=item.atom_index,
            angular_momentum=item.angular_momentum,
            spin_channel=item.spin_channel,
            lower_ev=lower_ev,
            upper_ev=upper_ev,
            integrated_local_states=integrate_pdos_window(
                data.energies_ev, item.local_dos_states_per_ev, lower_ev, upper_ev
            ),
            integrated_component_states=integrate_pdos_window(
                data.energies_ev, item.component_dos_states_per_ev, lower_ev, upper_ev
            ),
        ))
    return tuple(summaries)


@dataclass(frozen=True, slots=True)
class BandEdgePDOSAnalysis:
    """PDOS contributions adjacent to sampled VBM and CBM energies."""
    band_gap: BandGapAssessment
    pdos_energy_shift_ev: float
    vbm_window: tuple[PDOSWindowSummary, ...]
    cbm_window: tuple[PDOSWindowSummary, ...]

    def __post_init__(self):
        if self.band_gap.status is not BandGapStatus.GAP:
            raise ValueError("band-edge PDOS analysis requires a resolved band gap")
        if not isfinite(self.pdos_energy_shift_ev):
            raise ValueError("PDOS energy shift must be finite")
        object.__setattr__(self, 'vbm_window', tuple(self.vbm_window))
        object.__setattr__(self, 'cbm_window', tuple(self.cbm_window))


def analyze_band_edge_pdos(
    bands: BandStructureData,
    pdos: ProjectedDensityOfStatesData,
    *,
    window_ev: float = 0.25,
    pdos_energy_shift_ev: float = 0.0,
    group_by: str = 'element',
    gap_tolerance_ev: float = 1e-5,
) -> BandEdgePDOSAnalysis:
    """Correlate sampled VBM/CBM energies with element/orbital PDOS windows.

    The shift is added to every PDOS energy to express it on the band-energy
    reference. It is explicit by default: callers must ensure both datasets
    use compatible structures, methods and energy references. The VBM window
    lies below VBM; the CBM window lies above CBM. This is not k-resolved
    orbital character.
    """
    if not isfinite(window_ev) or window_ev <= 0:
        raise ValueError("window_ev must be finite and positive")
    if not isfinite(pdos_energy_shift_ev):
        raise ValueError("pdos_energy_shift_ev must be finite")
    assessment = assess_band_gap(bands, gap_tolerance_ev=gap_tolerance_ev)
    if assessment.status is not BandGapStatus.GAP:
        raise ValueError(f'band-edge PDOS analysis requires a resolved gap; got {assessment.status.value}')
    vbm = assessment.valence_band_maximum_ev
    cbm = assessment.conduction_band_minimum_ev
    if vbm is None or cbm is None:
        raise ValueError("band assessment does not contain both band edges")
    aligned = replace(pdos, energies_ev=tuple(e + pdos_energy_shift_ev for e in pdos.energies_ev))
    return BandEdgePDOSAnalysis(
        band_gap=assessment,
        pdos_energy_shift_ev=pdos_energy_shift_ev,
        vbm_window=summarize_pdos_window(aligned, vbm - window_ev, vbm, group_by=group_by),
        cbm_window=summarize_pdos_window(aligned, cbm, cbm + window_ev, group_by=group_by),
    )