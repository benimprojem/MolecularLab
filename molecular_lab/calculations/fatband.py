"""Backend-neutral fat-band records and QE projection-to-band mapping."""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from molecular_lab.calculations.qe_projection import QEAtomicProjectionData
from molecular_lab.materials.electronic_structure import BandStructureData, KPoint


@dataclass(frozen=True, slots=True)
class FatBandContribution:
    """Projection weight grouped by atom and angular-momentum projector."""

    atom_index: int
    element: str
    angular_momentum: int
    magnetic_number: int | None
    total_angular_momentum: float | None
    orbital_label: str
    source_wavefunction_indices: tuple[int, ...]
    raw_weight: float

    def __post_init__(self) -> None:
        if self.atom_index < 1:
            raise ValueError("atom_index must be one-based and positive")
        if not self.element or not self.orbital_label:
            raise ValueError("element and orbital_label cannot be empty")
        if not self.source_wavefunction_indices:
            raise ValueError("at least one source wavefunction index is required")
        if not isfinite(self.raw_weight) or self.raw_weight < 0:
            raise ValueError("raw_weight must be finite and non-negative")


@dataclass(frozen=True, slots=True)
class FatBandPoint:
    """One spin-resolved band eigenstate at one sampled k-point."""

    kpoint_index: int
    band_index: int
    spin_channel: str
    energy_ev: float
    contributions: tuple[FatBandContribution, ...]

    def __post_init__(self) -> None:
        if self.kpoint_index < 0 or self.band_index < 1:
            raise ValueError("kpoint_index is zero-based; band_index is one-based")
        if not self.spin_channel or not isfinite(self.energy_ev):
            raise ValueError("spin_channel and finite energy are required")
        object.__setattr__(self, "contributions", tuple(self.contributions))

    @property
    def total_projected_weight(self) -> float:
        """Raw sum of retained projector weights; not forced to one."""
        return sum(item.raw_weight for item in self.contributions)


@dataclass(frozen=True, slots=True)
class FatBandDataset:
    """Band energies and atom/orbital weights on the same sampled k-point grid."""

    kpoints: tuple[KPoint, ...]
    spin_channels: tuple[str, ...]
    number_of_bands: int
    points: tuple[FatBandPoint, ...]
    method: str = "Quantum ESPRESSO projwfc.x atomic projections"
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        kpoints = tuple(self.kpoints)
        spins = tuple(self.spin_channels)
        points = tuple(self.points)
        if not kpoints or not spins or self.number_of_bands < 1:
            raise ValueError("fat-band dataset requires k-points, spins, and bands")
        if len(spins) != len(set(spins)) or any(not spin.strip() for spin in spins):
            raise ValueError("fat-band spin channels must be non-empty and unique")
        expected = {
            (k, spin, band)
            for k in range(len(kpoints))
            for spin in spins
            for band in range(1, self.number_of_bands + 1)
        }
        observed = [(p.kpoint_index, p.spin_channel, p.band_index) for p in points]
        if len(observed) != len(set(observed)):
            raise ValueError("fat-band dataset contains duplicate state keys")
        if set(observed) != expected:
            raise ValueError("fat-band dataset state grid is incomplete or out of range")
        atom_elements: dict[int, str] = {}
        for point in points:
            seen_contributions: set[tuple] = set()
            for contribution in point.contributions:
                previous = atom_elements.setdefault(contribution.atom_index, contribution.element)
                if previous != contribution.element:
                    raise ValueError(
                        f"atom {contribution.atom_index} has conflicting element labels"
                    )
                key = (
                    contribution.atom_index, contribution.angular_momentum,
                    contribution.magnetic_number, contribution.total_angular_momentum,
                )
                if key in seen_contributions:
                    raise ValueError("duplicate atom/orbital contribution in one fat-band state")
                seen_contributions.add(key)
                source_indices = contribution.source_wavefunction_indices
                if (not source_indices or len(source_indices) != len(set(source_indices))
                        or any(isinstance(index, bool) or not isinstance(index, int) or index < 1
                               for index in source_indices)):
                    raise ValueError("source wavefunction indices must be unique positive integers")
        object.__setattr__(self, "kpoints", kpoints)
        object.__setattr__(self, "spin_channels", spins)
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "notes", tuple(self.notes))


def _orbital_label(angular_momentum: int, magnetic_number: int | None) -> str:
    if magnetic_number is None:
        return {0: "s", 1: "p", 2: "d", 3: "f"}.get(
            angular_momentum, f"l={angular_momentum}"
        )
    return f"l={angular_momentum},m={magnetic_number}"


def build_fatband_dataset(
    bands: BandStructureData,
    projections: QEAtomicProjectionData,
) -> FatBandDataset:
    """Join projections to energies by k-point order, band index and spin.

    Both inputs must originate from the same saved pw.x calculation and have
    identical k-point ordering. This function validates dimensions and labels;
    it cannot prove coordinate identity when the projection file omits a
    reciprocal-space grid.
    """
    if len(bands.kpoints) != projections.number_of_kpoints:
        raise ValueError("Band structure and atomic projections have different k-point counts")
    if len(bands.channels) != len(projections.spin_channels):
        raise ValueError("Band structure and atomic projections have different spin-channel counts")

    band_by_spin = {channel.spin_channel: channel for channel in bands.channels}
    aliases = {
        "total": ("total", "non_spin_polarized"),
        "spin_up": ("spin_up", "up"),
        "spin_down": ("spin_down", "down"),
    }
    resolved: dict[str, str] = {}
    for projection_spin in projections.spin_channels:
        candidates = aliases.get(projection_spin, (projection_spin,))
        match = next((name for name in candidates if name in band_by_spin), None)
        if match is None:
            raise ValueError(f"No matching band spin channel for projection channel {projection_spin!r}")
        resolved[projection_spin] = match
    if len(set(resolved.values())) != len(resolved):
        raise ValueError("Spin-channel aliases do not map one-to-one")

    number_of_bands = projections.number_of_bands
    for projection_spin, band_spin in resolved.items():
        channel = band_by_spin[band_spin]
        if channel.band_count != number_of_bands:
            raise ValueError(f"Band count mismatch for spin channel {projection_spin}")
        if channel.kpoint_count != projections.number_of_kpoints:
            raise ValueError(f"K-point count mismatch for spin channel {projection_spin}")

    wavefunctions = {item.index: item for item in projections.atomic_wavefunctions}
    grouped: dict[tuple, dict[str, object]] = {}
    for coefficient in projections.coefficients:
        wfc = wavefunctions.get(coefficient.atomic_wfc_index)
        if wfc is None:
            raise ValueError(
                f"Projection references unknown atomic wavefunction {coefficient.atomic_wfc_index}"
            )
        label = _orbital_label(wfc.angular_momentum, wfc.magnetic_number)
        key = (
            coefficient.kpoint_index, coefficient.spin_channel,
            coefficient.band_index, wfc.atom_index, wfc.element,
            wfc.angular_momentum, wfc.magnetic_number, wfc.total_angular_momentum, label,
        )
        aggregate = grouped.setdefault(key, {"weight": 0.0, "wfcs": set()})
        aggregate["weight"] = float(aggregate["weight"]) + coefficient.weight
        aggregate["wfcs"].add(coefficient.atomic_wfc_index)

    points: list[FatBandPoint] = []
    for spin in projections.spin_channels:
        channel = band_by_spin[resolved[spin]]
        for band_index in range(1, number_of_bands + 1):
            for kpoint_index in range(projections.number_of_kpoints):
                contributions = []
                for key, aggregate in grouped.items():
                    (key_k, key_spin, key_band, atom, element,
                     angular, magnetic, total_j, label) = key
                    if (key_k, key_spin, key_band) != (kpoint_index, spin, band_index):
                        continue
                    contributions.append(FatBandContribution(
                        atom_index=atom,
                        element=element,
                        angular_momentum=angular,
                        magnetic_number=magnetic,
                        total_angular_momentum=total_j,
                        orbital_label=label,
                        source_wavefunction_indices=tuple(sorted(aggregate["wfcs"])),
                        raw_weight=float(aggregate["weight"]),
                    ))
                contributions.sort(key=lambda item: (
                    item.atom_index, item.angular_momentum,
                    -1 if item.magnetic_number is None else item.magnetic_number,
                ))
                points.append(FatBandPoint(
                    kpoint_index=kpoint_index,
                    band_index=band_index,
                    spin_channel=spin,
                    energy_ev=channel.energies_ev[band_index - 1][kpoint_index],
                    contributions=tuple(contributions),
                ))

    return FatBandDataset(
        kpoints=bands.kpoints,
        spin_channels=projections.spin_channels,
        number_of_bands=number_of_bands,
        points=tuple(points),
        notes=(
            "Projection weights are raw |c|^2 values summed over matching atomic wavefunctions.",
            "Weights are not normalized; completeness depends on the pseudopotential projector basis.",
            "The join assumes identical k-point ordering from the same saved pw.x calculation.",
            "Band energies and k-point coordinates come from BandStructureData.",
        ),
    )
