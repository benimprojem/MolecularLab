"""Plot-ready transformations for atom/orbital-resolved fat-band datasets.

This module prepares data only. It does not depend on Matplotlib or a GUI and
never changes the raw QE projection weights held by FatBandDataset.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt
from typing import Iterable

from molecular_lab.calculations.dft_settings import BandPath
from molecular_lab.calculations.fatband import FatBandDataset


@dataclass(frozen=True, slots=True)
class FatBandPlotPoint:
    kpoint_index: int
    path_distance_inverse_angstrom: float
    energy_ev: float
    raw_weight: float
    display_weight: float

    def __post_init__(self) -> None:
        if self.kpoint_index < 0:
            raise ValueError("kpoint_index must be zero-based and non-negative")
        values = (
            self.path_distance_inverse_angstrom,
            self.energy_ev,
            self.raw_weight,
            self.display_weight,
        )
        if any(not isfinite(value) for value in values):
            raise ValueError("fat-band plot point values must be finite")
        if self.path_distance_inverse_angstrom < 0 or self.raw_weight < 0 or self.display_weight < 0:
            raise ValueError("path distance and projection weights cannot be negative")


@dataclass(frozen=True, slots=True)
class FatBandPlotSeries:
    group_key: str
    label: str
    band_index: int
    spin_channel: str
    points: tuple[FatBandPlotPoint, ...]

    def __post_init__(self) -> None:
        if not self.group_key or not self.label or not self.spin_channel:
            raise ValueError("fat-band plot series requires group, label and spin identifiers")
        if self.band_index < 1 or not self.points:
            raise ValueError("fat-band plot series requires a positive band and plot points")
        object.__setattr__(self, "points", tuple(self.points))


@dataclass(frozen=True, slots=True)
class FatBandSymmetryTick:
    label: str
    path_distance_inverse_angstrom: float
    kpoint_index: int

    def __post_init__(self) -> None:
        if not self.label or self.kpoint_index < 0:
            raise ValueError("symmetry tick requires a label and non-negative k-point index")
        if not isfinite(self.path_distance_inverse_angstrom) or self.path_distance_inverse_angstrom < 0:
            raise ValueError("symmetry tick distance must be finite and non-negative")


@dataclass(frozen=True, slots=True)
class FatBandPlotData:
    series: tuple[FatBandPlotSeries, ...]
    symmetry_ticks: tuple[FatBandSymmetryTick, ...]
    kpoint_distances_inverse_angstrom: tuple[float, ...]
    group_by: str
    normalization: str
    minimum_raw_weight: float
    x_axis_label: str = "Wave-vector path (Å⁻¹)"
    y_axis_label: str = "Energy (eV)"
    weight_axis_label: str = "Displayed projection weight"
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        distances = tuple(self.kpoint_distances_inverse_angstrom)
        if not distances or any(not isfinite(v) or v < 0 for v in distances):
            raise ValueError("plot data requires finite non-negative k-point path distances")
        if any(b < a for a, b in zip(distances, distances[1:])):
            raise ValueError("k-point path distances must be non-decreasing")
        if self.group_by not in {"atom_orbital", "element_orbital", "element", "atom", "orbital"}:
            raise ValueError("unsupported fat-band grouping")
        if self.normalization not in {"raw", "global_max", "per_band_max"}:
            raise ValueError("unsupported fat-band normalization")
        if not isfinite(self.minimum_raw_weight) or self.minimum_raw_weight < 0:
            raise ValueError("minimum_raw_weight must be finite and non-negative")
        object.__setattr__(self, "series", tuple(self.series))
        object.__setattr__(self, "symmetry_ticks", tuple(self.symmetry_ticks))
        object.__setattr__(self, "kpoint_distances_inverse_angstrom", distances)
        object.__setattr__(self, "notes", tuple(self.notes))


def _cross(a: tuple[float, float, float], b: tuple[float, float, float]):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def _reciprocal_lattice(
    lattice_vectors_angstrom: tuple[tuple[float, float, float], ...],
) -> tuple[tuple[float, float, float], ...]:
    if len(lattice_vectors_angstrom) != 3 or any(len(row) != 3 for row in lattice_vectors_angstrom):
        raise ValueError("lattice_vectors_angstrom must be a 3x3 matrix")
    a1, a2, a3 = tuple(tuple(float(v) for v in row) for row in lattice_vectors_angstrom)
    volume = _dot(a1, _cross(a2, a3))
    if not isfinite(volume) or abs(volume) < 1e-12:
        raise ValueError("direct lattice vectors must span a non-zero cell volume")
    factor = 2.0 * 3.141592653589793 / volume
    return tuple(
        tuple(factor * value for value in vector)
        for vector in (_cross(a2, a3), _cross(a3, a1), _cross(a1, a2))
    )


def _interpolate_path(path: BandPath) -> tuple[tuple[float, float, float], ...]:
    coordinates: list[tuple[float, float, float]] = []
    for index, start in enumerate(path.points[:-1]):
        end = path.points[index + 1]
        subdivisions = start.points_to_next
        for step in range(subdivisions):
            fraction = step / subdivisions
            coordinates.append(tuple(
                start.coordinates[axis]
                + fraction * (end.coordinates[axis] - start.coordinates[axis])
                for axis in range(3)
            ))
    coordinates.append(path.points[-1].coordinates)
    return tuple(coordinates)


def _group(contribution, group_by: str) -> tuple[str, str]:
    atom = contribution.atom_index
    element = contribution.element
    orbital = contribution.orbital_label
    if group_by == "atom_orbital":
        return f"atom:{atom}:{element}:{orbital}", f"{element}{atom} · {orbital}"
    if group_by == "element_orbital":
        return f"element:{element}:{orbital}", f"{element} · {orbital}"
    if group_by == "element":
        return f"element:{element}", element
    if group_by == "atom":
        return f"atom:{atom}:{element}", f"{element}{atom}"
    if group_by == "orbital":
        return f"orbital:{orbital}", orbital
    raise ValueError(f"unsupported fat-band grouping: {group_by}")


def prepare_fatband_plot_data(
    dataset: FatBandDataset,
    *,
    band_path: BandPath,
    lattice_vectors_angstrom: tuple[tuple[float, float, float], ...],
    group_by: str = "element_orbital",
    normalization: str = "global_max",
    minimum_raw_weight: float = 0.0,
    elements: Iterable[str] | None = None,
    atom_indices: Iterable[int] | None = None,
    orbital_labels: Iterable[str] | None = None,
    spin_channels: Iterable[str] | None = None,
    band_indices: Iterable[int] | None = None,
) -> FatBandPlotData:
    """Convert a FatBandDataset into grouped, distance-aware plot series.

    The raw mode keeps |c|^2 unchanged for display. global_max divides display
    weights by the largest selected raw weight. per_band_max scales each
    band/spin set independently. Raw values are always retained alongside
    display values; no normalization is written back to the source dataset.

    Distances are cumulative Cartesian reciprocal-space distances in inverse
    angstrom, calculated from the supplied direct lattice vectors and
    fractional coordinates. The supplied BandPath must match dataset order.
    """
    if group_by not in {"atom_orbital", "element_orbital", "element", "atom", "orbital"}:
        raise ValueError("unsupported fat-band grouping")
    if normalization not in {"raw", "global_max", "per_band_max"}:
        raise ValueError("unsupported fat-band normalization")
    if not isfinite(minimum_raw_weight) or minimum_raw_weight < 0:
        raise ValueError("minimum_raw_weight must be finite and non-negative")

    expected_path = _interpolate_path(band_path)
    if len(expected_path) != len(dataset.kpoints):
        raise ValueError(
            f"BandPath defines {len(expected_path)} k-points but fat-band data contains {len(dataset.kpoints)}."
        )
    for index, (expected, actual) in enumerate(zip(expected_path, dataset.kpoints)):
        if actual.coordinate_system != "fractional":
            raise ValueError("fat-band plot distances require fractional k-point coordinates")
        if any(abs(a - b) > 1e-8 for a, b in zip(expected, actual.coordinates)):
            raise ValueError(f"Fat-band k-point {index} does not match the supplied BandPath")

    reciprocal = _reciprocal_lattice(lattice_vectors_angstrom)
    cartesian = tuple(
        tuple(sum(point[axis] * reciprocal[axis][component] for axis in range(3))
              for component in range(3))
        for point in expected_path
    )
    distances = [0.0]
    for previous, current in zip(cartesian, cartesian[1:]):
        delta = tuple(current[i] - previous[i] for i in range(3))
        distances.append(distances[-1] + sqrt(_dot(delta, delta)))

    selected_elements = set(elements) if elements is not None else None
    selected_atoms = set(atom_indices) if atom_indices is not None else None
    selected_orbitals = set(orbital_labels) if orbital_labels is not None else None
    selected_spins = set(spin_channels) if spin_channels is not None else set(dataset.spin_channels)
    selected_bands = set(band_indices) if band_indices is not None else set(range(1, dataset.number_of_bands + 1))
    unknown_spins = selected_spins.difference(dataset.spin_channels)
    if unknown_spins:
        raise ValueError(f"unknown spin channels requested: {sorted(unknown_spins)}")
    invalid_bands = [band for band in selected_bands if band < 1 or band > dataset.number_of_bands]
    if invalid_bands:
        raise ValueError(f"band indices out of range: {sorted(invalid_bands)}")

    aggregated: dict[tuple[str, int, str, int], float] = {}
    labels: dict[str, str] = {}
    energies: dict[tuple[str, int, int], float] = {}
    for point in dataset.points:
        if point.spin_channel not in selected_spins or point.band_index not in selected_bands:
            continue
        energies[(point.spin_channel, point.band_index, point.kpoint_index)] = point.energy_ev
        for contribution in point.contributions:
            if selected_elements is not None and contribution.element not in selected_elements:
                continue
            if selected_atoms is not None and contribution.atom_index not in selected_atoms:
                continue
            if selected_orbitals is not None and contribution.orbital_label not in selected_orbitals:
                continue
            key, label = _group(contribution, group_by)
            labels[key] = label
            aggregate_key = (key, point.band_index, point.spin_channel, point.kpoint_index)
            aggregated[aggregate_key] = aggregated.get(aggregate_key, 0.0) + contribution.raw_weight

    raw_values = tuple(aggregated.values())
    global_max = max(raw_values, default=0.0)
    per_band_max: dict[tuple[str, int], float] = {}
    for (key, band, spin, kpoint), weight in aggregated.items():
        scope = (spin, band)
        per_band_max[scope] = max(per_band_max.get(scope, 0.0), weight)

    series: list[FatBandPlotSeries] = []
    for group_key in sorted(labels):
        for spin in dataset.spin_channels:
            if spin not in selected_spins:
                continue
            for band in sorted(selected_bands):
                if (spin, band, 0) not in energies:
                    continue
                points = []
                for kpoint_index, distance in enumerate(distances):
                    raw = aggregated.get((group_key, band, spin, kpoint_index), 0.0)
                    if raw < minimum_raw_weight:
                        display = 0.0
                    elif normalization == "raw":
                        display = raw
                    elif normalization == "global_max":
                        display = raw / global_max if global_max > 0 else 0.0
                    else:
                        maximum = per_band_max.get((spin, band), 0.0)
                        display = raw / maximum if maximum > 0 else 0.0
                    points.append(FatBandPlotPoint(
                        kpoint_index=kpoint_index,
                        path_distance_inverse_angstrom=distance,
                        energy_ev=energies[(spin, band, kpoint_index)],
                        raw_weight=raw,
                        display_weight=display,
                    ))
                series.append(FatBandPlotSeries(
                    group_key=group_key,
                    label=labels[group_key],
                    band_index=band,
                    spin_channel=spin,
                    points=tuple(points),
                ))

    ticks: list[FatBandSymmetryTick] = []
    path_index = 0
    for index, point in enumerate(band_path.points):
        if index:
            path_index += band_path.points[index - 1].points_to_next
        ticks.append(FatBandSymmetryTick(
            label=point.label,
            path_distance_inverse_angstrom=distances[path_index],
            kpoint_index=path_index,
        ))

    return FatBandPlotData(
        series=tuple(series),
        symmetry_ticks=tuple(ticks),
        kpoint_distances_inverse_angstrom=tuple(distances),
        group_by=group_by,
        normalization=normalization,
        minimum_raw_weight=minimum_raw_weight,
        notes=(
            "Projection weights are unchanged in the source FatBandDataset.",
            "display_weight is a visualization scale only; raw_weight retains the original |c|^2 sum.",
            "Path distances use the supplied direct lattice and fractional BandPath coordinates.",
            "Series are separated by group, band and spin channel.",
        ),
    )
