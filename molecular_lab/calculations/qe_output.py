"""Quantum ESPRESSO band and DOS output parsing."""
from __future__ import annotations

from pathlib import Path
import re
from collections.abc import Mapping

from molecular_lab.materials.electronic_structure import (
    BandChannelData, BandStructureData, DensityOfStatesChannel,
    DensityOfStatesData, KPoint, ProjectedDOSComponent,
    ProjectedDOSRecord, ProjectedDensityOfStatesData,
)


class QEOutputParseError(ValueError):
    pass


def _float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_qe_fermi_level(output: str) -> float | None:
    pattern = r"the\s+Fermi\s+energy\s+is\s+([-+0-9.EeDd]+)\s*ev"
    matches = re.findall(pattern, output, re.IGNORECASE)
    return _float(matches[-1]) if matches else None


def parse_qe_bands_gnu(
    text: str, *, spin_channel: str = "non_spin_polarized",
    fermi_level_ev: float | None = None, backend_version: str = "",
    exchange_correlation: str = "", notes: tuple[str, ...] = (),
) -> BandStructureData:
    """Parse blank-separated bands.x .gnu blocks (path distance, energy in eV)."""
    blocks: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            if current:
                blocks.append(current)
                current = []
            continue
        if line.startswith(("#", "!", "&")):
            continue
        fields = line.split()
        if len(fields) < 2:
            continue
        try:
            current.append((_float(fields[0]), _float(fields[1])))
        except ValueError:
            continue
    if current:
        blocks.append(current)
    if not blocks or len(blocks[0]) < 2:
        raise QEOutputParseError("No usable band data found in QE .gnu output")
    reference = [row[0] for row in blocks[0]]
    for index, block in enumerate(blocks):
        if len(block) != len(reference):
            raise QEOutputParseError(f"Band {index} has an inconsistent sample count")
        if any(abs(row[0] - reference[i]) > 1e-8 for i, row in enumerate(block)):
            raise QEOutputParseError("Band blocks do not share the same k-path grid")
    kpoints = tuple(KPoint((x, 0.0, 0.0), 1.0, "path_distance") for x in reference)
    channel = BandChannelData(spin_channel, tuple(
        tuple(row[1] for row in block) for block in blocks
    ))
    return BandStructureData(
        kpoints=kpoints, channels=(channel,), fermi_level_ev=fermi_level_ev,
        method="Quantum ESPRESSO bands.x (.gnu)", backend="quantum_espresso",
        backend_version=backend_version, exchange_correlation=exchange_correlation,
        notes=tuple(notes) + (
            "K-point coordinate is cumulative path distance, not a 3D reciprocal vector.",
            "QE .gnu band energies are in eV; occupations are not included.",
        ),
    )


def read_qe_bands_gnu(path: str | Path, **kwargs) -> BandStructureData:
    return parse_qe_bands_gnu(Path(path).read_text(encoding="utf-8"), **kwargs)


def parse_qe_total_dos(
    text: str, *, spin_polarized: bool = False,
    fermi_level_ev: float | None = None, backend_version: str = "",
    notes: tuple[str, ...] = (),
) -> DensityOfStatesData:
    """Parse QE dos.x total DOS columns in eV and states/eV."""
    rows: list[list[float]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith(("#", "!", "&")):
            continue
        try:
            values = [_float(value) for value in line.split()]
        except ValueError:
            continue
        if len(values) >= 2:
            rows.append(values)
    required = 3 if spin_polarized else 2
    if len(rows) < 2 or any(len(row) < required for row in rows):
        raise QEOutputParseError("QE DOS table is incomplete for the requested spin mode")
    energies = tuple(row[0] for row in rows)
    if any(b <= a for a, b in zip(energies, energies[1:])):
        raise QEOutputParseError("QE DOS energy grid must be strictly increasing")
    if spin_polarized:
        channels = (
            DensityOfStatesChannel("spin_up", tuple(row[1] for row in rows)),
            DensityOfStatesChannel("spin_down", tuple(row[2] for row in rows)),
        )
    else:
        channels = (DensityOfStatesChannel("total", tuple(row[1] for row in rows)),)
    return DensityOfStatesData(
        energies_ev=energies, channels=channels, fermi_level_ev=fermi_level_ev,
        projection="total", normalization="per_cell",
        method="Quantum ESPRESSO dos.x", backend="quantum_espresso",
        backend_version=backend_version,
        notes=tuple(notes) + (
            "QE DOS energies are in eV and DOS values in states/eV.",
            "Additional integrated-DOS columns are not stored by the common model.",
        ),
    )


def read_qe_total_dos(path: str | Path, **kwargs) -> DensityOfStatesData:
    return parse_qe_total_dos(Path(path).read_text(encoding="utf-8"), **kwargs)


_QE_PDOS_FILE_RE = re.compile(
    r"^.+\\.pdos_atm#(?P<atom>\\d+)\\((?P<element>[^)]+)\\)"
    r"_wfc#(?P<wfc>\\d+)\\((?P<orbital>[spdf])\\)$",
    re.IGNORECASE,
)


def _qe_orbital_component_labels(angular_momentum: int) -> tuple[str, ...]:
    labels = ["m=0"]
    for magnetic_number in range(1, angular_momentum + 1):
        labels.extend((f"m={magnetic_number}, cos(phi)", f"m={magnetic_number}, sin(phi)"))
    return tuple(labels)


def parse_qe_projected_dos(
    files: Mapping[str, str],
    *,
    spin_polarized: bool = False,
    fermi_level_ev: float | None = None,
    backend_version: str = "",
    notes: tuple[str, ...] = (),
) -> ProjectedDensityOfStatesData:
    """Parse projwfc.x atomic-wavefunction PDOS files keyed by basename.

    Supports collinear scalar-relativistic s/p/d/f output. Fully relativistic
    j-resolved and k-resolved files are deliberately rejected/not inferred.
    """
    if not isinstance(files, Mapping) or not files:
        raise QEOutputParseError("No QE projected-DOS files were supplied")
    parsed: list[ProjectedDOSRecord] = []
    reference_energies: tuple[float, ...] | None = None
    for filename, content in sorted(files.items()):
        basename = Path(filename).name
        if ".pdos_atm#" not in basename:
            continue
        match = _QE_PDOS_FILE_RE.fullmatch(basename)
        if match is None:
            raise QEOutputParseError(
                f"Unsupported QE PDOS filename (j-resolved or unknown format): {basename}"
            )
        atom_index = int(match.group("atom"))
        element = match.group("element")
        wavefunction_index = int(match.group("wfc"))
        orbital = match.group("orbital").lower()
        angular_momentum = {"s": 0, "p": 1, "d": 2, "f": 3}[orbital]
        component_count = 2 * angular_momentum + 1
        required = 2 + component_count * (2 if spin_polarized else 1)
        rows: list[list[float]] = []
        for raw in content.splitlines():
            line = raw.strip()
            if not line or line.startswith(("#", "!", "&")):
                continue
            try:
                values = [_float(value) for value in line.split()]
            except ValueError:
                continue
            if len(values) < required:
                raise QEOutputParseError(
                    f"{basename} has {len(values)} columns; expected at least {required}"
                )
            rows.append(values)
        if len(rows) < 2:
            raise QEOutputParseError(f"{basename} contains fewer than two usable energy samples")
        energies = tuple(row[0] for row in rows)
        if any(b <= a for a, b in zip(energies, energies[1:])):
            raise QEOutputParseError(f"{basename} energy grid must be strictly increasing")
        if reference_energies is None:
            reference_energies = energies
        elif len(energies) != len(reference_energies) or any(
            abs(a - b) > 1e-8 for a, b in zip(energies, reference_energies)
        ):
            raise QEOutputParseError("QE PDOS files do not share a common energy grid")
        labels = _qe_orbital_component_labels(angular_momentum)
        spin_names = ("spin_up", "spin_down") if spin_polarized else ("total",)
        for spin_index, spin_name in enumerate(spin_names):
            if spin_polarized:
                local_column = 1 + spin_index
                component_start = 3 + spin_index
            else:
                local_column = 1
                component_start = 2
            local_values = tuple(row[local_column] for row in rows)
            components = tuple(
                ProjectedDOSComponent(
                    component_index=index,
                    component_label=labels[index],
                    density_states_per_ev=tuple(row[component_start + index * (2 if spin_polarized else 1)] for row in rows),
                )
                for index in range(component_count)
            )
            parsed.append(ProjectedDOSRecord(
                atom_index=atom_index, element=element,
                wavefunction_index=wavefunction_index,
                angular_momentum=angular_momentum, spin_channel=spin_name,
                local_dos_states_per_ev=local_values, components=components,
            ))
    if reference_energies is None or not parsed:
        raise QEOutputParseError("No supported QE atomic PDOS files were found")
    return ProjectedDensityOfStatesData(
        energies_ev=reference_energies, records=tuple(parsed),
        fermi_level_ev=fermi_level_ev,
        method="Quantum ESPRESSO projwfc.x atomic PDOS",
        backend="quantum_espresso", backend_version=backend_version,
        notes=tuple(notes) + (
            "Values are in states/eV; atomic projections are defined by the supplied pseudopotentials.",
            "Orbital components retain QE's angular-component ordering; labels are not assumed to be global Cartesian d-orbital names.",
            "Only collinear scalar-relativistic s/p/d/f files are parsed; j-resolved and k-resolved formats are unsupported.",
        ),
    )


def read_qe_projected_dos(
    paths: list[str | Path] | tuple[str | Path, ...], **kwargs
) -> ProjectedDensityOfStatesData:
    files = {Path(path).name: Path(path).read_text(encoding="utf-8") for path in paths}
    return parse_qe_projected_dos(files, **kwargs)
