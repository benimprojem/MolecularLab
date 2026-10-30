"""Quantum ESPRESSO band and DOS output parsing."""
from __future__ import annotations

from pathlib import Path
import re

from molecular_lab.materials.electronic_structure import (
    BandChannelData, BandStructureData, DensityOfStatesChannel,
    DensityOfStatesData, KPoint,
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
