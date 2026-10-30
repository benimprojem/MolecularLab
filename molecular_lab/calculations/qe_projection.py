"""Parse k-resolved atomic projections written by Quantum ESPRESSO projwfc.x.

Supports the legacy PROJECTIONS layout and the newer EIGENSTATES/PROJS layout.
Complex amplitudes and raw |c|**2 weights are both retained; no normalization
is applied.
"""
from __future__ import annotations

from dataclasses import dataclass
import re
import math
import xml.etree.ElementTree as ET

from molecular_lab.calculations.dft_settings import BandPath
from molecular_lab.materials.electronic_structure import (
    BandChannelData, BandStructureData, KPoint,
)


class QEProjectionParseError(ValueError):
    """Raised when a QE atomic-projection file cannot be interpreted safely."""


@dataclass(frozen=True, slots=True)
class AtomicWavefunction:
    """One pseudopotential atomic-wavefunction projector (one-based indices)."""
    index: int
    atom_index: int
    element: str
    wavefunction_index: int
    angular_momentum: int
    magnetic_number: int | None = None
    total_angular_momentum: float | None = None
    source_label: str = ""

    def __post_init__(self) -> None:
        if self.index < 1 or self.atom_index < 1 or self.wavefunction_index < 1:
            raise ValueError("QE atomic-wavefunction indices must be positive")
        if not self.element:
            raise ValueError("element cannot be empty")
        if self.angular_momentum < 0:
            raise ValueError("angular_momentum cannot be negative")


@dataclass(frozen=True, slots=True)
class QEAtomicProjectionCoefficient:
    """Complex projection coefficient for one spin, k-point, band and projector."""
    kpoint_index: int
    spin_channel: str
    band_index: int
    atomic_wfc_index: int
    real: float
    imaginary: float
    weight: float

    def __post_init__(self) -> None:
        if isinstance(self.kpoint_index, bool) or not isinstance(self.kpoint_index, int) or self.kpoint_index < 0:
            raise ValueError("kpoint_index must be a non-negative integer")
        if isinstance(self.band_index, bool) or not isinstance(self.band_index, int) or self.band_index < 1:
            raise ValueError("band_index must be a positive one-based integer")
        if isinstance(self.atomic_wfc_index, bool) or not isinstance(self.atomic_wfc_index, int) or self.atomic_wfc_index < 1:
            raise ValueError("atomic_wfc_index must be a positive one-based integer")
        if not isinstance(self.spin_channel, str) or not self.spin_channel.strip():
            raise ValueError("spin_channel cannot be empty")
        if not all(isinstance(value, (int, float)) and not isinstance(value, bool)
                   and math.isfinite(value)
                   for value in (self.real, self.imaginary, self.weight)):
            raise ValueError("projection coefficient values must be finite numbers")
        if self.weight < 0:
            raise ValueError("projection weight cannot be negative")
        expected = self.real * self.real + self.imaginary * self.imaginary
        if not math.isclose(self.weight, expected, rel_tol=1e-9, abs_tol=1e-12):
            raise ValueError("projection weight must equal real² + imaginary²")

    @property
    def amplitude(self) -> complex:
        return complex(self.real, self.imaginary)


@dataclass(frozen=True, slots=True)
class QEAtomicProjectionData:
    """Normalized k-resolved QE atomic-projection data."""
    number_of_bands: int
    number_of_kpoints: int
    spin_channels: tuple[str, ...]
    atomic_wavefunctions: tuple[AtomicWavefunction, ...]
    coefficients: tuple[QEAtomicProjectionCoefficient, ...]
    format_generation: str

    def __post_init__(self) -> None:
        if (isinstance(self.number_of_bands, bool) or not isinstance(self.number_of_bands, int)
                or isinstance(self.number_of_kpoints, bool) or not isinstance(self.number_of_kpoints, int)
                or self.number_of_bands < 1 or self.number_of_kpoints < 1):
            raise ValueError("band and k-point counts must be positive integers")
        spins = tuple(self.spin_channels)
        wavefunctions = tuple(self.atomic_wavefunctions)
        coefficients = tuple(self.coefficients)
        if not spins or len(spins) != len(set(spins)) or any(not spin.strip() for spin in spins):
            raise ValueError("spin channels must be non-empty and unique")
        indices = [item.index for item in wavefunctions]
        if indices != list(range(1, len(wavefunctions) + 1)):
            raise ValueError("atomic wavefunction indices must be contiguous and one-based")
        atom_elements: dict[int, str] = {}
        atom_wfcs: set[tuple[int, int]] = set()
        for item in wavefunctions:
            previous = atom_elements.setdefault(item.atom_index, item.element)
            if previous != item.element:
                raise ValueError(f"atom {item.atom_index} has conflicting element metadata")
            atom_wfc = (item.atom_index, item.wavefunction_index)
            if atom_wfc in atom_wfcs:
                raise ValueError(f"duplicate wavefunction index for atom {item.atom_index}")
            atom_wfcs.add(atom_wfc)
        if sorted(atom_elements) != list(range(1, max(atom_elements) + 1)):
            raise ValueError("atom indices in projector metadata must be contiguous and one-based")
        expected = {
            (k, spin, band, wfc)
            for k in range(self.number_of_kpoints)
            for spin in spins
            for band in range(1, self.number_of_bands + 1)
            for wfc in indices
        }
        observed = [
            (item.kpoint_index, item.spin_channel, item.band_index, item.atomic_wfc_index)
            for item in coefficients
        ]
        if len(observed) != len(set(observed)):
            raise ValueError("duplicate projection coefficient for a state/projector key")
        if set(observed) != expected:
            missing = len(expected - set(observed))
            extra = len(set(observed) - expected)
            raise ValueError(
                f"projection coefficient grid is incomplete or out of range "
                f"({missing} missing, {extra} unexpected)"
            )
        object.__setattr__(self, "spin_channels", spins)
        object.__setattr__(self, "atomic_wavefunctions", wavefunctions)
        object.__setattr__(self, "coefficients", coefficients)


_STATE_RE = re.compile(
    r"state\s*#\s*(?P<state>\d+)\s*:\s*atom\s+(?P<atom>\d+)\s*"
    r"\((?P<element>[^)]+)\)\s*,\s*wfc\s+(?P<wfc>\d+)\s*"
    r"\((?P<details>[^)]*)\)", re.IGNORECASE,
)
_NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")


def _tag(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def _child(parent: ET.Element, name: str) -> ET.Element | None:
    return next((item for item in parent if _tag(item) == name), None)


def _parse_wavefunction_metadata(output: str) -> tuple[AtomicWavefunction, ...]:
    records: list[AtomicWavefunction] = []
    for match in _STATE_RE.finditer(output):
        details = match.group("details")
        l_match = re.search(r"\bl\s*=\s*(\d+)", details, re.IGNORECASE)
        m_match = re.search(r"\bm\s*=\s*(-?\d+)", details, re.IGNORECASE)
        j_match = re.search(r"\bj\s*=\s*([-+]?\d+(?:\.\d+)?)", details, re.IGNORECASE)
        if l_match is None:
            raise QEProjectionParseError(
                f"Could not read angular momentum from state line: {match.group(0)}"
            )
        records.append(AtomicWavefunction(
            index=int(match.group("state")),
            atom_index=int(match.group("atom")),
            element=match.group("element").strip(),
            wavefunction_index=int(match.group("wfc")),
            angular_momentum=int(l_match.group(1)),
            magnetic_number=int(m_match.group(1)) if m_match else None,
            total_angular_momentum=float(j_match.group(1)) if j_match else None,
            source_label=details.strip(),
        ))
    if not records:
        raise QEProjectionParseError(
            "No state # ...: atom ... (Element), wfc ... (l=...) records "
            "were found in projwfc.x standard output."
        )
    indices = [item.index for item in records]
    if len(indices) != len(set(indices)):
        raise QEProjectionParseError("Duplicate state indices in projwfc.x output")
    return tuple(sorted(records, key=lambda item: item.index))


def _header_value(header: ET.Element, key: str) -> str | None:
    value = header.attrib.get(key)
    if value is not None:
        return value
    child = _child(header, key)
    return child.text.strip() if child is not None and child.text else None


def _parse_complex_pairs(text: str | None, expected_bands: int, context: str) -> tuple[complex, ...]:
    if not text:
        raise QEProjectionParseError(f"Missing projection coefficients in {context}")
    normalized = text.replace(",", " ").replace("D", "E").replace("d", "e")
    tokens = _NUMBER_RE.findall(normalized)
    if len(tokens) != expected_bands * 2:
        raise QEProjectionParseError(
            f"{context} contains {len(tokens)} numeric values; expected {expected_bands * 2}."
        )
    values = [float(token) for token in tokens]
    return tuple(complex(values[i], values[i + 1]) for i in range(0, len(values), 2))


def _spin_names(count: int) -> tuple[str, ...]:
    if count == 1:
        return ("total",)
    if count == 2:
        return ("spin_up", "spin_down")
    if count == 4:
        raise QEProjectionParseError(
            "QE reports four spin components (non-collinear/spinor data); "
            "spinor projection parsing is not implemented."
        )
    raise QEProjectionParseError(f"Unsupported NUMBER_OF_SPIN_COMPONENTS={count}")


def parse_qe_atomic_proj_xml(xml_text: str, projwfc_output: str) -> QEAtomicProjectionData:
    """Parse atomic_proj.xml and projwfc.x stdout; energies are supplied separately.

    The XML supplies complex coefficients. The projwfc.x output maps each
    atomic-wavefunction index to atom, element and angular-momentum metadata.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        raise QEProjectionParseError(f"Invalid atomic_proj.xml: {exc}") from exc

    header = _child(root, "HEADER")
    if header is None:
        raise QEProjectionParseError("atomic_proj.xml has no HEADER element")

    keys = ("NUMBER_OF_BANDS", "NUMBER_OF_K-POINTS",
            "NUMBER_OF_SPIN_COMPONENTS", "NUMBER_OF_ATOMIC_WFC")
    dimensions: dict[str, int] = {}
    for key in keys:
        raw = _header_value(header, key)
        if raw is None:
            raise QEProjectionParseError(f"Missing {key} in atomic_proj.xml HEADER")
        try:
            dimensions[key] = int(raw)
        except ValueError as exc:
            raise QEProjectionParseError(f"Invalid {key} value: {raw!r}") from exc

    nbands = dimensions["NUMBER_OF_BANDS"]
    nkpoints = dimensions["NUMBER_OF_K-POINTS"]
    nspin = dimensions["NUMBER_OF_SPIN_COMPONENTS"]
    nwfc = dimensions["NUMBER_OF_ATOMIC_WFC"]
    if min(nbands, nkpoints, nwfc) < 1:
        raise QEProjectionParseError("QE XML header contains a non-positive dimension")
    spins = _spin_names(nspin)
    wavefunctions = _parse_wavefunction_metadata(projwfc_output)
    if len(wavefunctions) != nwfc:
        raise QEProjectionParseError(
            f"XML declares {nwfc} atomic wavefunctions but projwfc.x output "
            f"defines {len(wavefunctions)}."
        )

    coefficients: list[QEAtomicProjectionCoefficient] = []
    eigenstates = _child(root, "EIGENSTATES")
    if eigenstates is not None:
        # Newer QE: PROJS blocks are spin-major, k-point-minor.
        blocks = [item for item in eigenstates if _tag(item) == "PROJS"]
        expected_blocks = nkpoints * len(spins)
        if len(blocks) != expected_blocks:
            raise QEProjectionParseError(
                f"EIGENSTATES contains {len(blocks)} PROJS blocks; expected {expected_blocks}."
            )
        for block_index, block in enumerate(blocks):
            spin_index, kpoint_index = divmod(block_index, nkpoints)
            seen: set[int] = set()
            for projector in block:
                if _tag(projector) != "ATOMIC_WFC":
                    continue
                try:
                    wfc_index = int(projector.attrib["index"])
                except (KeyError, ValueError) as exc:
                    raise QEProjectionParseError(
                        "ATOMIC_WFC is missing a valid index attribute"
                    ) from exc
                if wfc_index < 1 or wfc_index > nwfc or wfc_index in seen:
                    raise QEProjectionParseError(f"Invalid or duplicate atomic wavefunction index {wfc_index}")
                seen.add(wfc_index)
                amplitudes = _parse_complex_pairs(
                    projector.text, nbands,
                    f"spin={spins[spin_index]}, k-point={kpoint_index}, wfc={wfc_index}",
                )
                for band_index, amplitude in enumerate(amplitudes, start=1):
                    coefficients.append(QEAtomicProjectionCoefficient(
                        kpoint_index, spins[spin_index], band_index, wfc_index,
                        amplitude.real, amplitude.imag,
                        amplitude.real ** 2 + amplitude.imag ** 2,
                    ))
            if seen != set(range(1, nwfc + 1)):
                raise QEProjectionParseError(
                    f"Incomplete projector set at spin={spins[spin_index]}, k-point={kpoint_index}."
                )
        generation = "eigenstates_projs"
    else:
        # Legacy QE: PROJECTIONS/K-POINT.N/(SPIN.N)/ATMWFC.N.
        projection_root = _child(root, "PROJECTIONS")
        if projection_root is None:
            raise QEProjectionParseError("atomic_proj.xml contains neither EIGENSTATES nor PROJECTIONS")
        kpoint_blocks = list(projection_root)
        if len(kpoint_blocks) != nkpoints:
            raise QEProjectionParseError(
                f"PROJECTIONS contains {len(kpoint_blocks)} k-point blocks; expected {nkpoints}."
            )
        for kpoint_index, block in enumerate(kpoint_blocks):
            spin_blocks = [item for item in block if _tag(item).startswith("SPIN.")]
            if spin_blocks:
                if len(spin_blocks) != len(spins):
                    raise QEProjectionParseError("Legacy XML spin-block count does not match HEADER")
                selected = list(enumerate(spin_blocks))
            else:
                if len(spins) != 1:
                    raise QEProjectionParseError("Legacy spin-polarized XML has no SPIN.N blocks")
                selected = [(0, block)]
            for spin_index, spin_block in selected:
                seen: set[int] = set()
                for projector in spin_block:
                    tag = _tag(projector)
                    if not tag.startswith("ATMWFC."):
                        continue
                    try:
                        wfc_index = int(tag.split(".", 1)[1])
                    except ValueError as exc:
                        raise QEProjectionParseError(f"Invalid legacy projector tag: {tag}") from exc
                    if wfc_index < 1 or wfc_index > nwfc or wfc_index in seen:
                        raise QEProjectionParseError(f"Invalid or duplicate projector {tag}")
                    seen.add(wfc_index)
                    amplitudes = _parse_complex_pairs(
                        projector.text, nbands,
                        f"spin={spins[spin_index]}, k-point={kpoint_index}, wfc={wfc_index}",
                    )
                    for band_index, amplitude in enumerate(amplitudes, start=1):
                        coefficients.append(QEAtomicProjectionCoefficient(
                            kpoint_index, spins[spin_index], band_index, wfc_index,
                            amplitude.real, amplitude.imag,
                            amplitude.real ** 2 + amplitude.imag ** 2,
                        ))
                if seen != set(range(1, nwfc + 1)):
                    raise QEProjectionParseError(
                        f"Incomplete legacy projector set at spin={spins[spin_index]}, k-point={kpoint_index}."
                    )
        generation = "legacy_projections"

    expected = nbands * nkpoints * len(spins) * nwfc
    if len(coefficients) != expected:
        raise QEProjectionParseError(f"Parsed {len(coefficients)} coefficients; expected {expected}.")
    return QEAtomicProjectionData(
        nbands, nkpoints, spins, wavefunctions, tuple(coefficients), generation,
    )



def parse_qe_projwfc_band_output(
    output: str,
    *,
    fermi_level_ev: float | None = None,
    backend_version: str = "",
    exchange_correlation: str = "",
    band_path: BandPath | None = None,
) -> BandStructureData:
    """Extract original k-point coordinates and eigenvalues from projwfc.x stdout.

    This intentionally avoids bands.x reordered eigenvalue curves: atomic
    projection coefficients are indexed by the original pw.x eigenstate order.
    Keeping energies from projwfc.x stdout preserves that state identity at
    crossings. When a BandPath is supplied, exact fractional path coordinates
    are reconstructed from its vertices/subdivisions rather than assuming a
    unit convention for QE's printed k-vector. Only scalar, non-spin-polarized
    output is supported here.
    """
    number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
    kpoint_re = re.compile(
        rf"^\s*k\s*=\s*({number})\s+({number})\s+({number})\s*$",
        re.IGNORECASE | re.MULTILINE,
    )
    energy_re = re.compile(
        rf"====\s*e\(\s*(\d+)\s*\)\s*=\s*({number})\s*eV\s*====",
        re.IGNORECASE,
    )
    k_matches = list(kpoint_re.finditer(output))
    if not k_matches:
        raise QEProjectionParseError(
            "No reciprocal-space k-point records were found in projwfc.x output."
        )

    raw_kpoints: list[tuple[float, float, float]] = []
    band_rows: list[dict[int, float]] = []
    for index, match in enumerate(k_matches):
        coordinates = tuple(
            float(value.replace("D", "E").replace("d", "e"))
            for value in match.groups()
        )
        raw_kpoints.append(coordinates)
        end = k_matches[index + 1].start() if index + 1 < len(k_matches) else len(output)
        block = output[match.end():end]
        energies: dict[int, float] = {}
        for energy_match in energy_re.finditer(block):
            band_index = int(energy_match.group(1))
            if band_index in energies:
                raise QEProjectionParseError(
                    f"Duplicate eigenvalue for band {band_index} at k-point {index}."
                )
            energies[band_index] = float(
                energy_match.group(2).replace("D", "E").replace("d", "e")
            )
        if not energies:
            raise QEProjectionParseError(
                f"No eigenvalues were found for k-point {index} in projwfc.x output."
            )
        if sorted(energies) != list(range(1, len(energies) + 1)):
            raise QEProjectionParseError(
                f"Band indices at k-point {index} are not a contiguous one-based sequence."
            )
        band_rows.append(energies)

    if band_path is not None:
        path_coordinates: list[tuple[float, float, float]] = []
        for index, start in enumerate(band_path.points[:-1]):
            end = band_path.points[index + 1]
            subdivisions = start.points_to_next
            for step in range(subdivisions):
                fraction = step / subdivisions
                path_coordinates.append(tuple(
                    start.coordinates[axis]
                    + fraction * (end.coordinates[axis] - start.coordinates[axis])
                    for axis in range(3)
                ))
        path_coordinates.append(band_path.points[-1].coordinates)
        if len(path_coordinates) != len(raw_kpoints):
            raise QEProjectionParseError(
                f"BandPath defines {len(path_coordinates)} k-points but projwfc.x reports {len(raw_kpoints)}."
            )
        kpoints = tuple(KPoint(point, 1.0, "fractional") for point in path_coordinates)
    else:
        kpoints = tuple(KPoint(point, 1.0, "qe_output") for point in raw_kpoints)

    band_count = len(band_rows[0])
    if any(len(row) != band_count for row in band_rows):
        raise QEProjectionParseError(
            "projwfc.x output has inconsistent band counts between k-points."
        )
    energies_by_band = tuple(
        tuple(row[band] for row in band_rows)
        for band in range(1, band_count + 1)
    )
    return BandStructureData(
        kpoints=kpoints,
        channels=(BandChannelData("total", energies_by_band),),
        fermi_level_ev=fermi_level_ev,
        method="Quantum ESPRESSO projwfc.x original eigenstate energies",
        backend="quantum_espresso",
        backend_version=backend_version,
        exchange_correlation=exchange_correlation,
        notes=(
            "Eigenvalues are taken from projwfc.x stdout to preserve the original pw.x state ordering.",
            (
                "K-point coordinates are reconstructed in fractional crystal coordinates from the supplied BandPath."
                if band_path is not None
                else "K-point coordinates retain QE's printed representation; its unit convention is not inferred here."
            ),
            "The conventional bands.x output may reorder bands; do not use its band indices for raw projection joins.",
        ),
    )

def read_qe_atomic_proj_xml(xml_path: str, projwfc_output_path: str) -> QEAtomicProjectionData:
    """Read atomic_proj.xml and projwfc.x stdout from filesystem paths."""
    from pathlib import Path

    return parse_qe_atomic_proj_xml(
        Path(xml_path).read_text(encoding="utf-8"),
        Path(projwfc_output_path).read_text(encoding="utf-8"),
    )
