"""Optional ASE bridge for periodic crystal structures.

The bridge is deliberately strict: it refuses disorder and defect annotations
that cannot be represented faithfully as a plain ASE Atoms object.
"""
from __future__ import annotations

from importlib.util import find_spec
from typing import Any

from molecular_lab.structures.crystal import CrystalStructure


class ASEBridgeError(RuntimeError):
    """Base error raised while converting a crystal to ASE."""


class ASEUnavailableError(ASEBridgeError):
    """ASE is not installed or cannot be imported."""


class ASEUnsupportedStructureError(ASEBridgeError):
    """The current ASE bridge cannot represent the supplied structure."""


def ase_availability() -> tuple[bool, str]:
    """Return whether ASE appears importable without importing it eagerly."""
    try:
        available = find_spec("ase") is not None
    except (ImportError, ValueError):
        available = False
    return (available, "" if available else "Optional dependency 'ase' is not installed")


def to_ase_atoms(structure: CrystalStructure) -> Any:
    """Convert an ordered, fully occupied crystal structure to ase.Atoms.

    Lattice vectors use the same row-vector convention as CrystalStructure.
    Fractional coordinates are passed unchanged as scaled positions.
    """
    if not isinstance(structure, CrystalStructure):
        raise TypeError("structure must be a CrystalStructure")
    if structure.defects:
        raise ASEUnsupportedStructureError(
            "defect annotations must be materialized into an explicit structure first"
        )
    if any(site.occupancy != 1.0 for site in structure.sites):
        raise ASEUnsupportedStructureError(
            "partial occupancies/disorder require an explicit ordered supercell model"
        )

    available, reason = ase_availability()
    if not available:
        raise ASEUnavailableError(reason)

    try:
        from ase import Atoms
    except ImportError as exc:
        raise ASEUnavailableError("ASE was detected but could not be imported") from exc

    try:
        atoms = Atoms(
            symbols=[site.element for site in structure.sites],
            scaled_positions=[site.fractional_coordinates for site in structure.sites],
            cell=structure.lattice_vectors_angstrom,
            pbc=structure.periodic_boundary_conditions,
        )
    except (ValueError, TypeError) as exc:
        raise ASEUnsupportedStructureError(
            f"ASE rejected the crystal structure: {exc}"
        ) from exc

    atoms.info.update({
        "structure_id": structure.structure_id,
        "source": structure.source,
        "space_group": structure.space_group or "",
        "magnetic_ordering": structure.magnetic_ordering or "",
        "total_charge": structure.total_charge,
        "notes": list(structure.notes),
        "site_ids": [site.site_id for site in structure.sites],
        "oxidation_states": [site.oxidation_state for site in structure.sites],
    })
    if any(site.magnetic_moment_mu_b is not None for site in structure.sites):
        atoms.set_initial_magnetic_moments([
            0.0 if site.magnetic_moment_mu_b is None else site.magnetic_moment_mu_b
            for site in structure.sites
        ])
    return atoms
