"""Periodic crystal structure records for backend-independent atomistic workflows.

Coordinates are fractional and lattice vectors are stored as row vectors in angstrom.
This module represents structures; it does not generate or validate physical crystals.
Version: 0.1 | Rev.: 00
"""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from math import isfinite

class DefectKind(str, Enum):
    VACANCY = "vacancy"
    SUBSTITUTION = "substitution"
    INTERSTITIAL = "interstitial"
    ANTISITE = "antisite"

@dataclass(frozen=True, slots=True)
class CrystalSite:
    """One crystallographic site; occupancy is the fraction assigned to this species."""
    site_id: str
    element: str
    fractional_coordinates: tuple[float, float, float]
    occupancy: float = 1.0
    oxidation_state: int | None = None
    magnetic_moment_mu_b: float | None = None
    label: str = ""
    def __post_init__(self):
        coordinates = tuple(self.fractional_coordinates)
        if not self.site_id.strip():
            raise ValueError("site_id cannot be empty")
        if not self.element or not self.element[0].isupper() or not self.element.isalpha():
            raise ValueError("element must be a canonical element-like symbol")
        if len(coordinates) != 3 or any(not isfinite(v) for v in coordinates):
            raise ValueError("fractional coordinates must contain three finite values")
        if not isfinite(self.occupancy) or not 0.0 < self.occupancy <= 1.0:
            raise ValueError("site occupancy must be in (0, 1]")
        if self.oxidation_state is not None and (isinstance(self.oxidation_state, bool) or not isinstance(self.oxidation_state, int)):
            raise ValueError("oxidation_state must be an integer or None")
        if self.magnetic_moment_mu_b is not None and not isfinite(self.magnetic_moment_mu_b):
            raise ValueError("magnetic moment must be finite")
        object.__setattr__(self, "fractional_coordinates", coordinates)

@dataclass(frozen=True, slots=True)
class CrystalDefect:
    """Defect annotation relative to a declared host crystal structure."""
    defect_id: str
    kind: DefectKind
    host_site_id: str | None = None
    inserted_element: str | None = None
    fractional_coordinates: tuple[float, float, float] | None = None
    charge_state: int | None = None
    notes: str = ""
    def __post_init__(self):
        if not self.defect_id.strip():
            raise ValueError("defect_id cannot be empty")
        if not isinstance(self.kind, DefectKind):
            object.__setattr__(self, "kind", DefectKind(self.kind))
        if self.charge_state is not None and (isinstance(self.charge_state, bool) or not isinstance(self.charge_state, int)):
            raise ValueError("defect charge_state must be an integer or None")
        if self.kind is DefectKind.INTERSTITIAL:
            if self.host_site_id is not None:
                raise ValueError("interstitial defect must not reference a host site")
            if not self.inserted_element or not self.inserted_element[0].isupper() or not self.inserted_element.isalpha():
                raise ValueError("interstitial defect requires an inserted element")
            if self.fractional_coordinates is None:
                raise ValueError("interstitial defect requires fractional coordinates")
            coordinates = tuple(self.fractional_coordinates)
            if len(coordinates) != 3 or any(not isfinite(v) for v in coordinates):
                raise ValueError("interstitial coordinates must contain three finite values")
            object.__setattr__(self, "fractional_coordinates", coordinates)
        else:
            if not self.host_site_id or not self.host_site_id.strip():
                raise ValueError(f"{self.kind.value} defect requires host_site_id")
            if self.fractional_coordinates is not None:
                raise ValueError("site-based defect must not supply separate fractional coordinates")
            if self.kind in (DefectKind.SUBSTITUTION, DefectKind.ANTISITE):
                if not self.inserted_element or not self.inserted_element[0].isupper() or not self.inserted_element.isalpha():
                    raise ValueError(f"{self.kind.value} defect requires inserted_element")
            elif self.inserted_element is not None:
                raise ValueError("vacancy defect must not specify inserted_element")

@dataclass(frozen=True, slots=True)
class CrystalStructure:
    """Periodic unit cell; lattice vectors are rows and coordinates are fractional."""
    structure_id: str
    lattice_vectors_angstrom: tuple[tuple[float, float, float], ...]
    sites: tuple[CrystalSite, ...]
    periodic_boundary_conditions: tuple[bool, bool, bool] = (True, True, True)
    defects: tuple[CrystalDefect, ...] = ()
    space_group: str | None = None
    magnetic_ordering: str | None = None
    total_charge: int | None = 0
    source: str = ""
    notes: tuple[str, ...] = ()
    def __post_init__(self):
        lattice = tuple(tuple(row) for row in self.lattice_vectors_angstrom)
        sites = tuple(self.sites)
        pbc = tuple(self.periodic_boundary_conditions)
        defects = tuple(self.defects)
        if not self.structure_id.strip():
            raise ValueError("structure_id cannot be empty")
        if len(lattice) != 3 or any(len(row) != 3 for row in lattice):
            raise ValueError("lattice_vectors_angstrom must be a 3x3 matrix")
        if any(not isfinite(v) for row in lattice for v in row):
            raise ValueError("lattice vectors must contain finite values")
        a, b, c = lattice
        volume = (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))
        if volume <= 1e-10:
            raise ValueError("lattice cell must have a positive, non-zero volume")
        if not sites:
            raise ValueError("crystal structure must contain at least one site")
        if len(pbc) != 3 or any(type(flag) is not bool for flag in pbc):
            raise ValueError("periodic_boundary_conditions must contain three booleans")
        if not any(pbc):
            raise ValueError("at least one periodic direction is required for a crystal structure")
        ids = [site.site_id for site in sites]
        if len(ids) != len(set(ids)):
            raise ValueError("crystal site IDs must be unique")
        known_ids = set(ids)
        occupancy_by_position = {}
        for site in sites:
            key = tuple(
                round(value % 1.0, 10) if pbc[index] else round(value, 10)
                for index, value in enumerate(site.fractional_coordinates)
            )
            occupancy_by_position[key] = occupancy_by_position.get(key, 0.0) + site.occupancy
        if any(total > 1.0 + 1e-8 for total in occupancy_by_position.values()):
            raise ValueError("occupancies for species sharing one periodic site cannot sum above one")
        defect_ids = [defect.defect_id for defect in defects]
        if len(defect_ids) != len(set(defect_ids)):
            raise ValueError("defect IDs must be unique")
        for defect in defects:
            if defect.host_site_id is not None and defect.host_site_id not in known_ids:
                raise ValueError(f"defect references unknown host site: {defect.host_site_id}")
        if self.total_charge is not None and (isinstance(self.total_charge, bool) or not isinstance(self.total_charge, int)):
            raise ValueError("total_charge must be an integer or None")
        object.__setattr__(self, "lattice_vectors_angstrom", lattice)
        object.__setattr__(self, "sites", sites)
        object.__setattr__(self, "periodic_boundary_conditions", pbc)
        object.__setattr__(self, "defects", defects)
        object.__setattr__(self, "notes", tuple(self.notes))
    @property
    def cell_volume_angstrom3(self):
        a, b, c = self.lattice_vectors_angstrom
        return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))
    @property
    def cartesian_coordinates_angstrom(self):
        """Cartesian coordinates in site order, using r = fractional @ lattice."""
        lattice = self.lattice_vectors_angstrom
        return tuple(
            tuple(sum(site.fractional_coordinates[i] * lattice[i][j] for i in range(3)) for j in range(3))
            for site in self.sites
        )

    @property
    def composition_by_site(self):
        """Occupancy-weighted site counts, not a reduced empirical formula."""
        counts = {}
        for site in self.sites:
            counts[site.element] = counts.get(site.element, 0.0) + site.occupancy
        return dict(sorted(counts.items()))
