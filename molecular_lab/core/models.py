"""
Molecular Lab domain models.

Version: 0.1
Rev.: 01
Status: Foundational composition, evidence, and molecular graph records.

These models intentionally use the Python standard library only. They describe
input and result records; they do not establish chemical stability.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from math import isfinite
from types import MappingProxyType
from typing import Mapping


class CompositionKind(str, Enum):
    MOLECULAR = "molecular"
    EMPIRICAL = "empirical"
    FORMULA_UNIT = "formula_unit"
    CELL = "cell"


class SearchDomain(str, Enum):
    MOLECULE = "molecule"
    PERIODIC = "periodic"
    MIXED = "mixed"


class EvidenceCategory(str, Enum):
    EXPERIMENTAL = "experimental"
    DATABASE = "database"
    THEORETICAL = "theoretical"
    HEURISTIC = "heuristic"
    HYPOTHESIS = "hypothesis"


@dataclass(frozen=True, slots=True)
class ElementRef:
    """Identity of a chemical element; properties are intentionally minimal."""
    atomic_number: int
    symbol: str

    def __post_init__(self) -> None:
        if self.atomic_number < 1:
            raise ValueError("atomic_number must be a positive integer")
        if not self.symbol or not self.symbol[0].isupper() or not self.symbol.isalpha():
            raise ValueError("symbol must be a canonical element-like symbol")


@dataclass(frozen=True, slots=True)
class IsotopeRef:
    """Specific isotope identified by atomic number and mass number."""
    atomic_number: int
    mass_number: int

    def __post_init__(self) -> None:
        if self.atomic_number < 1:
            raise ValueError("atomic_number must be positive")
        if self.mass_number < self.atomic_number:
            raise ValueError("mass_number cannot be smaller than atomic_number")


@dataclass(frozen=True, slots=True)
class Composition:
    """Element-count representation; isotope resolution is a separate future model."""
    counts: Mapping[str, int]
    net_charge: int | None = 0
    kind: CompositionKind = CompositionKind.MOLECULAR
    input_formula: str | None = None

    def __post_init__(self) -> None:
        if not self.counts:
            raise ValueError("composition must contain at least one element")
        normalized: dict[str, int] = {}
        for symbol, count in self.counts.items():
            if not symbol or not symbol[0].isupper() or not symbol.isalpha():
                raise ValueError(f"invalid element symbol: {symbol!r}")
            if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
                raise ValueError(f"count for {symbol} must be a positive integer")
            normalized[symbol] = count
        object.__setattr__(self, "counts", MappingProxyType(dict(sorted(normalized.items()))))
        if self.net_charge is not None and (
            isinstance(self.net_charge, bool) or not isinstance(self.net_charge, int)
        ):
            raise ValueError("net_charge must be an integer or None")
        if not isinstance(self.kind, CompositionKind):
            object.__setattr__(self, "kind", CompositionKind(self.kind))

    @property
    def total_atoms(self) -> int:
        return sum(self.counts.values())


@dataclass(frozen=True, slots=True)
class MolecularAtom:
    """One explicitly represented atom in a proposed molecular graph."""
    atom_id: int
    element: str

    def __post_init__(self) -> None:
        if isinstance(self.atom_id, bool) or self.atom_id < 0:
            raise ValueError("atom_id must be a non-negative integer")
        if not self.element or not self.element[0].isupper() or not self.element.isalpha():
            raise ValueError("element must be a canonical element-like symbol")


@dataclass(frozen=True, slots=True)
class MolecularBond:
    """An undirected bond between two distinct atom IDs."""
    atom_a: int
    atom_b: int
    order: int = 1

    def __post_init__(self) -> None:
        if self.atom_a == self.atom_b:
            raise ValueError("a bond cannot connect an atom to itself")
        if min(self.atom_a, self.atom_b) < 0:
            raise ValueError("atom IDs must be non-negative")
        if self.order not in (1, 2, 3):
            raise ValueError("bond order must be 1, 2, or 3")
        if self.atom_a > self.atom_b:
            a, b = self.atom_b, self.atom_a
            object.__setattr__(self, "atom_a", a)
            object.__setattr__(self, "atom_b", b)


@dataclass(frozen=True, slots=True)
class MolecularGraph:
    """Explicit connectivity proposal; geometry and stability are not implied."""
    atoms: tuple[MolecularAtom, ...]
    bonds: tuple[MolecularBond, ...]
    method: str
    limitations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        atoms = tuple(self.atoms)
        bonds = tuple(self.bonds)
        object.__setattr__(self, "atoms", atoms)
        object.__setattr__(self, "bonds", bonds)
        if not self.method.strip():
            raise ValueError("graph method cannot be empty")
        ids = [atom.atom_id for atom in atoms]
        if len(ids) != len(set(ids)):
            raise ValueError("atom IDs must be unique")
        known = set(ids)
        seen: set[tuple[int, int]] = set()
        for bond in bonds:
            if bond.atom_a not in known or bond.atom_b not in known:
                raise ValueError("bond references an unknown atom ID")
            pair = (bond.atom_a, bond.atom_b)
            if pair in seen:
                raise ValueError("duplicate bond between the same atoms")
            seen.add(pair)


@dataclass(frozen=True, slots=True)
class SearchConstraints:
    """Explicit bounds for a composition search; these are not chemistry rules."""
    allowed_elements: frozenset[str]
    max_total_atoms: int
    domain: SearchDomain = SearchDomain.MOLECULE
    min_total_atoms: int = 1
    per_element_max: Mapping[str, int] = field(default_factory=dict)
    allowed_charges: tuple[int, ...] | None = None

    def __post_init__(self) -> None:
        elements = frozenset(self.allowed_elements)
        if not elements:
            raise ValueError("allowed_elements cannot be empty")
        for symbol in elements:
            if not symbol or not symbol[0].isupper() or not symbol.isalpha():
                raise ValueError(f"invalid element symbol: {symbol!r}")
        object.__setattr__(self, "allowed_elements", elements)
        if self.min_total_atoms < 1:
            raise ValueError("min_total_atoms must be at least 1")
        if self.max_total_atoms < self.min_total_atoms:
            raise ValueError("max_total_atoms must be >= min_total_atoms")
        if not isinstance(self.domain, SearchDomain):
            object.__setattr__(self, "domain", SearchDomain(self.domain))
        limits: dict[str, int] = {}
        for symbol, maximum in self.per_element_max.items():
            if symbol not in elements:
                raise ValueError(f"per-element limit supplied for disallowed element: {symbol}")
            if isinstance(maximum, bool) or not isinstance(maximum, int) or maximum < 1:
                raise ValueError(f"maximum count for {symbol} must be a positive integer")
            limits[symbol] = maximum
        object.__setattr__(self, "per_element_max", MappingProxyType(limits))
        if self.allowed_charges is not None:
            charges = tuple(self.allowed_charges)
            if any(isinstance(q, bool) or not isinstance(q, int) for q in charges):
                raise ValueError("allowed_charges must contain integers only")
            object.__setattr__(self, "allowed_charges", charges)

    def accepts_composition(self, composition: Composition) -> bool:
        """Check only configured search bounds, not chemical plausibility."""
        if not set(composition.counts).issubset(self.allowed_elements):
            return False
        if not self.min_total_atoms <= composition.total_atoms <= self.max_total_atoms:
            return False
        if any(
            count > self.per_element_max[symbol]
            for symbol, count in composition.counts.items()
            if symbol in self.per_element_max
        ):
            return False
        if self.allowed_charges is not None and composition.net_charge not in self.allowed_charges:
            return False
        return True


@dataclass(frozen=True, slots=True)
class PropertyValue:
    """A numeric property with an explicit unit and optional conditions."""
    name: str
    value: float
    unit: str
    conditions: Mapping[str, str | float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("property name cannot be empty")
        if not self.unit.strip():
            raise ValueError("property unit cannot be empty")
        if not isfinite(self.value):
            raise ValueError("property value must be finite")
        object.__setattr__(self, "conditions", MappingProxyType(dict(self.conditions)))


@dataclass(frozen=True, slots=True)
class EvidenceRecord:
    """Provenance record; category does not itself encode confidence."""
    category: EvidenceCategory
    claim: str
    source: str

    def __post_init__(self) -> None:
        if not self.claim.strip():
            raise ValueError("evidence claim cannot be empty")
        if not self.source.strip():
            raise ValueError("evidence source cannot be empty")
        if not isinstance(self.category, EvidenceCategory):
            object.__setattr__(self, "category", EvidenceCategory(self.category))


@dataclass(frozen=True, slots=True)
class Candidate:
    """A composition-level discovery candidate with optional evidence."""
    candidate_id: str
    composition: Composition
    evidence: tuple[EvidenceRecord, ...] = ()

    def __post_init__(self) -> None:
        if not self.candidate_id.strip():
            raise ValueError("candidate_id cannot be empty")
        object.__setattr__(self, "evidence", tuple(self.evidence))


@dataclass(frozen=True, slots=True)
class FragmentPort:
    """A declared attachment site on a molecular fragment."""
    port_id: str
    atom_id: int
    allowed_bond_orders: tuple[int, ...] = (1,)
    label: str = ""

    def __post_init__(self) -> None:
        if not self.port_id.strip():
            raise ValueError("port_id cannot be empty")
        if isinstance(self.atom_id, bool) or not isinstance(self.atom_id, int) or self.atom_id < 0:
            raise ValueError("atom_id must be a non-negative integer")
        orders = tuple(self.allowed_bond_orders)
        if not orders or any(
            isinstance(order, bool) or not isinstance(order, int) or order not in (1, 2, 3)
            for order in orders
        ):
            raise ValueError("allowed_bond_orders must contain 1, 2, or 3")
        if len(orders) != len(set(orders)):
            raise ValueError("allowed_bond_orders cannot contain duplicates")
        object.__setattr__(self, "allowed_bond_orders", orders)


@dataclass(frozen=True, slots=True)
class MolecularFragment:
    """Reusable molecular subgraph with explicit, provenance-aware ports."""
    fragment_id: str
    name: str
    version: str
    graph: MolecularGraph
    ports: tuple[FragmentPort, ...]
    source: str
    description: str = ""

    def __post_init__(self) -> None:
        for field_name in ("fragment_id", "name", "version", "source"):
            if not getattr(self, field_name).strip():
                raise ValueError(f"{field_name} cannot be empty")
        ports = tuple(self.ports)
        object.__setattr__(self, "ports", ports)
        port_ids = [port.port_id for port in ports]
        if len(port_ids) != len(set(port_ids)):
            raise ValueError("fragment port IDs must be unique")
        atom_ids = {atom.atom_id for atom in self.graph.atoms}
        if any(port.atom_id not in atom_ids for port in ports):
            raise ValueError("fragment port references an unknown atom ID")
