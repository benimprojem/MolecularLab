"""Backend-neutral settings for periodic DFT calculation requests.

These records describe requested methods and inputs. They do not run DFT,
select scientifically adequate parameters, or assert backend availability.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from math import isfinite
from types import MappingProxyType
from typing import Mapping

from molecular_lab.structures.crystal import CrystalStructure


class EnergyUnit(str, Enum):
    EV = "eV"
    RYDBERG = "Ry"
    HARTREE = "Ha"


class KPointScheme(str, Enum):
    MONKHORST_PACK = "monkhorst_pack"
    GAMMA_CENTERED = "gamma_centered"


class SmearingScheme(str, Enum):
    FERMI_DIRAC = "fermi_dirac"
    GAUSSIAN = "gaussian"
    METHFESSEL_PAXTON = "methfessel_paxton"
    COLD = "cold"


class DFTTask(str, Enum):
    SCF = "scf"
    RELAX = "relax"
    VC_RELAX = "vc_relax"
    BANDS = "bands"
    DOS = "dos"


class DFTOutput(str, Enum):
    TOTAL_ENERGY = "total_energy"
    FORCES = "forces"
    STRESS = "stress"
    BAND_STRUCTURE = "band_structure"
    DENSITY_OF_STATES = "density_of_states"
    CHARGE_DENSITY = "charge_density"
    MAGNETIC_MOMENTS = "magnetic_moments"


@dataclass(frozen=True, slots=True)
class EnergyCutoff:
    value: float
    unit: EnergyUnit = EnergyUnit.EV

    def __post_init__(self) -> None:
        if not isfinite(self.value) or self.value <= 0:
            raise ValueError("energy cutoff must be finite and positive")
        if not isinstance(self.unit, EnergyUnit):
            object.__setattr__(self, "unit", EnergyUnit(self.unit))

    @property
    def value_ev(self) -> float:
        return self.value * {
            EnergyUnit.EV: 1.0,
            EnergyUnit.RYDBERG: 13.605693122994,
            EnergyUnit.HARTREE: 27.211386245988,
        }[self.unit]

    def to(self, unit: EnergyUnit) -> "EnergyCutoff":
        if not isinstance(unit, EnergyUnit):
            unit = EnergyUnit(unit)
        factor = {
            EnergyUnit.EV: 1.0,
            EnergyUnit.RYDBERG: 13.605693122994,
            EnergyUnit.HARTREE: 27.211386245988,
        }
        return EnergyCutoff(self.value_ev / factor[unit], unit)


@dataclass(frozen=True, slots=True)
class KPointMesh:
    grid: tuple[int, int, int]
    scheme: KPointScheme = KPointScheme.GAMMA_CENTERED
    shift: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self) -> None:
        grid = tuple(self.grid)
        shift = tuple(float(value) for value in self.shift)
        if len(grid) != 3 or any(isinstance(n, bool) or not isinstance(n, int) or n < 1 for n in grid):
            raise ValueError("k-point grid must contain three positive integers")
        if len(shift) != 3 or any(not isfinite(value) for value in shift):
            raise ValueError("k-point shift must contain three finite values")
        if not isinstance(self.scheme, KPointScheme):
            object.__setattr__(self, "scheme", KPointScheme(self.scheme))
        object.__setattr__(self, "grid", grid)
        object.__setattr__(self, "shift", shift)

    def validate_for_structure(self, structure: CrystalStructure) -> None:
        for axis, periodic in enumerate(structure.periodic_boundary_conditions):
            if not periodic and self.grid[axis] != 1:
                raise ValueError(
                    f"non-periodic axis {axis} must use a k-point grid count of 1"
                )


@dataclass(frozen=True, slots=True)
class BandPathPoint:
    """A reciprocal-space path vertex and number of samples to the next vertex."""
    coordinates: tuple[float, float, float]
    label: str
    points_to_next: int = 20

    def __post_init__(self) -> None:
        coordinates = tuple(float(value) for value in self.coordinates)
        if len(coordinates) != 3 or any(not isfinite(value) for value in coordinates):
            raise ValueError("band-path coordinates must contain three finite values")
        if not self.label.strip():
            raise ValueError("band-path point label cannot be empty")
        if isinstance(self.points_to_next, bool) or not isinstance(self.points_to_next, int) or self.points_to_next < 0:
            raise ValueError("points_to_next must be a non-negative integer")
        object.__setattr__(self, "coordinates", coordinates)


@dataclass(frozen=True, slots=True)
class BandPath:
    """Ordered high-symmetry vertices in fractional reciprocal coordinates."""
    points: tuple[BandPathPoint, ...]
    coordinate_system: str = "crystal"

    def __post_init__(self) -> None:
        points = tuple(self.points)
        if len(points) < 2 or any(not isinstance(point, BandPathPoint) for point in points):
            raise ValueError("band path requires at least two BandPathPoint values")
        if self.coordinate_system != "crystal":
            raise ValueError("only fractional crystal reciprocal coordinates are supported")
        if any(point.points_to_next < 1 for point in points[:-1]):
            raise ValueError("every non-final band-path point needs a positive segment subdivision count")
        if points[-1].points_to_next != 0:
            raise ValueError("the final band-path point must use points_to_next=0")
        object.__setattr__(self, "points", points)

    @property
    def sample_count(self) -> int:
        return 1 + sum(point.points_to_next for point in self.points[:-1])


@dataclass(frozen=True, slots=True)
class DOSSettings:
    """Explicit energy window and sampling for total DOS post-processing."""
    energy_min_ev: float
    energy_max_ev: float
    points: int = 2001
    broadening_ev: float = 0.05

    def __post_init__(self) -> None:
        if not isfinite(self.energy_min_ev) or not isfinite(self.energy_max_ev):
            raise ValueError("DOS energy bounds must be finite")
        if self.energy_max_ev <= self.energy_min_ev:
            raise ValueError("DOS energy_max_ev must exceed energy_min_ev")
        if isinstance(self.points, bool) or not isinstance(self.points, int) or self.points < 2:
            raise ValueError("DOS points must be an integer of at least 2")
        if not isfinite(self.broadening_ev) or self.broadening_ev <= 0:
            raise ValueError("DOS broadening must be finite and positive")


@dataclass(frozen=True, slots=True)
class SmearingSettings:
    scheme: SmearingScheme
    width: EnergyCutoff

    def __post_init__(self) -> None:
        if not isinstance(self.scheme, SmearingScheme):
            object.__setattr__(self, "scheme", SmearingScheme(self.scheme))


@dataclass(frozen=True, slots=True)
class ConvergenceSettings:
    energy_tolerance_ev: float = 1e-6
    force_tolerance_ev_per_angstrom: float | None = None
    max_scf_iterations: int = 100
    density_tolerance: float | None = None
    stress_tolerance_gpa: float | None = None

    def __post_init__(self) -> None:
        values = (self.energy_tolerance_ev, self.force_tolerance_ev_per_angstrom,
                  self.stress_tolerance_gpa, self.density_tolerance)
        if any(value is not None and (not isfinite(value) or value <= 0) for value in values):
            raise ValueError("convergence tolerances must be finite and positive")
        if isinstance(self.max_scf_iterations, bool) or self.max_scf_iterations < 1:
            raise ValueError("max_scf_iterations must be a positive integer")
        if not isinstance(self.max_scf_iterations, int):
            raise TypeError("max_scf_iterations must be an integer")


@dataclass(frozen=True, slots=True)
class PseudopotentialReference:
    element: str
    file_name: str
    family: str
    source: str
    valence_electrons: float | None = None
    checksum: str = ""
    license: str = ""

    def __post_init__(self) -> None:
        for name in ("element", "file_name", "family", "source"):
            if not getattr(self, name).strip():
                raise ValueError(f"{name} cannot be empty")
        if self.valence_electrons is not None and (
            not isfinite(self.valence_electrons) or self.valence_electrons <= 0
        ):
            raise ValueError("valence_electrons must be finite and positive")


@dataclass(frozen=True, slots=True)
class DFTCalculationSettings:
    task: DFTTask
    exchange_correlation: str
    energy_cutoff: EnergyCutoff
    kpoint_mesh: KPointMesh
    pseudopotentials: Mapping[str, PseudopotentialReference]
    spin_polarized: bool = False
    spin_orbit_coupling: bool = False
    smearing: SmearingSettings | None = None
    convergence: ConvergenceSettings = field(default_factory=ConvergenceSettings)
    requested_outputs: tuple[DFTOutput, ...] = (DFTOutput.TOTAL_ENERGY,)
    temperature_kelvin: float | None = None
    pressure_gpa: float | None = None
    band_path_ref: str = ""
    band_path: BandPath | None = None
    dos: DOSSettings | None = None
    notes: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.task, DFTTask):
            object.__setattr__(self, "task", DFTTask(self.task))
        if not self.exchange_correlation.strip():
            raise ValueError("exchange_correlation cannot be empty")
        if not isinstance(self.energy_cutoff, EnergyCutoff):
            raise TypeError("energy_cutoff must be an EnergyCutoff")
        if not isinstance(self.kpoint_mesh, KPointMesh):
            raise TypeError("kpoint_mesh must be a KPointMesh")
        if not isinstance(self.convergence, ConvergenceSettings):
            raise TypeError("convergence must be ConvergenceSettings")
        pseudo = dict(self.pseudopotentials)
        if any(not isinstance(value, PseudopotentialReference) for value in pseudo.values()):
            raise TypeError("pseudopotentials must contain PseudopotentialReference values")
        if any(key != value.element for key, value in pseudo.items()):
            raise ValueError("pseudopotential mapping keys must match their element")
        outputs = tuple(
            item if isinstance(item, DFTOutput) else DFTOutput(item)
            for item in self.requested_outputs
        )
        if not outputs or len(outputs) != len(set(outputs)):
            raise ValueError("requested_outputs must be non-empty and unique")
        if self.temperature_kelvin is not None and (
            not isfinite(self.temperature_kelvin) or self.temperature_kelvin < 0
        ):
            raise ValueError("temperature_kelvin must be finite and non-negative")
        if self.pressure_gpa is not None and not isfinite(self.pressure_gpa):
            raise ValueError("pressure_gpa must be finite")
        if self.task is DFTTask.BANDS:
            if self.band_path is None and not self.band_path_ref.strip():
                raise ValueError("BANDS task requires an explicit BandPath or band_path_ref")
            if self.band_path is not None and not isinstance(self.band_path, BandPath):
                raise TypeError("band_path must be a BandPath")
        elif self.band_path_ref or self.band_path is not None:
            raise ValueError("band_path is only valid for a BANDS task")
        if self.task is DFTTask.DOS:
            if self.dos is None or not isinstance(self.dos, DOSSettings):
                raise ValueError("DOS task requires explicit DOSSettings")
        elif self.dos is not None:
            raise ValueError("dos settings are only valid for a DOS task")
        object.__setattr__(self, "pseudopotentials", MappingProxyType(pseudo))
        object.__setattr__(self, "requested_outputs", outputs)

    def validate_for_structure(self, structure: CrystalStructure) -> None:
        self.kpoint_mesh.validate_for_structure(structure)
        required = set(structure.composition_by_site)
        missing = required.difference(self.pseudopotentials)
        if missing:
            raise ValueError(
                "missing pseudopotential references for: " + ", ".join(sorted(missing))
            )
        extra = set(self.pseudopotentials).difference(required)
        if extra:
            raise ValueError(
                "pseudopotential references contain elements absent from structure: "
                + ", ".join(sorted(extra))
            )
