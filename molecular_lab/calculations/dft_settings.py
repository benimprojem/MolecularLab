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

    def __post_init__(self) -> None:
        values = (self.energy_tolerance_ev, self.force_tolerance_ev_per_angstrom,
                  self.density_tolerance)
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
        if any(key != value.element for key, value in pseudo.items()):
            raise ValueError("pseudopotential mapping keys must match their element")
        if any(not isinstance(value, PseudopotentialReference) for value in pseudo.values()):
            raise TypeError("pseudopotentials must contain PseudopotentialReference values")
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
        if self.spin_orbit_coupling and not self.spin_polarized:
            raise ValueError("spin-orbit coupling requires spin-polarized/non-collinear treatment")
        if self.task == DFTTask.BANDS and not self.band_path_ref.strip():
            raise ValueError("BANDS task requires an explicit band_path_ref")
        if self.task != DFTTask.BANDS and self.band_path_ref:
            raise ValueError("band_path_ref is only valid for a BANDS task")
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
