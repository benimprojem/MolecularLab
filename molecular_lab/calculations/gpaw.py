"""GPAW calculator preparation for periodic DFT (no calculation execution).

This adapter builds ASE/GPAW objects lazily. It never starts SCF or geometry
optimization; callers must explicitly attach/run the returned calculator.
"""
from __future__ import annotations

from dataclasses import dataclass
import re
from types import MappingProxyType
from typing import Any, Mapping

from molecular_lab.calculations.ase_bridge import to_ase_atoms
from molecular_lab.calculations.dft_settings import (
    DFTCalculationSettings,
    DFTOutput,
    DFTTask,
    EnergyUnit,
    KPointScheme,
    SmearingScheme,
)
from molecular_lab.structures.crystal import CrystalStructure


class GPAWInputError(ValueError):
    """The requested settings cannot be represented by this GPAW adapter."""


class GPAWUnavailableError(RuntimeError):
    """The optional GPAW package is not importable."""


@dataclass(frozen=True, slots=True)
class PreparedGPAWInput:
    atoms: Any
    calculator_parameters: Mapping[str, Any]
    plane_wave_cutoff_ev: float
    task: DFTTask
    optimizer_fmax_ev_per_angstrom: float | None
    cell_filter_parameters: Mapping[str, float]
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "calculator_parameters",
            MappingProxyType(dict(self.calculator_parameters)),
        )
        object.__setattr__(
            self, "cell_filter_parameters",
            MappingProxyType(dict(self.cell_filter_parameters)),
        )


_TOKEN_RE = re.compile(r"^[A-Za-z0-9_.+\-]+$")
_GPAW_GPA_TO_EV_PER_ANG3 = 0.006241509074


def _validate_token(value: str, field: str) -> str:
    if not isinstance(value, str) or not _TOKEN_RE.fullmatch(value):
        raise GPAWInputError(
            f"{field} must contain only letters, digits, underscore, dot, plus or hyphen"
        )
    return value


def prepare_gpaw_input(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
) -> PreparedGPAWInput:
    """Prepare ASE Atoms and GPAW keyword arguments without running GPAW.

    GPAW uses its installed PAW setups rather than the QE-style external
    pseudopotential filenames in DFTCalculationSettings. Those references are
    therefore not consumed by this adapter.
    """
    if not isinstance(structure, CrystalStructure):
        raise TypeError("structure must be a CrystalStructure")
    if not isinstance(settings, DFTCalculationSettings):
        raise TypeError("settings must be DFTCalculationSettings")

    if settings.task not in (DFTTask.SCF, DFTTask.RELAX, DFTTask.VC_RELAX):
        raise GPAWInputError(
            "GPAW adapter currently supports SCF, RELAX and VC_RELAX; "
            "BANDS and DOS require separate post-SCF workflows"
        )
    if settings.spin_orbit_coupling:
        raise GPAWInputError(
            "SOC is not enabled in this initial adapter; a dedicated non-collinear workflow is required"
        )
    if any(site.magnetic_moment_mu_b is not None for site in structure.sites):
        raise GPAWInputError(
            "site-resolved initial magnetic moments are not yet mapped to GPAW magmoms"
        )
    if not all(structure.periodic_boundary_conditions):
        raise GPAWInputError(
            "GPAW plane-wave mode uses periodic boundary conditions; this adapter requires 3D PBC"
        )
    if settings.pressure_gpa is not None and settings.task is not DFTTask.VC_RELAX:
        raise GPAWInputError("pressure_gpa is supported only for VC_RELAX")

    supported_outputs = {
        DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES, DFTOutput.STRESS,
    }
    unsupported = set(settings.requested_outputs).difference(supported_outputs)
    if unsupported:
        names = ", ".join(sorted(item.value for item in unsupported))
        raise GPAWInputError(f"requested outputs are not implemented by this adapter: {names}")

    settings.kpoint_mesh.validate_for_structure(structure)
    if any(value != 0.0 for value in settings.kpoint_mesh.shift):
        raise GPAWInputError(
            "GPAW adapter currently supports only zero explicit k-point shifts"
        )
    xc = _validate_token(settings.exchange_correlation, "exchange_correlation")

    occupations = None
    if settings.smearing is not None:
        if settings.smearing.scheme is not SmearingScheme.FERMI_DIRAC:
            raise GPAWInputError(
                "GPAW adapter currently supports Fermi-Dirac smearing only"
            )
        occupations = {
            "name": "fermi-dirac",
            "width": settings.smearing.width.to(EnergyUnit.EV).value,
        }

    convergence: dict[str, float] = {
        "energy": settings.convergence.energy_tolerance_ev,
    }
    if settings.convergence.density_tolerance is not None:
        convergence["density"] = settings.convergence.density_tolerance
    if settings.convergence.force_tolerance_ev_per_angstrom is not None:
        convergence["forces"] = settings.convergence.force_tolerance_ev_per_angstrom

    parameters: dict[str, Any] = {
        "xc": xc,
        "kpts": {
            "size": tuple(settings.kpoint_mesh.grid),
            "gamma": settings.kpoint_mesh.scheme is KPointScheme.GAMMA_CENTERED,
        },
        "spinpol": settings.spin_polarized,
        "maxiter": settings.convergence.max_scf_iterations,
        "convergence": convergence,
    }
    if structure.total_charge is not None:
        parameters["charge"] = float(structure.total_charge)
    if occupations is not None:
        parameters["occupations"] = occupations

    atoms = to_ase_atoms(structure)
    cell_filter_parameters: dict[str, float] = {}
    if settings.pressure_gpa is not None:
        cell_filter_parameters["scalar_pressure"] = (
            settings.pressure_gpa * _GPAW_GPA_TO_EV_PER_ANG3
        )

    warnings = [
        "GPAW PAW setups are selected from the local GPAW installation; DFTCalculationSettings.pseudopotentials are not used.",
        "This function prepares calculator parameters only; it does not execute SCF or geometry optimization.",
        "The XC identifier is passed to GPAW but its availability in the installed build is not checked.",
    ]
    if settings.convergence.energy_tolerance_ev is not None:
        warnings.append(
            "The energy tolerance is mapped to GPAW convergence['energy'], whose documented criterion is expressed per valence electron; this is not numerically identical to every backend's energy threshold."
        )
    if settings.temperature_kelvin is not None:
        warnings.append(
            "temperature_kelvin is retained as metadata only; this is not a finite-temperature DFT workflow."
        )
    if settings.spin_polarized:
        warnings.append(
            "Spin-polarized mode is enabled without guessed initial magnetic moments."
        )
    if settings.task in (DFTTask.RELAX, DFTTask.VC_RELAX):
        warnings.append(
            "An ASE optimizer must be run explicitly; optimizer_fmax is supplied separately from GPAW SCF convergence."
        )
    if settings.task is DFTTask.VC_RELAX:
        warnings.append(
            "Variable-cell relaxation requires ASE FrechetCellFilter and GPAW stress support; GPAW plane-wave mode is selected for this reason."
        )

    return PreparedGPAWInput(
        atoms=atoms,
        calculator_parameters=parameters,
        plane_wave_cutoff_ev=settings.energy_cutoff.to(EnergyUnit.EV).value,
        task=settings.task,
        optimizer_fmax_ev_per_angstrom=(
            settings.convergence.force_tolerance_ev_per_angstrom
            if settings.task in (DFTTask.RELAX, DFTTask.VC_RELAX)
            else None
        ),
        cell_filter_parameters=cell_filter_parameters,
        warnings=tuple(warnings),
    )


def create_gpaw_calculator(
    prepared: PreparedGPAWInput,
    *,
    txt: str | None = None,
) -> Any:
    """Instantiate a GPAW calculator but do not attach it or run a calculation."""
    if not isinstance(prepared, PreparedGPAWInput):
        raise TypeError("prepared must be PreparedGPAWInput")
    try:
        from gpaw import GPAW, PW
    except ImportError as exc:
        raise GPAWUnavailableError(
            "Optional dependency 'gpaw' is not installed or could not be imported"
        ) from exc

    return GPAW(
        mode=PW(prepared.plane_wave_cutoff_ev),
        txt=txt,
        **dict(prepared.calculator_parameters),
    )


def build_optimization_target(prepared: PreparedGPAWInput) -> Any:
    """Return ASE's optimization target; does not instantiate or run an optimizer."""
    if not isinstance(prepared, PreparedGPAWInput):
        raise TypeError("prepared must be PreparedGPAWInput")
    if prepared.task is DFTTask.SCF:
        raise GPAWInputError("SCF has no geometry-optimization target")
    if prepared.task is DFTTask.RELAX:
        return prepared.atoms
    try:
        from ase.filters import FrechetCellFilter
    except ImportError as exc:
        raise GPAWUnavailableError(
            "ASE FrechetCellFilter is required for variable-cell relaxation"
        ) from exc
    return FrechetCellFilter(
        prepared.atoms,
        **dict(prepared.cell_filter_parameters),
    )
