"""Quantum ESPRESSO pw.x input-deck generation (no solver execution).

The adapter intentionally supports a conservative first subset of periodic
calculations. Unsupported settings are rejected rather than silently dropped.
"""
from __future__ import annotations

from dataclasses import dataclass
import re

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


class QuantumEspressoInputError(ValueError):
    """The requested settings cannot be represented safely in a pw.x input."""


@dataclass(frozen=True, slots=True)
class PreparedQEInput:
    input_text: str
    required_pseudopotentials: tuple[str, ...]
    warnings: tuple[str, ...] = ()


_TOKEN_RE = re.compile(r"^[A-Za-z0-9_.+\-]+$")
_SMEARING_NAMES = {
    SmearingScheme.FERMI_DIRAC: "fermi-dirac",
    SmearingScheme.GAUSSIAN: "gaussian",
    SmearingScheme.METHFESSEL_PAXTON: "methfessel-paxton",
    SmearingScheme.COLD: "marzari-vanderbilt",
}
_RY_EV = 13.605693122994
_BOHR_ANGSTROM = 0.529177210903


def _token(value: str, field: str) -> str:
    if not isinstance(value, str) or not _TOKEN_RE.fullmatch(value):
        raise QuantumEspressoInputError(
            f"{field} must contain only letters, digits, underscore, dot, plus or hyphen"
        )
    return value


def _number(value: float) -> str:
    return f"{value:.12g}"


def _quoted_path(value: str, field: str) -> str:
    if not isinstance(value, str) or not value or any(c in value for c in ("\n", "\r", "'", "!")):
        raise QuantumEspressoInputError(f"{field} must be a non-empty safe path")
    return "'" + value.replace("\\", "/") + "'"


def prepare_pw_input(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    *,
    prefix: str = "molecular_lab",
    outdir: str = "./tmp",
    pseudo_dir: str = "./pseudo",
) -> PreparedQEInput:
    """Build a pw.x input deck for SCF, RELAX, or VC-RELAX.

    ASE is required for this bridge because QE's ATOMIC_SPECIES card needs
    standard atomic masses. No pseudopotential file is opened or verified.
    """
    if not isinstance(structure, CrystalStructure):
        raise TypeError("structure must be a CrystalStructure")
    if not isinstance(settings, DFTCalculationSettings):
        raise TypeError("settings must be DFTCalculationSettings")

    if settings.task not in (DFTTask.SCF, DFTTask.RELAX, DFTTask.VC_RELAX):
        raise QuantumEspressoInputError(
            "pw.x adapter currently supports SCF, RELAX and VC_RELAX only; "
            "BANDS and DOS need separate workflow inputs"
        )
    if settings.spin_orbit_coupling:
        raise QuantumEspressoInputError(
            "SOC input is not enabled until fully relativistic pseudopotential metadata is modeled"
        )
    if any(site.magnetic_moment_mu_b is not None for site in structure.sites):
        raise QuantumEspressoInputError(
            "site-resolved initial magnetic moments cannot yet be mapped safely to QE species"
        )
    if settings.pressure_gpa is not None and settings.task is not DFTTask.VC_RELAX:
        raise QuantumEspressoInputError("pressure_gpa is supported only for VC_RELAX")
    if settings.task is DFTTask.VC_RELAX and not all(structure.periodic_boundary_conditions):
        raise QuantumEspressoInputError(
            "VC_RELAX currently requires three-dimensional periodic boundary conditions"
        )
    supported_outputs = {DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES, DFTOutput.STRESS}
    unsupported = set(settings.requested_outputs).difference(supported_outputs)
    if unsupported:
        names = ", ".join(sorted(item.value for item in unsupported))
        raise QuantumEspressoInputError(f"requested outputs are not implemented by this adapter: {names}")

    settings.validate_for_structure(structure)
    if settings.kpoint_mesh.scheme is not KPointScheme.MONKHORST_PACK:
        raise QuantumEspressoInputError(
            "this adapter currently accepts Monkhorst-Pack meshes only"
        )
    shifts = []
    for value in settings.kpoint_mesh.shift:
        if value not in (0.0, 1.0):
            raise QuantumEspressoInputError(
                "QE automatic k-point shifts must be exactly 0 or 1"
            )
        shifts.append(int(value))

    prefix = _token(prefix, "prefix")
    xc = _token(settings.exchange_correlation, "exchange_correlation")
    for element, pseudo in settings.pseudopotentials.items():
        _token(element, "pseudopotential element")
        _token(pseudo.file_name, f"{element} pseudopotential file_name")

    # The ASE bridge also enforces ordered, fully occupied sites and materialized defects.
    atoms = to_ase_atoms(structure)
    symbols = atoms.get_chemical_symbols()
    masses = atoms.get_masses()
    scaled = atoms.get_scaled_positions(wrap=False)
    cell = atoms.cell.array
    species = sorted(set(symbols))
    species_lines = []
    for element in species:
        index = symbols.index(element)
        species_lines.append(
            f"{element:4s} {_number(float(masses[index])):>14s} "
            f"{settings.pseudopotentials[element].file_name}"
        )

    calculation = {
        DFTTask.SCF: "scf",
        DFTTask.RELAX: "relax",
        DFTTask.VC_RELAX: "vc-relax",
    }[settings.task]
    cutoff_ry = settings.energy_cutoff.to(EnergyUnit.RYDBERG).value
    conv_thr_ry = settings.convergence.energy_tolerance_ev / _RY_EV
    force_thr_ry_bohr = None
    if settings.convergence.force_tolerance_ev_per_angstrom is not None:
        force_thr_ry_bohr = (
            settings.convergence.force_tolerance_ev_per_angstrom
            / _RY_EV
            * _BOHR_ANGSTROM
        )

    control = [
        "&CONTROL",
        f"  calculation = '{calculation}'",
        f"  prefix = '{prefix}'",
        f"  outdir = {_quoted_path(outdir, 'outdir')}",
        f"  pseudo_dir = {_quoted_path(pseudo_dir, 'pseudo_dir')}",
        f"  tprnfor = {'.true.' if (DFTOutput.FORCES in settings.requested_outputs or settings.task is not DFTTask.SCF) else '.false.'}",
        f"  tstress = {'.true.' if (DFTOutput.STRESS in settings.requested_outputs or settings.task is DFTTask.VC_RELAX) else '.false.'}",
    ]
    if force_thr_ry_bohr is not None:
        if settings.task is DFTTask.SCF:
            raise QuantumEspressoInputError(
                "force convergence tolerance is meaningful here only for RELAX or VC_RELAX"
            )
        control.append(f"  forc_conv_thr = {_number(force_thr_ry_bohr)}")
    control.extend(["/", ""])

    system = [
        "&SYSTEM",
        "  ibrav = 0",
        f"  nat = {len(symbols)}",
        f"  ntyp = {len(species)}",
        f"  ecutwfc = {_number(cutoff_ry)}",
        f"  input_dft = '{xc}'",
    ]
    if structure.total_charge is not None:
        system.append(f"  tot_charge = {_number(float(structure.total_charge))}")
    if settings.spin_polarized:
        system.append("  nspin = 2")
    if settings.smearing is not None:
        system.extend([
            "  occupations = 'smearing'",
            f"  smearing = '{_SMEARING_NAMES[settings.smearing.scheme]}'",
            f"  degauss = {_number(settings.smearing.width.to(EnergyUnit.RYDBERG).value)}",
        ])
    system.extend(["/", ""])

    electrons = [
        "&ELECTRONS",
        f"  conv_thr = {_number(conv_thr_ry)}",
        f"  electron_maxstep = {settings.convergence.max_scf_iterations}",
        "/",
        "",
    ]
    sections = control + system + electrons

    if settings.task in (DFTTask.RELAX, DFTTask.VC_RELAX):
        sections.extend(["&IONS", "  ion_dynamics = 'bfgs'", "/", ""])
    if settings.task is DFTTask.VC_RELAX:
        if settings.convergence.stress_tolerance_gpa is None:
            raise QuantumEspressoInputError(
                "VC_RELAX requires an explicit stress_tolerance_gpa"
            )
        sections.extend([
            "&CELL",
            "  cell_dynamics = 'bfgs'",
            f"  press_conv_thr = {_number(settings.convergence.stress_tolerance_gpa * 10.0)}",
        ])
        if settings.pressure_gpa is not None:
            sections.append(f"  press = {_number(settings.pressure_gpa * 10.0)}")
        sections.extend(["/", ""])

    sections.extend(["ATOMIC_SPECIES", *species_lines, "", "CELL_PARAMETERS angstrom"])
    sections.extend("  " + " ".join(_number(float(v)) for v in vector) for vector in cell)
    sections.extend(["", "ATOMIC_POSITIONS crystal"])
    sections.extend(
        f"  {symbol:4s} " + " ".join(_number(float(v)) for v in position)
        for symbol, position in zip(symbols, scaled)
    )
    sections.extend([
        "",
        "K_POINTS automatic",
        "  " + " ".join(str(v) for v in settings.kpoint_mesh.grid)
        + " " + " ".join(str(v) for v in shifts),
        "",
    ])

    warnings = [
        "Pseudopotential files are referenced by name but are not downloaded, opened or validated.",
        "The selected exchange-correlation identifier is passed to QE; availability in the installed QE build is not checked.",
    ]
    if settings.spin_polarized:
        warnings.append(
            "Spin-polarized mode is enabled, but no species-specific starting magnetization is generated."
        )
    if settings.temperature_kelvin is not None:
        warnings.append(
            "temperature_kelvin is retained as metadata only; this deck is not a finite-temperature DFT workflow."
        )
    if settings.task in (DFTTask.RELAX, DFTTask.VC_RELAX) and force_thr_ry_bohr is None:
        warnings.append(
            "No ionic force-convergence tolerance was supplied; QE's built-in default will be used."
        )
    if settings.task in (DFTTask.RELAX, DFTTask.VC_RELAX):
        warnings.append(
            "The backend-neutral energy tolerance maps to QE SCF conv_thr; QE ionic etot_conv_thr remains at its built-in default."
        )

    return PreparedQEInput(
        input_text="\n".join(sections),
        required_pseudopotentials=tuple(
            settings.pseudopotentials[element].file_name for element in species
        ),
        warnings=tuple(warnings),
    )
