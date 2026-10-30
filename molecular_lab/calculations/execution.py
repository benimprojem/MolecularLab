"""Explicit execution and result normalization for periodic DFT backends.

Input preparation is separate from execution. The runner supports SCF, fixed-cell
ionic RELAX and variable-cell VC_RELAX for Quantum ESPRESSO and GPAW, plus
sequential QE BANDS and total-DOS workflows.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from pathlib import Path
import re
import shutil
import subprocess
from typing import Sequence

from molecular_lab.calculations.ase_bridge import ASEUnavailableError
from molecular_lab.calculations.dft_settings import (
    DFTCalculationSettings, DFTOutput, DFTTask,
)
from molecular_lab.calculations.quantum_espresso import (
    prepare_pw_input, prepare_qe_bands_postprocess_input,
    prepare_qe_dos_postprocess_input, prepare_qe_projected_dos_postprocess_input,
    prepare_qe_fatband_projection_input,
)
from molecular_lab.calculations.qe_output import (
    QEOutputParseError, parse_qe_bands_gnu, parse_qe_fermi_level,
    parse_qe_total_dos, parse_qe_projected_dos,
)
from molecular_lab.materials.electronic_structure import (
    BandStructureData, DensityOfStatesData, ProjectedDensityOfStatesData,
)
from molecular_lab.calculations.qe_projection import (
    parse_qe_atomic_proj_xml, parse_qe_projwfc_band_output,
)
from molecular_lab.calculations.fatband import FatBandDataset, build_fatband_dataset
from molecular_lab.structures.crystal import CrystalSite, CrystalStructure


class DFTExecutionStatus(str, Enum):
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    BACKEND_UNAVAILABLE = "backend_unavailable"
    TIMED_OUT = "timed_out"
    NOT_CONVERGED = "not_converged"


@dataclass(frozen=True, slots=True)
class DFTExecutionResult:
    backend: str
    status: DFTExecutionStatus
    message: str
    total_energy_ev: float | None = None
    converged: bool = False
    return_code: int | None = None
    backend_version: str = ""
    stdout: str = ""
    stderr: str = ""
    optimized_structure: CrystalStructure | None = None
    forces_ev_per_angstrom: tuple[tuple[float, float, float], ...] = ()
    stress_gpa: tuple[tuple[float, float, float], ...] | None = None
    band_structure: BandStructureData | None = None
    density_of_states: DensityOfStatesData | None = None
    projected_density_of_states: ProjectedDensityOfStatesData | None = None
    fat_band_structure: FatBandDataset | None = None

    @property
    def succeeded(self) -> bool:
        return self.status is DFTExecutionStatus.SUCCEEDED


_QE_ENERGY_RE = re.compile(
    r"^\s*!\s+total energy\s*=\s*([-+0-9.EeDd]+)\s+Ry\b",
    re.MULTILINE,
)
_QE_VERSION_RE = re.compile(r"Program PWSCF\s+v\.?([0-9][^\s]*)", re.IGNORECASE)
_QE_FORCE_RE = re.compile(
    r"atom\s+\d+\s+type\s+\d+\s+force\s*=\s*"
    r"([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)",
    re.IGNORECASE,
)
_RY_TO_EV = 13.605693122994
_BOHR_ANGSTROM = 0.529177210903
_EV_PER_ANG3_TO_GPA = 160.21766208
_RY_PER_BOHR_TO_EV_PER_ANGSTROM = _RY_TO_EV / _BOHR_ANGSTROM


def _float_token(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _validate_request(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
) -> None:
    if not isinstance(structure, CrystalStructure):
        raise TypeError("structure must be a CrystalStructure")
    if not isinstance(settings, DFTCalculationSettings):
        raise TypeError("settings must be DFTCalculationSettings")
    if settings.task not in (DFTTask.SCF, DFTTask.RELAX, DFTTask.VC_RELAX):
        raise ValueError("execution currently supports SCF, RELAX and VC_RELAX only")
    allowed = {
        DFTTask.SCF: {DFTOutput.TOTAL_ENERGY, DFTOutput.STRESS},
        DFTTask.RELAX: {DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES, DFTOutput.STRESS},
        DFTTask.VC_RELAX: {DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES, DFTOutput.STRESS},
    }.get(settings.task, set())
    unsupported = set(settings.requested_outputs).difference(allowed)
    if unsupported:
        if settings.task is DFTTask.SCF:
            raise ValueError("SCF execution currently supports total energy only, with optional stress")
        raise ValueError(
            f"{settings.task.value} execution does not support requested outputs: "
            + ", ".join(sorted(item.value for item in unsupported))
        )
    if settings.task in (DFTTask.RELAX, DFTTask.VC_RELAX):
        if not all(structure.periodic_boundary_conditions):
            raise ValueError(f"{settings.task.value} currently requires 3D periodic boundary conditions")
        if settings.convergence.force_tolerance_ev_per_angstrom is None:
            raise ValueError(f"{settings.task.value} requires an explicit force_tolerance_ev_per_angstrom")
    if settings.task is DFTTask.VC_RELAX:
        required = {DFTOutput.TOTAL_ENERGY, DFTOutput.FORCES, DFTOutput.STRESS}
        if not required.issubset(settings.requested_outputs):
            raise ValueError("VC_RELAX requires TOTAL_ENERGY, FORCES and STRESS outputs")
        if settings.convergence.stress_tolerance_gpa is None:
            raise ValueError("VC_RELAX requires an explicit stress_tolerance_gpa")


def _updated_structure(
    original: CrystalStructure,
    lattice: Sequence[Sequence[float]],
    scaled_positions: Sequence[Sequence[float]],
) -> CrystalStructure:
    if len(scaled_positions) != len(original.sites):
        raise ValueError("backend returned a different atom count from the input structure")
    sites = tuple(
        CrystalSite(
            site_id=site.site_id,
            element=site.element,
            fractional_coordinates=tuple(float(v) for v in position),
            occupancy=site.occupancy,
            oxidation_state=site.oxidation_state,
            magnetic_moment_mu_b=site.magnetic_moment_mu_b,
            label=site.label,
        )
        for site, position in zip(original.sites, scaled_positions)
    )
    return CrystalStructure(
        structure_id=original.structure_id,
        lattice_vectors_angstrom=tuple(
            tuple(float(v) for v in row) for row in lattice
        ),
        sites=sites,
        periodic_boundary_conditions=original.periodic_boundary_conditions,
        defects=original.defects,
        space_group=original.space_group,
        magnetic_ordering=original.magnetic_ordering,
        total_charge=original.total_charge,
        source=original.source,
        notes=original.notes,
    )


def _parse_qe_final_structure(
    output: str,
    original: CrystalStructure,
) -> CrystalStructure | None:
    match = re.search(
        r"Begin final coordinates(.*?)End final coordinates",
        output, re.IGNORECASE | re.DOTALL,
    )
    if match is None:
        return None
    block = match.group(1)
    lines = [line.strip() for line in block.splitlines() if line.strip()]
    lattice = original.lattice_vectors_angstrom
    positions: list[tuple[float, float, float]] | None = None
    for index, line in enumerate(lines):
        upper = line.upper()
        if upper.startswith("CELL_PARAMETERS"):
            if index + 3 >= len(lines):
                return None
            try:
                lattice = tuple(
                    tuple(_float_token(value) for value in lines[index + offset].split()[:3])
                    for offset in (1, 2, 3)
                )
            except (ValueError, IndexError):
                return None
        if upper.startswith("ATOMIC_POSITIONS"):
            parsed = []
            for atom_line in lines[index + 1:index + 1 + len(original.sites)]:
                fields = atom_line.split()
                if len(fields) < 4:
                    return None
                try:
                    parsed.append(tuple(_float_token(value) for value in fields[1:4]))
                except ValueError:
                    return None
            if len(parsed) != len(original.sites):
                return None
            positions = parsed
            break
    if positions is None:
        return None
    return _updated_structure(original, lattice, positions)


def _parse_qe_final_forces(
    output: str,
    atom_count: int,
) -> tuple[tuple[float, float, float], ...]:
    matches = list(_QE_FORCE_RE.finditer(output))
    if len(matches) < atom_count:
        return ()
    final = matches[-atom_count:]
    return tuple(
        tuple(_float_token(value) * _RY_PER_BOHR_TO_EV_PER_ANGSTROM for value in item.groups())
        for item in final
    )

def _parse_qe_stress(output: str) -> tuple[tuple[float, float, float], ...] | None:
    """Return the final QE stress tensor in GPa, tension-positive (ASE convention).

    QE prints the stress components in kbar with positive hydrostatic pressure
    corresponding to compression; negate and divide by ten for GPa.
    """
    lines = output.splitlines()
    blocks = []
    for index, line in enumerate(lines):
        lowered = line.lower()
        if "total" not in lowered or "stress" not in lowered or "(kbar)" not in lowered:
            continue
        rows = []
        for candidate in lines[index + 1:index + 4]:
            values = candidate.split()
            if len(values) < 6:
                rows = []
                break
            try:
                numeric = [_float_token(value) for value in values[:6]]
            except ValueError:
                rows = []
                break
            rows.append(numeric)
        if len(rows) == 3:
            blocks.append(tuple(
                tuple(-row[column] / 10.0 for column in (3, 4, 5))
                for row in rows
            ))
    return blocks[-1] if blocks else None


def execute_quantum_espresso(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    *,
    working_directory: str | Path,
    executable: str = "pw.x",
    timeout_seconds: float = 3600.0,
    command_prefix: Sequence[str] = (),
) -> DFTExecutionResult:
    """Run QE pw.x for SCF, fixed-cell RELAX or variable-cell VC_RELAX."""
    _validate_request(structure, settings)
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")
    workdir = Path(working_directory).expanduser().resolve()
    if not workdir.is_dir():
        raise ValueError("working_directory must be an existing directory")
    executable_path = shutil.which(executable)
    if executable_path is None and Path(executable).is_file():
        executable_path = str(Path(executable).resolve())
    if executable_path is None:
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.BACKEND_UNAVAILABLE,
            f"QE executable was not found: {executable}",
        )

    pseudo_dir = workdir / "pseudo"
    missing = [
        settings.pseudopotentials[element].file_name
        for element in sorted(structure.composition_by_site)
        if not (pseudo_dir / settings.pseudopotentials[element].file_name).is_file()
    ]
    if missing:
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.BACKEND_UNAVAILABLE,
            "Required pseudopotential files are missing from ./pseudo: "
            + ", ".join(missing),
        )

    try:
        prepared = prepare_pw_input(
            structure, settings, outdir="./tmp", pseudo_dir="./pseudo",
        )
    except ASEUnavailableError as exc:
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.BACKEND_UNAVAILABLE,
            f"ASE is required to prepare the QE input: {exc}",
        )
    (workdir / "tmp").mkdir(exist_ok=True)
    try:
        completed = subprocess.run(
            [*command_prefix, executable_path],
            input=prepared.input_text,
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.TIMED_OUT,
            f"QE exceeded the {timeout_seconds:g} second time limit.",
            stdout=exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or ""),
            stderr=exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or ""),
        )
    except OSError as exc:
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.FAILED,
            f"Could not start QE: {type(exc).__name__}: {exc}",
        )

    output = completed.stdout or ""
    version_match = _QE_VERSION_RE.search(output)
    energies = list(_QE_ENERGY_RE.finditer(output))
    converged = (
        "convergence has been achieved" in output.lower()
        if settings.task is DFTTask.SCF
        else bool(re.search(
            r"(bfgs converged in|end of bfgs geometry optimization)",
            output, re.IGNORECASE,
        ))
    )
    job_done = "JOB DONE." in output
    structural_task = settings.task in (DFTTask.RELAX, DFTTask.VC_RELAX)
    optimized = _parse_qe_final_structure(output, structure) if structural_task else None
    forces = _parse_qe_final_forces(output, len(structure.sites)) if structural_task else ()
    stress = (
        _parse_qe_stress(output)
        if DFTOutput.STRESS in settings.requested_outputs or settings.task is DFTTask.VC_RELAX
        else None
    )

    if completed.returncode != 0:
        status = DFTExecutionStatus.FAILED
        message = f"QE exited with return code {completed.returncode}."
    elif not converged:
        status = DFTExecutionStatus.NOT_CONVERGED
        message = "QE finished without reporting the requested convergence."
    elif not job_done or not energies:
        status = DFTExecutionStatus.FAILED
        message = "QE output is incomplete: normal completion or total energy was not found."
    elif structural_task and optimized is None:
        status = DFTExecutionStatus.FAILED
        message = "QE reported convergence but its final atomic coordinates/cell could not be parsed."
    elif structural_task and len(forces) != len(structure.sites):
        status = DFTExecutionStatus.FAILED
        message = "QE final atomic forces could not be parsed for every site."
    elif (DFTOutput.STRESS in settings.requested_outputs or settings.task is DFTTask.VC_RELAX) and stress is None:
        status = DFTExecutionStatus.FAILED
        message = "QE final stress tensor could not be parsed."
    else:
        status = DFTExecutionStatus.SUCCEEDED
        message = (
            "QE SCF converged and reported a total energy."
            if settings.task is DFTTask.SCF
            else (
                "QE variable-cell relaxation converged; final cell, atomic positions, forces and stress were parsed."
                if settings.task is DFTTask.VC_RELAX
                else "QE fixed-cell ionic relaxation converged; final structure and forces were parsed."
            )
        )

    energy_ev = (
        _float_token(energies[-1].group(1)) * _RY_TO_EV
        if status is DFTExecutionStatus.SUCCEEDED else None
    )
    return DFTExecutionResult(
        "quantum_espresso", status, message,
        total_energy_ev=energy_ev,
        converged=converged,
        return_code=completed.returncode,
        backend_version=version_match.group(1) if version_match else "",
        stdout=output,
        stderr=completed.stderr or "",
        optimized_structure=optimized if status is DFTExecutionStatus.SUCCEEDED else None,
        forces_ev_per_angstrom=forces if structural_task and status is DFTExecutionStatus.SUCCEEDED else (),
        stress_gpa=stress if status is DFTExecutionStatus.SUCCEEDED else None,
    )


def execute_quantum_espresso_scf(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    **kwargs,
) -> DFTExecutionResult:
    """Backward-compatible SCF-only entry point."""
    if settings.task is not DFTTask.SCF:
        raise ValueError("execute_quantum_espresso_scf accepts SCF tasks only")
    return execute_quantum_espresso(structure, settings, **kwargs)


def execute_gpaw(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    *,
    log_file: str | Path | None = None,
) -> DFTExecutionResult:
    """Run GPAW SCF, fixed-cell RELAX or variable-cell VC_RELAX through ASE."""
    _validate_request(structure, settings)
    try:
        from molecular_lab.calculations.gpaw import (
            build_optimization_target, create_gpaw_calculator, prepare_gpaw_input,
        )
        from gpaw import __version__ as gpaw_version
        from ase.optimize import BFGS
    except ImportError as exc:
        return DFTExecutionResult(
            "gpaw", DFTExecutionStatus.BACKEND_UNAVAILABLE,
            f"GPAW, ASE or its dependencies are unavailable: {exc}",
        )

    prepared = prepare_gpaw_input(structure, settings)
    calculator = create_gpaw_calculator(
        prepared, txt=str(log_file) if log_file is not None else None,
    )
    prepared.atoms.calc = calculator
    try:
        if settings.task in (DFTTask.RELAX, DFTTask.VC_RELAX):
            target = (
                build_optimization_target(prepared)
                if settings.task is DFTTask.VC_RELAX
                else prepared.atoms
            )
            optimizer = BFGS(
                target,
                logfile=(str(log_file) + ".optimizer") if log_file is not None else "-",
            )
            optimizer_converged = optimizer.run(
                fmax=settings.convergence.force_tolerance_ev_per_angstrom,
            )
            forces_array = prepared.atoms.get_forces()
            forces = tuple(tuple(float(v) for v in row) for row in forces_array)
            max_force = max(
                (sum(component * component for component in row) ** 0.5 for row in forces),
                default=0.0,
            )
            force_ok = bool(optimizer_converged) or (
                max_force <= settings.convergence.force_tolerance_ev_per_angstrom
            )
            stress = None
            stress_ok = True
            if DFTOutput.STRESS in settings.requested_outputs or settings.task is DFTTask.VC_RELAX:
                stress_array = prepared.atoms.get_stress(voigt=False)
                stress = tuple(tuple(float(v) * _EV_PER_ANG3_TO_GPA for v in row) for row in stress_array)
            if settings.task is DFTTask.VC_RELAX:
                target_pressure = settings.pressure_gpa or 0.0
                tolerance = settings.convergence.stress_tolerance_gpa
                residual = [
                    abs(stress[i][j] + (target_pressure if i == j else 0.0))
                    for i in range(3) for j in range(3)
                ]
                stress_ok = max(residual, default=0.0) <= tolerance
            optimized = _updated_structure(
                structure,
                prepared.atoms.cell.array,
                prepared.atoms.get_scaled_positions(wrap=False),
            )
            energy = float(prepared.atoms.get_potential_energy())
            if not force_ok or not stress_ok:
                reason = (
                    "force tolerance" if not force_ok
                    else "stress tolerance (post-checked against the requested external pressure)"
                )
                return DFTExecutionResult(
                    "gpaw", DFTExecutionStatus.NOT_CONVERGED,
                    f"GPAW optimization stopped before satisfying the requested {reason}.",
                    total_energy_ev=energy,
                    converged=False,
                    backend_version=str(gpaw_version),
                    optimized_structure=optimized,
                    forces_ev_per_angstrom=forces,
                    stress_gpa=stress,
                )
        else:
            energy = float(prepared.atoms.get_potential_energy())
            forces = ()
            optimized = None
            stress = None
            if DFTOutput.STRESS in settings.requested_outputs:
                stress_array = prepared.atoms.get_stress(voigt=False)
                stress = tuple(tuple(float(v) * _EV_PER_ANG3_TO_GPA for v in row) for row in stress_array)
            scf_state = getattr(getattr(calculator, "scf", None), "converged", None)
            if scf_state is False:
                return DFTExecutionResult(
                    "gpaw", DFTExecutionStatus.NOT_CONVERGED,
                    "GPAW returned an energy but its SCF state is not converged.",
                    backend_version=str(gpaw_version),
                    stress_gpa=stress,
                )
            return DFTExecutionResult(
                "gpaw", DFTExecutionStatus.SUCCEEDED,
                "GPAW returned a total energy."
                + (" SCF convergence was confirmed." if scf_state is True else " SCF convergence flag was not exposed by this GPAW build."),
                total_energy_ev=energy,
                converged=scf_state is True,
                backend_version=str(gpaw_version),
                stress_gpa=stress,
            )
    except Exception as exc:
        return DFTExecutionResult(
            "gpaw", DFTExecutionStatus.FAILED,
            f"GPAW calculation failed: {type(exc).__name__}: {exc}",
            backend_version=str(gpaw_version),
        )

    return DFTExecutionResult(
        "gpaw", DFTExecutionStatus.SUCCEEDED,
        (
            "GPAW variable-cell relaxation converged; final cell, atomic positions, forces and stress were collected."
            if settings.task is DFTTask.VC_RELAX
            else "GPAW fixed-cell ionic relaxation converged; final structure and forces were collected."
        ),
        total_energy_ev=energy,
        converged=True,
        backend_version=str(gpaw_version),
        optimized_structure=optimized,
        forces_ev_per_angstrom=forces,
        stress_gpa=stress,
    )


def execute_gpaw_scf(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    *,
    log_file: str | Path | None = None,
) -> DFTExecutionResult:
    """Backward-compatible SCF-only entry point."""
    if settings.task is not DFTTask.SCF:
        raise ValueError("execute_gpaw_scf accepts SCF tasks only")
    return execute_gpaw(structure, settings, log_file=log_file)



def _resolve_qe_executable(executable: str) -> str | None:
    resolved = shutil.which(executable)
    if resolved is None and Path(executable).is_file():
        resolved = str(Path(executable).resolve())
    return resolved


def _run_qe_input(
    executable: str,
    input_text: str,
    *,
    workdir: Path,
    timeout_seconds: float,
    command_prefix: Sequence[str],
) -> tuple[subprocess.CompletedProcess | None, DFTExecutionResult | None]:
    executable_path = _resolve_qe_executable(executable)
    if executable_path is None:
        return None, DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.BACKEND_UNAVAILABLE,
            f"QE executable was not found: {executable}",
        )
    try:
        completed = subprocess.run(
            [*command_prefix, executable_path],
            input=input_text, cwd=workdir, capture_output=True,
            text=True, timeout=timeout_seconds, check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return None, DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.TIMED_OUT,
            f"{executable} exceeded the {timeout_seconds:g} second time limit.",
            stdout=exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or ""),
            stderr=exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or ""),
        )
    except OSError as exc:
        return None, DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.FAILED,
            f"Could not start {executable}: {type(exc).__name__}: {exc}",
        )
    if completed.returncode != 0:
        return completed, DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.FAILED,
            f"{executable} exited with return code {completed.returncode}.",
            return_code=completed.returncode,
            stdout=completed.stdout or "", stderr=completed.stderr or "",
        )
    if "JOB DONE." not in (completed.stdout or ""):
        return completed, DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.FAILED,
            f"{executable} did not report normal completion.",
            return_code=completed.returncode,
            stdout=completed.stdout or "", stderr=completed.stderr or "",
        )
    return completed, None


def execute_qe_bands(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    *,
    working_directory: str | Path,
    executable: str = "pw.x",
    bands_executable: str = "bands.x",
    projwfc_executable: str = "projwfc.x",
    timeout_seconds: float = 3600.0,
    command_prefix: Sequence[str] = (),
    output_file: str = "molecular_lab.bands",
) -> DFTExecutionResult:
    """Execute SCF -> pw.x BANDS -> bands.x and optional projwfc.x fat bands."""
    if settings.task is not DFTTask.BANDS or settings.band_path is None:
        raise ValueError("execute_qe_bands requires BANDS settings with a structured BandPath")
    if settings.spin_polarized:
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.BACKEND_UNAVAILABLE,
            "Spin-polarized BANDS execution is not yet supported by the common QE workflow.",
        )
    wants_fatband = DFTOutput.FAT_BAND_STRUCTURE in settings.requested_outputs
    workdir = Path(working_directory).expanduser().resolve()
    if not workdir.is_dir():
        raise ValueError("working_directory must be an existing directory")
    scf_settings = replace(
        settings, task=DFTTask.SCF, band_path=None, band_path_ref="",
        dos=None, requested_outputs=(DFTOutput.TOTAL_ENERGY,),
    )
    scf_result = execute_quantum_espresso(
        structure, scf_settings, working_directory=workdir,
        executable=executable, timeout_seconds=timeout_seconds,
        command_prefix=command_prefix,
    )
    if not scf_result.succeeded:
        return replace(scf_result, message="BANDS workflow stopped at prerequisite SCF: " + scf_result.message)
    try:
        bands_input = prepare_pw_input(structure, settings, outdir="./tmp", pseudo_dir="./pseudo")
        post_input = prepare_qe_bands_postprocess_input(
            settings, output_file=output_file,
        )
        projection_input = (
            prepare_qe_fatband_projection_input(settings)
            if wants_fatband else None
        )
    except Exception as exc:
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.FAILED,
            f"Could not prepare BANDS stages: {type(exc).__name__}: {exc}",
            total_energy_ev=scf_result.total_energy_ev,
            backend_version=scf_result.backend_version,
        )
    (workdir / "tmp").mkdir(exist_ok=True)
    bands_run, error = _run_qe_input(
        executable, bands_input.input_text, workdir=workdir,
        timeout_seconds=timeout_seconds, command_prefix=command_prefix,
    )
    if error is not None:
        return replace(error, message="BANDS workflow stopped at pw.x BANDS stage: " + error.message)
    post_run, error = _run_qe_input(
        bands_executable, post_input.input_text, workdir=workdir,
        timeout_seconds=timeout_seconds, command_prefix=command_prefix,
    )
    if error is not None:
        return replace(error, message="BANDS workflow stopped at bands.x stage: " + error.message)
    gnu_path = workdir / (output_file + ".gnu")
    if not gnu_path.is_file():
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.FAILED,
            f"bands.x completed but expected output was not found: {gnu_path.name}",
            total_energy_ev=scf_result.total_energy_ev,
            backend_version=scf_result.backend_version,
            stdout=(bands_run.stdout or "") + (post_run.stdout or ""),
            stderr=(bands_run.stderr or "") + (post_run.stderr or ""),
        )
    try:
        fermi_level = parse_qe_fermi_level(scf_result.stdout)
        band_data = parse_qe_bands_gnu(
            gnu_path.read_text(encoding="utf-8"),
            fermi_level_ev=fermi_level,
            backend_version=scf_result.backend_version,
            exchange_correlation=settings.exchange_correlation,
        )
    except (OSError, QEOutputParseError) as exc:
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.FAILED,
            f"BANDS workflow output parsing failed: {exc}",
            total_energy_ev=scf_result.total_energy_ev,
            backend_version=scf_result.backend_version,
        )

    fat_band_data = None
    if wants_fatband:
        projection_run, error = _run_qe_input(
            projwfc_executable, projection_input.input_text, workdir=workdir,
            timeout_seconds=timeout_seconds, command_prefix=command_prefix,
        )
        if error is not None:
            return replace(error, message="BANDS workflow stopped at projwfc.x stage: " + error.message)
        projection_xml = workdir / "tmp" / "molecular_lab.save" / "atomic_proj.xml"
        if not projection_xml.is_file():
            return DFTExecutionResult(
                "quantum_espresso", DFTExecutionStatus.FAILED,
                "projwfc.x completed but atomic_proj.xml was not found at tmp/molecular_lab.save/atomic_proj.xml.",
                total_energy_ev=scf_result.total_energy_ev,
                backend_version=scf_result.backend_version,
                stdout=(scf_result.stdout or "") + (bands_run.stdout or "")
                    + (post_run.stdout or "") + (projection_run.stdout or ""),
                stderr=(bands_run.stderr or "") + (post_run.stderr or "")
                    + (projection_run.stderr or ""),
            )
        try:
            projection_data = parse_qe_atomic_proj_xml(
                projection_xml.read_text(encoding="utf-8"),
                projection_run.stdout or "",
            )
            original_state_bands = parse_qe_projwfc_band_output(
                projection_run.stdout or "",
                fermi_level_ev=fermi_level,
                backend_version=scf_result.backend_version,
                exchange_correlation=settings.exchange_correlation,
            )
            fat_band_data = build_fatband_dataset(original_state_bands, projection_data)
        except (OSError, ValueError) as exc:
            return DFTExecutionResult(
                "quantum_espresso", DFTExecutionStatus.FAILED,
                f"Fat-band projection parsing/join failed: {exc}",
                total_energy_ev=scf_result.total_energy_ev,
                backend_version=scf_result.backend_version,
                stdout=(scf_result.stdout or "") + (bands_run.stdout or "")
                    + (post_run.stdout or "") + (projection_run.stdout or ""),
                stderr=(bands_run.stderr or "") + (post_run.stderr or "")
                    + (projection_run.stderr or ""),
            )

    return DFTExecutionResult(
        "quantum_espresso", DFTExecutionStatus.SUCCEEDED,
        "QE SCF, BANDS and bands.x stages completed; band data were parsed."
            + (" Original-state projwfc.x projections were joined into fat-band data." if wants_fatband else ""),
        total_energy_ev=scf_result.total_energy_ev, converged=True,
        return_code=post_run.returncode,
        backend_version=scf_result.backend_version,
        stdout=(scf_result.stdout or "") + (bands_run.stdout or "") + (post_run.stdout or ""),
        stderr=(bands_run.stderr or "") + (post_run.stderr or ""),
        band_structure=band_data,
        fat_band_structure=fat_band_data,
    )


def execute_qe_dos(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    *,
    working_directory: str | Path,
    executable: str = "pw.x",
    dos_executable: str = "dos.x",
    projwfc_executable: str = "projwfc.x",
    timeout_seconds: float = 3600.0,
    command_prefix: Sequence[str] = (),
    output_file: str = "molecular_lab.dos",
    pdos_prefix: str = "molecular_lab.pdos",
) -> DFTExecutionResult:
    """Execute SCF -> NSCF -> optional dos.x/projwfc.x post-processing."""
    if settings.task is not DFTTask.DOS or settings.dos is None:
        raise ValueError("execute_qe_dos requires DOS settings")
    wants_total = DFTOutput.DENSITY_OF_STATES in settings.requested_outputs
    wants_projected = DFTOutput.PROJECTED_DOS in settings.requested_outputs
    if not (wants_total or wants_projected):
        raise ValueError("DOS task must request total DOS, projected DOS, or both")
    workdir = Path(working_directory).expanduser().resolve()
    if not workdir.is_dir():
        raise ValueError("working_directory must be an existing directory")
    scf_settings = replace(
        settings, task=DFTTask.SCF, band_path=None, band_path_ref="",
        dos=None, requested_outputs=(DFTOutput.TOTAL_ENERGY,),
    )
    scf_result = execute_quantum_espresso(
        structure, scf_settings, working_directory=workdir,
        executable=executable, timeout_seconds=timeout_seconds,
        command_prefix=command_prefix,
    )
    if not scf_result.succeeded:
        return replace(scf_result, message="DOS workflow stopped at prerequisite SCF: " + scf_result.message)
    try:
        nscf_input = prepare_pw_input(structure, settings, outdir="./tmp", pseudo_dir="./pseudo")
        dos_input = prepare_qe_dos_postprocess_input(settings, output_file=output_file) if wants_total else None
        pdos_input = prepare_qe_projected_dos_postprocess_input(
            settings, pdos_prefix=pdos_prefix,
        ) if wants_projected else None
    except Exception as exc:
        return DFTExecutionResult(
            "quantum_espresso", DFTExecutionStatus.FAILED,
            f"Could not prepare DOS stages: {type(exc).__name__}: {exc}",
            total_energy_ev=scf_result.total_energy_ev,
            backend_version=scf_result.backend_version,
        )
    (workdir / "tmp").mkdir(exist_ok=True)
    nscf_run, error = _run_qe_input(
        executable, nscf_input.input_text, workdir=workdir,
        timeout_seconds=timeout_seconds, command_prefix=command_prefix,
    )
    if error is not None:
        return replace(error, message="DOS workflow stopped at pw.x NSCF stage: " + error.message)
    stdout_parts = [scf_result.stdout or "", nscf_run.stdout or ""]
    stderr_parts = [nscf_run.stderr or ""]
    last_run = nscf_run
    dos_data = None
    projected_data = None
    if wants_total:
        dos_run, error = _run_qe_input(
            dos_executable, dos_input.input_text, workdir=workdir,
            timeout_seconds=timeout_seconds, command_prefix=command_prefix,
        )
        if error is not None:
            return replace(error, message="DOS workflow stopped at dos.x stage: " + error.message)
        last_run = dos_run
        stdout_parts.append(dos_run.stdout or "")
        stderr_parts.append(dos_run.stderr or "")
        dos_path = workdir / output_file
        if not dos_path.is_file():
            return DFTExecutionResult(
                "quantum_espresso", DFTExecutionStatus.FAILED,
                f"dos.x completed but expected output was not found: {dos_path.name}",
                total_energy_ev=scf_result.total_energy_ev,
                backend_version=scf_result.backend_version,
                stdout="".join(stdout_parts), stderr="".join(stderr_parts),
            )
        try:
            dos_data = parse_qe_total_dos(
                dos_path.read_text(encoding="utf-8"),
                spin_polarized=settings.spin_polarized,
                fermi_level_ev=parse_qe_fermi_level(scf_result.stdout),
                backend_version=scf_result.backend_version,
            )
        except (OSError, QEOutputParseError) as exc:
            return DFTExecutionResult(
                "quantum_espresso", DFTExecutionStatus.FAILED,
                f"Total DOS output parsing failed: {exc}",
                total_energy_ev=scf_result.total_energy_ev,
                backend_version=scf_result.backend_version,
            )
    if wants_projected:
        pdos_run, error = _run_qe_input(
            projwfc_executable, pdos_input.input_text, workdir=workdir,
            timeout_seconds=timeout_seconds, command_prefix=command_prefix,
        )
        if error is not None:
            return replace(error, message="DOS workflow stopped at projwfc.x stage: " + error.message)
        last_run = pdos_run
        stdout_parts.append(pdos_run.stdout or "")
        stderr_parts.append(pdos_run.stderr or "")
        pdos_paths = sorted(workdir.glob(pdos_prefix + ".pdos_atm#*"))
        if not pdos_paths:
            return DFTExecutionResult(
                "quantum_espresso", DFTExecutionStatus.FAILED,
                "projwfc.x completed but no atomic PDOS files were found.",
                total_energy_ev=scf_result.total_energy_ev,
                backend_version=scf_result.backend_version,
                stdout="".join(stdout_parts), stderr="".join(stderr_parts),
            )
        try:
            projected_data = parse_qe_projected_dos(
                {path.name: path.read_text(encoding="utf-8") for path in pdos_paths},
                spin_polarized=settings.spin_polarized,
                fermi_level_ev=parse_qe_fermi_level(scf_result.stdout),
                backend_version=scf_result.backend_version,
            )
        except (OSError, QEOutputParseError) as exc:
            return DFTExecutionResult(
                "quantum_espresso", DFTExecutionStatus.FAILED,
                f"Projected DOS output parsing failed: {exc}",
                total_energy_ev=scf_result.total_energy_ev,
                backend_version=scf_result.backend_version,
            )
    requested = "total DOS" if wants_total else ""
    if wants_projected:
        requested += (" and " if requested else "") + "atom/orbital PDOS"
    return DFTExecutionResult(
        "quantum_espresso", DFTExecutionStatus.SUCCEEDED,
        f"QE SCF, NSCF and requested {requested} post-processing completed.",
        total_energy_ev=scf_result.total_energy_ev, converged=True,
        return_code=last_run.returncode,
        backend_version=scf_result.backend_version,
        stdout="".join(stdout_parts), stderr="".join(stderr_parts),
        density_of_states=dos_data,
        projected_density_of_states=projected_data,
    )


def execute_dft_calculation(
    backend: str,
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    *,
    working_directory: str | Path | None = None,
    executable: str = "pw.x",
    timeout_seconds: float = 3600.0,
    command_prefix: Sequence[str] = (),
    log_file: str | Path | None = None,
    bands_executable: str = "bands.x",
    projwfc_executable: str = "projwfc.x",
) -> DFTExecutionResult:
    """Common user-facing entry point dispatching to a selected DFT backend."""
    backend_key = backend.strip().lower()
    if backend_key in ("qe", "quantum_espresso", "quantum-espresso"):
        if settings.task is DFTTask.BANDS:
            if working_directory is None:
                raise ValueError("working_directory is required for Quantum ESPRESSO")
            return execute_qe_bands(
                structure, settings, working_directory=working_directory,
                executable=executable, bands_executable=bands_executable,
                projwfc_executable=projwfc_executable,
                timeout_seconds=timeout_seconds, command_prefix=command_prefix,
            )
        if settings.task is DFTTask.DOS:
            if working_directory is None:
                raise ValueError("working_directory is required for Quantum ESPRESSO")
            return execute_qe_dos(
                structure, settings, working_directory=working_directory,
                executable=executable, timeout_seconds=timeout_seconds,
                command_prefix=command_prefix,
            )
        if working_directory is None:
            raise ValueError("working_directory is required for Quantum ESPRESSO")
        return execute_quantum_espresso(
            structure, settings, working_directory=working_directory,
            executable=executable, timeout_seconds=timeout_seconds,
            command_prefix=command_prefix,
        )
    if backend_key == "gpaw":
        if settings.task in (DFTTask.BANDS, DFTTask.DOS):
            return DFTExecutionResult(
                "gpaw", DFTExecutionStatus.BACKEND_UNAVAILABLE,
                "GPAW BANDS/DOS workflows are not implemented yet.",
            )
        return execute_gpaw(structure, settings, log_file=log_file)
    return DFTExecutionResult(
        backend_key or "unknown", DFTExecutionStatus.BACKEND_UNAVAILABLE,
        f"Unknown DFT backend: {backend}",
    )
