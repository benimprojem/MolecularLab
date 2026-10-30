"""Explicit execution entry points for periodic DFT backends.

The first execution slice deliberately supports only Quantum ESPRESSO SCF.
Input preparation remains separate from execution, and success requires both
QE's convergence marker and its normal job-completion marker.
"""
from __future__ import annotations

from dataclasses import dataclass
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
from molecular_lab.calculations.quantum_espresso import prepare_pw_input
from molecular_lab.structures.crystal import CrystalStructure


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

    @property
    def succeeded(self) -> bool:
        return self.status is DFTExecutionStatus.SUCCEEDED


_QE_ENERGY_RE = re.compile(
    r"^\s*!\s+total energy\s*=\s*([-+0-9.Ee]+)\s+Ry\b",
    re.MULTILINE,
)
_QE_VERSION_RE = re.compile(r"Program PWSCF\s+v\.?([0-9][^\s]*)", re.IGNORECASE)
_RY_TO_EV = 13.605693122994


def _validate_scf_request(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
) -> None:
    if not isinstance(structure, CrystalStructure):
        raise TypeError("structure must be a CrystalStructure")
    if not isinstance(settings, DFTCalculationSettings):
        raise TypeError("settings must be DFTCalculationSettings")
    if settings.task is not DFTTask.SCF:
        raise ValueError("the initial execution layer supports SCF tasks only")
    if set(settings.requested_outputs) != {DFTOutput.TOTAL_ENERGY}:
        raise ValueError("the initial execution layer currently returns total energy only")


def execute_quantum_espresso_scf(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    *,
    working_directory: str | Path,
    executable: str = "pw.x",
    timeout_seconds: float = 3600.0,
    command_prefix: Sequence[str] = (),
) -> DFTExecutionResult:
    """Run one QE pw.x SCF calculation in an existing, caller-owned directory.

    The caller is responsible for providing the requested pseudopotentials in
    ./pseudo (or the configured pseudo_dir) and for choosing scientifically
    converged settings. No network download or implicit installation occurs.
    """
    _validate_scf_request(structure, settings)
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
        name for name in (
            settings.pseudopotentials[element].file_name
            for element in sorted(structure.composition_by_site)
        )
        if not (pseudo_dir / name).is_file()
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
    energy_match = _QE_ENERGY_RE.search(output)
    converged = "convergence has been achieved" in output.lower()
    job_done = "JOB DONE." in output
    if completed.returncode != 0:
        status = DFTExecutionStatus.FAILED
        message = f"QE exited with return code {completed.returncode}."
    elif not converged:
        status = DFTExecutionStatus.NOT_CONVERGED
        message = "QE finished without reporting SCF convergence."
    elif not job_done or energy_match is None:
        status = DFTExecutionStatus.FAILED
        message = "QE output is incomplete: normal completion or total energy was not found."
    else:
        status = DFTExecutionStatus.SUCCEEDED
        message = "QE SCF converged and reported a total energy."

    energy_ev = (
        float(energy_match.group(1)) * _RY_TO_EV
        if status is DFTExecutionStatus.SUCCEEDED and energy_match is not None
        else None
    )
    return DFTExecutionResult(
        "quantum_espresso", status, message,
        total_energy_ev=energy_ev,
        converged=converged,
        return_code=completed.returncode,
        backend_version=version_match.group(1) if version_match else "",
        stdout=output,
        stderr=completed.stderr or "",
    )


def execute_gpaw_scf(
    structure: CrystalStructure,
    settings: DFTCalculationSettings,
    *,
    log_file: str | Path | None = None,
) -> DFTExecutionResult:
    """Run a GPAW SCF calculation directly through ASE's calculator API."""
    _validate_scf_request(structure, settings)
    try:
        from molecular_lab.calculations.gpaw import (
            create_gpaw_calculator, prepare_gpaw_input,
        )
        from gpaw import __version__ as gpaw_version
    except ImportError as exc:
        return DFTExecutionResult(
            "gpaw", DFTExecutionStatus.BACKEND_UNAVAILABLE,
            f"GPAW or its dependencies are unavailable: {exc}",
        )

    prepared = prepare_gpaw_input(structure, settings)
    calculator = create_gpaw_calculator(
        prepared, txt=str(log_file) if log_file is not None else None,
    )
    prepared.atoms.calc = calculator
    try:
        energy = float(prepared.atoms.get_potential_energy())
    except Exception as exc:
        return DFTExecutionResult(
            "gpaw", DFTExecutionStatus.FAILED,
            f"GPAW calculation failed: {type(exc).__name__}: {exc}",
            backend_version=str(gpaw_version),
        )

    scf_state = getattr(getattr(calculator, "scf", None), "converged", None)
    if scf_state is False:
        return DFTExecutionResult(
            "gpaw", DFTExecutionStatus.NOT_CONVERGED,
            "GPAW returned an energy but its SCF state is not converged.",
            backend_version=str(gpaw_version),
        )
    return DFTExecutionResult(
        "gpaw", DFTExecutionStatus.SUCCEEDED,
        "GPAW returned a total energy."
        + (" SCF convergence was confirmed." if scf_state is True else " SCF convergence flag was not exposed by this GPAW build."),
        total_energy_ev=energy,
        converged=scf_state is True,
        backend_version=str(gpaw_version),
    )
