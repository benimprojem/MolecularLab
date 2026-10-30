# DFT Execution Layer

**Status:** Initial execution prototype  
**Scope:** Explicit SCF execution for Quantum ESPRESSO and GPAW

## Separation of responsibilities

- `dft_settings.py` describes backend-neutral calculation settings.
- Backend preparation modules translate settings into QE input text or GPAW/ASE objects.
- `execution.py` is the explicit execution entry point and normalizes the first result fields.
- No backend is installed, downloaded, or launched implicitly by importing MolecularLab.

## Current execution scope

The first execution slice supports only periodic SCF calculations requesting total energy.

### Quantum ESPRESSO

`execute_quantum_espresso_scf`:
- Requires an existing working directory and QE executable.
- Requires every referenced UPF file to exist under `working_directory/pseudo`.
- Creates the local `tmp` output directory and runs `pw.x` with the prepared input on standard input.
- Enforces a caller-configurable process timeout.
- Reports success only when the process exits with code 0, the SCF convergence marker and `JOB DONE.` are present, and a total energy is parsed.
- Converts the reported Rydberg energy to eV.
- Does not download pseudopotentials or select convergence parameters.

### GPAW

`execute_gpaw_scf`:
- Imports GPAW only when explicitly called.
- Prepares the ASE structure and GPAW calculator, attaches it, and requests the potential energy.
- Reports calculation exceptions as failures.
- Reports SCF convergence as confirmed only when the calculator exposes a positive convergence flag; otherwise the result clearly marks convergence as unconfirmed.

## Deliberate limitations

- Geometry relaxation, variable-cell relaxation, bands, DOS, forces, stress, charge density, and magnetic moments are not yet part of the execution API.
- No automatic MPI launch, scheduler integration, cancellation, resource estimation, or checkpoint/restart is implemented.
- QE execution assumes the conventional local `pseudo` directory used by this entry point.
- A returned total energy is not, by itself, proof of thermodynamic stability, synthesizability, or agreement with experiment.
- Scientific adequacy still depends on convergence testing for cutoff, k-point mesh, pseudopotential choice, cell size, and the chosen functional.

