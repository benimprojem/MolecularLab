# DFT Execution Layer

**Status:** Initial execution prototype  
**Scope:** SCF and fixed-cell ionic RELAX for Quantum ESPRESSO and GPAW

## Separation of responsibilities

- `dft_settings.py` describes backend-neutral calculation settings.
- Backend preparation modules translate settings into QE input text or GPAW/ASE objects.
- `execution.py` is the explicit execution entry point and normalizes result fields.
- No backend is installed, downloaded, or launched implicitly by importing MolecularLab.

## Common API

Use `execute_dft_calculation(backend, structure, settings, ...)` to select a backend through one entry point. Supported backend identifiers are `qe`, `quantum_espresso`, and `gpaw`.

The returned `DFTExecutionResult` has a consistent status, energy in eV, backend/version, convergence flag, and—after a successful RELAX—an updated `CrystalStructure` plus per-site forces in eV/angstrom. The original structure is not mutated.

## Supported tasks

### SCF

- Requests total energy.
- QE success requires a successful process exit, the SCF convergence marker, normal job completion, and a parsed total energy.
- GPAW executes through ASE's calculator API. Its convergence flag is marked confirmed only when the calculator exposes a positive convergence state.

### Fixed-cell ionic RELAX

- Requires 3D periodic boundary conditions and an explicit `force_tolerance_ev_per_angstrom`.
- QE uses `pw.x` ionic BFGS, requires the geometry-convergence marker, normal job completion, final energy, final atomic coordinates, and a parsed force for every atom.
- QE final forces are converted from Ry/Bohr to eV/angstrom.
- GPAW uses ASE BFGS with the requested maximum force threshold, then collects the final potential energy, forces, cell and scaled positions.
- The optimized structure preserves site IDs and site metadata, and is returned separately from the input structure.

## Quantum ESPRESSO requirements

`execute_quantum_espresso` requires an existing working directory, a QE executable, and all referenced UPF files under `working_directory/pseudo`. It creates a local `tmp` directory and enforces a caller-configurable process timeout. It does not download pseudopotentials or select convergence parameters.

## GPAW requirements

`execute_gpaw` imports GPAW and ASE only when called. GPAW uses local PAW setups; the QE-style pseudopotential references are not consumed. An optional log path is accepted; the optimizer writes to a separate `.optimizer` log.

## Deliberate limitations

- Variable-cell relaxation (`VC_RELAX`), bands, DOS, stress normalization, charge density, magnetic moments, and post-SCF workflows are not part of the common execution API yet.
- No automatic MPI launch, scheduler integration, cancellation, resource estimation, or checkpoint/restart is implemented.
- QE execution assumes the conventional local `pseudo` directory.
- GPAW optimizer execution requires an explicit force threshold; the application does not invent one.
- A returned total energy or optimized local structure is not, by itself, proof of thermodynamic stability, synthesizability, or agreement with experiment.
- Scientific adequacy still depends on convergence testing for cutoff, k-point mesh, pseudopotential/PAW setup, cell size, and the chosen functional.
