# DFT Execution Layer

**Status:** QE electronic post-processing execution prototype  
**Scope:** SCF/RELAX/VC_RELAX for QE/GPAW; BANDS, total DOS and atom/orbital PDOS workflows for QE

## Separation of responsibilities

- `dft_settings.py` describes backend-neutral calculation settings.
- Backend preparation modules translate settings into QE input text or GPAW/ASE objects.
- `execution.py` is the explicit execution entry point and normalizes result fields.
- No backend is installed, downloaded, or launched implicitly by importing MolecularLab.

## Common API

Use `execute_dft_calculation(backend, structure, settings, ...)` to select a backend through one entry point. Supported backend identifiers are `qe`, `quantum_espresso`, and `gpaw`.

The returned `DFTExecutionResult` has a consistent status, energy in eV, backend/version, convergence flag, and—after a successful structural optimization—an updated `CrystalStructure` plus per-site forces in eV/angstrom. When requested, the stress tensor is normalized to GPa using a tension-positive convention. The original structure is not mutated.

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

### Variable-cell ionic VC_RELAX

- Requires 3D periodic boundary conditions and explicit force and stress tolerances.
- The request must include total energy, forces and stress.
- QE uses `vc-relax`, ionic and cell BFGS, and maps `stress_tolerance_gpa` to `press_conv_thr` in kbar (1 GPa = 10 kbar). The target external pressure is mapped to QE's `press` in kbar.
- QE final cell vectors and fractional positions are returned as a new structure. The final stress tensor is parsed from the kbar columns and converted to GPa with the sign changed to the common tension-positive convention.
- GPAW uses ASE `FrechetCellFilter` with BFGS. GPAW stress is collected from ASE in eV/angstrom^3 and converted to GPa.
- GPAW's optimizer stopping threshold combines atomic and cell-filter generalized forces; it is not numerically identical to a pure stress threshold. MolecularLab therefore independently checks the final stress residual against the requested external pressure and reports `NOT_CONVERGED` if that check fails. This conservative post-check prevents a false success, though it can report non-convergence after the optimizer has stopped.

## Quantum ESPRESSO requirements

`execute_quantum_espresso` requires an existing working directory, a QE executable, and all referenced UPF files under `working_directory/pseudo`. It creates a local `tmp` directory and enforces a caller-configurable process timeout. It does not download pseudopotentials or select convergence parameters.

## GPAW requirements

`execute_gpaw` imports GPAW and ASE only when called. GPAW uses local PAW setups; the QE-style pseudopotential references are not consumed. An optional log path is accepted; the optimizer writes to a separate `.optimizer` log.

## Deliberate limitations

- QE BANDS, total DOS and atom/orbital PDOS execute as sequential workflows and populate common result records. GPAW BANDS/DOS, spin-polarized BANDS, j-resolved/k-resolved PDOS, charge density, magnetic moments, and other post-SCF workflows remain unsupported.
- No automatic MPI launch, scheduler integration, cancellation, resource estimation, or checkpoint/restart is implemented.
- QE execution assumes the conventional local `pseudo` directory.
- GPAW optimizer execution requires an explicit force threshold; the application does not invent one.
- A returned total energy or optimized local structure is not, by itself, proof of thermodynamic stability, synthesizability, or agreement with experiment.
- Scientific adequacy still depends on convergence testing for cutoff, k-point mesh, pseudopotential/PAW setup, cell size, and the chosen functional.


## Stress convention and provenance

The common `stress_gpa` tensor is expressed in GPa with tension positive. QE's printed hydrostatic pressure convention is compression positive, so its parsed kbar tensor is negated before conversion to GPa. GPAW values returned through ASE are already in ASE's stress convention and are converted from eV/angstrom^3 using 1 eV/angstrom^3 = 160.21766208 GPa.

The stress tolerance is an explicit request parameter, not an application-selected default. QE maps the requested tolerance to its native `press_conv_thr` pressure criterion; this is not represented as an independent component-by-component stress residual test in the common layer. GPAW's final stress residual is checked component-by-component against the requested external pressure after ASE optimization. The output convention is shared, but backend stopping criteria are not numerically interchangeable.


## QE BANDS workflow

execute_dft_calculation("qe", structure, bands_settings, working_directory=...) runs three gated stages: prerequisite SCF, non-self-consistent pw.x BANDS along the explicit BandPath, and bands.x. The workflow expects the generated filband.gnu file and parses it into BandStructureData. QE documents the .gnu band output in eV and path coordinates; MolecularLab labels the first coordinate as path_distance, not as a 3D reciprocal vector. [QE bands.x reference](https://www.quantum-espresso.org/Doc/INPUT_BANDS.html).

The current workflow rejects spin-polarized BANDS rather than silently merging spin channels. The parsed .gnu data does not contain occupations; a Fermi level is attached only when QE explicitly reports one. A gap assessment can therefore remain undetermined for insulating calculations whose output does not provide a usable Fermi level.

## QE total-DOS workflow

execute_dft_calculation("qe", structure, dos_settings, working_directory=...) runs prerequisite SCF, uniform-grid NSCF, and dos.x in sequence. The expected total-DOS file is parsed into DensityOfStatesData; energy is stored in eV and DOS in states/eV as specified by QE. [QE dos.x reference](https://www.quantum-espresso.org/Doc/INPUT_DOS.html).

For spin-polarized DOS, the parser stores the up/down columns separately. Extra integrated-DOS columns are not currently retained. Projected DOS can be requested with DFTOutput.PROJECTED_DOS and is generated by projwfc.x from the same NSCF wavefunctions.

Both workflows stop when any prerequisite process fails, times out, does not report normal completion, or fails to produce a parseable result file. The automated test suite mocks the external executables; a real QE installation, compatible pseudopotentials, and converged scientific settings are still required for actual calculations.


## QE atom/orbital projected-DOS workflow

Request `DFTTask.DOS` with explicit `DOSSettings` and include `DFTOutput.PROJECTED_DOS` in `requested_outputs`. The common dispatcher runs prerequisite SCF, uniform-grid NSCF, then `projwfc.x`. If `DFTOutput.DENSITY_OF_STATES` is also requested, it additionally runs `dos.x`. The returned `DFTExecutionResult.projected_density_of_states` contains `ProjectedDensityOfStatesData`.

The PDOS model retains:
- QE's one-based atom index and element symbol
- pseudopotential atomic-wavefunction index and angular momentum (l)
- separate spin-up and spin-down records for collinear spin-polarized calculations
- the local DOS for that atomic wavefunction and its (2l+1) angular components
- the common energy grid in eV and projected DOS in states/eV

QE defines these projections using the atomic wavefunctions supplied by the pseudopotentials. They are basis/projector-dependent decompositions, not unique experimentally observable quantities. The p and d component labels follow QE's documented ordering; f components remain generic angular labels. Fully relativistic (j)-resolved, k-resolved and non-collinear PDOS are not supported. [QE projwfc.x documentation](https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html).

The automated workflow tests mock `pw.x`, `dos.x` and `projwfc.x`; no real solver run or scientific convergence is implied.
