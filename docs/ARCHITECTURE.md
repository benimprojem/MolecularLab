<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 12 | Status: Initial SCF execution layer --> 

# System Architecture

## Purpose
Define the boundaries and responsibilities of MolecularLab modules. This document distinguishes implemented software contracts from scientific calculations that still require external backends.

## High-level layers

1. **Presentation layer:** desktop UI, input validation feedback, result visualization.
2. **Application layer:** discovery-session orchestration, job scheduling, configuration.
3. **Domain layer:** chemical entities, structures, candidate records, material properties.
4. **Calculation layer:** backend adapters, execution settings, convergence and result normalization.
5. **Persistence layer:** local database, calculation metadata, cached structures, provenance.

## Planned package layout

- `molecular_lab/core/`: shared domain models, configuration, orchestration, and task lifecycle.
- `molecular_lab/chemistry/`: element definitions, formula parsing, charge, spin, and chemical constraints.
- `molecular_lab/structures/`: molecular graphs, conformers, 3D coordinates, and crystal structures.
- `molecular_lab/calculations/`: calculation interfaces and backend adapters.
- `molecular_lab/materials/`: phases, defects, optical, electrical, mechanical, and thermal properties.
- `molecular_lab/engines/`: forward/inverse discovery workflows and solver orchestration.
- `molecular_lab/database/`: local persistence, reference-molecule catalogue, external database adapters, and provenance.
- `ui/`: desktop application.
- `tests/`: unit, integration, and scientific validation tests.

## Domain and calculation contracts

The shared domain contract is specified in [DATA_MODEL.md](DATA_MODEL.md). Domain records remain independent of scientific backends. Adapters convert between Molecular Lab structures/settings and solver-specific inputs while preserving atom identity, units, conditions, and provenance.

The common solver orchestration contract is implemented in `molecular_lab/engines/calculation_protocol.py`. Backend-neutral periodic DFT inputs are modeled in `molecular_lab/calculations/dft_settings.py`. The optional `ase_bridge.py` converts supported ordered crystal records to ASE Atoms without making ASE a core dependency.

- `CalculationRequest` carries a request ID, physics domain, named method, scalar parameters, operating conditions, requested outputs, and optional structure reference.
- `SolverCapability` declares a solver's supported domains, methods, and optional dependencies.
- `SolverAvailability` distinguishes a declared adapter from an installed/usable backend.
- `CalculationResult` normalizes success/failure/unsupported/unavailable states and carries scalar outputs, frequency-dependent response datasets, spin-resolved band structures, and total/projected DOS records.
- `DFTCalculationSettings` records the task, exchange-correlation identifier, plane-wave cutoff with units, k-point mesh, spin/SOC flags, optional smearing, convergence criteria, pseudopotential references, conditions, and requested outputs.
- The optional ASE bridge preserves cell vectors, fractional positions, PBC flags, and site metadata. Partial occupancy and defect annotations are rejected until explicit ordered structures are materialized.
- The Quantum ESPRESSO adapter generates `pw.x` input decks for SCF/RELAX/VC-RELAX. It does not launch QE during input preparation and rejects SOC, site-specific magnetic initialization, unsupported outputs, BANDS/DOS workflows, and unsupported k-point conventions.
- The GPAW adapter prepares ASE Atoms and GPAW calculator parameters for SCF/RELAX/VC-RELAX, and can instantiate a calculator or ASE optimization target. GPAW uses local PAW setups; QE-style pseudopotential file references are deliberately not consumed.
- The initial explicit execution API lives in `molecular_lab/calculations/execution.py`. It currently runs SCF total-energy requests only: QE is launched as a subprocess with a timeout and preflight checks for its executable and referenced UPF files; GPAW is invoked lazily through ASE's calculator API. No job is launched during module import or input preparation.
- QE execution is marked successful only when the process exits normally, its output contains the SCF convergence marker and `JOB DONE.`, and a total energy is parsed. GPAW's convergence is marked confirmed only if the installed calculator exposes a positive SCF convergence flag; an energy with an unavailable convergence flag is explicitly distinguished.
- Both execution paths normalize total energy to eV and retain backend identity/version where available. The current result record is a first execution slice, not yet the full scientific provenance record required for persisted publication-quality results.

## Periodic structure contract

Backend-independent crystal, site, and defect records are implemented in `molecular_lab/structures/crystal.py`. They preserve row-vector lattice geometry in angstrom, fractional atomic coordinates, per-site occupancy, per-axis periodic boundary conditions, optional magnetic/oxidation metadata, and explicit vacancy/substitution/interstitial/antisite annotations. These are structural records only: defect annotations do not build supercells or calculate formation energies. See [CRYSTAL_STRUCTURES.md](CRYSTAL_STRUCTURES.md).

## Electronic-structure result contract

The initial electronic-structure records live in `molecular_lab/materials/electronic_structure.py`. They retain k-point coordinates/weights, spin-resolved band energies, optional occupation arrays, Fermi level, direct and indirect band-edge estimates, and total/projected DOS grids. The band-gap helper reports metallic or undetermined states explicitly when occupations, Fermi-level information, or sampled edges do not support a gap estimate. It does not run DFT, interpolate bands, or establish experimental conductivity. See [ELECTRONIC_STRUCTURE.md](ELECTRONIC_STRUCTURE.md).

## Adapter boundary and execution rules

Each optional adapter imports its scientific dependencies lazily. The common solver protocol uses the Python standard library and does not require optional scientific software.

A backend exception is converted to an explicit failed calculation result at the adapter boundary. A result with a mismatched request ID or solver identity is rejected. Scientific outputs must include units, conditions, method/backend provenance, evidence category, limitations, and convergence information where relevant.

The initial DFT execution entry points are synchronous and deliberately limited. A persistent job queue, cancellation, progress reporting, resource scheduling, retry policy, MPI orchestration, checkpoint/restart, and distributed execution remain application-layer work and must not be simulated by this registry.

## Design requirements

- Keep calculation backends behind explicit interfaces.
- Store units and calculation conditions with every numerical property.
- Keep candidate generation separate from candidate validation and ranking.
- Make long-running calculations cancellable and report their status.
- Preserve provenance from input through final result.
- Check generated molecular candidates against local and optional external reference databases before presenting them as unregistered candidates.
- Do not silently substitute one computational method for another.
- Keep composition, structural identity, and calculation result as distinct records.
- Version persisted data schemas independently from application release versions.

## Open decisions

The UI framework, persistence technology, job scheduler, internal unit convention, and scientific validation policy will be selected during implementation and documented here when decided.

## Implementation status

Implemented: backend-neutral DFT settings; optional ASE structure conversion; QE `pw.x` input-deck generation and initial SCF runner; GPAW calculator preparation and initial SCF runner; common solver orchestration; molecular and material data models; bounded forward discovery; electronic-structure and periodic crystal records.

Not implemented: geometry optimization execution; BANDS/DOS multi-stage workflows; automatic pseudopotential/PAW setup acquisition and validation; SOC and species/sublattice magnetic initialization; broad crystal structure generation; stability assessment; persistent job queue/cancellation; and actual material-property calculations.
