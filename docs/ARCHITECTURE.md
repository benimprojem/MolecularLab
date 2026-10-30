<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 10 | Status: Initial architecture refinement -->

# System Architecture

## Purpose

Define the boundaries and responsibilities of Molecular Lab modules. This document distinguishes implemented software contracts from scientific calculations that still require external backends.

## High-level layers

1. **Presentation layer:** desktop UI, input validation feedback, result visualization.
2. **Application layer:** discovery-session orchestration, job scheduling, configuration.
3. **Domain layer:** chemical entities, structures, candidate records, material properties.
4. **Calculation layer:** backend adapters, execution settings, convergence and result normalization.
5. **Persistence layer:** local database, calculation metadata, cached structures, provenance.

## Planned package layout

- `molecular_lab/core/`: shared models, configuration, orchestration, and task lifecycle.
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

The common solver orchestration contract is implemented in `molecular_lab/engines/calculation_protocol.py`. Backend-neutral periodic DFT inputs are modeled in `molecular_lab/calculations/dft_settings.py`, independently of any particular solver. The optional `molecular_lab/calculations/ase_bridge.py` converts supported ordered crystal records to ASE Atoms without making ASE a core dependency. The first Quantum ESPRESSO adapter in `molecular_lab/calculations/quantum_espresso.py` generates `pw.x` input decks for a conservative subset of tasks; it does not launch the solver.

- `CalculationRequest` carries a request ID, physics domain, named method, scalar parameters, operating conditions, requested outputs, and optional structure reference.
- `SolverCapability` declares a solver's supported domains, methods, and optional dependencies.
- `SolverAvailability` distinguishes a declared adapter from an installed/usable backend.
- `SolverRegistry` registers adapters, lists availability, matches requests, and dispatches only to a compatible available adapter.
- `CalculationResult` normalizes success/failure/unsupported/unavailable states and carries scalar outputs, frequency-dependent response datasets, spin-resolved band structures, and total/projected DOS records.
- `DFTCalculationSettings` records the task, exchange-correlation identifier, plane-wave cutoff with units, k-point mesh, spin/SOC flags, optional smearing, convergence criteria, pseudopotential references, conditions, and requested outputs.
- The optional ASE bridge preserves cell vectors, fractional positions, PBC flags, and site metadata. Partial occupancy and defect annotations are rejected until explicit ordered structures are materialized.
- The QE input writer supports SCF, fixed-cell RELAX, and fully periodic VC-RELAX. It explicitly rejects SOC, site-specific magnetic initialization, unsupported outputs, BANDS/DOS workflows, and unsupported k-point conventions rather than silently dropping them.
- Structure-aware DFT validation checks non-periodic k-point directions and exact element coverage by pseudopotential references. Neither the settings model nor the QE writer verifies pseudopotential files or solver build capabilities.

## Periodic structure contract

Backend-independent crystal, site, and defect records are implemented in `molecular_lab/structures/crystal.py`. They preserve row-vector lattice geometry in angstrom, fractional atomic coordinates, per-site occupancy, per-axis periodic boundary conditions, optional magnetic/oxidation metadata, and explicit vacancy/substitution/interstitial/antisite annotations. These are structural records only: defect annotations do not build supercells or calculate formation energies. See [CRYSTAL_STRUCTURES.md](CRYSTAL_STRUCTURES.md).

## Electronic-structure result contract

The initial electronic-structure records live in `molecular_lab/materials/electronic_structure.py`. They retain k-point coordinates/weights, spin-resolved band energies, optional occupation arrays, Fermi level, direct and indirect band-edge estimates, and total/projected DOS grids. The band-gap helper reports metallic or undetermined states explicitly when occupations, Fermi-level information, or sampled edges do not support a gap estimate. It does not run DFT, interpolate bands, or establish experimental conductivity. See [ELECTRONIC_STRUCTURE.md](ELECTRONIC_STRUCTURE.md) for schema and limitations.

## Adapter boundary and execution rules

Each optional adapter implements a capability declaration, an availability check, and a run method. Adapter-specific packages are imported only by their adapter modules; the common protocol uses the Python standard library and does not require optional scientific software.

A backend exception is converted to an explicit failed calculation result at the adapter boundary. A result with a mismatched request ID or solver identity is rejected. Scientific outputs must include units, conditions, method/backend provenance, evidence category, limitations, and convergence information where relevant.

The current protocol is synchronous and intentionally small. A persistent job queue, cancellation, progress reporting, resource scheduling, retry policy, and distributed execution remain application-layer work and must not be simulated by this registry.

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

Implemented: backend-neutral DFT calculation input settings and validation; optional ASE structure conversion; a Quantum ESPRESSO `pw.x` input-deck writer for SCF/RELAX/VC-RELAX; initial offline reference-molecule catalogue and lookup functions; optional RDKit graph-identification adapter; bounded forward-discovery pipeline; material design/target and frequency-dependent response records; criterion-bound electrical behavior screening; electronic-structure band/DOS result records with band-edge screening; periodic crystal/site/defect records; and the common solver request/capability/availability/dispatch protocol with unit tests.

Not implemented: automatic pseudopotential acquisition/validation; QE executable detection or execution; GPAW adapter; BANDS/DOS multi-stage workflow; SOC and species/sublattice magnetic initialization; geometry optimization execution; broad crystal structure generation; stability assessment; job queue/cancellation; and actual material-property calculations.
