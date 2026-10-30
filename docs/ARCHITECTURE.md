<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 01 | Status: Initial architecture refinement -->

# System Architecture

## Purpose

Define the planned boundaries and responsibilities of Molecular Lab modules. This document is an architectural outline; it does not describe implemented scientific behavior.

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
- `molecular_lab/engines/`: forward and inverse discovery workflows.
- `molecular_lab/database/`: local persistence, reference-molecule catalogue, external database adapters, and provenance-preserving imports.
- `ui/`: desktop application.
- `tests/`: unit, integration, and scientific validation tests.

## Domain data contract

The shared domain contract is specified in [DATA_MODEL.md](DATA_MODEL.md). It defines the initial records and relationships for elements, isotopes, compositions, search constraints, molecular graphs, geometries, crystal structures, candidates, calculations, properties, and evidence.

Domain records should remain independent of specific scientific backends. Adapters are responsible for conversion and must preserve atom identity, units, conditions, and provenance.

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

The UI framework, persistence technology, job scheduler, exact backend interfaces, internal unit convention, and validation library will be selected during implementation and documented here when decided.

## Implementation status

The initial offline reference-molecule catalogue and lookup functions are implemented under `molecular_lab/database/`. General persistence, external database synchronization, application orchestration, and scientific backend integration remain unimplemented.
