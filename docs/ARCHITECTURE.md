<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 | Status: Initial outline -->

# System Architecture

## Purpose

Define the planned boundaries and responsibilities of Molecular Lab modules. This document is an architectural outline; it does not describe implemented behavior.

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
- `molecular_lab/database/`: persistence and reference-data access.
- `ui/`: desktop application.
- `tests/`: unit, integration, and scientific validation tests.

## Design requirements

- Keep calculation backends behind explicit interfaces.
- Store units and calculation conditions with every numerical property.
- Keep candidate generation separate from candidate validation and ranking.
- Make long-running calculations cancellable and report their status.
- Preserve provenance from input through final result.
- Do not silently substitute one computational method for another.

## Open decisions

The UI framework, persistence technology, job scheduler, and exact backend interfaces will be selected during implementation and documented here when decided.
