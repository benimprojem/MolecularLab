<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 | Status: Initial outline -->

# Calculation Methods

## Purpose

Define how calculations will be selected, recorded, and explained. No calculation method is implemented by this document.

## Calculation record requirements

Every calculation record should include:

- Input structure and composition identifier.
- Requested property and units.
- Method, backend, backend version, and relevant settings.
- Basis set, pseudopotential, functional, or model identifier where applicable.
- Charge, spin, temperature, pressure, and boundary conditions when relevant.
- Convergence criteria and final status.
- Output values, uncertainty information, warnings, and provenance.

## Planned workflow tiers

1. **Reference and heuristic screening:** fast checks and database lookups.
2. **Approximate calculations:** lower-cost methods for candidate prioritization.
3. **Detailed calculations:** suitable molecular quantum chemistry or periodic electronic-structure methods.
4. **Specialized calculations:** phonons, molecular dynamics, optical response, defects, and other targeted analyses when supported.

## Formula documentation standard

Whenever a calculation formula is introduced, its documentation must state:

1. Mathematical expression.
2. Meaning and units of every variable.
3. Physical or chemical property being calculated.
4. Why the calculation is needed.
5. Where it appears in the workflow.
6. Validity range, assumptions, and limitations.
7. Which backend and numerical method implement it.

## Important distinctions

A geometry optimizer finding a stationary structure does not establish global stability. A calculated band gap depends on the electronic-structure method. Optical response may require frequency-dependent dielectric calculations, while photonic-crystal behavior requires electromagnetic modeling at a different length scale.

## Implementation status

No scientific calculations are currently implemented. Backend selection and method validation will be documented as integrations are added.
