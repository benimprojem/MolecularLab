<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 | Status: Initial outline -->

# Scientific Limitations and Responsible Reporting

## Purpose

Record limitations that must remain visible throughout development.

## Search-space limitations

The number of possible compositions, connectivities, conformers, crystal structures, and defect configurations can grow rapidly. Search bounds and pruning rules are necessary. A bounded search is not proof that all possible candidates have been considered.

## Model limitations

Every computational method has a domain of applicability. Approximate methods, machine-learning models, force fields, and quantum-chemistry methods may fail for particular elements, oxidation states, bonding patterns, phases, or physical conditions.

## Stability limitations

- Geometry convergence does not prove thermodynamic stability.
- A local minimum does not prove kinetic persistence.
- A low calculated energy does not by itself prove experimental existence.
- Stability may depend on temperature, pressure, environment, and competing phases.
- Synthesis feasibility requires evidence beyond a computed structure.

## Property limitations

Optical, electrical, mechanical, and thermal properties may require different methods and length scales. Molecular calculations cannot automatically predict bulk-material behavior. Photonic structures require electromagnetic modeling of spatially patterned media and should not be confused with intrinsic molecular optical response.

## Data and provenance

Reference databases may be incomplete, inconsistent, or method-dependent. Each imported datum should retain its source, conditions, units, and access date where available.

## Reporting rule

When a backend is unavailable, a calculation fails, or evidence is insufficient, Molecular Lab must state that clearly. It must not fabricate numerical results, imply that a calculation was run, or label an unverified candidate as experimentally confirmed.
