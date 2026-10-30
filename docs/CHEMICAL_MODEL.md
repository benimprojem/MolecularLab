<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 | Status: Initial outline -->

# Chemical Model

## Purpose

Describe how Molecular Lab will represent chemical elements, compositions, molecular structures, and crystal structures.

## Planned domain concepts

- **Element:** atomic number, symbol, isotope information, and supported reference data.
- **Composition:** element-count mapping, total atom count, and optional net charge.
- **Molecular graph:** atoms as vertices and bonds as edges, with bond-order and stereochemical information where applicable.
- **3D structure:** atomic coordinates, units, coordinate frame, and structural constraints.
- **Crystal structure:** lattice vectors, atomic sites, periodicity, symmetry information, and occupancy.
- **Candidate:** a composition and one or more proposed structures, together with validation and evidence records.

## Composition generation

Forward discovery must accept an element set without requiring exact atom counts. It will enumerate only within explicit search bounds, such as maximum total atoms and per-element limits.

Composition enumeration is a proposal stage, not a declaration that a compound exists. Chemical plausibility filters should consider relevant valence patterns, electron count, charge, spin, bonding possibilities, and structural constraints. These filters must be method-aware and should not incorrectly rule out unusual but chemically meaningful systems.

## Isomers and structures

A single formula may correspond to multiple connectivities and 3D arrangements. Molecular graph generation, conformer generation, geometry optimization, and structure validation are separate steps.

For solids, composition alone is insufficient: candidate crystal prototypes, lattice parameters, symmetry, defects, and competing phases may need to be considered.

## Stability terminology

The software must distinguish local optimization convergence, local-minimum evidence, thermodynamic stability, kinetic stability, database or experimental confirmation, and synthesis feasibility. These are not interchangeable labels.

## Validation

The eventual implementation should include tests for formula parsing, atom counts, charge handling, graph validity, coordinate units, periodic boundary conditions, and round-trip serialization.
