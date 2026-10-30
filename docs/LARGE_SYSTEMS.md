<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 -->

# Large Molecular Systems and Scalability

## Purpose

This document defines the planned path from the current small-molecule prototypes toward large organic molecules, macromolecules, polymers, biomolecules, supramolecular assemblies, and periodic materials.

Large-system support is a planned capability, not a claim about the current implementation.

## Three distinct scaling problems

Molecular Lab must treat these as separate capabilities:

1. **Representation:** storing large structures, atoms, bonds, fragments, stereochemistry, charges, coordinates, and metadata.
2. **Structure generation:** proposing connectivities, conformations, substitutions, and configurations without unbounded combinatorial search.
3. **Scientific calculation:** optimizing structures and calculating energies or properties with a method suitable for the system size and chemistry.

A model capable of storing 100,000 atoms does not imply that the application can enumerate all possible structures or perform a high-level quantum calculation on all of them.

## Current implementation limits

The present graph generators are intentionally bounded:

- General acyclic and unicyclic graph prototypes: default maximum of 8 heavy atoms.
- Registered fused-benzenoid templates: separate default maximum of 16 heavy atoms.
- Supported graph chemistry: neutral H/C/N/O/F/Cl with fixed common valences.
- Fused benzenoid generation uses registered layouts, not arbitrary polycyclic topology generation.

These are limits of the current prototypes, not fundamental limits of the overall application.

## Planned architecture

### Scalable representation

The structure model should evolve to support:

- Stable atom and bond identifiers.
- Formal charge, isotope, radical state, aromaticity, and stereochemistry.
- Optional 3D coordinates and periodic cell data.
- Fragment, residue, and repeat-unit annotations.
- Efficient sparse neighbor indexing.
- Streaming serialization and chunked loading.

### Hierarchical search

For large molecules, exhaustive atom-by-atom graph enumeration is generally impractical. The planned workflow is:

1. Search known structures and reference databases.
2. Identify reusable fragments, functional groups, ring systems, and repeat units.
3. Assemble candidates from validated fragments under valence and charge constraints.
4. Apply inexpensive graph and chemistry filters.
5. Deduplicate candidates using canonical identifiers where available.
6. Generate conformers for shortlisted candidates.
7. Refine a limited candidate set with progressively more expensive methods.

Search limits should include atom, bond, fragment, ring, conformer, candidate, runtime, and memory bounds. The application must report when a search was truncated.

### Specialized structure domains

Different domains need dedicated generators and validators:

- **Large organic molecules:** fragment assembly, functional-group templates, ring libraries, and stereochemistry-aware construction.
- **Polymers:** repeat units, chain length, tacticity, branching, cross-linking, and chain conformations.
- **Biomolecules:** residue-based models for peptides, proteins, nucleic acids, and cofactors.
- **Supramolecular systems:** multiple components, non-covalent contacts, host–guest arrangements, and assembly states.
- **Periodic materials:** crystal cells, symmetry, atomic sites, defects, substitutions, and surfaces. These require a periodic structure model, not an ordinary molecular graph.

These should be separate modules and adapters rather than extensions that force every domain through the small-molecule graph enumerator.

## Calculation methods and limits

The method must match the system size and target property:

- **Molecular mechanics / force fields:** useful for large-system conformational exploration and pre-optimization when the parameter set covers the chemistry. They do not generally predict electronic excitation or optical spectra.
- **Semiempirical methods such as xTB:** faster approximate electronic-structure calculations for supported elements and chemistries; applicability must be checked.
- **DFT:** useful for selected higher-accuracy calculations, but cost grows with system size and settings. Large systems may require fragmentation, embedding, linear-scaling methods, or HPC resources.
- **Fragmentation / QM-MM / embedding:** treats a chemically important region at higher accuracy while approximating its environment; boundary assumptions must be reported.
- **Molecular dynamics:** requires a suitable force field or electronic method, ensemble, timestep, and adequate sampling.

No method is universally suitable. Results must retain method, parameters, supported chemistry, system size, conditions, convergence, approximations, and uncertainty.

## Resource-aware execution

Large-system workflows should support runtime and memory budgets, job queues, cancellation, checkpoint/restart, parallel resource limits, reproducible random seeds, and separate screening/refinement stages.

Resource exhaustion, unsupported chemistry, or a truncated search must never be reported as proof that a candidate is chemically impossible.

## Evidence and reporting

Each candidate should retain its input structure, normalized representation, generator and version, search bounds, truncation notes, identifiers, calculation backend and settings, convergence status, conditions, assumptions, evidence category, and uncertainty.

A generated structure or computational minimum is not by itself proof of experimental existence, kinetic persistence, thermodynamic stability, or synthesizability.

## Roadmap

1. Extend the structure representation without breaking current small-molecule APIs.
2. Add a validated fragment library and fragment connection rules.
3. Implement bounded fragment-based large-molecule assembly and deduplication.
4. Add conformer generation and geometry validation adapters.
5. Integrate suitable force-field and semiempirical backends, followed by higher-level methods.
6. Add polymer, biomolecular, supramolecular, and periodic-material workflows as separate modules.
7. Add scaling tests for large graphs, memory use, indexing, serialization, and deterministic bounds.

## Current status

This is an architectural roadmap. Large-molecule fragment generation, polymer and biomolecular workflows, conformer generation, geometry optimization, and large-system calculations are not implemented yet. The current graph prototype limits remain in force until the corresponding modules are implemented and tested.
