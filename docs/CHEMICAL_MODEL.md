<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 01 | Status: Domain model specification -->

# Chemical Model

## 1. Purpose and scope

This document defines the first domain-level model for representing chemical systems in Molecular Lab. It is a specification, not a claim that the corresponding classes or algorithms are already implemented.

The model must support both finite molecular systems and periodic materials. These share elements and compositions, but their structures and physical descriptions are not interchangeable.

## 2. Core concepts

### 2.1 Element and isotope

An **Element** is identified by its atomic number and chemical symbol. Element reference data may include standard atomic weight, allowed isotope records, covalent and van der Waals radii, electronegativity, and other sourced values.

An **Isotope** is identified by atomic number and mass number. Isotopic mass and abundance must be stored with source and provenance where known. An isotope choice is optional for ordinary composition searches.

Element data must not silently imply a fixed valence or oxidation state. Chemical behavior depends on bonding environment, charge, spin, and physical conditions.

### 2.2 Composition

A **Composition** is a normalized mapping from element/isotope identity to a non-negative integer count, plus a net charge when defined.

Examples:
- H₂O: H = 2, O = 1
- H₂O₂: H = 2, O = 2

The composition is not a structure. It does not specify which atoms are bonded, the 3D geometry, the phase, or whether the system is stable.

The model must distinguish:
- **Formula as entered:** the user's original formula text.
- **Normalized composition:** canonical element-count representation.
- **Empirical formula:** the simplest integer ratio, where relevant.
- **Molecular formula:** actual atom counts for a discrete molecule.
- **Cell composition:** atom counts associated with a periodic unit cell.

Do not reduce a molecular formula to an empirical formula when that would discard meaningful atom counts.

### 2.3 Charge and spin

Net charge is an integer property of a modeled system. Spin information should be represented explicitly, preferably through a well-defined multiplicity or spin-state model.

Charge and multiplicity must not be guessed silently. If they are unknown, represent them as unspecified and allow the workflow to enumerate supported charge/spin alternatives within user-defined limits.

Electron-count consistency checks are useful validation steps, but they do not prove that a proposed structure is chemically valid or stable.

### 2.4 Molecular graph

A **MolecularGraph** represents a finite molecular connectivity model:
- Atoms are vertices with element/isotope identity and optional atom-level attributes.
- Bonds are edges with bond order and optional aromatic, stereochemical, or coordination metadata.
- Formal charge and radical/electron annotations may be associated with atoms where the chosen representation supports them.

The graph represents connectivity, not final 3D coordinates. Bond order and valence are model-dependent; unusual bonding, delocalization, multicenter bonding, and coordination chemistry may require representations richer than a simple single/double/triple bond graph.

### 2.5 Molecular geometry

A **MolecularGeometry** contains atom coordinates, coordinate units, atom ordering/mapping to a molecular graph, and optional constraints. Coordinate arrays must identify their unit explicitly; internal canonical storage should use a single documented unit.

A geometry may be:
- an initial generated guess,
- an optimized stationary structure,
- a transition-state candidate,
- or a geometry imported from an external source.

These statuses must not be conflated.

### 2.6 Periodic crystal structure

A **CrystalStructure** contains:
- three lattice vectors and their units,
- atomic sites and fractional coordinates,
- element/isotope identity and occupancy,
- periodic boundary conditions,
- optional symmetry and magnetic information.

A crystal structure is not merely a large molecule. Periodicity, unit-cell choice, symmetry, defects, and magnetic ordering may affect both representation and calculations.

### 2.7 Candidate and structural alternatives

A **Candidate** is a proposed chemical system with one composition and one or more associated structural alternatives. Each alternative has its own identity, geometry or crystal structure, charge/spin state, validation results, and calculation history.

Candidates with the same composition must not be merged if they represent different isomers, conformers, polymorphs, magnetic states, or defect configurations.

## 3. Forward composition search

The forward engine accepts an allowed element set without requiring the user to supply exact atom counts.

Search must be bounded using explicit parameters, for example:
- maximum total atom count,
- minimum and maximum count per selected element,
- allowed net-charge range,
- allowed spin/multiplicity choices,
- molecular versus periodic search mode,
- optional element exclusions and required-element rules.

A generated composition is only a search hypothesis. It must pass through further structure generation and validation before it can become a ranked candidate.

### 3.1 Chemical plausibility filters

Filters may use:
- electron count and charge consistency,
- candidate oxidation-state or valence patterns,
- bond-order and coordination constraints,
- known chemical motifs and reference databases,
- geometric feasibility,
- method-specific model applicability.

Filters should be staged and explainable. A cheap filter may mark a candidate as low priority, but it should not permanently reject unusual chemistry unless the rejection rule is justified for the selected search domain.

A valence-rule match is not proof of stability. A valence-rule mismatch is not universally proof of impossibility.

### 3.2 Composition enumeration is not molecule enumeration

For a given composition, multiple molecular graphs or crystal structures may exist. Conversely, some compositions may not correspond to a bound molecule under the selected conditions. The software must keep these stages separate:

1. Composition proposal.
2. Structure/connectivity proposal.
3. Geometry or lattice construction.
4. Optimization and physical validation.
5. Stability and evidence assessment.

## 4. Formula-based exploration

When a user supplies a formula, Molecular Lab should preserve the exact atom counts and explore supported structural alternatives. It must not silently change the formula during isomer enumeration.

For molecular systems, alternatives may include connectivity graphs, stereoisomers, conformers, charge states, and spin states. For solids, alternatives may include crystal prototypes, polymorphs, unit cells, magnetic orderings, and defect configurations.

## 5. Stability and evidence vocabulary

Use separate fields for separate scientific questions:

- **Optimization status:** whether a numerical optimization converged.
- **Local-minimum evidence:** whether available analysis supports a local minimum.
- **Thermodynamic stability:** stability relative to specified competing states under stated conditions.
- **Kinetic persistence:** barriers or other evidence relevant to persistence over time.
- **Experimental evidence:** whether a matching material or species has been experimentally reported.
- **Synthesis evidence:** whether a preparation route has been reported or assessed.

A converged optimization is not proof of thermodynamic stability. A local minimum is not proof of experimental existence. Synthesis feasibility is not implied by a low computed energy.

## 6. Units, identity, and provenance

Every physical quantity must carry a unit. Every external reference datum must retain its source and, where possible, source version or access date. Atom indices must be stable within a structure record and mappings must be retained when structures are reordered by a backend.

## 7. Validation requirements

Future tests should cover:
- element and isotope identity,
- formula parsing and canonicalization,
- preservation of molecular atom counts,
- charge and electron-count bookkeeping,
- spin/multiplicity validation,
- graph connectivity and bond metadata,
- coordinate shape, atom mapping, and units,
- lattice-vector validity and fractional coordinates,
- serialization/deserialization without loss of chemical identity.

## 8. Implementation status

This specification is being established before implementation. No chemistry validation, composition enumeration, structure generation, or stability calculation is claimed to exist yet.
