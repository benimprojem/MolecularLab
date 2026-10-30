<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 02 -->

# Discovery Engine

## 1. Purpose

Specify the forward and inverse candidate-search workflows. This document also records the first implemented forward-search stage: bounded composition enumeration.

## 2. Forward discovery

**Input:** a set of allowed elements, plus optional bounds and chemical or structural constraints.

### 2.1 Planned workflow

1. Validate element symbols and user constraints.
2. Enumerate bounded compositions.
3. Apply chemically informed plausibility filters.
4. Generate candidate molecular graphs or crystal prototypes.
5. Enumerate relevant isomers, conformers, or structural variants.
6. Build initial geometries.
7. Optimize and validate structures using available backends.
8. Evaluate requested properties and stability evidence.
9. Group, rank, and report candidates with confidence and provenance.

The engine must not treat all integer element-count combinations as chemically meaningful. Search limits must be visible and configurable. Unusual compositions should not be rejected solely by intuition; they should be assessed using documented rules and suitable methods.

### 2.2 Implemented composition enumeration

The first stage is implemented in molecular_lab/engines/composition_enumerator.py as the lazy iterator iter_compositions.

It enumerates element-count combinations within the shared SearchConstraints record:

- allowed_elements: permitted element symbols.
- min_total_atoms and max_total_atoms: inclusive total-atom bounds.
- per_element_max: optional per-element upper bounds. An omitted limit defaults to max_total_atoms.
- allowed_charges: optional net-charge choices. If omitted, neutral charge is used.
- domain: molecular, periodic, or mixed interpretation.

Selected elements are optional by default, so a result may contain any non-empty subset. The caller can pass required_elements to require one or more selected elements in every result. Required symbols must belong to allowed_elements.

The iterator is deterministic and lazy. Results are emitted by increasing total atom count, then sorted element symbols and ascending counts. This avoids holding the full result set in memory, although large search bounds can still produce a very large number of results. Cancellation, progress reporting, and a user-facing result cap remain future work.

Molecular searches emit CompositionKind.MOLECULAR. Periodic searches emit CompositionKind.FORMULA_UNIT. Mixed searches emit both interpretations for each count vector. These are search-domain labels only; no molecular graph or crystal structure is created.

For an H/O search with total atom counts bounded from 1 through 4, the iterator includes the count vectors for H₂O and H₂O₂. It also emits other mathematically permitted compositions. That is intentional: chemical plausibility is a separate stage and has not yet been implemented.

**Important:** enumeration is not a chemical validation step. A generated composition is not automatically a molecule, crystal, stable compound, experimentally known substance, or synthesizable material.

### 2.3 Initial composition diagnostics — implemented

The first advisory diagnostic module is molecular_lab/chemistry/plausibility.py. It reports:

- Total electron count and even/odd parity when net charge is specified.
- Whether a formal oxidation-state assignment can balance the supplied charge using a limited curated list of common states.
- Whether that check is inconclusive because one or more elements are not covered.

These checks are deliberately conservative. Electron parity is useful later for spin-state enumeration but is not a stability test. A matching oxidation-state assignment is not proof of bonding or existence. Failure to find a match is not a universal impossibility rule and must only trigger review. The current reference table is incomplete and is not a definitive list of allowed oxidation states.

### 2.4 Further chemical plausibility filters — planned

Potential filters may consider electron count, charge consistency, oxidation-state and valence patterns, bond-order and coordination constraints, known motifs, geometric feasibility, reference databases, and method applicability.

Filters must be staged and explainable. A cheap filter may lower a candidate's priority, but unusual chemistry must not be permanently rejected unless the rule is justified for the selected domain. A valence-rule match is not proof of stability, and a mismatch is not universally proof of impossibility.

## 3. Formula-based exploration

When a fixed formula is supplied, the engine should focus on possible structures, connectivity, isomers, conformers, and geometries consistent with that formula and the specified conditions. It must preserve the requested atom counts.

## 4. Inverse discovery

**Input:** target properties, operating conditions, and restrictions.

Potential search dimensions include element combinations, stoichiometry, crystal prototypes, substitutions, vacancies, interstitials, dopants, and structural parameters.

The engine should compare calculated or database properties with target ranges, rank candidates by a transparent objective, and report trade-offs. A target combination may be physically incompatible or practically unsynthesizable; the software must be allowed to return no suitable candidate.

## 5. Candidate ranking

Ranking should be multi-criteria and explainable. It may include target-property fit, stability evidence, uncertainty, element availability, toxicity, cost, and synthesis considerations. Weights and hard constraints must be visible to the user.

## 6. Evidence categories

Candidate records should distinguish experimental validation, database-supported evidence, theoretical prediction, heuristic screening, exploratory hypothesis, and unresolved status. A score must never be presented as experimental certainty.

## 7. Implementation status

Implemented:
- Bounded, deterministic, lazy composition enumeration.
- Per-element count bounds, required-element selection, charge alternatives, and molecular/periodic/mixed domain labels.

Not implemented:
- Chemical plausibility scoring or rejection.
- Molecular graph or crystal structure generation.
- Geometry optimization, stability assessment, or material-property calculations.

Tests are authored in tests/test_composition_enumerator.py. Execution verification is pending.
