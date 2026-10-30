<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 | Status: Initial outline -->

# Discovery Engine

## Purpose

Specify the planned forward and inverse candidate-search workflows.

## Forward discovery

**Input:** a set of allowed elements, plus optional bounds and chemical or structural constraints.

**Planned stages:**

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

## Formula-based exploration

When a fixed formula is supplied, the engine should focus on possible structures, connectivity, isomers, conformers, and geometries consistent with that formula and the specified conditions.

## Inverse discovery

**Input:** target properties, operating conditions, and restrictions.

Potential search dimensions include element combinations, stoichiometry, crystal prototypes, substitutions, vacancies, interstitials, dopants, and structural parameters.

The engine should compare calculated or database properties with target ranges, rank candidates by a transparent objective, and report trade-offs. A target combination may be physically incompatible or practically unsynthesizable; the software must be allowed to return no suitable candidate.

## Candidate ranking

Ranking should be multi-criteria and explainable. It may include target-property fit, stability evidence, uncertainty, element availability, toxicity, cost, and synthesis considerations. Weights and hard constraints must be visible to the user.

## Evidence categories

Candidate records should distinguish experimental validation, database-supported evidence, theoretical prediction, exploratory hypothesis, and unresolved status. A score must never be presented as experimental certainty.
