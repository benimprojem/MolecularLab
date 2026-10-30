<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 04 -->

# Discovery Engine

## 1. Purpose

Specify the forward and inverse candidate-search workflows and record the implemented early stages.

## 2. Forward discovery

**Input:** a set of allowed elements, plus optional bounds and chemical or structural constraints.

### 2.1 Workflow

1. Validate element symbols and user constraints.
2. Enumerate bounded compositions.
3. Apply explainable composition-level diagnostics.
4. Query the local reference catalogue and, when enabled, external reference databases.
5. Propose molecular graphs or crystal prototypes within explicitly stated model coverage.
5. Enumerate relevant isomers, conformers, or structural variants.
6. Build initial geometries.
7. Optimize and validate structures using available backends.
8. Evaluate requested properties and stability evidence.
9. Group, rank, and report candidates with confidence and provenance.

The engine must not treat all integer element-count combinations as chemically meaningful. Search limits must be visible and configurable. Unusual compositions should not be rejected solely by intuition; they should be assessed using documented rules and suitable methods.

### 2.2 Bounded composition enumeration — implemented

The first stage is implemented in `molecular_lab/engines/composition_enumerator.py` as the lazy iterator `iter_compositions`.

It enumerates element-count combinations within the shared `SearchConstraints` record:

- `allowed_elements`: permitted element symbols.
- `min_total_atoms` and `max_total_atoms`: inclusive total-atom bounds.
- `per_element_max`: optional per-element upper bounds.
- `allowed_charges`: optional net-charge choices. If omitted, neutral charge is used.
- `domain`: molecular, periodic, or mixed interpretation.

Selected elements are optional by default. The caller can pass `required_elements` to require one or more selected elements in every result. The iterator is deterministic and lazy. Molecular searches emit molecular compositions; periodic searches emit formula units; mixed searches emit both interpretations. These labels do not create structures.

Enumeration is not chemical validation. A composition is not automatically a molecule, crystal, stable compound, experimentally known substance, or synthesizable material.

### 2.3 Initial composition diagnostics — implemented

`molecular_lab/chemistry/plausibility.py` reports total electron count and parity when charge is known, and checks formal charge balance against a limited curated set of common oxidation states.

Electron parity is useful for later spin-state enumeration but is not a stability test. A matching formal oxidation-state assignment is not proof of bonding or existence. Failure to find a match is not a universal impossibility rule and triggers review rather than automatic rejection.

### 2.4 Reference-molecule lookup — implemented

The initial offline catalogue is stored in `molecular_lab/database/reference_molecules.json` and accessed through `molecular_lab/database/reference_molecules.py`.

The seed catalogue contains five PubChem-derived records: water, hydrogen peroxide, methane, ethanol, and dimethyl ether. Lookup is available by formula and preferred name/synonym. Each record retains its PubChem CID, source URL, retrieval date, formula, canonical SMILES, InChIKey, and selected external identifiers.

The lookup distinguishes three outcomes:

- `EXACT_STRUCTURE`: an input InChIKey matches a reference record.
- `FORMULA_ONLY`: one or more reference records have the same element counts, but exact structure identity is not established.
- `NOT_FOUND`: the local catalogue contains no formula match.

This distinction matters for isomers. Ethanol and dimethyl ether both have formula C2H6O, so formula alone cannot identify which structure was generated. A local `NOT_FOUND` result is not evidence of novelty; it only describes the current local catalogue. External PubChem synchronization/import is planned, not yet implemented.

### 2.5 Molecular graph records — implemented

The core domain layer now contains:

- `MolecularAtom`: explicit atom ID and element symbol.
- `MolecularBond`: an undirected bond between two atom IDs with bond order 1–3.
- `MolecularGraph`: atoms, bonds, proposal method, and explicit limitations.

The graph record validates unique atom IDs, valid bond endpoints, and duplicate bonds. It does not validate chemistry, geometry, or stability by itself.

### 2.6 Restricted molecular graph proposals — prototype implemented

`molecular_lab/engines/molecular_graphs.py` provides `iter_saturated_acyclic_graphs`.

The current prototype is intentionally restricted:

- Neutral compositions only.
- Supported elements: H, C, N, O, F, and Cl.
- Fixed common valences: H(1), C(4), N(3), O(2), F(1), Cl(1).
- Single bonds only.
- Connected, acyclic heavy-atom skeletons; hydrogen atoms are added explicitly to fill remaining valences.
- A configurable heavy-atom bound and output candidate cap.

The heavy-atom tree topologies are generated from Prüfer sequences. Equivalent assignments among identical element symbols are deduplicated within the canonicalization search bound. The prototype can propose familiar saturated structures such as H2O, H2O2, CH4, and C2H6O.

This model does **not** cover double/triple bonds, aromatic systems, rings, charged species, radicals, hypervalent chemistry, coordination compounds, or elements outside the listed set. For example, C2H4 requires unsaturation and therefore produces no result here. That means “not generated within this prototype's scope,” not “chemically impossible.”

A proposed graph is only a connectivity hypothesis. It is not a 3D geometry, optimized structure, energy minimum, stability prediction, experimental observation, or synthesis guarantee.

## 3. Formula-based exploration

When a fixed formula is supplied, the engine should preserve its atom counts and explore only structures consistent with those counts. As coverage expands, this will include connectivity, isomers, conformers, and geometries. Current graph proposals cover only the restricted saturated acyclic domain above.

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
- Offline reference-molecule seed catalogue with formula/name lookup and InChIKey identity matching.
- Initial electron-count and formal oxidation-state diagnostics.
- Backend-neutral molecular atom, bond, and graph records.
- Restricted neutral saturated acyclic graph proposal prototype.

Not implemented:
- General-purpose molecular graph generation beyond the stated prototype scope.
- Double/triple bond, ring, aromatic, radical, charged, hypervalent, or coordination graph generation.
- Crystal structure generation.
- 3D geometry generation/optimization, stability assessment, or material-property calculations.

Tests are authored in `tests/`. Execution verification remains pending.
