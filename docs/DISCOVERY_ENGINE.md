<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 11 -->

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
6. Enumerate relevant isomers, conformers, or structural variants.
7. Build initial geometries.
8. Optimize and validate structures using available backends.
9. Evaluate requested properties and stability evidence.
10. Group, rank, and report candidates with confidence and provenance.

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

The molecular lookup distinguishes these outcomes:

- `EXACT_STRUCTURE`: an input InChIKey matches a reference record.
- `FORMULA_ONLY`: one or more reference records have the same element counts, but exact structure identity is not established.
- `NOT_FOUND`: the local catalogue contains no formula match.
- `NOT_APPLICABLE`: used by the integrated workflow when a periodic formula unit is presented to the molecule-only catalogue; a separate crystal/material reference source is required.

This distinction matters for isomers. Ethanol and dimethyl ether both have formula C2H6O, so formula alone cannot identify which structure was generated. A local `NOT_FOUND` result is not evidence of novelty; it only describes the current local catalogue. External PubChem synchronization/import is planned, not yet implemented.

### 2.5 Molecular graph records — implemented

The core domain layer now contains:

- `MolecularAtom`: explicit atom ID and element symbol.
- `MolecularBond`: an undirected bond between two atom IDs with bond order 1–3.
- `MolecularGraph`: atoms, bonds, proposal method, and explicit limitations.

The graph record validates unique atom IDs, valid bond endpoints, and duplicate bonds. It does not validate chemistry, geometry, or stability by itself.

### 2.6 Restricted molecular graph proposals — prototype implemented

`molecular_lab/engines/molecular_graphs.py` provides both `iter_saturated_acyclic_graphs` (the original compatibility prototype) and `iter_acyclic_graphs` (the expanded bond-order-aware generator), and `iter_unicyclic_graphs` (a bounded single-cycle generator).

The current prototype is intentionally restricted:

- Neutral compositions only.
- Supported elements: H, C, N, O, F, and Cl.
- Fixed common valences: H(1), C(4), N(3), O(2), F(1), Cl(1).
- Single, double, and triple bonds in the expanded generator.
- Connected acyclic or unicyclic heavy-atom skeletons, generated by separate functions; hydrogen atoms are added explicitly to fill remaining valences.
- Configurable candidate, heavy-atom, and bond-order-assignment bounds.

The heavy-atom tree topologies are generated from Prüfer sequences. Equivalent assignments among identical element symbols are deduplicated within the canonicalization search bound. The expanded generator can propose saturated structures such as H2O, H2O2, CH4, and C2H6O, as well as simple unsaturated structures such as C2H4, C2H2, CH2O, and CO2 when they fit the fixed-valence tree model. Bond-order assignments are bounded by `max_assignments` to prevent uncontrolled combinatorial growth.

The unicyclic generator adds one non-tree edge to a connected tree, producing a connected heavy-atom graph with exactly one independent cycle. It can represent simple three-, four-, five-, and six-membered rings as well as a single ring with acyclic substituents, subject to atom count, valence, formula, and search bounds. Simple cycles are enumerated before general ring-plus-substituent topologies so small ring systems are not crowded out by the larger topology search. For example, C6H6 can be proposed as a six-membered ring with alternating single/double (Kekule) bonds. When RDKit is available, the structure-identification adapter reports whether RDKit perceives aromatic atoms and bonds after sanitization. This is an aromaticity *recognition* step on a Kekule graph, not a general resonance-structure generator or a quantum-chemical aromaticity calculation.

It does **not** yet cover fused rings, bridged or spiro systems, multiple independent rings, charged species, radicals, hypervalent chemistry, coordination compounds, or elements outside the listed set. A composition that produces no graph means “not generated within this prototype's scope,” not “chemically impossible.”

A proposed graph is only a connectivity hypothesis. It is not a 3D geometry, optimized structure, energy minimum, stability prediction, experimental observation, or synthesis guarantee.

### 2.7 Graph-to-catalogue structure identification — implemented

`molecular_lab/engines/structure_identification.py` connects the restricted graph proposals to the offline reference catalogue.

For each graph, the adapter derives its elemental formula and attempts to construct a molecular representation with optional RDKit. When RDKit is available and sanitization succeeds, it calculates canonical SMILES and an InChIKey, then uses the InChIKey for exact catalogue matching. This allows the two generated C2H6O connectivity isomers to be distinguished as ethanol and dimethyl ether.

RDKit is a lazy optional dependency. If it is not installed, or if the graph cannot be converted, the adapter falls back to formula-level lookup and clearly reports that exact identity is unresolved. A conversion failure is not silently treated as a novel structure. RDKit's documented Python interface provides molecular construction, sanitization, canonical SMILES, and optional InChI support; installation is environment-dependent.

This adapter inherits the graph generator's narrow chemical coverage. It does not make the graph generator universal, and a successful identifier match does not validate 3D geometry, stability, or synthesis feasibility.

### 2.8 Integrated forward discovery pipeline — implemented

`molecular_lab/engines/forward_discovery.py` composes the existing stages behind `iter_forward_discovery` and `run_forward_discovery`.

For each bounded composition, the workflow:
1. Runs advisory electron-count and formal oxidation-state diagnostics.
2. Checks the local reference catalogue at formula level.
3. For neutral molecular compositions, asks the restricted graph prototype for connectivity proposals.
4. Runs each proposed graph through the optional structure-identification adapter.
5. Returns one `ForwardDiscoveryResult` containing the composition, diagnostics, reference lookup, graphs, identifiers, and explicit graph-generation status.

The workflow is bounded by `max_compositions`, `max_graphs_per_composition`, `max_heavy_atoms`, and `max_graph_assignments` (the last bounds bond-order combinations examined per composition in each graph family). It is lazy through `iter_forward_discovery`; `run_forward_discovery` materializes the bounded results as a tuple.

Diagnostics are advisory and do not automatically reject compositions. A composition outside the graph prototype's scope remains visible with `no_graph_within_prototype_scope`; charged and formula-unit records are also retained with explicit statuses rather than being misreported as impossible. Periodic-only searches are rejected because crystal structure generation is not implemented.

The pipeline does not yet create 3D coordinates, optimize geometry, calculate energies, assess thermodynamic/kinetic stability, or calculate material properties. Formula-level catalogue matches remain distinct from exact graph identity, which depends on successful structure-identifier generation. Formula-unit compositions are explicitly marked `NOT_APPLICABLE` for the molecule-only catalogue rather than being treated as unknown molecules.

A minimal H/O search can be run as follows:

```python
from molecular_lab.core.models import SearchConstraints, SearchDomain
from molecular_lab.engines.forward_discovery import iter_forward_discovery

constraints = SearchConstraints(
    allowed_elements=frozenset({"H", "O"}),
    min_total_atoms=1,
    max_total_atoms=4,
    domain=SearchDomain.MOLECULE,
)

for result in iter_forward_discovery(
    constraints,
    required_elements={"H", "O"},
    max_compositions=100,
):
    formula = result.composition.counts
    print(
        formula,
        result.graph_generation_status,
        result.reference_lookup.match_level.value,
        len(result.graphs),
    )
```

The output includes compositions that have no graph proposal under the current prototype. Consumers should inspect `graph_generation_note` and diagnostic findings rather than treating an empty graph list as a chemical impossibility. If RDKit is unavailable, graph identity remains formula-level.

## 3. Formula-based exploration

When a fixed formula is supplied, the engine should preserve its atom counts and explore only structures consistent with those counts. As coverage expands, this will include connectivity, isomers, conformers, and geometries. Current graph proposals cover only the restricted bond-order-aware acyclic and unicyclic domains above.

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
- Restricted neutral saturated and bond-order-aware acyclic graph proposal prototypes.
- Bounded neutral unicyclic graph proposal prototype with simple-cycle-first enumeration and single, double, and triple bonds; optional RDKit aromaticity reporting after graph sanitization.
- Graph-to-catalogue identification adapter with optional RDKit canonical SMILES/InChIKey generation and formula-only fallback.
- Integrated bounded forward molecular discovery pipeline combining enumeration, advisory diagnostics, catalogue lookup, graph proposals, and structure identification.

Not implemented:
- General-purpose molecular graph generation beyond the stated prototype scope.
- RDKit installation/configuration in a packaged application environment; the adapter is optional and lazy-loaded.
- Periodic/crystal structure generation and periodic-only forward discovery.
- Fused, bridged, spiro, or multi-ring graph generation; general resonance enumeration; radical, charged, hypervalent, or coordination graph generation. Aromaticity is only recognized opportunistically by RDKit for supported Kekule graphs.
- Crystal structure generation.
- 3D geometry generation/optimization, stability assessment, or material-property calculations.

Tests are authored in `tests/`. Execution verification remains pending.
