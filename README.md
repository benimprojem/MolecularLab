<!--
Project: Molecular Lab
Version: 0.1
Rev.: 03
Status: Initial project foundation
-->

# Molecular Lab

**Molecular Lab** is a planned Windows desktop environment for exploring molecules, compounds, crystals, and functional materials through computational chemistry and materials-science workflows.

> **Project status:** Early foundation. The repository includes the initial domain records, 118-element identity data, a syntax-level chemical formula parser, a bounded composition enumerator, advisory composition diagnostics, a restricted molecular graph prototype, a small offline reference-molecule catalogue, and an optional RDKit adapter that identifies generated graphs against that catalogue. General structure generation, geometry optimization, stability assessment, and material-property calculations are not implemented yet.

## 1. Project goal

Molecular Lab is intended to support two complementary discovery directions:

- **Forward discovery:** start with a set of available chemical elements and search for chemically plausible compositions, structures, and candidate compounds.
- **Inverse discovery:** start with desired material properties and constraints, then search for compositions, crystal structures, defects, and dopants that may satisfy them.

The scope is intended to include molecules, solids, crystals, and eventually liquid, gas, and phase-related investigations.

## 2. Discovery modes

### 2.1 Forward discovery — elements to candidates

The user selects the elements that may participate in a candidate system. Exact atom counts do not have to be supplied.

For example, an H/O search may explore bounded, chemically plausible compositions such as H₂O and H₂O₂, then consider possible structures and isomers. The search engine must not interpret every mathematical combination of atom counts as a valid molecule. Chemical constraints, charge, electron count, bonding, geometry, and computational screening must guide candidate generation.

The search space will be bounded by configurable limits such as maximum total atoms and maximum count per element. Candidates should be classified by evidence level, rather than presented as equally reliable discoveries.

### 2.2 Formula-based structure exploration

A user may also provide a fixed formula, such as H₂O₂, and ask the application to explore possible connectivity patterns, isomers, conformers, and geometries.

### 2.3 Inverse discovery — target properties to candidates

The user specifies desired properties and constraints. Possible targets include:

- Optical transmission, reflection, absorption, or response over selected wavelength ranges.
- Electrical behavior, conductivity class, and band gap.
- Mechanical properties such as hardness, brittleness, and elastic moduli.
- Thermal and phase stability under specified temperature and pressure.
- Defects, vacancies, substitutions, interstitials, dopants, and their concentrations.
- Element availability, cost, toxicity, and avoidance constraints.

The search may include molecular compounds, bulk materials, crystal prototypes, and photonic structures where appropriate.

## 3. Scientific workflow

The intended balanced workflow is staged:

1. Validate inputs and define search bounds.
2. Generate candidate compositions and structures.
3. Apply chemical and structural plausibility filters.
4. Use databases, heuristics, or suitable fast models for initial screening.
5. Run approximate calculations on promising candidates.
6. Use higher-accuracy methods, including suitable quantum-chemistry or periodic DFT backends, for selected candidates.
7. Evaluate stability, competing phases, decomposition pathways, and synthesis considerations where the available methods permit.
8. Calculate requested properties, rank candidates, and report uncertainty and provenance.

A local energy minimum is not, by itself, proof of thermodynamic stability, kinetic persistence, experimental existence, or practical synthesizability. Results must clearly distinguish experimental evidence, database evidence, theoretical prediction, and unverified hypotheses.

## 4. Planned software architecture

The project is designed around modular components:

- **UI:** Windows desktop interface and user workflows.
- **Core:** configuration, orchestration, task management, and shared data models.
- **Chemistry:** element data, formula handling, valence and electron-count checks, and chemical constraints.
- **Structures:** molecular graphs, isomer enumeration, 3D structures, crystal prototypes, and geometry handling.
- **Calculations:** interfaces to computational backends and standardized calculation records.
- **Materials:** solid-state properties, phases, defects, dopants, and optical/material analysis.
- **Engines:** forward and inverse discovery workflows.
- **Database:** local data, reference datasets, candidate records, and provenance.

The architecture is intended to support optional external engines rather than assuming every dependency is installed.

## 5. Planned scientific software integrations

Potential integrations include:

| Tool | Intended role | Important limitation |
|---|---|---|
| RDKit | Molecular representations and cheminformatics | Not a universal chemistry engine; element and bonding coverage has limits |
| xTB | Approximate molecular calculations | Accuracy and supported systems depend on the method and elements |
| Psi4 | Molecular quantum chemistry | Primarily molecular workflows; method and basis-set choices matter |
| ASE | Atomistic workflow management | Orchestrates calculations; does not itself provide every calculator |
| Quantum ESPRESSO | Periodic electronic-structure calculations | Requires suitable pseudopotentials and careful convergence |
| GPAW / CP2K | Additional electronic-structure workflows | Capabilities, accuracy, and computational cost vary by setup |

These are candidate integrations, not claims that they are already installed or operational. RDKit is currently an optional, lazily loaded dependency for graph identification; install it with `pip install rdkit` in the project environment when this feature is needed. Its availability and InChI support should be checked in the target environment. Other backends may be added as the project develops.

## 6. Result quality and evidence

Every result should record, where available:

- Formula, charge, spin, structure, and calculation conditions.
- Backend, method, version, basis set or pseudopotential, and relevant settings.
- Energy and calculated properties with units.
- Convergence and validation status.
- Stability category and the evidence supporting it.
- Uncertainty, known limitations, and data provenance.

The interface should distinguish at least experimentally validated results, database-supported results, high-confidence theoretical candidates, exploratory predictions, and unresolved candidates. These categories must not be conflated.

## 7. Repository map

```text
MolecularLab/
├── README.md
├── docs/
│   ├── ARCHITECTURE.md
│   ├── CHEMICAL_MODEL.md
│   ├── DATA_MODEL.md
│   ├── DISCOVERY_ENGINE.md
│   ├── CALCULATIONS.md
│   ├── DISCOVERY_ENGINE.md
│   └── LIMITATIONS.md
├── molecular_lab/
│   ├── core/
│   ├── chemistry/
│   ├── engines/
│   ├── structures/
│   ├── calculations/
│   ├── materials/
│   ├── engines/
│   └── database/
├── ui/
└── tests/
```

## 8. Development principles

1. Scientific assumptions and search bounds must be explicit.
2. A generated candidate must not automatically be called a stable compound.
3. Experimental observations and theoretical predictions must be labeled separately.
4. Every numerical result must identify its method, conditions, units, and provenance.
5. Missing calculation backends must be reported honestly; the application must never fabricate completed calculations.
6. Fast screening and detailed verification should be separate, traceable stages.
7. Modules should remain replaceable so new algorithms and calculation backends can be integrated without redesigning the whole application.

## 9. Current status and roadmap

### Version 0.1 — project foundation

- [x] Establish repository structure.
- [x] Define the two discovery directions.
- [x] Document the scientific scope and development principles.
- [x] Define the initial shared domain data contract and foundational Python records.
- [x] Add the 118-element identity table and initial formula parser (syntax tests added; execution verification pending).
- [x] Implement bounded, deterministic composition enumeration (chemical plausibility filtering remains pending).
- [x] Add initial composition-level electron-count and limited oxidation-state diagnostics (advisory only).
- [x] Add a small offline reference-molecule catalogue with formula, synonym, and InChIKey lookup.
- [ ] Expand structure generation and validation beyond the restricted prototype.
- [ ] Add optional PubChem synchronization/import with rate limits, caching, and provenance retention.
- [ ] Integrate the first calculation backend.
- [ ] Build the first usable desktop workflow.

Version **1.0** will be reserved for the first complete, usable release. Documentation-only updates and small corrections do not require a project version increase; they are tracked through per-file revision numbers.

## 10. License and citation

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for the complete license text. The copyright notice currently uses the GitHub repository owner identifier `benimprojem`.
