<!--
Project: Molecular Lab
Version: 0.1
Rev.: 19
Status: Early foundation with ASE bridge and Quantum ESPRESSO input-deck generation
-->

# Molecular Lab

**Molecular Lab** is a planned Windows desktop environment for exploring molecules, compounds, crystals, and functional materials through computational chemistry and materials-science workflows.

> **Project status:** Early foundation. The repository includes the initial domain records, 118-element identity data, a syntax-level chemical formula parser, a bounded composition enumerator, advisory composition diagnostics, restricted bond-order-aware acyclic and unicyclic molecular graph prototypes, registered fused-benzenoid graph prototypes including anthracene and phenanthrene, a small offline reference-molecule catalogue, an optional RDKit graph-identification adapter, and an integrated bounded forward-discovery pipeline. Material design/target records, frequency-dependent complex response records, explicit-criterion electrical behavior screening, spin-resolved band-structure and DOS result models with direct/indirect gap screening, periodic crystal/site/defect records, backend-neutral periodic DFT input settings, an optional ASE crystal-structure bridge, a Quantum ESPRESSO `pw.x` input-deck writer, and a common solver orchestration protocol carrying structured electronic outputs are also implemented. These are data and software contracts and input-generation tools, not evidence that external scientific solvers have run. General structure generation beyond the bounded prototypes, geometry optimization execution, stability assessment, production solver execution adapters, and actual material-property calculations are not implemented yet.

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

- Optical transmission, reflection, absorption, refractive index, phase, polarization, and response over selected wavelength ranges.
- Magnetic response, permeability, susceptibility, anisotropy, resonance, and magnetic losses.
- Electrical conductivity/resistivity, dielectric response, impedance, and breakdown behavior.
- Electronic structure, band gap, density of states, carrier concentration, and semiconductor/insulator/conductor behavior (with conditions and classification criteria stated).
- Mechanical properties such as hardness, brittleness, and elastic moduli.
- Thermal and phase stability under specified temperature and pressure.
- Defects, vacancies, substitutions, dopants, and compositional constraints such as cost, abundance, toxicity, and rare-element avoidance.

## 3. Scientific integrity

Molecular Lab must distinguish experimental validation, database support, theoretical calculation, heuristic screening, and hypothesis. A proposed structure is not a confirmed substance; a local minimum is not automatically thermodynamic or kinetic stability; and a prediction must not be presented as an experimental fact.

Every numerical result should preserve units, operating conditions, method, software/backend version, evidence, uncertainty interpretation, convergence status, and limitations where applicable.

## 4. Calculation architecture

The initial common solver orchestration contract lives in `molecular_lab/engines/calculation_protocol.py`. It separates four calculation domains:

1. Electronic structure — band structures, density of states, band gaps, and Fermi levels.
2. Electrical transport — conductivity, resistivity, and carrier transport.
3. Magnetic response — susceptibility, permeability, anisotropy, and magnetic dynamics.
4. Electromagnetic response — optical spectra, scattering, fields, resonances, and photonic bands.

Each backend must declare the domains and methods it supports and report whether its software is actually available. The registry dispatches only compatible requests and returns explicit unsupported/unavailable/failed states rather than fabricating outputs.

Backend-neutral DFT settings are described in [docs/DFT_SETTINGS.md](docs/DFT_SETTINGS.md). The optional ASE conversion contract is described in [docs/ASE_BRIDGE.md](docs/ASE_BRIDGE.md). The initial Quantum ESPRESSO input writer and its limitations are documented in [docs/QUANTUM_ESPRESSO.md](docs/QUANTUM_ESPRESSO.md). Electronic-structure result schemas and their limits are described in [docs/ELECTRONIC_STRUCTURE.md](docs/ELECTRONIC_STRUCTURE.md) and [docs/CRYSTAL_STRUCTURES.md](docs/CRYSTAL_STRUCTURES.md). See also [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and [docs/METAMATERIALS.md](docs/METAMATERIALS.md).

## 5. Large systems and metamaterials

Large molecular systems require hierarchical fragment construction rather than exhaustive atom-by-atom enumeration. See [docs/LARGE_SYSTEMS.md](docs/LARGE_SYSTEMS.md) and [docs/FRAGMENTS.md](docs/FRAGMENTS.md).

Metamaterials are modeled as periodic geometries linked to constituent materials. Their optical, magnetic, electrical, and electronic responses are treated as distinct but potentially coupled physical domains. See [docs/METAMATERIALS.md](docs/METAMATERIALS.md).

## 6. Development principles

- Prefer bounded, deterministic, inspectable search procedures.
- Keep optional scientific packages behind adapters.
- Do not silently replace a requested method with a less accurate one.
- Keep units and provenance attached to every result.
- Distinguish search truncation from chemical rejection.
- Record model limitations and uncertainty instead of hiding them.
- Add tests for domain validation and backend boundaries before adding complex integrations.

## 7. Repository structure

- `molecular_lab/core/` — foundational domain models.
- `molecular_lab/chemistry/` — element data, formula parsing, and chemical diagnostics.
- `molecular_lab/engines/` — discovery engines and calculation orchestration.
- `molecular_lab/calculations/` — backend-neutral calculation settings, ASE conversion, and solver input adapters.
- `molecular_lab/materials/` — metamaterial, electronic-structure, response, and classification records.
- `docs/` — architecture and scientific design notes.
- `tests/` — automated tests.

## 8. License

This project is distributed under the MIT License. See [LICENSE](LICENSE).
