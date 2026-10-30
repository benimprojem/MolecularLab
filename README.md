<!--
Project: Molecular Lab
Version: 0.1
Rev.: 25
Status: QE BANDS and total DOS workflow execution prototype; SCF/RELAX/VC_RELAX on QE and GPAW
-->

# Molecular Lab

**Molecular Lab** is a planned Windows desktop environment for exploring molecules, compounds, crystals, and functional materials through computational chemistry and materials-science workflows.

> **Project status:** Early foundation. The repository includes the initial domain records, 118-element identity data, a syntax-level chemical formula parser, bounded composition enumeration, advisory composition diagnostics, restricted molecular graph prototypes, fused-benzenoid graph prototypes, a small offline reference-molecule catalogue, an optional RDKit graph-identification adapter, and a bounded forward-discovery pipeline. Material design/target records, frequency-dependent complex response records, explicit-criterion electrical behavior screening, spin-resolved band-structure and DOS result models, periodic crystal/site/defect records, backend-neutral periodic DFT settings, an optional ASE crystal-structure bridge, Quantum ESPRESSO input generation, GPAW calculator preparation, and a common explicit SCF/fixed-cell RELAX/variable-cell VC_RELAX execution API are also implemented. Structured reciprocal-space band paths, sequential QE SCF to BANDS to bands.x and SCF to NSCF to dos.x workflows, and parsers that normalize QE band/DOS files into shared electronic-structure records are implemented. GPAW BANDS/DOS and spin-polarized BANDS remain unsupported. Structural runs return a separate optimized crystal structure, per-site forces, and—when requested—a normalized stress tensor in GPa. VC_RELAX requires explicit force and stress tolerances. This is not yet a production calculation manager: BANDS/DOS workflows, stability assessment, broad structure generation, job scheduling/cancellation, and actual experimental validation remain future work. No real QE/GPAW scientific run is implied by the automated tests.

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

The common solver orchestration contract lives in `molecular_lab/engines/calculation_protocol.py`. It separates electronic structure, electrical transport, magnetic response, and electromagnetic response.

Backend-neutral DFT settings are described in [docs/DFT_SETTINGS.md](docs/DFT_SETTINGS.md). The optional ASE conversion contract is described in [docs/ASE_BRIDGE.md](docs/ASE_BRIDGE.md). The Quantum ESPRESSO input writer is documented in [docs/QUANTUM_ESPRESSO.md](docs/QUANTUM_ESPRESSO.md), and the GPAW calculator preparation adapter is documented in [docs/GPAW.md](docs/GPAW.md). The common SCF and fixed-cell RELAX execution API is described in [docs/DFT_EXECUTION.md](docs/DFT_EXECUTION.md). Electronic-structure result schemas and their limits are described in [docs/ELECTRONIC_STRUCTURE.md](docs/ELECTRONIC_STRUCTURE.md) and [docs/CRYSTAL_STRUCTURES.md](docs/CRYSTAL_STRUCTURES.md). See also [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and [docs/METAMATERIALS.md](docs/METAMATERIALS.md).

The common execution entry point is `execute_dft_calculation`, with explicit backend selection (`qe` or `gpaw`). SCF returns total energy. A successful fixed-cell RELAX returns the final total energy, a separate optimized `CrystalStructure`, and per-site forces in eV/angstrom. The original structure is not mutated. QE requires a local `pw.x` executable and all referenced UPF files in the caller's `pseudo` directory; its final coordinates and forces must be present in the output before a RELAX result is accepted. GPAW is imported lazily and uses ASE BFGS for ionic optimization. Both backends require an explicit force tolerance for RELAX. Variable-cell optimization is not yet exposed through the common execution API.

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
- `molecular_lab/calculations/` — backend-neutral settings, ASE conversion, solver adapters, and explicit execution entry points.
- `molecular_lab/materials/` — metamaterial, electronic-structure, response, and classification records.
- `docs/` — architecture and scientific design notes.
- `tests/` — automated tests.

## 8. License

This project is distributed under the MIT License. See [LICENSE](LICENSE).


### DFT execution update (Rev. 23)

The common DFT execution layer now includes SCF, fixed-cell RELAX and variable-cell VC_RELAX for Quantum ESPRESSO and GPAW. Structural results preserve the original input structure and return optimized structures separately. Requested stress tensors are normalized to GPa with a documented tension-positive sign convention. VC_RELAX requires explicit atomic-force and stress tolerances; GPAW's final stress residual is independently checked after ASE cell-filter optimization. This is an execution prototype, not evidence that a real backend calculation has been run or scientifically converged for any particular material.
