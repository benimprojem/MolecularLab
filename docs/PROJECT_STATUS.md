# Molecular Lab — Project Status & Development Roadmap

**Project:** Molecular Lab  
**Repository:** [benimprojem/MolecularLab](https://github.com/benimprojem/MolecularLab)  
**Status date:** 2026-10-30  
**Document role:** Canonical cross-chat progress record  
**Current workstream:** Periodic electronic structure → Quantum ESPRESSO fat bands

> This file is the hand-off point between separate Molecular Lab conversation areas. Before starting work in a new conversation, read this document and inspect the latest repository state. Update this file whenever a meaningful feature, limitation, test result, or next milestone changes.

## 1. Current project overview

Molecular Lab is an early-stage computational chemistry and materials-science application. It is being developed in independent but connected workstreams:

1. Chemical identity, formula parsing, molecular graphs and candidate generation.
2. Molecular geometry and future molecular quantum-chemistry workflows.
3. Crystal structures, defects and periodic material models.
4. DFT execution through Quantum ESPRESSO (QE) and GPAW.
5. Electronic structure: bands, band gaps, DOS, PDOS and fat bands.
6. Material response and screening: optical, magnetic and electrical properties.
7. Metamaterials and multi-physics material design.
8. Future user interface, calculation management and visualization.

The repository is a research/prototype codebase. A feature being implemented or passing mocked tests does not imply that a real scientific calculation has been run or experimentally validated.

## 2. Overall status

| Workstream | Current state | Notes |
|---|---|---|
| Element identity and chemical formula parsing | Implemented foundation | 118-element identity data and syntax-level formula parser |
| Composition enumeration and diagnostics | Prototype implemented | Bounded enumeration; advisory plausibility checks |
| Molecular graph generation | Prototype implemented | Restricted acyclic/unicyclic graph paths and fused benzenoid examples |
| General molecular 3D geometry and conformers | Not implemented | Needs geometry engine and validation |
| Crystal structures and sites | Implemented foundation | Periodic cell, sites, defects and related records |
| Shared DFT settings and ASE bridge | Implemented foundation | Backend-neutral settings and structure conversion |
| QE SCF / RELAX / VC_RELAX | Implemented prototype | QE execution and parsed energy/structure/forces/stress |
| GPAW SCF / RELAX / VC_RELAX | Implemented prototype | GPAW execution through ASE |
| QE band structure | Implemented prototype | SCF → pw.x BANDS → bands.x → .gnu parser |
| QE total DOS | Implemented prototype | SCF → NSCF → dos.x |
| QE atom/orbital PDOS | Implemented prototype | SCF → NSCF → projwfc.x; collinear scalar-relativistic support |
| Band-edge / PDOS correlation | Implemented analysis helper | Energy-window summaries; not k-resolved band character |
| QE k-resolved projection parser, fat-band data model and QE execution | Implemented prototype | Opt-in `DFTOutput.FAT_BAND_STRUCTURE`; synthetic workflow tests pass; real QE run still pending |
| Fat-band plot-data preparation | Implemented prototype | Grouped series, reciprocal path distance, symmetry ticks, filters and display-only normalization |
| GPAW BANDS / DOS | Not implemented | Explicitly reported as unavailable |
| Spin-polarized QE BANDS | Not implemented | Current BANDS runner rejects it |
| Molecular electronic structure and spectra | Not implemented | Future work |
| Job scheduling, cancellation and production calculation manager | Not implemented | Current runners are direct sequential execution prototypes |
| Experimental validation | Not performed | No material property claims should be presented as measured |

## 3. Latest completed work

### QE fat-band foundation

Files added:

- `molecular_lab/calculations/qe_projection.py`
  - Parses newer `EIGENSTATES/PROJS/ATOMIC_WFC` and legacy `PROJECTIONS/K-POINT.N/(SPIN.N)/ATMWFC.N` XML layouts.
  - Retains complex projection coefficients and raw `|c|²` weights.
  - Maps QE atomic-wavefunction indices to atom, element and angular-momentum metadata from `projwfc.x` standard output.
  - Rejects unsupported spinor/four-component input explicitly.

- `molecular_lab/calculations/fatband.py`
  - Defines atom/orbital-resolved contributions and spin-resolved fat-band points/datasets.
  - Joins projections to band energies by k-point order, band index and spin.
  - Retains raw weights; it does not silently normalize them.

- `tests/test_qe_projection.py`
  - Covers new and legacy XML layouts, complex coefficients, spin channels, malformed inputs and band/projection joining.

- `docs/FAT_BANDS.md`
  - Describes data sources, mathematical weights, index conventions and scientific limitations.

Also updated `molecular_lab/calculations/__init__.py`, `docs/ELECTRONIC_STRUCTURE.md` and `README.md`.

**Validation status:** Latest GitHub Actions run passed on commit `bb38b9fde288bbbbebcc28ad5bdefaf76d594b88` (274 passed, 1 skipped with RDKit; 269 passed, 6 skipped without RDKit). The run includes a regression test proving stale `atomic_proj.xml` is not reused. No real QE calculation has been run.

## 4. Important technical decisions and cautions

- PDOS and fat bands are different data products. PDOS is energy-resolved and summed/broadened over k-points; fat bands retain the individual band and k-point identity.
- Projection weights are pseudopotential/projector-dependent quantities, not unique experimental observables.
- Preserve raw projection weights. Any display normalization must be a separate derived value.
- Projection and energy records must originate from the same saved QE calculation and identical k-point ordering.
- Band crossings require care: `bands.x` may reorder eigenvalue curves. Do not blindly attach raw `projwfc.x` state weights to a reordered band index. Prefer the original QE eigenstate identity/energy sequence for the fat-band join, and keep the conventional `.gnu` band plot as a separate representation.
- QE's `atomic_proj.xml` contains the per-state projection coefficients. The XML is written under the calculation save directory; `lwrite_overlaps` concerns optional overlap-matrix output and is not the switch that enables the projection coefficients.
- The current parser handles scalar/collinear spin channels but explicitly does not implement non-collinear spinor projections or fully relativistic j-resolved metadata.
- No automated test result should be described as a real QE run.

## 5. Completed milestone: QE fat-band execution pipeline

**Goal achieved at prototype level:** A single BANDS request can return both conventional band data and atom/orbital-resolved fat-band data.

Implemented sequence:

1. Reuse the existing QE SCF stage.
2. Run `pw.x` BANDS on the structured high-symmetry path.
3. Run `bands.x` and parse the conventional `.gnu` result.
4. Run `projwfc.x` on the same saved BANDS wavefunctions when `DFTOutput.FAT_BAND_STRUCTURE` is requested.
5. Read `atomic_proj.xml` from `tmp/molecular_lab.save`.
6. Parse the original QE eigenstate energies from `projwfc.x` output.
7. Reconstruct fractional k-point coordinates from the exact input `BandPath`; do not infer units from QE's printed k-vector.
8. Join raw projection weights to original eigenstate energies and return `DFTExecutionResult.fat_band_structure`.
9. Reject mismatched band/k-point dimensions and missing projection XML.
10. Clear stale `.gnu` and `atomic_proj.xml` outputs before rerunning stages, preventing accidental reuse of old calculation data.

The conventional `band_structure` remains based on `bands.x`. Fat-band energies retain original `pw.x` state ordering so projection weights are not silently attached to reordered band indices at crossings.

**Remaining validation:** The implementation has synthetic/mocked tests. Real QE execution, version-specific output validation and convergence checks remain outstanding.

## 6. Next milestones after fat-band execution

### Electronic-structure sequence

- [x] Complete QE fat-band execution integration with synthetic/mocked workflow tests.
- [ ] Verify parser and workflow against a real QE version and a small reference crystal.
- [x] Add plot-ready fat-band data export, grouping, filters and display-only normalization.
- [x] Add a dependency-free static SVG renderer for fat-band plots.
- [ ] Add an interactive desktop visualization component after selecting the UI framework.
- [ ] Add atom, element and orbital filters to the visualization.
- [ ] Add explicit display-only normalization options without modifying raw weights.
- [ ] Revisit band-crossing continuity and state tracking with reference cases.
- [ ] Add spin-polarized BANDS and spin-resolved fat bands.
- [ ] Add GPAW BANDS / DOS / projection capabilities where the backend supports them.
- [ ] Add k-resolved and fully relativistic projection support as separate, explicitly tested extensions.
- [ ] Add effective-mass and Fermi-surface analysis only after sampling/interpolation contracts are defined.

### Broader project sequence

- [ ] General molecular 3D geometry and conformer generation.
- [ ] Molecular optimization and force-field / quantum-chemistry backend integration.
- [ ] Molecular orbital, vibrational and optical spectra workflows.
- [ ] Expand crystal generation, defect/dopant handling and structural validation.
- [ ] Add calculation provenance, persistent job records, scheduling and cancellation.
- [ ] Build user-facing result visualization and material comparison tools.
- [ ] Keep optical, magnetic and electrical response models tied to explicit methods, units and evidence levels.
- [ ] Add real backend benchmark/reference calculations and record convergence settings.

## 7. Test and validation record

Latest completed GitHub Actions run for the parser, execution integration, coordinate reconstruction and stale-output regression:

- RDKit-enabled suite: **274 passed, 1 skipped**.
- Without RDKit: **269 passed, 6 skipped**.
- Verified code commit: `bb38b9fde288bbbbebcc28ad5bdefaf76d594b88`.
- Run: https://github.com/benimprojem/MolecularLab/actions/runs/41091007057

The QE fat-band parser/execution and plot-data tests are green. The SVG renderer test run passed on commit `1e3c97ff68f80f0276757816300b16a2fa021f24` (280 passed, 1 skipped with RDKit). This validates synthetic/mocked data only, not a real QE solver execution.

Validation levels must be recorded separately:

1. Static/code review.
2. Unit tests with synthetic or mocked QE output.
3. GitHub Actions CI.
4. Real QE/GPAW execution with identified software version and input files.
5. Scientific convergence checks.
6. Experimental comparison, if available.

## 8. Cross-chat hand-off protocol

At the beginning of each new Molecular Lab conversation:

1. Confirm the repository is `benimprojem/MolecularLab`.
2. Read this `docs/PROJECT_STATUS.md` file from the current default branch.
3. Check the latest commits and relevant source/tests before editing.
4. Continue only the active milestone unless the user explicitly changes workstream.
5. Do not recreate files that already exist.
6. After a meaningful change, update this report with:
   - completed work,
   - files changed,
   - test/validation result,
   - remaining limitations,
   - exact next step.
7. Update the README only when the user-facing project capability/status materially changes.

**Active hand-off point:** Design the interactive fat-band UI and validate the complete flow against a real QE reference calculation. The plot-data layer and static SVG renderer are implemented and synthetic tests pass. Continue from `DFTExecutionResult.fat_band_structure`, `prepare_fatband_plot_data` and `render_fatband_svg`. Do not restart the parser, execution pipeline, chemical graph, crystal foundation or PDOS implementation from scratch.
