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
| Fat-band SVG renderer | Implemented prototype | Self-contained SVG; static output only, interactive UI not yet selected |
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
- [x] Add atom, element and orbital filters to plot-data preparation.
- [x] Add explicit display-only normalization options without modifying raw weights.
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

Latest completed GitHub Actions run including the fat-band plot-data and SVG renderer tests:

- RDKit-enabled suite: **280 passed, 1 skipped**.
- Without RDKit: test job **passed** (the latest no-RDKit count was not captured in this hand-off).
- Verified code commit: `1e3c97ff68f80f0276757816300b16a2fa021f24`.
- Run: https://github.com/benimprojem/MolecularLab/actions/runs/41091518752

The QE fat-band parser/execution, reciprocal path-distance calculation, grouping/filtering, display-only scaling and SVG renderer all pass synthetic/mocked tests. This does not constitute a real QE solver execution or scientific convergence validation.

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

**Documentation update (2026-10-30):** `docs/FAT_BANDS.md` now includes an end-to-end Python example showing the opt-in QE fat-band execution request, retrieval of `result.fat_band_structure`, plot-data preparation using the same `BandPath` and direct lattice vectors, and writing the self-contained SVG output. The example is illustrative and is not a validated material input deck. The latest status-report commit `dc0d8566a6734c66b2dc540c907ba3f892b0974e` also has a successful GitHub Actions run: https://github.com/benimprojem/MolecularLab/actions/runs/41091772699. The latest full code-test run remains the successful 280 passed / 1 skipped RDKit-enabled run recorded above.

**Example hardening (2026-10-30):** The end-to-end usage snippet now creates its QE working directory before calling the execution API, matching the runner's existing-directory contract. **Active hand-off point:** Validate the complete flow against a real QE reference calculation, then design the interactive fat-band UI. The plot-data layer, static SVG renderer and usage documentation are implemented; automated tests remain synthetic/mocked. Continue from `DFTExecutionResult.fat_band_structure`, `prepare_fatband_plot_data` and `render_fatband_svg`. Do not restart the parser, execution pipeline, chemical graph, crystal foundation or PDOS implementation from scratch.


## 9. Real QE validation preparation (2026-10-30)

- Added optional `DFTCalculationSettings.charge_density_cutoff` with type and
  lower-bound validation against `energy_cutoff`.
- QE input generation now emits `ecutrho` in Ry when explicitly supplied;
  it does not invent a pseudopotential-dependent default.
- Added tests for `ecutrho` emission and invalid lower values.
- Added `examples/qe_silicon_fatbands.py`: a real-solver smoke-test entry point
  using a primitive diamond silicon cell, caller-supplied PBE UPF/cutoffs,
  SCF → BANDS → `bands.x` → `projwfc.x` → SVG.
- Added `docs/QE_REFERENCE_VALIDATION.md` with setup, run command, provenance
  capture and convergence requirements.
- The smoke-test path is deliberately short (Γ-X-L-Γ); it verifies the
  workflow, not a complete publication-quality band path.
- No actual QE solver run has been performed by this change. Real execution,
  output inspection and convergence remain pending.

**Next exact step:** Run the silicon example in an environment with QE/ASE and
a licensed, identified PBE silicon UPF; capture the solver version, UPF hash,
cutoffs, k-mesh, output files and convergence evidence before marking the
real-QE validation milestone complete.


**CI update (2026-10-30):** The QE charge-density cutoff changes and the
silicon example are covered by the latest successful GitHub Actions run:
[run 41096358524](https://github.com/benimprojem/MolecularLab/actions/runs/41096358524),
commit `436ad12b69ca5412b64c812d3a4ba875f26755f1`. Both RDKit-enabled and
no-RDKit jobs passed. The new example test imports the script and verifies the
primitive diamond silicon cell builder. The real QE solver itself remains
unrun.


## 10. QE smoke-test provenance hardening (2026-10-30)

- The silicon example now writes `results/run_manifest.json` with QE result
  status/version/return code, energy, cell and atomic coordinates, pseudopotential
  source/license/family/SHA-256, cutoff values, SCF mesh and exact band-path
  vertices.
- Captures the execution result's stdout and stderr into separate files beside
  the SVG, including when the solver workflow returns a failure status.
- Added an optional pseudopotential license argument so the run record can
  preserve provenance/usage terms.
- Extracted the silicon smoke-test path builder and added tests for its
  fractional reciprocal coordinates, labels, sample count and invalid sampling.
- Fixed JSON manifest newline serialization.
- These changes improve reproducibility and diagnostics; they do not substitute
  for a real QE run. The solver execution and scientific convergence remain
  pending.


## 11. Silicon example input validation (2026-10-30)

- Tightened `build_silicon_band_path()` so its public helper rejects booleans
  and non-integer subdivision values, not only zero/negative values.
- Expanded the unit test to cover zero, negative, floating-point and boolean
  inputs. This closes a small validation gap: the CLI already parsed this
  argument as an integer, but direct Python callers could previously pass a
  float despite the helper's documented positive-integer contract.
- Commits: `65d9899586eae7f144ff5ad9d4f11e573fc7e1c2` (implementation) and
  `640e8203186469b56bca9992d4dc418086e5ba1c` (tests).
- CI for the earlier provenance/path/doc commits was still queued at the
  beginning of this follow-up. The two new commits also require CI confirmation;
  do not mark this change validated until their Actions runs complete.
- The QE result object and example were cross-checked: `status`,
  `message`, `backend_version`, `return_code`, `converged`,
  `total_energy_ev`, `stdout`, `stderr` and `fat_band_structure`
  match the actual `DFTExecutionResult` fields used by the example.
- Real QE execution remains pending until QE/ASE and the exact PBE Si UPF
  with its recommended cutoffs are available.


## 12. QE BANDS preflight checks (2026-10-30)

- Added a preflight before the prerequisite SCF stage in the QE BANDS workflow.
- Checks structure/settings pseudopotential references, required UPF files in
  the work directory, ASE discoverability, and executable resolution for
  `pw.x` and `bands.x`; `projwfc.x` is checked only for fat-band requests.
- Missing prerequisites return `BACKEND_UNAVAILABLE` with a list of issues.
  SCF is not started when preflight fails.
- Added tests for missing `projwfc.x` preventing SCF, missing UPF reporting,
  and plain BANDS not requiring `projwfc.x`.
- Updated `docs/QE_REFERENCE_VALIDATION.md` with the preflight contract.
- Commits: `4f3aaa5` (preflight implementation),
  `609c707` (ASE availability import),
  `f5cdd0d` (tests), `edac31c` (documentation).
- CI has not yet been confirmed for these commits. This preflight validates
  discoverability and required files only; it does not validate runtime
  compatibility, UPF contents, or scientific convergence.


## 13. Fat-band internal consistency checks (2026-10-30)

- Strengthened QE atomic projection validation: contiguous one-based projector metadata, contiguous atom indices, consistent element labels per atom, unique spin labels, and a complete unique spin × k-point × band × projector coefficient grid.
- Each complex coefficient now checks finite real/imaginary/weight values, non-negative weight, and weight = real² + imaginary² within numerical tolerance.
- Strengthened FatBandDataset validation: complete unique spin × k-point × band state grid, consistent atom/element identity, unique atom/orbital contribution keys and valid source-projector indices.
- Added regression tests for inconsistent weights, NaN values, missing or duplicate coefficient keys, out-of-range state keys, duplicate fat-band states and duplicate orbital contributions.
- Deliberately do not require projected weights to sum to one; incomplete projector bases and pseudopotential choices make that an invalid universal constraint. No automatic renormalization is introduced.
- This is structural/data-integrity validation, not proof that QE physical projection basis is complete or that a real calculation has converged.


**CI confirmation for milestone 13 (2026-10-30):** The final hardening commit
`aeadeb08285774b46a304e2eb1f944449ba13142` passed both GitHub Actions jobs:
- RDKit enabled: **294 passed, 1 skipped**.
- Without RDKit: **289 passed, 6 skipped**.
- Run: https://github.com/benimprojem/MolecularLab/actions/runs/41106841789

The final validation covers synthetic projection records and dataset invariants.
The next step is to add cross-checks against the crystal structure itself
(atom count and per-index element identity) at the QE workflow boundary, then
prepare a small reference-data comparison when a real QE environment and UPF
are available. Real QE execution and scientific convergence remain pending.


**Structure identity cross-check (2026-10-30):** The QE fat-band execution path now passes the input CrystalStructure into the dataset builder. The builder verifies that projection atom indices cover every input site and that each one-based QE atom index carries the same element as the corresponding input site. Added a regression test for matching and deliberately swapped Fe/O site order. This protects the projection join from a subtle but serious atom-order mismatch; it does not validate the physical correctness of the QE pseudopotential projectors.


**Atom-order validation CI (2026-10-30):** The QE workflow now checks projection atom indices and element labels against the exact input crystal-site order. An initial mocked workflow fixture intentionally/accidentally described a second O atom while the crystal fixture contained one Si site; the new validation correctly caught that inconsistency. The fixture was corrected to three projector states on the single Si site. Final CI passed on commit `c7160c80ed28b164a469b68a41b11c3e50d30f82`: RDKit enabled **295 passed, 1 skipped**; without RDKit **290 passed, 6 skipped**. Run: https://github.com/benimprojem/MolecularLab/actions/runs/41107182574.
