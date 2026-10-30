# Real QE reference validation: silicon fat bands

This procedure is the next validation level after the synthetic/mocked tests.
It is intentionally separate from CI: the GitHub test jobs do not install QE,
download pseudopotentials, or claim that a solver calculation has completed.

## 1. Requirements

- A working Quantum ESPRESSO installation with `pw.x`, `bands.x` and
  `projwfc.x` available on `PATH` (or pass explicit executable paths).
- Python dependencies installed for Molecular Lab, including ASE.
- One silicon **PBE** UPF pseudopotential compatible with the installed QE
  version. Keep its original license and source information.
- The exact `ecutwfc` and `ecutrho` recommendations for that pseudopotential.
  Do not copy cutoffs from another UPF, even if the filenames look similar.

## 2. Run the smoke test

From the repository root:

```bash
python examples/qe_silicon_fatbands.py \
  --pseudo /path/to/your/Si-PBE.UPF \
  --pseudo-family "your pseudo family and release" \
  --pseudo-source "URL or citation for the UPF" \
  --ecutwfc-ry YOUR_ECUTWFC_RY \
  --ecutrho-ry YOUR_ECUTRHO_RY
```

The example builds the two-atom primitive diamond cell, copies the selected UPF
into the runner's `pseudo/` directory, records its SHA-256 checksum, and runs:

1. SCF with the supplied Monkhorst-Pack mesh;
2. `pw.x` BANDS along a deliberately short fractional-coordinate
   `Γ-X-L-Γ` smoke-test path;
3. `bands.x` for the conventional band output;
4. `projwfc.x` for atomic/orbital projections;
5. the Molecular Lab fat-band parser, plot-data preparation and SVG renderer.

The short path verifies the integration and coordinate hand-off; it is not a
complete publication-quality silicon band path. The default lattice parameter
(5.431 Å) is a starting value, not a claim of a relaxed or converged structure.

## 3. Required evidence to record

For every real run, retain:

- QE version and executable paths;
- UPF source, license, family/release, filename and SHA-256 checksum;
- lattice parameter and cell vectors;
- `ecutwfc`, `ecutrho`, k-point mesh and band-path coordinates;
- full SCF, BANDS, `bands.x` and `projwfc.x` stdout/stderr;
- whether each stage ended normally and the resulting total energy;
- the generated `atomic_proj.xml`, `.gnu` output and SVG.

The runner currently expects the pseudopotential file under `<workdir>/pseudo`
and the saved calculation under `<workdir>/tmp/molecular_lab.save`.

## 4. Convergence and scientific checks (separate from smoke test)

A successful process exit only verifies that the workflow ran. Before treating
the band structure as a scientific result:

- converge `ecutwfc` and `ecutrho` using the selected UPF's recommendations;
- converge the SCF k-point mesh and check total-energy stability;
- verify the lattice parameter against a separately converged relaxation or a
  cited reference;
- use a complete, crystallographically verified silicon high-symmetry path;
- compare the resulting gap and dispersion with an identified reference method,
  while accounting for the known PBE band-gap underestimation;
- inspect projection weights and confirm that atom/orbital labels and state
  counts match the same saved BANDS calculation.

Do not interpret fat-band weights as unique observables: they depend on the
pseudopotential/projector construction. Do not claim experimental agreement
from this smoke test.

## 5. Current adapter detail

`DFTCalculationSettings.charge_density_cutoff` is optional. When supplied, the
QE adapter writes `ecutrho` in Ry after unit conversion and validates that it
is not below `ecutwfc`. The value is not guessed automatically because the
appropriate ratio depends on the selected pseudopotential family.
