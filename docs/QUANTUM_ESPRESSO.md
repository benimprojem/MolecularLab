# Quantum ESPRESSO Input Adapter

## Purpose

The first Quantum ESPRESSO adapter generates a `pw.x` input deck. It does not launch Quantum ESPRESSO, download pseudopotentials, inspect a QE installation, or execute a calculation.

Implementation: `molecular_lab/calculations/quantum_espresso.py`.

## Supported scope

- SCF (`calculation='scf'`)
- Fixed-cell ionic relaxation (`calculation='relax'`, BFGS)
- Variable-cell relaxation (`calculation='vc-relax'`, BFGS), only for fully three-dimensional periodic structures
- Plane-wave cutoff conversion to Ry
- Monkhorst-Pack automatic k-point mesh with QE-compatible integer shifts (0 or 1)
- Explicit exchange-correlation identifier, pseudopotential filenames, total charge, collinear spin-polarized mode, and supported smearing options
- Atomic species masses are taken from ASE; no atomic mass table is duplicated in this adapter

## Explicitly unsupported

- BANDS and DOS are now supported at the input-preparation level: BANDS uses an explicit reciprocal-space path and `bands.x`; DOS prepares an NSCF `pw.x` deck and a `dos.x` post-processing input. Execution orchestration and output parsing are not yet implemented.
- Spin-orbit coupling: the pseudopotential model does not yet state whether a file is fully relativistic
- Site-resolved initial magnetic moments: QE species-level magnetic initialization needs a deliberate species/sublattice mapping
- Gamma-centered meshes: only Monkhorst-Pack is currently translated
- Partial occupancy and defect annotations: the ASE bridge requires an explicit ordered, fully occupied structure
- Non-3D VC-RELAX cell constraints
- Charge-density and magnetic-moment output extraction
- Projected DOS and other post-SCF workflows

Unsupported requests raise `QuantumEspressoInputError`; they are not silently omitted.

## Units and mappings

- `ecutwfc`, `conv_thr`, `degauss`, and `forc_conv_thr` are emitted in QE's Rydberg-based units.
- Force tolerance is converted from eV/Å to Ry/Bohr.
- Pressure is converted from GPa to kbar (1 GPa = 10 kbar).
- The backend-neutral energy tolerance maps to QE's electronic SCF `conv_thr`. QE's separate ionic `etot_conv_thr` is not modeled and remains at the QE default.
- Optional temperature is metadata only; it does not create a finite-temperature DFT calculation.
- Spin-polarized mode emits `nspin=2`; no starting magnetization is guessed.

## File and runtime boundaries

The result reports the pseudopotential filenames expected in `pseudo_dir`. Their presence, contents, family compatibility, valence configuration, and license are not verified. The XC identifier is passed through, but support in the user's QE build is not checked. ASE is required when preparing the deck because QE's `ATOMIC_SPECIES` card requires atomic masses.

## Example

```python
from molecular_lab.calculations.quantum_espresso import prepare_pw_input

prepared = prepare_pw_input(crystal, settings, pseudo_dir="./pseudo")
print(prepared.input_text)
print(prepared.required_pseudopotentials)
```

## Scientific boundary

A generated input file is only a calculation request. It is not evidence that the selected pseudopotentials, functional, cutoff, k-point mesh, magnetic state, or convergence settings are scientifically adequate. Convergence testing and independent validation remain necessary.


## Band structure and DOS input preparation

For BANDS, supply a structured `BandPath` with fractional reciprocal coordinates, readable point labels, and a subdivision count for each segment. QE input generation emits `calculation='bands'` with `K_POINTS crystal_b`; a separate `bands.x` namelist can be prepared for the same prefix/outdir. A legacy `band_path_ref` is metadata only and is not enough to generate a QE path.

For total DOS, supply `DOSSettings` with an explicit energy window, grid size, and broadening in eV. The QE `pw.x` input is prepared as an NSCF stage using the configured uniform k-point mesh. A separate `dos.x` input is generated with energy bounds and grid spacing in eV; broadening is converted to Ry as required by QE. The caller must first complete the prerequisite SCF calculation and then the NSCF calculation with a consistent prefix/outdir.

These functions prepare input text only. They do not launch `pw.x`, `bands.x`, or `dos.x`, parse band/DOS output files, or assert convergence.


## Prepared BANDS workflow

Provide `DFTTask.BANDS`, `DFTOutput.BAND_STRUCTURE`, and a structured `BandPath`. Each non-final vertex contains its reciprocal fractional coordinates, a label, and a positive `points_to_next` subdivision count; the final vertex uses zero. The adapter emits a `pw.x` deck using `calculation='bands'` and `K_POINTS crystal_b`. After that run completes, `prepare_qe_bands_postprocess_input` emits a `bands.x` namelist. QE documents `bands.x` as the utility that regroups/reorders the computed eigenvalues and writes a band file. [QE bands.x input reference](https://www.quantum-espresso.org/Doc/INPUT_BANDS.html).

A named `band_path_ref` without actual vertices is not enough to generate an input deck. The application must supply or construct a real path for the chosen crystal; no universal high-symmetry path is guessed.

## Prepared total-DOS workflow

Provide `DFTTask.DOS`, `DFTOutput.DENSITY_OF_STATES`, and explicit `DOSSettings`. `prepare_pw_input` emits an NSCF deck using the configured uniform k-point mesh. After a compatible SCF and NSCF run in the same prefix/outdir, `prepare_qe_dos_postprocess_input` emits a `dos.x` namelist. The requested broadening is converted from eV to Ry; the energy window and grid spacing remain in eV, matching QE's input definition. [QE dos.x input reference](https://www.quantum-espresso.org/Doc/INPUT_DOS.html).

The current DOS path is total DOS only. Projected DOS through `projwfc.x` is not yet modeled. These helpers generate input text only: no process launch, multi-stage dependency handling, output parsing, or normalized `BandStructureData`/`DensityOfStatesData` construction is implemented yet.


## Sequential execution

The common execution entry point supports:

- BANDS: prerequisite SCF, pw.x BANDS path, bands.x, parse filband.gnu.
- Total DOS: prerequisite SCF, uniform-grid NSCF, dos.x, parse total DOS table.

Each stage uses the same prefix and outdir. A later stage is not launched if an earlier process fails. A successful process without the expected result file is not considered a successful workflow.

## Parsed output conventions

- bands.x writes a .gnu file in eV, with one band per blank-separated block. MolecularLab stores its horizontal coordinate as KPoint.coordinate_system=path_distance; it is not represented as a 3D reciprocal coordinate.
- dos.x total DOS is stored as energies in eV and DOS in states/eV. For spin-polarized DOS, up/down channels are kept separate. Additional integrated-DOS columns are not currently stored.
- A Fermi level is stored only when explicitly found in the QE output. It is not guessed from the middle of a gap.
- The current band parser does not receive occupation numbers from .gnu; consequently, band-gap assessment may remain undetermined without an explicit Fermi level or a richer occupation-bearing output source.

The automated workflow tests mock external executables and result files. They verify orchestration and parsing contracts, not the physical correctness of a QE calculation.
