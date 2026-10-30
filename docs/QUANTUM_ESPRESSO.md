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

- BANDS and DOS tasks: these require separate band-path and multi-stage SCF/NSCF/DOS workflows
- Spin-orbit coupling: the pseudopotential model does not yet state whether a file is fully relativistic
- Site-resolved initial magnetic moments: QE species-level magnetic initialization needs a deliberate species/sublattice mapping
- Gamma-centered meshes: only Monkhorst-Pack is currently translated
- Partial occupancy and defect annotations: the ASE bridge requires an explicit ordered, fully occupied structure
- Non-3D VC-RELAX cell constraints
- Charge-density and magnetic-moment output extraction

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
