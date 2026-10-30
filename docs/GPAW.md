# GPAW Adapter

## Purpose

The GPAW adapter prepares an ASE `Atoms` object and a GPAW calculator configuration. It does not run SCF, relax a structure, or execute a solver job.

Implementation: `molecular_lab/calculations/gpaw.py`.

## Supported initial scope

- SCF, fixed-cell ionic relaxation, and fully periodic variable-cell relaxation preparation
- GPAW plane-wave mode, with the shared energy cutoff converted to eV
- Gamma-centered or Monkhorst-Pack k-point mesh selection, with zero explicit fractional shifts
- Exchange-correlation identifier, total charge, spin-polarized flag, maximum SCF iterations, and convergence dictionary
- Fermi-Dirac smearing
- ASE `FrechetCellFilter` target preparation for variable-cell relaxation
- Optional scalar pressure for VC-RELAX, converted from GPa to eV/Å³

## GPAW-specific PAW setups

GPAW uses PAW setups supplied by the local GPAW installation. It does not consume the QE-style external pseudopotential filenames in `DFTCalculationSettings.pseudopotentials`. The adapter deliberately does not pretend those files are GPAW setup identifiers. Setup selection/availability is left to GPAW's installed setup library and is not verified here.

## Explicitly unsupported

- BANDS and DOS as primary tasks; these need separate post-SCF workflows
- Spin-orbit coupling in this initial adapter
- Site-resolved initial magnetic moments
- Non-zero explicit k-point shifts
- Gaussian, Methfessel-Paxton, and cold smearing (Fermi-Dirac only)
- Partial periodicity in the plane-wave mode
- Charge-density, band-structure, DOS, and magnetic-moment result extraction

Unsupported requests raise `GPAWInputError`.

## Units and interpretation

- Plane-wave cutoff is passed to GPAW in eV.
- K-point mesh is passed as ASE/GPAW `{'size': (...), 'gamma': bool}`.
- Smearing width remains in eV.
- The shared energy tolerance is passed to GPAW `convergence['energy']`. GPAW documents this criterion per valence electron, so it is not numerically identical to every backend's similarly named threshold.
- Force tolerance is supplied both as GPAW's SCF force convergence criterion and as the separate ASE optimizer `fmax` suggestion. The caller must explicitly configure and run the optimizer.
- Pressure for the ASE cell filter is converted to eV/Å³ (1 GPa ≈ 0.006241509074 eV/Å³).

## Usage

```python
from ase.optimize import BFGS
from molecular_lab.calculations.gpaw import (
    prepare_gpaw_input, create_gpaw_calculator, build_optimization_target,
)

prepared = prepare_gpaw_input(crystal, settings)
calculator = create_gpaw_calculator(prepared, txt="gpaw.log")
prepared.atoms.calc = calculator

if prepared.task.value == "scf":
    energy = prepared.atoms.get_potential_energy()
else:
    target = build_optimization_target(prepared)
    optimizer = BFGS(target, trajectory="relax.traj")
    optimizer.run(fmax=prepared.optimizer_fmax_ev_per_angstrom or 0.05)
```

The example is an execution recipe; the adapter functions themselves do not execute it. Variable-cell relaxation uses ASE's `FrechetCellFilter` and requires GPAW stress support, for which plane-wave mode is selected.

## Scientific boundary

Preparing a calculator is not a calculation result. The user must verify GPAW and ASE versions, PAW setup availability, XC support, cutoff/k-point/SCF convergence, and physical suitability of the model. No experimental property or stability claim follows from this adapter.
