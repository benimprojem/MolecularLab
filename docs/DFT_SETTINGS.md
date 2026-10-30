# DFT Calculation Settings

## Purpose

This module defines backend-neutral input settings for periodic density-functional-theory (DFT) calculations. It is a data contract, not a DFT engine. It does not execute calculations, choose scientifically adequate parameters, download pseudopotentials, or claim that a selected backend supports a requested option.

Implementation: `molecular_lab/calculations/dft_settings.py`.

## Settings represented

- **Task:** SCF, fixed-cell relaxation, variable-cell relaxation, band calculation, or DOS task.
- **Exchange-correlation:** explicit functional identifier (for example PBE, PBEsol, SCAN, r2SCAN, or HSE06). The identifier is not a backend support guarantee.
- **Plane-wave energy cutoff:** positive finite value with eV, Ry, or Ha unit and explicit conversion helpers.
- **Brillouin-zone sampling:** three-dimensional integer mesh, Monkhorst-Pack or Gamma-centered scheme, and explicit shift.
- **Spin:** collinear spin-polarized flag and a separate spin-orbit-coupling flag. SOC does not imply a magnetic state; non-collinear/spinor details remain backend-adapter concerns.
- **Smearing:** optional explicit scheme and energy width. No automatic smearing choice is made.
- **Convergence:** energy tolerance, optional force, density and stress tolerances, and maximum SCF iterations. Stress tolerance is expressed in GPa and is required for the common VC_RELAX execution path.
- **Pseudopotential references:** element, file name, family, source, optional valence electron count, checksum, and license metadata.
- **Outputs:** total energy, forces, stress, band structure, DOS, charge density, and magnetic moments.
- **Conditions:** optional temperature and pressure metadata.
- **Band path:** a separate explicit reference for band-path tasks; a uniform mesh is not treated as a high-symmetry path.

## Validation

- Values and units are validated; cutoffs convert between eV, Ry, and Ha.
- K-point grids require three positive integer counts.
- A structure-aware check requires one k-point along each non-periodic axis and checks exact pseudopotential element coverage.
- Pseudopotential references are immutable as a mapping and must match their element keys.
- Band tasks require an explicit band-path reference.
- SOC is represented independently from magnetic spin polarization.

## Scientific boundaries

Cutoff and k-point convergence are material-, property-, pseudopotential-, and method-dependent. The model intentionally supplies no universal recommended cutoff or mesh. A backend adapter must validate its own supported functionals, pseudopotential formats/families, smearing methods, relativistic treatment, and task combinations.

The common model does not yet represent a complete high-symmetry path, projected-orbital selection, Hubbard U, hybrid-functional exact-exchange controls, van der Waals correction parameters, electric-field boundary conditions, or finite-size corrections for charged defects. These should be added as explicit typed settings when a concrete workflow requires them, rather than as undocumented free-form strings.

A successful SCF convergence is not proof of structural, thermodynamic, or kinetic stability. A band gap is not by itself proof of measured conductivity or device suitability.

## Next integration

Implement backend-specific capability declarations and input translators (initial candidates: Quantum ESPRESSO and GPAW). The adapter must preserve units and provenance and return an explicit unsupported/unavailable status instead of silently substituting methods. Actual execution should only be enabled after dependency and pseudopotential availability checks.


## Variable-cell convergence

`ConvergenceSettings.stress_tolerance_gpa` is an explicit positive tolerance for variable-cell optimization. It is intentionally optional for other tasks and has no silently chosen default. The common VC_RELAX runner also requires an explicit atomic force tolerance. Backend adapters may implement different internal convergence algorithms; the common result reports the final stress tensor separately so that backend-specific stopping behavior is not hidden.


## Band-path and DOS request settings

- `BandPathPoint` stores fractional reciprocal coordinates, a display label, and the number of subdivisions to the next point. Every non-final point needs a positive subdivision count; the final point uses zero.
- `BandPath` is an ordered sequence of at least two vertices. Only fractional crystal reciprocal coordinates are currently accepted.
- `DOSSettings` requires an explicit minimum/maximum energy in eV, at least two grid points, and a positive broadening in eV.
- A BANDS task requires either a structured path or a reference identifier at the settings level. QE input generation specifically requires the structured path; a reference alone is not enough to emit k-points.
- A DOS task requires explicit `DOSSettings`. The QE input adapter maps the broadening width to Ry for `dos.x`.
