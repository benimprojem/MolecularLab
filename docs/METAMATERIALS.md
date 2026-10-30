<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 02 -->

# Metamaterial Discovery Architecture

## Purpose and scope

Metamaterials are a first-class research domain in Molecular Lab. Their response can be dominated by engineered subwavelength or wavelength-scale structure, not only by the intrinsic chemistry of a constituent material. The application therefore models material composition and repeating geometry as separate but linkable layers.

The project distinguishes a proposed design, a numerical result, and experimental validation. A design record is not a calculation; a computed local minimum or a single simulation is not proof of a stable, robust, or fabricable material.

## Multiphysics property domains

The metamaterial engine is intended to support coupled investigation across:

- **Optical:** complex refractive index, permittivity, reflectance, transmittance, absorption, phase delay, polarization conversion, optical activity, dispersion, and optical loss.
- **Magnetic:** complex permeability where applicable, susceptibility, anisotropy, saturation magnetization, coercivity, remanence, magnetic resonance, and frequency-dependent magnetic loss.
- **Electrical transport:** conductivity and resistivity, anisotropy, dielectric response, permittivity, impedance, charge transport, and breakdown strength.
- **Electronic structure:** band structure, band gap, density of states, Fermi level, carrier concentration, effective masses, and mobility where the selected method supports them.
- **Conductor / semiconductor / insulator behavior:** classifications derived from explicit electronic and transport evidence. No universal band-gap or conductivity threshold is applied across all materials, temperatures, frequencies, and measurement geometries.

DC conductivity, AC conductivity, dielectric loss, and optical absorption are related but not interchangeable. Each reported value must retain its unit, conditions, spectral range, method, uncertainty, and evidence category.

## Electrical behavior classification

The initial classification helper is deliberately criterion-driven. A caller must provide two cutoffs and identify their source:

- For a band-gap criterion (eV), values at/below the lower cutoff are labelled conductor-like, values at/above the upper cutoff insulator-like, and intermediate values semiconductor-like.
- For a conductivity criterion (S/m), values at/below the lower cutoff are labelled insulator-like, values at/above the upper cutoff conductor-like, and intermediate values semiconductor-like.

Every assessment records the selected metric, thresholds, criterion source, temperature, and frequency (0 Hz denotes DC). The label is a screening classification under that criterion, not a universal statement about the material. It does not infer mobility, carrier type, anisotropy, contacts, or transport mechanism. Conflicting evidence from band structure and transport must be retained and reviewed rather than silently forced into one class.

## Distinguish the physical mechanisms

- **Resonant metamaterials:** patterned inclusions or resonators produce frequency-dependent electric and/or magnetic response. Loss, bandwidth, coupling, and fabrication tolerances matter.
- **Geometric-phase structures:** spatially varying anisotropic elements control phase, polarization, or wavefront.
- **Photonic crystals:** periodic dielectric contrast creates band structures and potentially photonic band gaps. This is not automatically an effective-medium metamaterial.
- **Effective-medium designs:** homogenized descriptions are useful only when the unit cell is sufficiently subwavelength and homogenization assumptions hold.
- **Hybrid designs:** combine mechanisms and preserve their separate assumptions and validation methods.

Negative effective refractive index, permittivity, or permeability is a possible frequency-bounded target, not a universal definition or guaranteed outcome.

## Design and target representation

A design stores its geometry family, parameterized dimensions with units, 2D/3D periodicity, lattice periods, constituent materials, wavelength band, operating conditions, target properties, evidence, and limitations.

Targets are distinct from results. Each typed target records its domain, metric, relation (minimize, maximize, approach, or within range), unit, optional numerical bounds, and optional spectral interval.

## Computed response data model

The initial result schema stores frequency-dependent response without pretending to run a solver:

- **ComplexValue** stores real and imaginary components separately. The solver's time-harmonic sign convention must be recorded in provenance.
- **ResponseSample** stores one response tensor at a positive frequency, with optional free-space wavelength, temperature, propagation direction, polarization, and uncertainty.
- A 1×1 tensor represents a scalar/isotropic response; 2×2 and 3×3 tensors represent in-plane or full 3D response. Tensor components are not silently reduced to a scalar.
- **MaterialResponseDataset** groups samples for one domain and metric, sorts them by frequency, rejects duplicate frequency samples, and retains units, evidence, and limitations.
- **SolverProvenance** records method, backend, backend version, settings, and source. Settings are immutable after construction.
- **ConvergenceInfo** records whether convergence was checked, passed, failed, or remains unknown, together with error estimate, iteration count, and discretization notes.

Uncertainty values currently represent non-negative per-component magnitudes; the schema does not prescribe whether they are standard deviations, confidence bounds, or numerical error bars. That interpretation must be stated in dataset provenance or limitations until a richer uncertainty model is added. Wavelength and frequency are both retained when supplied, but their consistency is not inferred because wavelength may refer to vacuum or a medium-dependent wavelength.

## Forward and inverse discovery

### Forward mode

1. Select constituent materials and a parameterized unit-cell family.
2. Define lattice periodicity and geometry bounds.
3. Generate a bounded set of geometries and apply validity/manufacturability filters.
4. Simulate selected structures under explicit spectral, polarization, incidence, and boundary conditions.
5. Extract spectra, dispersion/band structures, and field distributions.
6. Retrieve effective parameters only when homogenization is justified, recording the retrieval method.
7. Rank candidates against targets and uncertainty.

### Inverse mode

1. Specify target response and operating spectral band.
2. Specify polarization, incidence, substrate, environment, thickness, and loss constraints.
3. Search geometry families and material combinations.
4. Use surrogate/reduced models only within their validated domain.
5. Verify shortlisted designs with a full-wave solver and convergence checks.
6. Evaluate robustness to dimensional variation, material dispersion, loss, and fabrication tolerances.

Inverse discovery must not assume a target spectrum uniquely determines a structure. Multiple geometries may be equivalent within tolerances, and some target combinations may be physically incompatible.

## Calculation backends

Potential modular backends include:

- **FDTD:** broadband time-domain response and field evolution, subject to grid and timestep convergence.
- **FEM:** frequency-domain and complex geometries, with mesh convergence and boundary-condition checks.
- **Frequency-domain / RCWA:** periodic layered or patterned structures where appropriate.
- **Plane-wave expansion:** photonic band structures for suitable periodic systems.
- **Eigenmode solvers:** resonances and dispersion.
- **Full-wave periodic solvers:** scattering, reflection/transmission, and field maps.
- **Electronic-structure and transport solvers:** DFT/band structure and suitable transport methods; these complement rather than replace electromagnetic full-wave solvers.

The project must not claim a solver is installed or a simulation completed unless the backend actually ran and its settings/results were captured. Molecular quantum chemistry or ordinary molecular force fields cannot replace full-wave electromagnetic simulation of a patterned structure.

## Validation and reporting

Candidate reports distinguish geometry proposal, analytical/reduced estimate, surrogate prediction, converged full-wave numerical result, independent solver validation, and experimental validation.

For effective parameters, record retrieval procedure, reference planes, branch choices, polarization, incidence angle, and frequency range. For photonic crystals, report band method and dimensionality. For full-wave results, record mesh/grid resolution, domain size, boundary conditions, material dispersion/loss model, and convergence tests.

## Current implementation status

Implemented: immutable design/unit-cell/typed-target records; frequency-dependent material-response schema with complex scalar/tensor values, temperature and propagation metadata, uncertainty magnitudes, convergence status, solver provenance, and evidence records; and explicit-criterion electrical behavior screening using band gap or conductivity.

Not implemented: actual electromagnetic/electronic/transport/magnetic solver adapters, geometry generators, material dispersion database, homogenization/retrieval, inverse optimizer, fabrication constraints, experimental-data ingestion, cross-metric conflict resolution, or a standardized statistical uncertainty interpretation.
