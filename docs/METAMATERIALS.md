<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 -->

# Metamaterial Discovery Architecture

## Purpose and scope

Metamaterials are a first-class research domain in Molecular Lab. Their response can be dominated by engineered subwavelength or wavelength-scale structure, not only by the intrinsic chemistry of a constituent material. The application must therefore model material composition and repeating geometry as separate but linkable layers.

This document defines the intended discovery and validation workflow. The initial Python records store a proposed unit-cell geometry, lattice periods, spectral band, operating conditions, and target properties. They do not calculate or validate electromagnetic response.

## Distinguish the physical mechanisms

The discovery engine must not treat every engineered optical material as the same phenomenon:

- **Resonant metamaterials:** patterned inclusions or resonators produce frequency-dependent electric and/or magnetic response. Loss, resonance bandwidth, coupling, and fabrication tolerances matter.
- **Geometric-phase structures:** spatially varying anisotropic elements can control phase, polarization, or wavefront. Orientation and element geometry are essential.
- **Photonic crystals:** periodic dielectric contrast creates band structures and, in suitable cases, photonic band gaps. The lattice period may be comparable to the wavelength in the medium; this is not automatically an effective-medium metamaterial.
- **Effective-medium designs:** a homogenized description may be useful when the unit cell is sufficiently subwavelength and the assumptions of homogenization are satisfied.
- **Hybrid designs:** combine mechanisms and must preserve the separate assumptions and validation methods.

A negative effective refractive index, epsilon, or mu is a possible target in some frequency ranges, not a universal definition or a guaranteed outcome.

## Design representation

Each design should retain:

- Unit-cell geometry family and parameterized dimensions with units.
- Two-dimensional or three-dimensional periodicity and lattice periods.
- Constituent materials and, when available, their frequency-dependent complex permittivity/permeability.
- Operating wavelength or frequency interval.
- Incidence angle, polarization, temperature, surrounding medium, and boundary conditions.
- Requested outputs such as reflection, transmission, absorption, phase shift, polarization conversion, dispersion, band gap, effective impedance, or retrieved effective parameters.
- Geometry/material source, model version, mesh and solver settings, convergence, uncertainty, and evidence category.

The initial model stores lengths in metres for lattice periods and wavelength-band bounds. Unit-cell dimensions use explicit units. Future schema work should standardize all physical quantities and support complex, frequency-dependent material response.

## Forward and inverse discovery

### Forward mode

1. Select constituent materials and a parameterized unit-cell family.
2. Define lattice periodicity and geometry bounds.
3. Generate a bounded set of geometries.
4. Apply geometric validity and manufacturability filters.
5. Simulate selected structures under explicit spectral, polarization, incidence, and boundary conditions.
6. Extract spectra, dispersion/band structures, and field distributions.
7. Where homogenization is justified, retrieve effective parameters and document the retrieval method.
8. Rank candidates against the requested behavior and uncertainty.

### Inverse mode

1. Specify a target response and operating spectral band.
2. Specify polarization, incidence, substrate, environment, thickness, and loss constraints.
3. Search geometry families and material combinations.
4. Use surrogate or reduced-order models only within their validated domain.
5. Verify shortlisted designs with a full-wave solver and convergence checks.
6. Evaluate robustness to dimensional variation, material dispersion, loss, and fabrication tolerances.

Inverse discovery must not assume that a target spectrum uniquely determines a structure. Multiple geometries may be equivalent within tolerances, and some target combinations may be physically incompatible.

## Calculation backends

Potential backends should be modular and selected by task:

- **FDTD:** broadband time-domain response and field evolution, subject to grid and timestep convergence.
- **FEM:** frequency-domain and complex geometries, with mesh convergence and boundary-condition checks.
- **Frequency-domain / RCWA:** periodic layered or patterned structures where the formulation is appropriate.
- **Plane-wave expansion:** photonic band structures for suitable periodic systems.
- **Eigenmode solvers:** resonances and dispersion.
- **Full-wave periodic solvers:** scattering, reflection/transmission, and field maps.

The project must not claim a solver is installed or a simulation completed unless the backend actually ran and its settings/results were captured. Molecular quantum chemistry or ordinary molecular force fields cannot replace full-wave electromagnetic simulation of a patterned structure.

## Validation and reporting

A candidate report should distinguish:

- Geometry proposal only.
- Analytical or reduced-model estimate.
- Surrogate prediction.
- Full-wave numerical result with convergence status.
- Independent solver or experimental validation.

For effective parameters, record the retrieval procedure, reference planes, branch choices, polarization, incidence angle, and frequency range. For photonic crystals, report the band calculation method and dimensionality. For all full-wave results, record mesh/grid resolution, domain size, boundary conditions, material dispersion/loss model, and convergence tests.

A visually plausible unit cell or a single simulation snapshot is not proof of a robust, fabricable metamaterial.

## Current implementation status

Implemented: immutable records for unit-cell dimensions, 2D/3D lattice periods, mechanism classification, wavelength bands, target-property labels, operating conditions, evidence, and limitations.

Not implemented: geometry generators, material dispersion database, electromagnetic solver adapters, homogenization/retrieval, inverse optimizer, fabrication constraints, and experimental-data ingestion.
