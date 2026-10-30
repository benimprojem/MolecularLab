<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 05 -->

# Electronic-Structure Result Model

## Scope

This module defines normalized records for periodic electronic-structure calculations and a conservative band-edge screening helper. It is a data/analysis layer only: no DFT calculation is performed by constructing these objects.

## Band structure

- **KPoint** stores three fractional reciprocal-space coordinates and a positive integration weight.
- **BandChannelData** stores band energies in eV as a bands-by-k-points array for one spin channel. Optional occupations retain the backend occupancy scale through maximum_occupation (commonly 0–2 for spin-degenerate data or 0–1 for a single spin channel).
- **BandStructureData** groups k-points and spin channels and can retain Fermi level, occupation tolerance, method, backend/version, exchange-correlation approximation, pseudopotential set, and notes.
- **BandGapAssessment** reports gap, metallic, or undetermined; indirect and direct gaps when available; VBM/CBM energies and sampled k-point indices.

The assess_band_gap helper uses occupation information when present. Partial occupations are not automatically labelled metallic because finite-temperature or smearing schemes can create fractional values in otherwise gapped systems. A partial state at the supplied Fermi level is labelled metallic; otherwise the result remains undetermined until smearing/temperature context is available. If occupations are absent, an explicit Fermi level is required to separate occupied and unoccupied states. If the global valence/conduction edges overlap within the configured tolerance, the result is metallic. Missing evidence yields undetermined.

The direct gap is the minimum occupied-to-unoccupied separation at the same sampled k-point; the indirect gap is global CBM minus global VBM. These are estimates over the supplied k-point set, not guarantees that the true extrema were sampled. A production adapter must record k-point convergence and whether supplied points form a path, integration mesh, or combined dataset. Spin-channel edge extrema are combined for the global indirect gap; the reported direct gap is evaluated within each spin channel at a common sampled k-point. This convention should be revisited for non-collinear or spin-mixing calculations.

## Density of states

- **DensityOfStatesChannel** stores non-negative total DOS values in states/eV for one spin channel.
- **ProjectedDOSRecord** stores atom index, element, pseudopotential wavefunction index, angular momentum, spin channel, local orbital DOS and (2l+1) angular components.
- **ProjectedDensityOfStatesData** groups atom/orbital records on a shared energy grid, with Fermi level and backend provenance.
- **DensityOfStatesData** stores a strictly increasing energy grid, spin channels, optional Fermi level, projection label, normalization description, and solver metadata.

The normalization must state whether DOS is per cell, atom, formula unit, or another convention. Projected DOS must preserve the projection definition (atom, orbital, spin, or other projector) in metadata/notes. The initial model does not integrate DOS, infer carrier density, or estimate mobility.

## Scientific limitations

- DFT band gaps depend on functional, pseudopotential, basis/cutoff, k-point sampling, smearing, spin treatment, and structural model. Standard semilocal DFT commonly underestimates many semiconductor gaps; the app must report the method rather than present a computed value as an experimental gap.
- A finite k-point path can miss a band extremum. Direct/indirect classification is only as reliable as the underlying sampling and band interpolation.
- A Fermi level alone is a weaker fallback than explicit occupations, particularly for smeared metallic systems. Adapters should supply occupations whenever available.
- SOC, non-collinear magnetism, hybrid functionals, GW, excitonic effects, finite temperature, disorder, defects, and strong correlation require explicit method support and cannot be assumed.
- DOS at the Fermi level can be useful evidence but must not silently override a band-edge assessment; transport classification is separate.
- A band gap does not by itself establish conductivity, mobility, semiconductor quality, or device behavior.

## Current status

Implemented: standard-library records, validation, and tests for k-points, spin-resolved bands, occupation arrays, direct/indirect gap screening, total DOS, and atom/orbital-resolved PDOS.

Partially implemented: QE BANDS, total DOS and collinear scalar-relativistic atom/orbital PDOS execution and output parsing. Not implemented: GPAW BANDS/DOS, spin-polarized BANDS, fully relativistic j-resolved and k-resolved PDOS, band interpolation, effective masses, Fermi-surface analysis, optical transitions, or actual experimental validation.


## Quantum ESPRESSO output ingestion

The QE adapter's parsed BANDS result currently uses the .gnu output produced by bands.x. Its first column is cumulative path distance, not a Cartesian or fractional reciprocal vector; the common KPoint record therefore includes a coordinate-system label. The energies are stored in eV.

The total-DOS parser reads the energy and DOS columns in eV and states/eV. Spin-polarized tables are represented as separate up/down channels. Integrated-DOS columns are ignored because the current common DOS model does not yet carry integrated counts. The parser attaches a Fermi level only when explicitly reported; it does not synthesize a mid-gap reference.


## Projected density of states

The QE PDOS adapter uses `projwfc.x` output files. Each record preserves the one-based atom index, element, pseudopotential wavefunction index, angular momentum (l), spin channel, local DOS and the (2l+1) angular components. Collinear spin-polarized output keeps spin-up and spin-down separate. QE's p and d component ordering is retained as documented; f components are kept as generic angular labels.

The decomposition depends on the pseudopotential atomic wavefunctions and projection procedure. It is not a unique observable assignment of electronic charge to atoms or orbitals. j-resolved relativistic and k-resolved output are not yet supported. See [QE projwfc.x documentation](https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html).\n\nThe analysis helpers in `molecular_lab.calculations.pdos_analysis` sum projected components by element or atom while retaining angular momentum and spin channels. They also integrate a supplied PDOS curve over an energy window using trapezoidal integration with interpolated boundaries. These are numerical summaries of the supplied projections, not independent physical observables. See [PDOS Aggregation and Energy-Window Analysis](PDOS_ANALYSIS.md).
