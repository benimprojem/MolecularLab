<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 04 -->

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

- **DensityOfStatesChannel** stores non-negative DOS values in states/eV for one spin channel.
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

Implemented: standard-library records, validation, and tests for k-points, spin-resolved bands, occupation arrays, direct/indirect gap screening, and total/projected DOS grids.

Partially implemented: Quantum ESPRESSO input preparation for a structured BANDS path and total-DOS NSCF/dos.x workflow. Not implemented: BANDS/DOS process orchestration and output parsers for QE/GPAW, construction of normalized result records from solver output, band interpolation, projected-DOS adapters, effective masses, Fermi-surface analysis, optical transitions, or actual DFT execution.


## Quantum ESPRESSO output ingestion

The QE adapter's parsed BANDS result currently uses the .gnu output produced by bands.x. Its first column is cumulative path distance, not a Cartesian or fractional reciprocal vector; the common KPoint record therefore includes a coordinate-system label. The energies are stored in eV.

The total-DOS parser reads the energy and DOS columns in eV and states/eV. Spin-polarized tables are represented as separate up/down channels. Integrated-DOS columns are ignored because the current common DOS model does not yet carry integrated counts. The parser attaches a Fermi level only when explicitly reported; it does not synthesize a mid-gap reference.
