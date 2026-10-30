# Fat-Band Data Model

## Purpose

The fat-band layer associates each sampled band eigenstate with atomic
wavefunction projection weights. It is distinct from projected density of
states (PDOS): PDOS is energy-resolved and integrated/broadened over k-points,
whereas fat bands retain the k-point and band indices.

## Quantum ESPRESSO sources

- `atomic_proj.xml` supplies complex projection coefficients produced by
  `projwfc.x`.
- `projwfc.x` standard output supplies the mapping from each atomic
  wavefunction index to its atom, element and angular-momentum metadata.
- `BandStructureData` supplies band energies and the k-point axis from the
  matching `pw.x`/`bands.x` workflow.

The parser supports the newer `EIGENSTATES/PROJS/ATOMIC_WFC` layout and the
legacy `PROJECTIONS/K-POINT.N/(SPIN.N)/ATMWFC.N` layout. Header dimensions
may be XML attributes or child elements. Four spin components (non-collinear
or spinor data) are rejected explicitly until a spinor-aware model exists.

## Projection mathematics

For a complex projection coefficient

\[
c_{n\mathbf{k},a} = \operatorname{Re}(c) + i\operatorname{Im}(c),
\]

the raw projector weight is

\[
w_{n\mathbf{k},a} = |c_{n\mathbf{k},a}|^2
= \operatorname{Re}(c)^2 + \operatorname{Im}(c)^2.
\]

Weights belonging to the same atom and angular-momentum projector are summed
across matching pseudopotential atomic wavefunctions. The parsed records retain
the original real and imaginary parts.

## Index and unit conventions

- k-point indices in Molecular Lab are zero-based.
- Band indices and QE atomic-wavefunction indices are one-based.
- Band energies come from `BandStructureData` and remain in eV.
- Projection weights are dimensionless.
- Atom indices follow QE's one-based atom ordering.
- Spin channels are `total`, `spin_up`, or `spin_down` for the currently
  supported scalar/collinear cases.

## Validation and scientific limits

The builder validates the number of bands, k-points and spin channels. It joins
by band index, k-point order and spin label; it does not match by approximate
energy. The caller must ensure that `atomic_proj.xml` and the supplied band
structure originate from the same saved `pw.x` calculation and exactly the
same k-point sequence. The projection XML coefficients alone do not provide
the plotted band-path distance.

Projection weights are not normalized automatically and their sum is not
required to equal one. Completeness depends on the pseudopotential's atomic
projector basis and the projection procedure. The stored `raw_weight` values
preserve the QE-derived quantities; visualization may derive a separate
normalized scale without overwriting them.

This first implementation supports collinear scalar-relativistic projections.
Fully relativistic j-resolved and non-collinear/spinor formats, band
interpolation, and automatic verification of k-point coordinate identity
remain unsupported.


## Execution workflow

A BANDS request opts into fat-band output by including
`DFTOutput.FAT_BAND_STRUCTURE` in `requested_outputs`. The QE runner then
executes:

1. SCF with `pw.x`.
2. The structured high-symmetry path with `pw.x` using `calculation='bands'`.
3. `bands.x` to produce the conventional `.gnu` band representation.
4. `projwfc.x` on the same saved BANDS wavefunctions.
5. Read `tmp/molecular_lab.save/atomic_proj.xml` and the state metadata/eigenvalue
   records printed by `projwfc.x`.
6. Build the `FatBandDataset` from the original QE eigenstate sequence.

The standard `band_structure` result is still populated from `bands.x`.
The fat-band energies are deliberately taken from the `projwfc.x` eigenstate
output, not copied from `bands.x` reordered curves. QE's `bands.x` can reorder
states to improve band continuity; raw projection coefficients are indexed by
the original `pw.x` state ordering. Fractional path coordinates are reconstructed
from the exact `BandPath` supplied to `pw.x`; the parser does not assume a unit
convention for the k-vector text printed by `projwfc.x`. Pairing those arrays solely by a reordered
band number can assign an orbital weight to the wrong energy at a crossing.
The two representations are therefore kept distinct until a state-tracking
layer can establish a verified mapping.

The current runner supports non-spin-polarized scalar projections only.
Spin-polarized BANDS, SOC/spinor data and fully relativistic j-resolved
projections remain unavailable. The automated workflow tests use synthetic XML
and mocked executables; a real QE reference calculation is still required.

## Projection output source

QE documents `atomic_proj.xml` as the per-state projection file written by
`projwfc.x`. The `lwrite_overlaps` option adds the atomic-orbital overlap
matrix; it is not required to enable the projection coefficients.
The `kresolveddos` option produces energy-resolved k-point DOS output and is
not required for the per-eigenstate fat-band data model.


## Plot-ready transformation

\`prepare_fatband_plot_data\` converts a \`FatBandDataset\` into grouped
\`FatBandPlotSeries\` records without depending on a plotting library.

Supported grouping levels:

- \`atom_orbital\`: one atom and orbital component per series.
- \`element_orbital\`: sum the same element/orbital over atoms.
- \`element\`: sum all selected orbital contributions for each element.
- \`atom\`: sum all selected orbital contributions for each atom.
- \`orbital\`: sum the selected orbital contribution over atoms.

The transformation supports filters for elements, atom indices, orbital labels,
spin channels and band indices. A minimum raw-weight threshold can suppress
small markers while preserving the original value in the plot record.

Display scaling is explicit:

- \`raw\`: display the unscaled projection weight.
- \`global_max\`: divide all selected display weights by the global maximum.
- \`per_band_max\`: scale each band/spin group independently.

The source dataset is never modified. Each plot point retains both
\`raw_weight\` and \`display_weight\`.

The x-axis is a cumulative reciprocal-space path distance in Å⁻¹. It is
calculated from the fractional k-points, the direct lattice vectors and the
reciprocal lattice (including the 2π convention). High-symmetry labels and
their exact path positions are returned as symmetry ticks. The supplied
\`BandPath\` must match the fat-band k-point sequence exactly.
