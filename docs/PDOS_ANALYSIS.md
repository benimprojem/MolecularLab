# PDOS Aggregation and Energy-Window Analysis

## Purpose

The aggregation layer converts atom/wavefunction-resolved Quantum ESPRESSO PDOS records into convenient element-level or atom-level orbital sums. It does not run a solver and does not alter the source PDOS records.

## Aggregation

`aggregate_projected_dos(data, group_by="element")` groups by element, angular momentum (s/p/d/f as l=0/1/2/3), and spin channel. It sums the local wavefunction DOS and the angular-component DOS independently. `group_by="atom"` retains individual atom indices while still summing wavefunctions with the same l and spin.

Spin-up and spin-down are deliberately not merged. Consumers can form a total spin-summed curve explicitly when scientifically appropriate. The angular-component sum is the primary orbital-resolved sum; the local DOS is retained separately as a cross-check against the component sum. Small numerical differences can occur due to output formatting or parser conventions.

## Energy-window integration

`integrate_pdos_window(energies_ev, density_states_per_ev, lower_ev, upper_ev)` applies trapezoidal integration, linearly interpolating the window boundaries. The result is an integrated projected-state count under the source curve's normalization. For example, integrate from the Fermi level to a selected upper energy to inspect unoccupied projected weight, provided the energy grid covers that range.

## Interpretation limits

- A projected sum is only as complete as the projectors and files supplied by the calculation.
- Atomic-orbital projections depend on pseudopotential wavefunctions and are not a unique observable partition.
- Do not assume projected DOS must exactly equal total DOS. Differences may reflect incomplete projectors, normalization, omitted channels, or parser/backend conventions.
- Integrated PDOS is not automatically a carrier concentration, charge transfer, conductivity, or magnetic moment.
- Compare spin channels only when their normalization and calculation settings are consistent.


## Integrated orbital contributions

`summarize_pdos_window(data, lower_ev, upper_ev, group_by="element")` first aggregates records by atom or element, angular momentum and spin, then integrates the local-DOS and angular-component curves over the requested energy interval. It returns one `PDOSWindowSummary` per element/atom, orbital angular momentum and spin channel. This makes s/p/d/f contributions and spin-up/down contributions directly comparable without merging distinct channels.

The interval must lie inside the sampled energy grid; boundaries are linearly interpolated. The local-DOS and angular-component integrals are reported separately as a consistency check. They are not automatically equal to a total-DOS integral, and their difference must not be interpreted as a measured charge transfer.
