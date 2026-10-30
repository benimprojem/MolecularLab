# ASE Structure Bridge

## Purpose

molecular_lab/calculations/ase_bridge.py provides an optional conversion from the project's backend-neutral CrystalStructure model to ASE's Atoms object. ASE is not a mandatory dependency and is imported only when conversion is requested.

## Coordinate convention

The project stores lattice vectors as three rows in angstrom and atomic positions as fractional coordinates. The bridge passes the row-vector cell and fractional positions to ASE as cell and scaled_positions. Periodic boundary flags are preserved axis by axis.

## Metadata preserved

The bridge places the structure identifier, source, space-group label, magnetic-ordering label, total charge, notes, site identifiers, and oxidation-state annotations in Atoms.info. If any site has a specified magnetic moment, initial magnetic moments are set for all sites; unspecified site moments are represented as zero in that ASE initial-moment array.

## Strict unsupported cases

- Partial occupancy/disorder: rejected. An ordered supercell or explicit disorder model is required; the bridge does not round occupancies or choose species randomly.
- Defect annotations: rejected. A vacancy/substitution/interstitial annotation is not silently assumed to have modified the atom list. A future defect materializer must produce a new explicit structure.
- Missing ASE: reported as ASEUnavailableError; importing the core package does not require ASE.

This bridge performs data conversion only. It does not calculate energy, forces, stress, electronic structure, or stability.

## Next step

Backend adapters can consume the ASE structure alongside DFTCalculationSettings. They must separately validate functional, pseudopotential, relativistic, smearing, and task support before writing solver input.
