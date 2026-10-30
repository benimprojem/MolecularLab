<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 -->

# Periodic Crystal Structure Model

## Purpose

A backend-independent representation of a periodic crystal unit cell. This is a structure schema, not a crystal predictor or a claim of physical stability.

## Coordinate convention

- Lattice vectors are three row vectors in angstrom.
- Atomic positions are fractional coordinates relative to those lattice vectors.
- Fractional coordinates are preserved as supplied; values outside [0, 1) are allowed and are not silently wrapped.
- Cell volume is the signed determinant of the lattice matrix and must be positive and non-zero.
- Periodic boundary conditions are explicit for each lattice direction. Slab and wire models can disable non-periodic directions, but at least one direction must be periodic.

For fractional coordinate f and row-vector lattice matrix A, Cartesian position is r = f A. Every backend adapter must preserve this convention.

## Atomic sites and occupancy

Each site has a stable ID, element symbol, fractional coordinates, occupancy in (0, 1], and optional oxidation state, magnetic moment in Bohr magnetons, and label. Disordered structures with multiple species at one crystallographic position are represented by separate records sharing the same coordinates. The model groups positions modulo periodic translations and rejects combined occupancies above one. This is a structural consistency check; it does not resolve occupational disorder or generate an ordered supercell.

The occupancy-weighted composition is a convenience summary, not an integer formula or proof of a chemically valid ordered supercell.

## Defect annotations

Defects can be vacancies, substitutions, interstitials, or antisites. Site-based defects reference a stable host site ID. Interstitial defects carry an inserted element and fractional position. Charge state is explicit when known.

These records annotate a proposed defect configuration; they do not automatically edit the host site list, create a supercell, calculate defect formation energy, or establish defect stability. A later defect builder should create an explicit defect-containing structure and retain a mapping to the pristine host.

## Provenance and physical interpretation

The structure may retain a source, space-group label, magnetic-ordering description, total charge, and notes. The space-group string is descriptive metadata only; symmetry is not verified by the initial model. Magnetic ordering is not inferred from per-site moments. Atom identity/order must be mapped explicitly when converting to external solver inputs.

## Current status

Implemented: standard-library crystal/site/defect records, lattice-volume validation, fractional-to-Cartesian conversion, fractional-coordinate preservation, periodic-position occupancy validation, occupancy-weighted composition summary, and unit tests.

Not implemented: CIF/POSCAR import/export, symmetry detection, primitive/conventional cell conversion, supercell generation, defect-structure construction, geometry optimization, or DFT execution.
