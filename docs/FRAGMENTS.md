<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 -->

# Molecular Fragment Model

## Purpose

Define reusable molecular subgraphs that can later be assembled into larger candidates. This first increment adds data records and validation only; it does not yet provide a fragment catalogue, connection algorithm, or large-molecule generator.

## Data model

- `MolecularFragment` wraps an existing `MolecularGraph`, a stable fragment ID, name, version, source, and description.
- `FragmentPort` identifies an atom that may participate in a connection and lists the allowed proposed bond orders.
- Port identifiers are unique within a fragment and must reference an atom in that fragment.
- Provenance is mandatory so curated, imported, and generated fragments can be distinguished.

A port is an assembly permission, not a guarantee that a connection is chemically valid. An assembler must check combined valence, formal charge, aromaticity, stereochemistry, atom mapping, and duplicate structures after connecting fragments.

## Planned assembly workflow

1. Select fragments from a versioned, provenance-aware library.
2. Match compatible ports and allowed bond orders.
3. Apply element-, valence-, charge-, and chemistry-specific connection rules.
4. Enforce atom, fragment, ring, and candidate budgets.
5. Canonicalize and deduplicate assembled graphs.
6. Preserve fragment identities and connection history in the resulting candidate.
7. Send selected candidates to conformer generation and suitable calculation backends.

Search truncation and unsupported connection chemistry must be reported explicitly. Neither should be treated as evidence that a molecule cannot exist.

## Initial exclusions

The first model does not encode stereochemistry, aromatic attachment semantics, formal charge, radicals, isotopes, partial bonds, polymer repeat-unit semantics, or non-covalent interactions. These require explicit extensions and validators before being used in scientific assembly.

## Current status

Implemented: immutable fragment and port records, basic identity checks, unique port checks, and validation that ports reference atoms in the wrapped molecular graph.

Not implemented: built-in curated fragment library, chemical connection validator, fragment assembler, canonical assembly deduplication, conformer generation, and large-molecule calculations.
