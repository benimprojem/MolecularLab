<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 -->

# Reference molecule catalogue

The initial offline catalogue is stored in `reference_molecules.json`. It is a
small verified seed set, not a complete chemical database.

Each record stores:
- Preferred name and known synonyms.
- Molecular formula.
- PubChem Compound ID (CID).
- Canonical SMILES and InChIKey as supplied by the source record.
- CAS Registry Number when included in the seed.
- Source name, source URL, and retrieval date.

## Identity rules

A formula match is not an exact molecular identity match. Multiple isomers can
share the same formula; for example, ethanol and dimethyl ether both have
formula C2H6O.

The current local lookup therefore reports three levels:
- `EXACT_STRUCTURE`: an input InChIKey matches a reference record.
- `FORMULA_ONLY`: one or more records share the element counts, but exact
  structure identity has not been established.
- `NOT_FOUND`: the local catalogue contains no formula match.

`NOT_FOUND` means only that this local catalogue has no record. It must never
be displayed as proof that a molecule is new, unknown, or undiscovered.

## Data source and expansion

The seed records were transcribed from PubChem Compound records. PubChem offers
PUG-REST for programmatic access and bulk data downloads. A later import adapter
can retrieve additional records, preserve source identifiers and provenance,
and update the local cache. Network access must remain optional; offline lookup
must continue to work.

PubChem aggregates contributions from many sources. Provenance and licensing
must be retained per imported record, and the source terms must be reviewed
before redistributing large imported datasets.
