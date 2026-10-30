<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 03 | Status: Initial data contract -->

# Domain Data Model

## 1. Purpose

This document defines the initial data contract shared by discovery engines, chemistry modules, structure builders, calculation backends, persistence, and the user interface.

The model is deliberately backend-neutral. External libraries may use different internal representations; adapters must translate to and from the Molecular Lab domain model without losing identity, units, conditions, or provenance.

## 2. Entity relationships

At a conceptual level:

```text
Element / Isotope
       |
       v
Composition ---- SearchConstraints
       |
       v
Candidate
   |       |
   |       +---- Evidence / Provenance
   |
   +---- StructuralAlternative (one or more)
              |
              +---- MolecularGraph (molecular systems)
              +---- MolecularGeometry (finite 3D structure)
              +---- CrystalStructure (periodic systems)
              |
              +---- CalculationRecord (zero or more)
                         |
                         +---- PropertyValue (zero or more)
```

A composition may have multiple candidates or structural alternatives. A calculation belongs to a specific structure/state and a specific set of conditions; it must not be attached only to a formula if multiple structures are being compared.

## 3. Core records

### 3.1 ElementRef

Required:
- `atomic_number`: positive integer.
- `symbol`: canonical chemical symbol.

Optional:
- `name`
- `standard_atomic_weight` with unit or defined convention.
- `data_sources`

Element records are reference data. User-defined or imported values must not overwrite trusted reference data without provenance.

### 3.2 IsotopeRef

Required:
- `atomic_number`
- `mass_number`

Optional:
- `isotopic_mass` with unit.
- `natural_abundance` as a fraction and source.

### 3.3 Composition

Required:
- `counts`: mapping of element or isotope keys to positive integer counts.
- `net_charge`: integer or explicitly unknown.
- `composition_kind`: molecular, empirical, formula_unit, or cell.

Optional:
- `input_formula`: original user text.
- `normalized_formula`: canonical display formula.

Invariant: counts must be positive integers for present species; zero-count entries are omitted. The model must preserve actual atom counts for molecular formulas and cell compositions.

### 3.4 SearchConstraints

Required:
- `allowed_elements`: non-empty set of element references.
- `search_domain`: molecule, periodic, or mixed.
- `max_total_atoms`: positive integer for bounded enumeration.

Optional:
- `min_total_atoms`
- `per_element_count_limits`
- `allowed_charge_range`
- `allowed_multiplicities`
- `required_elements`
- `excluded_elements`
- `temperature`, `pressure`, and environment conditions with units.

Constraints are search instructions, not scientific conclusions. Defaults must be visible to the user and stored with the search session.

### 3.5 MolecularAtom and MolecularBond

Atom required fields:
- stable `atom_id` within the structure,
- element/isotope identity.

Optional atom fields:
- formal charge,
- radical/electron annotation,
- stereochemical annotation,
- atom label.

Bond required fields:
- endpoint atom IDs,
- bond representation supported by the selected graph model.

Optional bond fields:
- bond order,
- aromatic flag,
- stereochemistry,
- coordination or special-bond metadata.

The data model must not assume that every chemical system can be represented completely by a simple classical bond-order graph.

### 3.6 MolecularGeometry

Required:
- ordered atom IDs,
- coordinates,
- coordinate unit,
- geometry role/status.

Coordinates must map unambiguously to atom IDs. A backend adapter must declare any coordinate-unit conversion it performs.

### 3.7 CrystalStructure

The initial runtime implementation is in `molecular_lab/structures/crystal.py`.

Required:
- three row-vector lattice vectors in angstrom,
- atomic sites with stable site IDs, element symbols, and fractional coordinates,
- periodic boundary conditions for each lattice direction.

Optional:
- site occupancy in (0, 1],
- oxidation state and magnetic moment,
- descriptive space-group and magnetic-ordering metadata,
- total cell charge,
- vacancy, substitution, interstitial, and antisite defect annotations.

The cell volume is validated as positive and non-zero. Fractional coordinates are preserved without automatic wrapping. A Cartesian-coordinate view is calculated using r = f A. Species sharing the same periodic site are checked so their occupancies do not sum above one. Defect annotations reference host site IDs where applicable; they do not themselves construct a defect supercell or calculate defect formation energies.

See [CRYSTAL_STRUCTURES.md](CRYSTAL_STRUCTURES.md) for the implemented conventions and limitations.

### 3.8 StructuralAlternative

Required:
- stable identifier,
- parent candidate identifier,
- structure representation,
- charge and spin state or explicit unknown state.

Optional:
- conformer/isomer/polymorph label,
- generation method,
- validation status.

### 3.9 CalculationRecord

Required:
- calculation identifier,
- target structural-alternative identifier,
- backend and method identifiers,
- calculation status,
- input conditions,
- creation timestamp.

Optional:
- backend version,
- basis set or pseudopotential,
- functional/model details,
- convergence settings,
- runtime and resource use,
- output artifact references,
- warnings and failure reason.

### 3.10 DFTCalculationSettings

The initial backend-neutral periodic DFT settings model is implemented in `molecular_lab/calculations/dft_settings.py`.

It records:
- calculation task (SCF, relaxation, variable-cell relaxation, bands, or DOS),
- exchange-correlation functional identifier,
- plane-wave energy cutoff with explicit eV/Ry/Ha unit,
- regular k-point mesh scheme, grid, and shift,
- spin-polarization and spin-orbit-coupling flags,
- optional smearing scheme and width,
- SCF energy, force, and density convergence criteria,
- element-keyed pseudopotential references and provenance fields,
- requested output types and optional temperature/pressure conditions,
- a separate explicit band-path reference for band tasks.

Structure-aware validation checks k-point counts along non-periodic axes and ensures that pseudopotential references cover exactly the structure's elements. It does not verify pseudopotential files, backend compatibility, or scientific convergence. No universal cutoff or k-point recommendation is embedded in the model.

See [DFT_SETTINGS.md](DFT_SETTINGS.md) for units, validation rules, and scientific limitations.

### 3.11 PropertyValue

Required:
- property identifier,
- numeric value or structured value,
- unit,
- calculation/evidence reference,
- conditions when property depends on conditions.

A value without a unit or a clear dimension must not be treated as a valid physical result. Vector, tensor, spectrum, and frequency-dependent values may require structured payloads rather than a single scalar.

### 3.11 ReferenceMolecule

A reference-molecule record represents a previously catalogued molecular identity.

Required fields:
- stable local `record_id`,
- preferred name and molecular formula,
- source database and source record identifier,
- source URL and retrieval date,
- structure identifier such as InChIKey when available.

Optional fields:
- synonyms,
- canonical SMILES or other source-provided structure encodings,
- CAS Registry Number,
- source-specific annotations and licensing metadata.

Formula is a composition key, not a unique structure key. Multiple isomers may share one formula. A formula match must therefore be reported separately from an exact structure match. Exact identity matching requires a standardized structure identifier (initially InChIKey); string equality of arbitrary SMILES is not sufficient unless both strings have been normalized by the same defined canonicalization procedure.

A missing local match means only that the local catalogue has no matching record. It does not establish that a molecule is new or undiscovered.

### 3.12 EvidenceRecord

Required:
- evidence category,
- source or method,
- claim being supported.

Evidence categories should distinguish experimental report, curated database, theoretical calculation, heuristic estimate, and unverified hypothesis. Evidence strength and uncertainty should be separate from the category.

## 4. Identity and relationships

- Use stable IDs for candidates, structural alternatives, calculations, and evidence records.
- Formula strings are display/input forms, not unique structure identifiers.
- Different isomers or polymorphs must have distinct structural-alternative IDs.
- A calculation result must identify the exact structure, state, method, and conditions used.
- Imported records must preserve external identifiers and source provenance.
- Reference-molecule matches must distinguish formula-only matches from exact structure identity.

## 5. Serialization and storage

The first implementation should favor explicit, versioned, human-readable records. JSON is a suitable interchange format for early development, provided schemas and units are explicit.

The persistent schema must include a schema version independent of the project release version. Schema migrations must be planned before changing persisted field meanings.

## 6. Example conceptual record

The following is illustrative only; it is not a calculation result:

```json
{
  "schema_version": "0.1",
  "candidate_id": "candidate-example-001",
  "composition": {
    "counts": {"H": 2, "O": 2},
    "net_charge": 0,
    "composition_kind": "molecular",
    "input_formula": "H2O2"
  },
  "structural_alternatives": [],
  "evidence": []
}
```

An empty structural-alternatives list means that no structure has yet been generated or attached. It does not imply that a valid molecule has been found.

## 7. Open decisions

The following decisions remain open until implementation and testing:
- Exact Python types and validation library.
- Canonical internal unit system and conversion boundary.
- Formula grammar and isotope notation.
- Spin representation and supported electronic states.
- Persistent database technology.
- Schema migration strategy.

## 8. Reference catalogue status

An initial offline reference-molecule catalogue is implemented under `molecular_lab/database/`. It contains a small PubChem-derived seed set and supports formula and synonym lookup. Exact structure matching currently requires an input InChIKey. This is a seed catalogue, not a comprehensive chemical database.

## 9. Implementation status

This document specifies the shared domain contract. Initial runtime models now exist for composition, molecular graphs, periodic crystal/site/defect records, material response, and electronic-structure results. A complete versioned persistence serializer, database schema, and broad molecular/crystal structure conversion pipeline remain unimplemented.
