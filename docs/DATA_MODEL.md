<!-- Project: Molecular Lab | Version: 0.1 | Rev.: 00 | Status: Initial data contract -->

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

Required:
- lattice vectors,
- lattice unit,
- atomic sites with fractional coordinates,
- periodicity definition.

Optional:
- occupancies,
- symmetry metadata,
- magnetic moments/order,
- supercell and defect annotations.

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

### 3.10 PropertyValue

Required:
- property identifier,
- numeric value or structured value,
- unit,
- calculation/evidence reference,
- conditions when property depends on conditions.

A value without a unit or a clear dimension must not be treated as a valid physical result. Vector, tensor, spectrum, and frequency-dependent values may require structured payloads rather than a single scalar.

### 3.11 EvidenceRecord

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

## 8. Implementation status

This is a domain contract only. No runtime models, serializers, database tables, or validation routines are implemented by this document.
