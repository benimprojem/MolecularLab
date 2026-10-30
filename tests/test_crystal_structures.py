"""Tests for periodic crystal structure records."""
import pytest
from molecular_lab.structures.crystal import CrystalDefect, CrystalSite, CrystalStructure, DefectKind

def make_structure(**kwargs):
    values = dict(
        structure_id="silicon-cell",
        lattice_vectors_angstrom=((5.43, 0.0, 0.0), (0.0, 5.43, 0.0), (0.0, 0.0, 5.43)),
        sites=(CrystalSite("Si-1", "Si", (0.0, 0.0, 0.0)), CrystalSite("Si-2", "Si", (0.25, 0.25, 0.25))),
    )
    values.update(kwargs)
    return CrystalStructure(**values)

def test_volume_and_site_composition():
    structure = make_structure()
    assert structure.cell_volume_angstrom3 == pytest.approx(5.43 ** 3)
    assert structure.composition_by_site == {"Si": 2.0}

def test_fractional_coordinates_are_not_wrapped():
    structure = make_structure(sites=(CrystalSite("O1", "O", (1.1, -0.2, 0.5), 0.75),))
    assert structure.sites[0].fractional_coordinates == (1.1, -0.2, 0.5)
    assert structure.composition_by_site == {"O": 0.75}

def test_rejects_singular_and_left_handed_cells():
    with pytest.raises(ValueError, match="positive, non-zero volume"):
        make_structure(lattice_vectors_angstrom=((1, 0, 0), (2, 0, 0), (0, 0, 1)))
    with pytest.raises(ValueError, match="positive, non-zero volume"):
        make_structure(lattice_vectors_angstrom=((-1, 0, 0), (0, 1, 0), (0, 0, 1)))

def test_rejects_invalid_occupancy_and_duplicate_site_ids():
    with pytest.raises(ValueError, match="occupancy"):
        CrystalSite("x", "Si", (0, 0, 0), 0.0)
    site = CrystalSite("Si1", "Si", (0, 0, 0))
    with pytest.raises(ValueError, match="unique"):
        make_structure(sites=(site, site))

def test_vacancy_and_substitution_annotations():
    defects = (
        CrystalDefect("vac1", DefectKind.VACANCY, host_site_id="Si-2"),
        CrystalDefect("dop1", DefectKind.SUBSTITUTION, host_site_id="Si-1", inserted_element="P"),
    )
    structure = make_structure(defects=defects)
    assert structure.defects[0].kind is DefectKind.VACANCY
    assert structure.defects[1].inserted_element == "P"

def test_interstitial_uses_fractional_coordinates():
    defect = CrystalDefect("int1", DefectKind.INTERSTITIAL, inserted_element="Li", fractional_coordinates=(0.5, 0.5, 0.5), charge_state=1)
    assert make_structure(defects=(defect,)).defects[0] == defect
    with pytest.raises(ValueError, match="requires fractional coordinates"):
        CrystalDefect("bad", DefectKind.INTERSTITIAL, inserted_element="Li")

def test_rejects_unknown_defect_site_and_invalid_pbc():
    defect = CrystalDefect("vac1", DefectKind.VACANCY, host_site_id="unknown")
    with pytest.raises(ValueError, match="unknown host site"):
        make_structure(defects=(defect,))
    with pytest.raises(ValueError, match="three booleans"):
        make_structure(periodic_boundary_conditions=(True, False))

def test_rejects_empty_or_fully_nonperiodic_cell():
    with pytest.raises(ValueError, match="at least one site"):
        make_structure(sites=())
    with pytest.raises(ValueError, match="at least one periodic direction"):
        make_structure(periodic_boundary_conditions=(False, False, False))

def test_partial_occupancy_is_retained():
    sites = (CrystalSite("Na", "Na", (0, 0, 0), 0.5), CrystalSite("K", "K", (0, 0, 0), 0.5))
    assert make_structure(sites=sites).composition_by_site == {"K": 0.5, "Na": 0.5}


def test_fractional_to_cartesian_coordinate_conversion():
    structure = make_structure(
        lattice_vectors_angstrom=((2.0, 0.0, 0.0), (1.0, 3.0, 0.0), (0.0, 0.0, 4.0)),
        sites=(CrystalSite("X", "Si", (0.5, 0.5, 0.25)),),
    )
    assert structure.cartesian_coordinates_angstrom == ((1.5, 1.5, 1.0),)


def test_disordered_species_occupancy_cannot_exceed_one():
    sites = (
        CrystalSite("Na", "Na", (0.0, 0.0, 0.0), 0.75),
        CrystalSite("K", "K", (1.0, 0.0, 0.0), 0.5),
    )
    with pytest.raises(ValueError, match="cannot sum above one"):
        make_structure(sites=sites)
