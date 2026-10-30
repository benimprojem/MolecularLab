import pytest

from molecular_lab.calculations.ase_bridge import (
    ASEUnavailableError, ASEUnsupportedStructureError, ase_availability,
    to_ase_atoms,
)
from molecular_lab.structures.crystal import (
    CrystalDefect, CrystalSite, CrystalStructure, DefectKind,
)


def make_crystal(**kwargs):
    values = dict(
        structure_id="Si",
        lattice_vectors_angstrom=((5.43, 0, 0), (0, 5.43, 0), (0, 0, 5.43)),
        sites=(
            CrystalSite("Si1", "Si", (0, 0, 0), magnetic_moment_mu_b=0.5),
            CrystalSite("Si2", "Si", (0.25, 0.25, 0.25)),
        ),
    )
    values.update(kwargs)
    return CrystalStructure(**values)


def test_ase_availability_returns_boolean_and_reason():
    available, reason = ase_availability()
    assert isinstance(available, bool)
    assert isinstance(reason, str)
    if available:
        assert reason == ""
    else:
        assert "ase" in reason.lower()


def test_rejects_non_crystal_input():
    with pytest.raises(TypeError):
        to_ase_atoms(object())


def test_partial_occupancy_is_not_silently_discarded():
    structure = make_crystal(sites=(
        CrystalSite("Si1", "Si", (0, 0, 0), occupancy=0.5),
    ))
    with pytest.raises(ASEUnsupportedStructureError, match="partial occupancies"):
        to_ase_atoms(structure)


def test_defect_annotations_are_not_silently_discarded():
    structure = make_crystal(defects=(
        CrystalDefect("vac1", DefectKind.VACANCY, host_site_id="Si1"),
    ))
    with pytest.raises(ASEUnsupportedStructureError, match="defect annotations"):
        to_ase_atoms(structure)


def test_conversion_when_ase_is_installed():
    pytest.importorskip("ase")
    atoms = to_ase_atoms(make_crystal())
    assert atoms.get_chemical_symbols() == ["Si", "Si"]
    assert atoms.get_pbc().tolist() == [True, True, True]
    assert atoms.info["structure_id"] == "Si"
    assert atoms.info["site_ids"] == ["Si1", "Si2"]
    assert atoms.get_initial_magnetic_moments().tolist() == [0.5, 0.0]
    assert atoms.get_scaled_positions(wrap=False)[1].tolist() == pytest.approx(
        [0.25, 0.25, 0.25]
    )


def test_missing_ase_is_reported_explicitly(monkeypatch):
    import molecular_lab.calculations.ase_bridge as bridge

    monkeypatch.setattr(bridge, "ase_availability", lambda: (
        False, "Optional dependency 'ase' is not installed"
    ))
    with pytest.raises(ASEUnavailableError, match="not installed"):
        bridge.to_ase_atoms(make_crystal())
