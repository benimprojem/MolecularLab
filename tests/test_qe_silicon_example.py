import importlib.util
from pathlib import Path

import pytest


EXAMPLE_PATH = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "qe_silicon_fatbands.py"
)
SPEC = importlib.util.spec_from_file_location("qe_silicon_fatbands_example", EXAMPLE_PATH)
assert SPEC is not None and SPEC.loader is not None
EXAMPLE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXAMPLE)


def test_silicon_smoke_test_example_builds_primitive_diamond_cell():
    lattice_parameter = 5.431
    structure = EXAMPLE.build_silicon(lattice_parameter)

    assert structure.structure_id == "silicon-diamond-primitive"
    assert structure.composition_by_site == {"Si": 2.0}
    assert len(structure.sites) == 2
    assert structure.sites[0].fractional_coordinates == (0.0, 0.0, 0.0)
    assert structure.sites[1].fractional_coordinates == (0.25, 0.25, 0.25)
    assert structure.cell_volume_angstrom3 == pytest.approx(lattice_parameter**3 / 4.0)
