"""Run a small silicon QE fat-band smoke test.

This is a real-solver example, not a convergence study. Supply a PBE silicon
pseudopotential and its recommended cutoffs for that exact file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

from molecular_lab.calculations import (
    BandPath,
    BandPathPoint,
    DFTCalculationSettings,
    DFTOutput,
    DFTTask,
    EnergyCutoff,
    EnergyUnit,
    KPointMesh,
    KPointScheme,
    PseudopotentialReference,
    execute_dft_calculation,
    prepare_fatband_plot_data,
    render_fatband_svg,
)
from molecular_lab.structures.crystal import CrystalSite, CrystalStructure


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudo", type=Path, required=True, help="PBE Si UPF file")
    parser.add_argument("--pseudo-family", required=True, help="Pseudo family/version label")
    parser.add_argument("--pseudo-source", required=True, help="Download URL or source citation")
    parser.add_argument("--pseudo-license", default="", help="Pseudopotential license or usage terms")
    parser.add_argument("--ecutwfc-ry", type=float, required=True, help="Recommended wavefunction cutoff")
    parser.add_argument("--ecutrho-ry", type=float, required=True, help="Recommended charge-density cutoff")
    parser.add_argument("--a", type=float, default=5.431, help="Starting cubic lattice parameter in angstrom")
    parser.add_argument("--kmesh", type=int, default=6, help="SCF Monkhorst-Pack grid count per axis")
    parser.add_argument("--points-per-segment", type=int, default=30)
    parser.add_argument("--workdir", type=Path, default=Path("calculations/qe-silicon-fatbands"))
    parser.add_argument("--pw", default="pw.x")
    parser.add_argument("--bands", default="bands.x")
    parser.add_argument("--projwfc", default="projwfc.x")
    return parser.parse_args()


def build_silicon(a: float) -> CrystalStructure:
    half = a / 2.0
    return CrystalStructure(
        structure_id="silicon-diamond-primitive",
        lattice_vectors_angstrom=(
            (0.0, half, half),
            (half, 0.0, half),
            (half, half, 0.0),
        ),
        sites=(
            CrystalSite("Si1", "Si", (0.0, 0.0, 0.0)),
            CrystalSite("Si2", "Si", (0.25, 0.25, 0.25)),
        ),
        source="Ideal diamond-cubic silicon primitive-cell smoke-test model",
        notes=(
            "Starting lattice parameter only; relax/converge before scientific interpretation.",
        ),
    )


def build_silicon_band_path(points_per_segment: int) -> BandPath:
    if (
        isinstance(points_per_segment, bool)
        or not isinstance(points_per_segment, int)
        or points_per_segment < 1
    ):
        raise ValueError("points_per_segment must be a positive integer")
    return BandPath((
        BandPathPoint((0.0, 0.0, 0.0), "Γ", points_per_segment),
        BandPathPoint((0.5, 0.0, 0.5), "X", points_per_segment),
        BandPathPoint((0.5, 0.5, 0.5), "L", points_per_segment),
        BandPathPoint((0.0, 0.0, 0.0), "Γ", 0),
    ))


def main() -> int:
    args = parse_args()
    pseudo_source = args.pseudo.expanduser().resolve()
    if not pseudo_source.is_file():
        raise FileNotFoundError(f"Silicon pseudopotential not found: {pseudo_source}")
    if args.kmesh < 1 or args.points_per_segment < 1:
        raise ValueError("kmesh and points-per-segment must be positive")
    if args.a <= 0:
        raise ValueError("lattice parameter must be positive")

    workdir = args.workdir.expanduser().resolve()
    pseudo_dir = workdir / "pseudo"
    pseudo_dir.mkdir(parents=True, exist_ok=True)
    pseudo_target = pseudo_dir / pseudo_source.name
    if pseudo_source != pseudo_target.resolve():
        shutil.copy2(pseudo_source, pseudo_target)
    checksum = hashlib.sha256(pseudo_target.read_bytes()).hexdigest()

    structure = build_silicon(args.a)
    # Fractional coordinates in the primitive FCC reciprocal basis.
    # This deliberately short Γ-X-L-Γ path is a smoke-test path, not a
    # complete high-symmetry path for publication-quality band analysis.
    band_path = build_silicon_band_path(args.points_per_segment)
    pseudo = PseudopotentialReference(
        element="Si",
        file_name=pseudo_target.name,
        family=args.pseudo_family,
        source=args.pseudo_source,
        checksum=checksum,
        license=args.pseudo_license,
    )
    settings = DFTCalculationSettings(
        task=DFTTask.BANDS,
        exchange_correlation="PBE",
        energy_cutoff=EnergyCutoff(args.ecutwfc_ry, EnergyUnit.RYDBERG),
        charge_density_cutoff=EnergyCutoff(args.ecutrho_ry, EnergyUnit.RYDBERG),
        kpoint_mesh=KPointMesh(
            (args.kmesh, args.kmesh, args.kmesh),
            KPointScheme.MONKHORST_PACK,
        ),
        pseudopotentials={"Si": pseudo},
        requested_outputs=(DFTOutput.FAT_BAND_STRUCTURE,),
        band_path=band_path,
        notes="Silicon QE fat-band smoke test; convergence must be assessed separately.",
    )

    result = execute_dft_calculation(
        backend="qe",
        structure=structure,
        settings=settings,
        working_directory=workdir,
        executable=args.pw,
        bands_executable=args.bands,
        projwfc_executable=args.projwfc,
    )
    results_dir = workdir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = results_dir / "qe_stdout.log"
    stderr_path = results_dir / "qe_stderr.log"
    stdout_path.write_text(result.stdout, encoding="utf-8")
    stderr_path.write_text(result.stderr, encoding="utf-8")

    manifest = {
        "workflow": "silicon_qe_fatband_smoke_test",
        "status": result.status.value,
        "message": result.message,
        "qe_version": result.backend_version or None,
        "return_code": result.return_code,
        "converged": result.converged,
        "total_energy_ev": result.total_energy_ev,
        "structure": {
            "id": structure.structure_id,
            "lattice_vectors_angstrom": structure.lattice_vectors_angstrom,
            "fractional_sites": [
                {"site_id": site.site_id, "element": site.element,
                 "fractional_coordinates": site.fractional_coordinates}
                for site in structure.sites
            ],
        },
        "pseudopotential": {
            "element": "Si",
            "file_name": pseudo_target.name,
            "family": args.pseudo_family,
            "source": args.pseudo_source,
            "license": args.pseudo_license,
            "sha256": checksum,
        },
        "settings": {
            "exchange_correlation": "PBE",
            "ecutwfc_ry": args.ecutwfc_ry,
            "ecutrho_ry": args.ecutrho_ry,
            "scf_kpoint_mesh": [args.kmesh, args.kmesh, args.kmesh],
            "band_path_coordinate_system": band_path.coordinate_system,
            "band_path": [
                {"coordinates": point.coordinates, "label": point.label,
                 "points_to_next": point.points_to_next}
                for point in band_path.points
            ],
        },
        "logs": {"stdout": str(stdout_path), "stderr": str(stderr_path)},
    }
    manifest_path = results_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if not result.succeeded:
        raise RuntimeError(
            f"QE fat-band workflow failed ({result.status.value}): {result.message}. "
            f"See {stdout_path} and {stderr_path}"
        )
    if result.fat_band_structure is None:
        raise RuntimeError(
            f"QE returned no fat-band dataset. See {manifest_path}"
        )

    plot_data = prepare_fatband_plot_data(
        result.fat_band_structure,
        band_path=band_path,
        lattice_vectors_angstrom=structure.lattice_vectors_angstrom,
        group_by="element_orbital",
        normalization="global_max",
    )
    svg = render_fatband_svg(
        plot_data,
        title="Silicon — QE fat-band smoke test",
        width=1280,
        height=800,
    )
    output = results_dir / "silicon-fat-bands.svg"
    output.write_text(svg, encoding="utf-8")
    manifest["outputs"] = {
        "fat_band_svg": str(output),
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"QE version: {result.backend_version or 'not reported'}")
    print(f"Total energy: {result.total_energy_ev} eV")
    print(f"Fat-band SVG: {output}")
    print(f"Run manifest: {manifest_path}")
    print("Reminder: this run is not a cutoff/k-mesh convergence study.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
