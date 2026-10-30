"""Calculation interfaces and computational backend adapters.

Version: 0.1 | Rev.: 04
"""

from molecular_lab.calculations.qe_projection import (
    AtomicWavefunction, QEAtomicProjectionCoefficient, QEAtomicProjectionData,
    QEProjectionParseError, parse_qe_atomic_proj_xml, read_qe_atomic_proj_xml,
    parse_qe_projwfc_band_output,
)
from molecular_lab.calculations.fatband import (
    FatBandContribution, FatBandDataset, FatBandPoint, build_fatband_dataset,
)
from molecular_lab.calculations.fatband_plot import (
    FatBandPlotData, FatBandPlotPoint, FatBandPlotSeries, FatBandSymmetryTick,
    prepare_fatband_plot_data,
)

from molecular_lab.calculations.pdos_analysis import (
    BandEdgePDOSAnalysis, PDOSWindowSummary, ProjectedDOSAggregate,
    aggregate_projected_dos, analyze_band_edge_pdos,
    integrate_pdos_window, summarize_pdos_window,
)
from molecular_lab.calculations.dft_settings import (
    BandPath, BandPathPoint, ConvergenceSettings, DOSSettings,
    DFTCalculationSettings, DFTOutput, DFTTask, EnergyCutoff, EnergyUnit,
    KPointMesh, KPointScheme, PseudopotentialReference,
    SmearingScheme, SmearingSettings,
)
from molecular_lab.calculations.ase_bridge import (
    ASEBridgeError, ASEUnavailableError, ASEUnsupportedStructureError,
    ase_availability, to_ase_atoms,
)
from molecular_lab.calculations.quantum_espresso import (
    PreparedQEInput, PreparedQEPostprocessInput, QuantumEspressoInputError,
    prepare_pw_input, prepare_qe_bands_postprocess_input,
    prepare_qe_dos_postprocess_input, prepare_qe_projected_dos_postprocess_input,
    prepare_qe_fatband_projection_input,
)
from molecular_lab.calculations.gpaw import (
    GPAWInputError, GPAWUnavailableError, PreparedGPAWInput,
    build_optimization_target, create_gpaw_calculator, prepare_gpaw_input,
)
from molecular_lab.calculations.qe_output import (
    QEOutputParseError, parse_qe_bands_gnu, parse_qe_fermi_level,
    parse_qe_total_dos, read_qe_bands_gnu, read_qe_total_dos,
    parse_qe_projected_dos, read_qe_projected_dos,
)
from molecular_lab.calculations.execution import (
    DFTExecutionResult, DFTExecutionStatus,
    execute_dft_calculation, execute_gpaw, execute_gpaw_scf,
    execute_qe_bands, execute_qe_dos,
    execute_quantum_espresso, execute_quantum_espresso_scf,
)

__all__ = [
    "AtomicWavefunction", "QEAtomicProjectionCoefficient", "QEAtomicProjectionData",
    "QEProjectionParseError", "parse_qe_atomic_proj_xml", "read_qe_atomic_proj_xml",
    "parse_qe_projwfc_band_output",
    "FatBandContribution", "FatBandDataset", "FatBandPoint", "build_fatband_dataset",
    "FatBandPlotData", "FatBandPlotPoint", "FatBandPlotSeries", "FatBandSymmetryTick",
    "prepare_fatband_plot_data",
    "BandEdgePDOSAnalysis", "PDOSWindowSummary", "ProjectedDOSAggregate",
    "aggregate_projected_dos", "analyze_band_edge_pdos",
    "integrate_pdos_window", "summarize_pdos_window",
    "ASEBridgeError", "ASEUnavailableError", "ASEUnsupportedStructureError",
    "ase_availability", "to_ase_atoms",
    "BandPath", "BandPathPoint", "ConvergenceSettings", "DOSSettings",
    "DFTCalculationSettings", "DFTOutput", "DFTTask",
    "EnergyCutoff", "EnergyUnit", "KPointMesh", "KPointScheme",
    "PseudopotentialReference", "SmearingScheme", "SmearingSettings",
    "PreparedQEInput", "PreparedQEPostprocessInput", "QuantumEspressoInputError",
    "prepare_pw_input", "prepare_qe_bands_postprocess_input",
    "prepare_qe_dos_postprocess_input", "prepare_qe_projected_dos_postprocess_input",
    "prepare_qe_fatband_projection_input",
    "QEOutputParseError", "parse_qe_bands_gnu", "parse_qe_fermi_level",
    "parse_qe_total_dos", "read_qe_bands_gnu", "read_qe_total_dos",
    "parse_qe_projected_dos", "read_qe_projected_dos",
    "GPAWInputError", "GPAWUnavailableError", "PreparedGPAWInput",
    "build_optimization_target", "create_gpaw_calculator", "prepare_gpaw_input",
    "DFTExecutionResult", "DFTExecutionStatus",
    "execute_dft_calculation", "execute_gpaw", "execute_gpaw_scf",
    "execute_qe_bands", "execute_qe_dos",
    "execute_quantum_espresso", "execute_quantum_espresso_scf",
]
