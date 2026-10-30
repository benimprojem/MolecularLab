"""Calculation interfaces and computational backend adapters.

Version: 0.1 | Rev.: 04
"""

from molecular_lab.calculations.dft_settings import (
    ConvergenceSettings, DFTCalculationSettings, DFTOutput, DFTTask,
    EnergyCutoff, EnergyUnit, KPointMesh, KPointScheme, PseudopotentialReference,
    SmearingScheme, SmearingSettings,
)
from molecular_lab.calculations.ase_bridge import (
    ASEBridgeError, ASEUnavailableError, ASEUnsupportedStructureError,
    ase_availability, to_ase_atoms,
)
from molecular_lab.calculations.quantum_espresso import (
    PreparedQEInput, QuantumEspressoInputError, prepare_pw_input,
)
from molecular_lab.calculations.gpaw import (
    GPAWInputError, GPAWUnavailableError, PreparedGPAWInput,
    build_optimization_target, create_gpaw_calculator, prepare_gpaw_input,
)
from molecular_lab.calculations.execution import (
    DFTExecutionResult, DFTExecutionStatus,
    execute_dft_calculation, execute_gpaw, execute_gpaw_scf,
    execute_quantum_espresso, execute_quantum_espresso_scf,
)

__all__ = [
    "ASEBridgeError", "ASEUnavailableError", "ASEUnsupportedStructureError",
    "ase_availability", "to_ase_atoms",
    "ConvergenceSettings", "DFTCalculationSettings", "DFTOutput", "DFTTask",
    "EnergyCutoff", "EnergyUnit", "KPointMesh", "KPointScheme",
    "PseudopotentialReference", "SmearingScheme", "SmearingSettings",
    "PreparedQEInput", "QuantumEspressoInputError", "prepare_pw_input",
    "GPAWInputError", "GPAWUnavailableError", "PreparedGPAWInput",
    "build_optimization_target", "create_gpaw_calculator", "prepare_gpaw_input",
    "DFTExecutionResult", "DFTExecutionStatus",
    "execute_dft_calculation", "execute_gpaw", "execute_gpaw_scf",
    "execute_quantum_espresso", "execute_quantum_espresso_scf",
]
