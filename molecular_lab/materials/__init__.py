"""Material phases, defects, and physical-property models.

Version: 0.1 | Rev.: 04
"""

from molecular_lab.materials.electronic_structure import (
    BandChannelData, BandGapAssessment, BandGapStatus, BandStructureData,
    DensityOfStatesChannel, DensityOfStatesData, KPoint, assess_band_gap,
    ProjectedDOSComponent, ProjectedDOSRecord, ProjectedDensityOfStatesData,
)
from molecular_lab.materials.electrical_classification import (
    ElectricalBehavior, ElectricalBehaviorAssessment, ElectricalClassificationCriteria,
    ElectricalMetric, assess_electrical_behavior,
)
from molecular_lab.materials.metamaterials import (
    MetamaterialDesign, MetamaterialMechanism, MetamaterialPropertyDomain,
    MetamaterialPropertyTarget, MetamaterialUnitCell, TargetRelation, UnitCellParameter,
)
from molecular_lab.materials.response import (
    ComplexValue, ConvergenceInfo, ConvergenceStatus, MaterialResponseDataset,
    ResponseSample, SolverProvenance,
)

__all__ = [
    "BandChannelData", "BandGapAssessment", "BandGapStatus", "BandStructureData",
    "DensityOfStatesChannel", "DensityOfStatesData", "KPoint", "assess_band_gap",
    "ProjectedDOSComponent", "ProjectedDOSRecord", "ProjectedDensityOfStatesData",
    "ComplexValue", "ConvergenceInfo", "ConvergenceStatus", "ElectricalBehavior",
    "ElectricalBehaviorAssessment", "ElectricalClassificationCriteria", "ElectricalMetric",
    "MaterialResponseDataset", "MetamaterialDesign", "MetamaterialMechanism",
    "MetamaterialPropertyDomain", "MetamaterialPropertyTarget", "MetamaterialUnitCell",
    "ResponseSample", "SolverProvenance", "TargetRelation", "UnitCellParameter",
    "assess_electrical_behavior",
]
