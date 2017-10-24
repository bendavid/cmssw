import FWCore.ParameterSet.Config as cms
from CommonTools.RecoAlgos.sortedPFPrimaryVertices_cfi import sortedPFPrimaryVertices

primaryVertexAssociation = sortedPFPrimaryVertices.clone(
  qualityForPrimary = cms.int32(2),
  produceSortedVertices = cms.bool(False),
  producePileUpCollection  = cms.bool(False),  
  produceNoPileUpCollection = cms.bool(False)
)

from Configuration.Eras.Modifier_phase2_common_cff import phase2_common
phase2_common.toModify(
    primaryVertexAssociation,
    assignment=dict(maxDzErrorForPrimaryAssignment = 999.,
                    maxDzSigForHighRankedAssignment = 4.0,
                    ),
)

from Configuration.Eras.Modifier_phase2_timing_layer_cff import phase2_timing_layer
phase2_timing_layer.toModify(
    primaryVertexAssociation,
    assignment=dict(maxDtSigForHighRankedAssignment = 3.0,
                    useTiming = True,
                    ),
)
