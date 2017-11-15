import FWCore.ParameterSet.Config as cms
from CommonTools.RecoAlgos.sortedPFPrimaryVertices_cfi import sortedPFPrimaryVertices

primaryVertexAssociation = sortedPFPrimaryVertices.clone(
  qualityForPrimary = cms.int32(2),
  produceSortedVertices = cms.bool(False),
  producePileUpCollection  = cms.bool(False),  
  produceNoPileUpCollection = cms.bool(False)
)

primaryVertexAssociation1D = sortedPFPrimaryVertices.clone(
  qualityForPrimary = cms.int32(2),
  produceSortedVertices = cms.bool(False),
  producePileUpCollection  = cms.bool(False),  
  produceNoPileUpCollection = cms.bool(False)
)

from Configuration.Eras.Modifier_phase2_common_cff import phase2_common
phase2_common.toModify(
    primaryVertexAssociation,
    assignment=dict(maxDzErrorForPrimaryAssignment = 999.,
                    maxDzForHighRankedAssignment = 0.1,
                    ),
)

from Configuration.Eras.Modifier_phase2_timing_layer_cff import phase2_timing_layer
phase2_timing_layer.toModify(
    primaryVertexAssociation,
    assignment=dict(useTiming = True,
                    maxDtSigForHighRankedAssignment = 3.0,
                    ),
)

phase2_timing_layer.toModify(
    primaryVertexAssociation1D,
    vertices= cms.InputTag("offlinePrimaryVertices1D"),
    assignment=dict(maxDzErrorForPrimaryAssignment = 999.,
                    maxDzForHighRankedAssignment = 0.1,
                    ),
)
