import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.pfPileUp_cfi  import *
from CommonTools.ParticleFlow.TopProjectors.pfNoPileUp_cfi import *
from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import *

pfPileUpJME = pfPileUp.clone( PFCandidates='particleFlowPtrs',
                              Vertices = 'goodOfflinePrimaryVertices',
                              Jets=cms.InputTag('ak4PFJets'),
                              checkClosestZVertex = False )
pfNoPileUpJME = pfNoPileUp.clone( topCollection = 'pfPileUpJME',
                                  bottomCollection = 'particleFlowPtrs' )



pfNoPileUpJMESequence = cms.Sequence(
    goodOfflinePrimaryVertices +
    pfPileUpJME +
    pfNoPileUpJME
    )

from Configuration.Eras.Modifier_phase2_common_cff import phase2_common
phase2_common.toModify(
    pfPileUpJME,
    Vertices = cms.InputTag("offlinePrimaryVertices"),
)
