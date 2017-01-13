import FWCore.ParameterSet.Config as cms

packedPFCandidates = cms.EDProducer("PATPackedCandidateProducer",
    inputCollection = cms.InputTag("particleFlow"),
    inputVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    originalVertices = cms.InputTag("offlinePrimaryVertices"),
    originalTracks = cms.InputTag("generalTracks"),
    vertexAssociator = cms.InputTag("primaryVertexAssociation","original"),
    PuppiSrc = cms.InputTag("puppi"),
    PuppiNoLepSrc = cms.InputTag("puppiNoLep"),    
    secondaryVerticesForWhiteList = cms.VInputTag(
      cms.InputTag("inclusiveCandidateSecondaryVertices"),
      cms.InputTag("inclusiveCandidateSecondaryVerticesCvsL"),
      ),      
    minPtForTrackProperties = cms.double(0.95),
    storeTiming = cms.bool(False),
)

from Configuration.Eras.Modifier_phase2_timing_cff import phase2_timing
phase2_timing.toModify(packedPFCandidates, 
    storeTiming = cms.bool(True),
    originalVertices = cms.InputTag("offlinePrimaryVertices4D"),
)
