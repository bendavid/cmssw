import FWCore.ParameterSet.Config as cms
sortedPFPrimaryVertices = cms.EDProducer("PFCandidatePrimaryVertexSorter",
    sorting = cms.PSet(),
    assignment = cms.PSet(
    #cuts to assign primary tracks not used in PV fit based on dZ compatibility
    maxDzSigForHighRankedAssignment = cms.double(-1.), # in AND with next
    maxDzForHighRankedAssignment = cms.double(-1.), # in AND with prev
    maxDtSigForHighRankedAssignment = cms.double(-1.), # in AND with prev
    maxDzSigForPrimaryAssignment = cms.double(5.0), # in AND with next
    maxDzForPrimaryAssignment = cms.double(0.1), # in AND with prev
    maxDzErrorForPrimaryAssignment = cms.double(0.05), # in AND with prev, tracks with uncertainty above 500um cannot tell us which pv they come from
    maxDtSigForPrimaryAssignment = cms.double(4.0),
    # cuts used to recover b-tracks if they are closed to jet axis
    maxJetDeltaR = cms.double(0.5),
    minJetPt = cms.double(25),
    maxDistanceToJetAxis = cms.double(0.07), # std cut in b-tag is 700um
    maxDzForJetAxisAssigment = cms.double(0.1), # 1mm, because b-track IP is boost invariant
    maxDxyForJetAxisAssigment = cms.double(0.1), # 1mm, because b-track IP is boost invariant
    maxDtSigForJetAxisAssigment = cms.double(4.0),

    #cuts used to identify primary tracks compatible with beamspot
    maxDxySigForNotReconstructedPrimary = cms.double(2), #in AND with next
    maxDxyForNotReconstructedPrimary = cms.double(0.01), #in AND with prev
    useTiming = cms.bool(False),
    useProbability = cms.bool(False),
    minProbability = cms.double(0.5),
    ),
  particles = cms.InputTag("particleFlow"),
  vertices= cms.InputTag("offlinePrimaryVertices"),
  jets= cms.InputTag("ak4PFJets"),
  qualityForPrimary = cms.int32(3),
  usePVMET = cms.bool(True),
  produceAssociationToOriginalVertices = cms.bool(True),
  produceSortedVertices = cms.bool(True),
  producePileUpCollection  = cms.bool(True),
  produceNoPileUpCollection = cms.bool(True),

)

