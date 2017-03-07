import FWCore.ParameterSet.Config as cms

# Reco Vertex
# initialize magnetic field #########################
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesWithBS_cfi import *
from RecoVertex.V0Producer.generalV0Candidates_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *

from CommonTools.RecoAlgos.TrackWithVertexRefSelector_cfi import *
from RecoJets.JetProducers.TracksForJets_cff import *
from CommonTools.RecoAlgos.sortedPrimaryVertices_cfi import *
from RecoJets.JetProducers.caloJetsForTrk_cff import *

#timing
from RecoVertex.PrimaryVertexProducer.TkClusParameters_cff import DA2DParameters
DA2DParameters.TkDAClusParameters.verbose = cms.untracked.bool(False)

unsortedOfflinePrimaryVertices1D = offlinePrimaryVertices.clone()
unsortedOfflinePrimaryVertices1D.TkFilterParameters.minPt = cms.double(0.7)
offlinePrimaryVertices1D=sortedPrimaryVertices.clone(vertices="unsortedOfflinePrimaryVertices1D", particles="trackRefsForJetsBeforeSorting1D")
offlinePrimaryVertices1DWithBS=sortedPrimaryVertices.clone(vertices="unsortedOfflinePrimaryVertices1D:WithBS", particles="trackRefsForJetsBeforeSorting1D")

unsortedOfflinePrimaryVertices4D = offlinePrimaryVertices.clone()
unsortedOfflinePrimaryVertices4D.verbose = cms.untracked.bool(False)
unsortedOfflinePrimaryVertices4D.TkClusParameters = DA2DParameters
unsortedOfflinePrimaryVertices4D.TkFilterParameters.minPt = cms.double(0.7)
unsortedOfflinePrimaryVertices4D.TrackTimesLabel = cms.InputTag("trackTimeValueMapProducer:generalTracksConfigurableFlatResolutionModel")
unsortedOfflinePrimaryVertices4D.TrackTimeResosLabel = cms.InputTag("trackTimeValueMapProducer:generalTracksConfigurableFlatResolutionModelResolution")
offlinePrimaryVertices4D=sortedPrimaryVertices.clone(vertices="unsortedOfflinePrimaryVertices4D", particles="trackRefsForJetsBeforeSorting4D")
offlinePrimaryVertices4D.trackTimeTag=cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel")
offlinePrimaryVertices4D.trackTimeResoTag=cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution")
offlinePrimaryVertices4D.assignment.useTiming = True
offlinePrimaryVertices4DWithBS=sortedPrimaryVertices.clone(vertices="unsortedOfflinePrimaryVertices4D:WithBS", particles="trackRefsForJetsBeforeSorting4D")
offlinePrimaryVertices4DWithBS.trackTimeTag=cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel")
offlinePrimaryVertices4DWithBS.trackTimeResoTag=cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution")
offlinePrimaryVertices4DWithBS.assignment.useTiming = True

unsortedOfflinePrimaryVerticesLegacy = offlinePrimaryVertices.clone()
offlinePrimaryVerticesLegacy=sortedPrimaryVertices.clone(vertices="unsortedOfflinePrimaryVerticesLegacy", particles="trackRefsForJetsBeforeSortingLegacy")
offlinePrimaryVerticesLegacyWithBS=sortedPrimaryVertices.clone(vertices="unsortedOfflinePrimaryVerticesLegacy:WithBS", particles="trackRefsForJetsBeforeSortingLegacy")


trackWithVertexRefSelectorBeforeSorting1D = trackWithVertexRefSelector.clone(vertexTag="unsortedOfflinePrimaryVertices1D")
trackWithVertexRefSelectorBeforeSorting1D.ptMax=9e99
trackWithVertexRefSelectorBeforeSorting1D.ptErrorCut=9e99
trackRefsForJetsBeforeSorting1D = trackRefsForJets.clone(src="trackWithVertexRefSelectorBeforeSorting1D")

trackWithVertexRefSelectorBeforeSorting4D = trackWithVertexRefSelector.clone(vertexTag="unsortedOfflinePrimaryVertices4D")
trackWithVertexRefSelectorBeforeSorting4D.ptMax=9e99
trackWithVertexRefSelectorBeforeSorting4D.ptErrorCut=9e99
trackRefsForJetsBeforeSorting4D = trackRefsForJets.clone(src="trackWithVertexRefSelectorBeforeSorting4D")

trackWithVertexRefSelectorBeforeSortingLegacy = trackWithVertexRefSelector.clone(vertexTag="unsortedOfflinePrimaryVerticesLegacy")
trackWithVertexRefSelectorBeforeSortingLegacy.ptMax=9e99
trackWithVertexRefSelectorBeforeSortingLegacy.ptErrorCut=9e99
trackRefsForJetsBeforeSortingLegacy = trackRefsForJets.clone(src="trackWithVertexRefSelectorBeforeSortingLegacy")

unsortedOfflinePrimaryVertices=offlinePrimaryVertices.clone()
offlinePrimaryVertices=sortedPrimaryVertices.clone(vertices="unsortedOfflinePrimaryVertices", particles="trackRefsForJetsBeforeSorting")
offlinePrimaryVerticesWithBS=sortedPrimaryVertices.clone(vertices=cms.InputTag("unsortedOfflinePrimaryVertices","WithBS"), particles="trackRefsForJetsBeforeSorting")
trackWithVertexRefSelectorBeforeSorting = trackWithVertexRefSelector.clone(vertexTag="unsortedOfflinePrimaryVertices")
trackWithVertexRefSelectorBeforeSorting.ptMax=9e99
trackWithVertexRefSelectorBeforeSorting.ptErrorCut=9e99
trackRefsForJetsBeforeSorting = trackRefsForJets.clone(src="trackWithVertexRefSelectorBeforeSorting")


vertexreco = cms.Sequence(unsortedOfflinePrimaryVertices*
                          trackWithVertexRefSelectorBeforeSorting*
                          trackRefsForJetsBeforeSorting*
                          caloJetsForTrk * 
                          offlinePrimaryVertices*
                          offlinePrimaryVerticesWithBS*
                          generalV0Candidates*
                          inclusiveVertexing
                          )

#timing

from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import tpClusterProducer
from SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi import quickTrackAssociatorByHits
from SimTracker.TrackAssociation.trackTimeValueMapProducer_cfi import trackTimeValueMapProducer
_phase2_tktiming_vertexreco = cms.Sequence( tpClusterProducer *
                                            quickTrackAssociatorByHits *
                                            trackTimeValueMapProducer *
                                            vertexreco.copy() *
                                            unsortedOfflinePrimaryVertices1D *
                                            trackWithVertexRefSelectorBeforeSorting1D *
                                            trackRefsForJetsBeforeSorting1D *
                                            offlinePrimaryVertices1D *
                                            offlinePrimaryVertices1DWithBS *
                                            unsortedOfflinePrimaryVertices4D *
                                            trackWithVertexRefSelectorBeforeSorting4D *
                                            trackRefsForJetsBeforeSorting4D *
                                            offlinePrimaryVertices4D *
                                            offlinePrimaryVertices4DWithBS *
                                            unsortedOfflinePrimaryVerticesLegacy *
                                            trackWithVertexRefSelectorBeforeSortingLegacy *
                                            trackRefsForJetsBeforeSortingLegacy *
                                            offlinePrimaryVerticesLegacy *
                                            offlinePrimaryVerticesLegacyWithBS                                            
                                            )

from Configuration.Eras.Modifier_phase2_timing_cff import phase2_timing
phase2_timing.toReplaceWith(vertexreco, _phase2_tktiming_vertexreco)

#phase2_timing.toModify(
    #unsortedOfflinePrimaryVertices,
    #verbose = cms.untracked.bool(False),
    #TkClusParameters = DA2DParameters,
    #TkFilterParameters=dict(minPt=0.7),
    #TrackTimesLabel = cms.InputTag("trackTimeValueMapProducer:generalTracksConfigurableFlatResolutionModel"),
    #TrackTimeResosLabel = cms.InputTag("trackTimeValueMapProducer:generalTracksConfigurableFlatResolutionModelResolution"),
#)

#phase2_timing.toModify(
    #offlinePrimaryVertices,
    #trackTimeTag=cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    #trackTimeResoTag=cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution"),
    #assignment=dict(useTiming=True)
#)

#phase2_timing.toModify(
    #offlinePrimaryVerticesWithBS,
    #trackTimeTag=cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    #trackTimeResoTag=cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution"),
    #assignment=dict(useTiming=True)
#)
