from RecoMuon.MuonIdentification.muonTrackExtraThinningProducer_cfi import muonTrackExtraThinningProducer
from RecoMuon.MuonIdentification.muonTrackingRecHitThinningProducer_cfi import muonTrackingRecHitThinningProducer
from RecoTracker.TrackProducer.trackingRecHitThinningProducer_cfi import trackingRecHitThinningProducer
from RecoTracker.TrackProducer.siPixelClusterThinningProducer_cfi import siPixelClusterThinningProducer
from RecoTracker.TrackProducer.siStripClusterThinningProducer_cfi import siStripClusterThinningProducer

import FWCore.ParameterSet.Config as cms

thinnedGeneralTrackExtras = muonTrackExtraThinningProducer.clone(inputTag = "generalTracks",
                                                                  cut = cms.string("pt > 4.5"),
                                                                  slimTrajParams = cms.bool(True),
                                                                  slimResiduals = cms.bool(True))

#standalone muons not needed here because full collection of both TrackExtras and TrackingRecHits are stored in AOD

#global muons and tev muons have full TrackExtra collection stored in AOD, so only thinned TrackingRecHits collection is needed

thinnedGeneralTrackHits = trackingRecHitThinningProducer.clone(inputTag ="generalTracks",
                                                               trackExtraTag = "thinnedGeneralTrackExtras")

#thinned TrackingRecHits collections here are starting directly from muons rather than trackExtras to avoid
#transient intermediate TrackExtra collections

thinnedGlobalMuonHits = muonTrackingRecHitThinningProducer.clone(inputTag = "globalMuons",
                                                                 cut = cms.string("pt > 4.5"))

thinnedTevMuonHitsFirstHit = thinnedGlobalMuonHits.clone(inputTag = "tevMuons:firstHit")

thinnedTevMuonHitsPicky = thinnedGlobalMuonHits.clone(inputTag = "tevMuons:picky")

thinnedTevMuonHitsDyt = thinnedGlobalMuonHits.clone(inputTag = "tevMuons:dyt")

thinnedSiPixelClusters = siPixelClusterThinningProducer.clone(inputTag = "siPixelClusters",
                                                              trackingRecHitsTags = ["thinnedGeneralTrackHits",
                                                                                    "thinnedGlobalMuonHits",
                                                                                    "thinnedTevMuonHitsFirstHit",
                                                                                    "thinnedTevMuonHitsPicky",
                                                                                    "thinnedTevMuonHitsDyt"])

thinnedSiStripClusters = siStripClusterThinningProducer.clone(inputTag = "siStripClusters",
                                                              trackingRecHitsTags = ["thinnedGeneralTrackHits",
                                                                                    "thinnedGlobalMuonHits",
                                                                                    "thinnedTevMuonHitsFirstHit",
                                                                                    "thinnedTevMuonHitsPicky",
                                                                                    "thinnedTevMuonHitsDyt"])

thinnedTrackExtrasTask = cms.Task(thinnedGeneralTrackExtras,
                                  thinnedGeneralTrackHits,
                                  thinnedGlobalMuonHits,
                                  thinnedTevMuonHitsFirstHit,
                                  thinnedTevMuonHitsPicky,
                                  thinnedTevMuonHitsDyt,
                                  thinnedSiPixelClusters,
                                  thinnedSiStripClusters,
                                  )
