from RecoMuon.MuonIdentification.muonTrackExtraThinningProducer_cfi import muonTrackExtraThinningProducer
from RecoTracker.TrackProducer.trackingRecHitThinningProducer_cfi import trackingRecHitThinningProducer
from RecoTracker.TrackProducer.siPixelClusterThinningProducer_cfi import siPixelClusterThinningProducer
from RecoTracker.TrackProducer.siStripClusterThinningProducer_cfi import siStripClusterThinningProducer

import FWCore.ParameterSet.Config as cms

slimmedGeneralTrackExtras = muonTrackExtraThinningProducer.clone(inputTag = "thinnedGeneralTrackExtras",
                                                  muonTag = "slimmedMuons",
                                                  cut = cms.string("pt > 4.5"),
                                                  slimTrajParams = cms.bool(True),
                                                  slimResiduals = cms.bool(True),
                                                  )

#these ones point to the original full collection of TrackExtras since they are available in the AOD
slimmedStandAloneMuonExtras = slimmedGeneralTrackExtras.clone(inputTag = "standAloneMuons")

slimmedGlobalMuonExtras = slimmedGeneralTrackExtras.clone(inputTag = "globalMuons")

slimmedTevMuonExtrasFirstHit = slimmedGeneralTrackExtras.clone(inputTag = "tevMuons:firstHit")

slimmedTevMuonExtrasPicky = slimmedGeneralTrackExtras.clone(inputTag = "tevMuons:picky")

slimmedTevMuonExtrasDyt = slimmedGeneralTrackExtras.clone(inputTag = "tevMuons:dyt")

slimmedGeneralTrackHits = trackingRecHitThinningProducer.clone(inputTag = "thinnedGeneralTrackHits",
                                                               trackExtraTag = "slimmedGeneralTrackExtras")

#this one points to the original full collection of TrackingRecHits since it is available in the AOD
slimmedStandAloneMuonHits = trackingRecHitThinningProducer.clone(inputTag = "standAloneMuons",
                                                               trackExtraTag = "slimmedStandAloneMuonExtras")

slimmedGlobalMuonHits = trackingRecHitThinningProducer.clone(inputTag = "thinnedGlobalMuonHits",
                                                               trackExtraTag = "slimmedGlobalMuonExtras")

slimmedTevMuonHitsFirstHit = trackingRecHitThinningProducer.clone(inputTag = "thinnedTevMuonHitsFirstHit",
                                                               trackExtraTag = "slimmedTevMuonExtrasFirstHit")

slimmedTevMuonsHitsPicky = trackingRecHitThinningProducer.clone(inputTag = "thinnedTevMuonHitsPicky",
                                                               trackExtraTag = "slimmedTevMuonExtrasPicky")

slimmedTevMuonHitsDyt = trackingRecHitThinningProducer.clone(inputTag = "thinnedTevMuonHitsDyt",
                                                               trackExtraTag = "slimmedTevMuonExtrasDyt")

slimmedSiPixelClusters = siPixelClusterThinningProducer.clone(inputTag = "thinnedSiPixelClusters",
                                                              trackingRecHitsTags = ["slimmedGeneralTrackHits",
                                                                                    "slimmedGlobalMuonHits",
                                                                                    "slimmedTevMuonHitsFirstHit",
                                                                                    "slimmedTevMuonsHitsPicky",
                                                                                    "slimmedTevMuonHitsDyt"])
                                                              
slimmedSiStripClusters = siStripClusterThinningProducer.clone(inputTag = cms.InputTag("thinnedSiStripClusters"),
                                                              trackingRecHitsTags = ["slimmedGeneralTrackHits",
                                                                                    "slimmedGlobalMuonHits",
                                                                                    "slimmedTevMuonHitsFirstHit",
                                                                                    "slimmedTevMuonsHitsPicky",
                                                                                    "slimmedTevMuonHitsDyt"])

slimmedTrackExtrasTask = cms.Task(slimmedGeneralTrackExtras,
                                  slimmedStandAloneMuonExtras,
                                  slimmedGlobalMuonExtras,
                                  slimmedTevMuonExtrasFirstHit,
                                  slimmedTevMuonExtrasPicky,
                                  slimmedTevMuonExtrasDyt,
                                  slimmedGeneralTrackHits,
                                  slimmedStandAloneMuonHits,
                                  slimmedGlobalMuonHits,
                                  slimmedTevMuonHitsFirstHit,
                                  slimmedTevMuonsHitsPicky,
                                  slimmedTevMuonHitsDyt,
                                  slimmedSiPixelClusters,
                                  slimmedSiStripClusters,
                                  )
