from RecoMuon.MuonIdentification.thinnedGeneralTracks_cfi import thinnedGeneralTracks

import FWCore.ParameterSet.Config as cms

thinnedStandAloneMuons = thinnedGeneralTracks.clone(inputTag = cms.InputTag("standAloneMuons"),
                                                  )

thinnedGlobalMuons = thinnedGeneralTracks.clone(inputTag = "globalMuons",
                                                  )

thinnedTevMuonsFirstHit = thinnedGeneralTracks.clone(inputTag = cms.InputTag("tevMuons","firstHit"),
                                                  )

thinnedTevMuonsPicky = thinnedGeneralTracks.clone(inputTag = cms.InputTag("tevMuons","picky"),
                                                  )

thinnedTevMuonsDyt = thinnedGeneralTracks.clone(inputTag = cms.InputTag("tevMuons","dyt"),
                                                  )

thinnedTrackExtrasTask = cms.Task(thinnedGeneralTracks,
                                  thinnedStandAloneMuons,
                                  thinnedGlobalMuons,
                                  thinnedTevMuonsFirstHit,
                                  thinnedTevMuonsPicky,
                                  thinnedTevMuonsDyt,
                                  )
