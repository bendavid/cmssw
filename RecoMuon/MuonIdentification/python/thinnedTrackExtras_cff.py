from RecoMuon.MuonIdentification.thinnedGeneralTracks_cfi import thinnedGeneralTracks

import FWCore.ParameterSet.Config as cms

thinnedStandAloneMuons = thinnedGeneralTracks.clone(inputTag = cms.InputTag("standAloneMuons"),
                                                  ptMin = 3.,
                                                  )

thinnedGlobalMuons = thinnedGeneralTracks.clone(inputTag = "globalMuons",
                                                  ptMin = 3.,
                                                  )

thinnedTevMuonsFirstHit = thinnedGeneralTracks.clone(inputTag = cms.InputTag("tevMuons","firstHit"),
                                                  ptMin = 3.,
                                                  )

thinnedTevMuonsPicky = thinnedGeneralTracks.clone(inputTag = cms.InputTag("tevMuons","picky"),
                                                  ptMin = 3.,
                                                  )

thinnedTevMuonsDyt = thinnedGeneralTracks.clone(inputTag = cms.InputTag("tevMuons","dyt"),
                                                  ptMin = 3.,
                                                  )

thinnedTrackExtrasTask = cms.Task(thinnedGeneralTracks,
                                  thinnedStandAloneMuons,
                                  thinnedGlobalMuons,
                                  thinnedTevMuonsFirstHit,
                                  thinnedTevMuonsPicky,
                                  thinnedTevMuonsDyt,
                                  )
