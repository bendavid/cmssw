from RecoMuon.MuonIdentification.thinnedGeneralTracks_cfi import thinnedGeneralTracks

import FWCore.ParameterSet.Config as cms


slimmedGeneralTracks = thinnedGeneralTracks.clone(inputTag = "thinnedGeneralTracks",
                                                  isInputThinned = True,
                                                  muonTag = "slimmedMuons",
                                                  cut = "",
                                                  )

slimmedStandAloneMuons = slimmedGeneralTracks.clone(inputTag = "thinnedStandAloneMuons",
                                                  )

slimmedGlobalMuons = slimmedGeneralTracks.clone(inputTag = "thinnedGlobalMuons",
                                                  )

slimmedTevMuonsFirstHit = slimmedGeneralTracks.clone(inputTag = "thinnedTevMuonsFirstHit",
                                                  )

slimmedTevMuonsPicky = slimmedGeneralTracks.clone(inputTag = "thinnedTevMuonsPicky",
                                                  )

slimmedTevMuonsDyt = slimmedGeneralTracks.clone(inputTag = "thinnedTevMuonsDyt",
                                                  )

slimmedTrackExtrasTask = cms.Task(slimmedGeneralTracks,
                                  slimmedStandAloneMuons,
                                  slimmedGlobalMuons,
                                  slimmedTevMuonsFirstHit,
                                  slimmedTevMuonsPicky,
                                  slimmedTevMuonsDyt,
                                  )
