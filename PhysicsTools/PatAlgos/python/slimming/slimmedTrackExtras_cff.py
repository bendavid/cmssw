from RecoMuon.MuonIdentification.thinnedGeneralTracks_cfi import thinnedGeneralTracks

import FWCore.ParameterSet.Config as cms


slimmedGeneralTracks = thinnedGeneralTracks.clone(inputTag = "thinnedGeneralTracks",
                                                  isInputThinned = True,
                                                  muonTag = "slimmedMuons",
                                                  cut = "",
                                                  )

slimmedStandAloneMuons = thinnedGeneralTracks.clone(inputTag = "thinnedStandAloneMuons",
                                                  isInputThinned = True,
                                                  muonTag = "slimmedMuons",
                                                  cut = "",
                                                  )

slimmedGlobalMuons = thinnedGeneralTracks.clone(inputTag = "thinnedGlobalMuons",
                                                  isInputThinned = True,
                                                  muonTag = "slimmedMuons",
                                                  cut = "",
                                                  )

slimmedTevMuonsFirstHit = thinnedGeneralTracks.clone(inputTag = "thinnedTevMuonsFirstHit",
                                                  isInputThinned = True,
                                                  muonTag = "slimmedMuons",
                                                  cut = "",
                                                  )

slimmedTevMuonsPicky = thinnedGeneralTracks.clone(inputTag = "thinnedTevMuonsPicky",
                                                  isInputThinned = True,
                                                  muonTag = "slimmedMuons",
                                                  cut = "",
                                                  )

slimmedTevMuonsDyt = thinnedGeneralTracks.clone(inputTag = "thinnedTevMuonsDyt",
                                                  isInputThinned = True,
                                                  muonTag = "slimmedMuons",
                                                  cut = "",
                                                  )

slimmedTrackExtrasTask = cms.Task(slimmedGeneralTracks,
                                  slimmedStandAloneMuons,
                                  slimmedGlobalMuons,
                                  slimmedTevMuonsFirstHit,
                                  slimmedTevMuonsPicky,
                                  slimmedTevMuonsDyt,
                                  )
