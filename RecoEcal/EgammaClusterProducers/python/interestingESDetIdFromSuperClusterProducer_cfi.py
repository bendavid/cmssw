import FWCore.ParameterSet.Config as cms

interestingESDetIdFromSuperClusterProducer = cms.EDProducer("InterestingESDetIdFromSuperClusterProducer",
                                                    superClustersLabel = cms.InputTag(''),
                                                    recHitsLabel = cms.InputTag('ecalPreshowerRecHit','EcalRecHitsES'),
                                                    interestingDetIdCollection = cms.string(''),
                                                    )
