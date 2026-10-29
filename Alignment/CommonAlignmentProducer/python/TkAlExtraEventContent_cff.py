import FWCore.ParameterSet.Config as cms

# Extra event content for the resonance TkAl ALCARECO streams.

# Trigger content needed to fill the HLT and L1 bits and the trigger objects
# (TrigObj, matching of offline muons to L1/HLT objects) downstream:
#   - hltTriggerSummaryAOD: HLT objects and their filter keys
#   - gtStage2Digis (algorithm and external blocks only, the GT input object
#     copies are redundant with gmtStage2Digis/caloStage2Digis): L1 algorithm
#     decisions (initial/final, prescale column)
#   - gmtStage2Digis:Muon, caloStage2Digis:{EGamma,EtSum,Jet,Tau}: L1 objects
#     read by the NanoAOD TriggerObjectTableProducer
#   - hltGtStage2ObjectMap: L1 seed object combinations (only if available)
# The legacy gtDigis record kept by the stock streams is empty for Stage-2.
TkAlTriggerOutputCommands = cms.untracked.vstring(
    'keep triggerTriggerEvent_hltTriggerSummaryAOD_*_*',
    'keep *_hltGtStage2ObjectMap_*_*',
    'keep *_gtStage2Digis__*',
    'keep *_gmtStage2Digis_Muon_*',
    'keep *_caloStage2Digis_EGamma_*',
    'keep *_caloStage2Digis_EtSum_*',
    'keep *_caloStage2Digis_Jet_*',
    'keep *_caloStage2Digis_Tau_*',
    'keep recoBeamSpot_offlineBeamSpot_*_*',
)

# MC truth: generator record, pileup and the sim content needed for the
# gen and sim-hit association of the track refit (ResidualGlobalCorrectionMaker*:
# genParticles incl. barcodes and xyz0, generator, pileup, SimTracks and the
# tracker PSimHits). Added for MC only, see CustomiseMCJpsiX.keepMCTruth.
TkAlMCTruthOutputCommands = cms.untracked.vstring(
    'keep GenRunInfoProduct_generator_*_*',
    'keep GenLumiInfoHeader_generator_*_*',
    'keep GenLumiInfoProduct_generator_*_*',
    'keep GenEventInfoProduct_generator_*_*',
    'keep *_genParticles_*_*',
    'keep PileupSummaryInfos_*_*_*',
    'keep int_*_bunchSpacing_*',
    'keep SimTracks_g4SimHits_*_*',
    'keep SimVertexs_g4SimHits_*_*',
    'keep PSimHits_g4SimHits_TrackerHits*_*',
)
