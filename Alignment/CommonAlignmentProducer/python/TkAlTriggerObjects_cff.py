import FWCore.ParameterSet.Config as cms

# MiniAOD trigger objects for the resonance TkAl ALCARECO streams: the packed
# pat::TriggerObjectStandAlone collection slimmedPatTrigger, configured as in
# MiniAOD (switchOnTriggerStandAlone + packTriggerPathNames), which is what the
# NanoAOD trigger object table reads. ~0.7 kB/event instead of ~4.7 kB/event
# for hltTriggerSummaryAOD in UL16 J/psi MC.
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger as _patTrigger
from PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi import selectedPatTrigger
from PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi import slimmedPatTrigger

patTrigger = _patTrigger.clone(
    onlyStandAlone = True,
    processName = 'HLT',
    packTriggerPathNames = cms.bool(True),
)

seqTkAlTriggerObjects = cms.Sequence(patTrigger + selectedPatTrigger + slimmedPatTrigger)
