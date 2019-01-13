import FWCore.ParameterSet.Config as cms

#process = cms.Process("Demo")

from Configuration.StandardSequences.Eras import eras

process = cms.Process('ANALYSIS',eras.Phase2C4_timing_layer_bar)


process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'file:testpipt0p7/step3.root',
            'file:step3_RAW2DIGI_RECO.root',
            #'file:testpipt0p7/step3nom.root',
            #'file:testmupt10/step3offset.root',
                )
                            )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.demo = cms.EDAnalyzer('BackProp'
                              )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('timeres.root')
                                   )

process.p = cms.Path(process.demo)

