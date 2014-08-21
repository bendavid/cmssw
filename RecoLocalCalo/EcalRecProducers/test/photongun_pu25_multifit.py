# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/SingleGammaFlatPt10MeVTo100Gev_cfi.py --conditions auto:upgradePLS1 -n 100 --eventcontent RECOSIM -s GEN,SIM,DIGI,DIGI2RAW,RAW2DIGI,RECO --datatier GEN-SIM-RECO --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --geometry ExtendedZeroMaterial --magField 38T_PostLS1 --io photongun_pu25.io --python photongun_pu25.py --no_exec --fileout file:photongun_pu25.root --pileup AVE_20_BX_25ns --pileup_input das:/RelValMinBias_13/CMSSW_7_1_0_pre5-POSTLS171_V1-v1/GEN-SIM
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
#process.load('Configuration.Geometry.GeometryExtendedZeroMaterial_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
#process.load('Configuration.StandardSequences.Generator_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
#process.load('GeneratorInterface.Core.genFilterSummary_cff')
#process.load('Configuration.StandardSequences.SimIdeal_cff')
#process.load('Configuration.StandardSequences.Digi_cff')
#process.load('Configuration.StandardSequences.DigiToRaw_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#import os
#filesperjob = 5
#jobnum = int(os.getenv("JOBNUM"))
#outfile = os.getenv("OUTFILE")
#filelist = [line.strip() for line in open("/afs/cern.ch/user/b/bendavid/work/CMSSWcluster3/test/lsf/filelist.txt", 'r')]

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    #fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_2_0_pre1/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RAW-HLTDEBUG/POSTLS172_V1-v1/00000/3859E6A1-48FE-E311-B451-0026189438E2.root', 
        #'/store/relval/CMSSW_7_2_0_pre1/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RAW-HLTDEBUG/POSTLS172_V1-v1/00000/4EE44DA8-48FE-E311-BB97-0026189438F2.root', 
        #'/store/relval/CMSSW_7_2_0_pre1/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RAW-HLTDEBUG/POSTLS172_V1-v1/00000/6C5D3ED0-48FE-E311-8EA8-0025905A60B8.root', 
        #'/store/relval/CMSSW_7_2_0_pre1/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RAW-HLTDEBUG/POSTLS172_V1-v1/00000/D05131AD-48FE-E311-B306-002618943821.root')
    #fileNames = cms.untracked.vstring('/store/cmst3/user/bendavid/photongun_nopu/photongun_nopu_963_1_Yga.root'),    
    fileNames = cms.untracked.vstring('/store/cmst3/user/bendavid/photongun_pu25/photongun_pu25_980_2_OoL.root'),
    #fileNames = cms.untracked.vstring(),
)

#for filename in filelist[jobnum*filesperjob:(jobnum+1)*filesperjob]:
#  process.source.fileNames.append(filename)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('Configuration/Generator/python/SingleGammaFlatPt10MeVTo100Gev_cfi.py nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    fastCloning = cms.untracked.bool(False),                                         
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string(outfile),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    ),
)

process.RECOSIMoutput.outputCommands.extend(cms.untracked.vstring('keep EBDigiCollection_ecalDigis_*_*',
        'keep EEDigiCollection_ecalDigis_*_*',
        'keep ESDigiCollection_ecalPreshowerDigis_*_*'))

# Additional output definition

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS1', '')


# Path and EndPath definitions
#process.generation_step = cms.Path(process.pgen)
#process.simulation_step = cms.Path(process.psim)
#process.digitisation_step = cms.Path(process.pdigi)
#process.digi2raw_step = cms.Path(process.DigiToRaw)
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.reconstruction_step = cms.Path(process.reconstruction)
#process.digitisation_step = cms.Path(cms.SequencePlaceholder("randomEngineStateProducer")*cms.SequencePlaceholder("mix")*process.ecalDigiSequence*process.addPileupInfo) 
#process.digi2raw_step = cms.Path(process.ecalPacker*process.esDigiToRaw*process.rawDataCollector)
#process.raw2digi_step = cms.Path(process.ecalDigis*process.ecalPreshowerDigis)
process.reconstruction_step = cms.Path(process.offlineBeamSpot*process.ecalLocalRecoSequence*process.pfClusteringPS*process.pfClusteringECAL*process.ecalClusters)

#process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.ecalRecHit.recoverEBFE = False
process.ecalRecHit.recoverEEFE = False
process.ecalRecHit.killDeadChannels = False
process.particleFlowSuperClusterECAL.useRegression = False

process.ecalLocalRecoMinimal = cms.Sequence(process.ecalMaxSampleUncalibRecHit*process.ecalRecHit*process.ecalPreshowerRecHit)

process.reconstruction_step = cms.Path(process.offlineBeamSpot*process.ecalLocalRecoMinimal*process.pfClusteringPS*process.pfClusteringECAL*process.ecalClusters)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)

# Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.digi2raw_step,process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
## filter all path with the production filter sequence
#for path in process.paths:
	#getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
