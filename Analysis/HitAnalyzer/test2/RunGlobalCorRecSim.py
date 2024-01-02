# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --mc --runUnscheduled --conditions auto:run2_design --step RECO --nThreads 32 --geometry DB:Extended --era Run2_2016 --dasquery file dataset=/MuonGunDesignFwd/bendavid-MuonGunDesignFwd-449be899dc2f7c04b70495711fe9dfd6/USER instance=prod/phys03 --processName RECO2 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016

process = cms.Process('RECO2',Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.MagneticField_0T_cff')
#process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.load('TrackPropagation.Geant4e.Geant4ePropagator_cfi')

#process.UniformMagneticFieldESProducer = cms.ESProducer("UniformMagneticFieldESProducer",
    #ZFieldInTesla = cms.double(3.8),
    #label = cms.untracked.string('ParabolicMf')
#)

#process.UniformMagneticFieldESProducer2 = cms.ESProducer("UniformMagneticFieldESProducer",
    #ZFieldInTesla = cms.double(3.8),
    #label = cms.untracked.string('')
#)

#process.VolumeBasedMagneticFieldESProducer.valueOverride = 9500
#process.ParabolicParametrizedMagneticFieldProducer.valueOverride = 9500

#0.2T
#process.VolumeBasedMagneticFieldESProducer.valueOverride = 950
#process.ParabolicParametrizedMagneticFieldProducer.valueOverride = 950

#process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(2)
#)

# Input source
#process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring(
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_1.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_10.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_11.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_12.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_13.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_14.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_15.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_16.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_17.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_18.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_19.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_2.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_20.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_3.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_4.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_5.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_6.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_7.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_8.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150/MuonGunDesignFwd150/200917_172653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_9.root'
    #),
    #secondaryFileNames = cms.untracked.vstring()
#)

#process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring(
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_1.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_10.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_11.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_12.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_13.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_14.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_15.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_16.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_17.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_18.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_19.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_2.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_20.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_3.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_4.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_5.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_6.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_7.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_8.root', 
        #'/store/group/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0T/MuonGunDesignFwd150_0T/200918_095655/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_9.root'
    #),
    #secondaryFileNames = cms.untracked.vstring()
#)


#process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring(
        ##'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_1.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_10.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_11.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_12.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_13.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_14.root', 
        ##'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_15.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_16.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_17.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_18.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_19.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_2.root', 
        ##'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_20.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_3.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_4.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_5.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_6.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_7.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_8.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_2T/MuonGunDesignFwd150_2T/200918_124828/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_9.root'
    #),
    #secondaryFileNames = cms.untracked.vstring()
#)


#process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring(
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_1.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_10.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_11.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_12.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_13.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_14.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_15.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_16.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_17.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_18.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_19.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_2.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_20.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_3.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_4.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_5.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_6.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_7.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_8.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_ZeroMaterial/MuonGunDesignFwd150_ZeroMaterial/200918_121653/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_9.root'
    #),
    #secondaryFileNames = cms.untracked.vstring()
#)

#process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring(
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_1.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_10.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_11.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_12.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_13.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_14.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_15.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_16.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_17.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_18.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_19.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_2.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_20.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_3.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_4.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_5.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_6.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_7.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_8.root', 
        #'/store/cmst3/group/wmass/bendavid/muoncal/MuonGunDesignFwd150_0p2T/MuonGunDesignFwd150_0p2T/200918_163610/0000/muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM_9.root'
    #),
    #secondaryFileNames = cms.untracked.vstring()
#)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #'file:/data/bendavid/cmsswdev/muonscale/CMSSW_10_6_17_patch1/work/muongundet/muongundet_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM.root',
       #'file:/data/bendavid/cmsswdev/muonscale/CMSSW_10_6_17_patch1/work/muongundetmoreprecise/muongundet_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM.root',
       #'file:/data/bendavid/cmsswdev/muonscale/CMSSW_10_6_17_patch1/work/muongundetmoreprecisenoms/muongundet_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM.root',
       #'file:/data/bendavid/cmsswdev/muonscale/CMSSW_10_6_17_patch1/work/muongundetnoms/muongundet_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM.root',
        #'file:/data/bendavid/cmsswdev/muonscale/CMSSW_10_6_17_patch1/work/muongundetmoremoreprecisenoms/muongundet_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM.root',
        'file:/data/bendavid/cmsswdev/muonscale/CMSSW_10_6_17_patch1/work/muongundet_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM.root',
        #'file:/data/bendavid/cmsswdev/muonscale/CMSSW_10_6_17_patch1/work/simtestdi1en7/muongundet_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM.root',
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('step1_RECO.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_design', '')

#process.geopro = cms.EDProducer("GeometryProducer",
     #UseMagneticField = cms.bool(True),
     #UseSensitiveDetectors = cms.bool(False),
     #MagneticField = cms.PSet(
         #UseLocalMagFieldManager = cms.bool(False),
         #Verbosity = cms.untracked.bool(False),
         #ConfGlobalMFM = cms.PSet(
             #Volume = cms.string('OCMS'),
             #OCMS = cms.PSet(
                 #Stepper = cms.string('G4ClassicalRK4'),
                 #Type = cms.string('CMSIMField'),
                    #StepperParam = cms.PSet(
                        #VacRegions = cms.vstring(),
        ##                   VacRegions = cms.vstring('DefaultRegionForTheWorld','BeamPipeVacuum','BeamPipeOutside'),
                        #MaximumEpsilonStep = cms.untracked.double(0.01),   ## in mm
                        #DeltaOneStep = cms.double(0.001),      ## in mm
                        #MaximumLoopCounts = cms.untracked.double(1000.0),
                        #DeltaChord = cms.double(0.001), ## in mm
                        #MinStep = cms.double(0.1),      ## in mm
                        #DeltaIntersectionAndOneStep = cms.untracked.double(-1.0),
                        #DeltaIntersection = cms.double(0.0001),## in mm
                        #MaxStep = cms.double(100000.),         ## in cm
                        #MinimumEpsilonStep = cms.untracked.double(1e-05), ## in mm
                        #EnergyThSimple = cms.double(0.002),               ## in GeV
                        #DeltaChordSimple = cms.double(0.1),    ## in mm
                        #DeltaOneStepSimple = cms.double(0.1),  ## in mm
                        #DeltaIntersectionSimple = cms.double(0.01),       ## in mm
                        #MaxStepSimple = cms.double(100000.),   ## in cm
                    #)
                ##StepperParam = cms.PSet(
                     ##MaximumEpsilonStep = cms.untracked.double(0.01), ## in mm
                     ##DeltaOneStep = cms.double(0.001), ## in mm
                     ##MaximumLoopCounts = cms.untracked.double(1000.0),
                     ##DeltaChord = cms.double(0.001), ## in mm
                     ##MinStep = cms.double(0.1), ## in mm
                     ##DeltaIntersectionAndOneStep = cms.untracked.double(-1.0),
                     ##DeltaIntersection = cms.double(0.0001), ## in mm
                     ##MinimumEpsilonStep = cms.untracked.double(1e-05) ## in mm
                 ##)
             #)
         #),
         #delta = cms.double(1.0)
     #)
   #)

#process.TTRHBuilderAngleAndTemplate.StripCPE = "stripCPE"

#process.TFileService = cms.Service("TFileService", fileName = cms.string("globalcor.root") )

process.globalCor = cms.EDAnalyzer('ResidualGlobalCorrectionMaker',
                                   src = cms.InputTag("generalTracks"),
                                   #fitFromGenParms = cms.bool(False),
                                   fitFromGenParms = cms.bool(True),
                                   fillTrackTree = cms.bool(True),
                                   doSim = cms.bool(True),
)

#process.globalCorOut = cms.EDAnalyzer('ResidualGlobalCorrectionOutput',
                                   #src = cms.InputTag("globalCor"),
                                   #fitFromGenParms = cms.bool(False),
                                   #fillTrackTree = cms.bool(True),
#)

#process.globalCorGen = cms.EDAnalyzer('ResidualGlobalCorrectionMaker',
                                   #src = cms.InputTag("generalTracks"),
                                   #fitFromGenParms = cms.bool(True),
                                   #fillTrackTree = cms.bool(False),
#)


# Path and EndPath definitions
#process.reconstruction_step = cms.Path(process.reconstruction_fromRECO)
#process.reconstruction_step = cms.Path(process.globalCor*process.globalCorOut)
process.reconstruction_step = cms.Path(process.globalCor)
#process.reconstruction_step = cms.Path(process.geopro*process.globalCor)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
process.schedule = cms.Schedule(process.reconstruction_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(32)
#process.options.numberOfStreams=cms.untracked.uint32(0)
#process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
