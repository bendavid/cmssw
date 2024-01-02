from CRABClient.UserUtilities import config

config = config()

#label = "MuonGunGlobalCorGen_v29"
#label = "MuonGunGlobalCorRec_GBLNo2dStrip"
#label = "MuonGunGlobalCorRec_GBLAnalyticPropSimpleMF"
#label = "MuonGunGlobalCorRec_GBLAnalyticProp"
#label = "MuonGunDesignFwd_GeantIntPrecise"
#label = "MuonGunDesignFwd_GeantIntIdealField"
#label = "MuonGunDesignFwd_GeantInt"
#label = "MuonGunFwd_GeantInt"
label = "MuonGunDesign_GeantIntZeroMaterial"
#label = "MuonGunDesignFwd150"
#label = "MuonGunDesignFwd150_0T"
#label = "MuonGunDesignFwd150_2T"
#label = "MuonGunDesignFwd150_0p2T"
#label = "MuonGunDesignFwd150_ZeroMaterial"
#label = "MuonGunGlobalCorRec_GBLchisq"
#label = "MuonGunGlobalCorGen_GBLchisq"
#label = "MuonGunGlobalCorRecRefit_v35"

config.General.requestName = label
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'PrivateMC'
#config.JobType.psetName = 'runGlobalCorGen.py'
#config.JobType.psetName = 'muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM.py'
#config.JobType.psetName = 'muongun150.py'
#config.JobType.psetName = 'muongun150zerobfield.py'
#config.JobType.psetName = 'muongun150zeromaterial.py'
#config.JobType.psetName = 'muongun150_2T.py'
#config.JobType.psetName = 'muongun150_0p2T.py'
#config.JobType.psetName = 'muongunfevtsim.py'
#config.JobType.psetName = 'muongunfevtsimnonideal.py'
config.JobType.psetName = 'muongunfevtsimfull.py'


#config.JobType.psetName = 'runGlobalCorRecRefit.py'
#config.JobType.outputFiles = ["trackTreeGrads.root"]
config.JobType.eventsPerLumi = 1000
config.JobType.maxMemoryMB = 16000
config.JobType.numCores = 8

config.Data.outputPrimaryDataset = label
config.Data.splitting = 'EventBased'
#config.Data.unitsPerJob = 10000
#config.Data.unitsPerJob = 32000
#config.Data.unitsPerJob = 20000
config.Data.unitsPerJob = 50000
NJOBS = 1800
#NJOBS = 250  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
#NJOBS = 20  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True
config.Data.outputDatasetTag = label

config.Site.storageSite = "T2_CH_CERN"
#config.Data.outLFNDirBase = "/store/group/cmst3/group/wmass/bendavid/muoncal"
config.Data.outLFNDirBase = "/store/group/phys_smp/bendavid"

#config.Data.ignoreLocality = True
#config.Data.ignoreLocality = False
#config.Site.blacklist = ["T2_IT_Pisa"]
#config.Site.blacklist = ["T2_US_Florida"]
#config.Site.whitelist = ["T2_US_Nebraska", "T2_IT_Pisa"]
#config.Site.whitelist = ["T2_US_Nebraska"]
#config.Site.whitelist = ["T2*"]
#config.Site.whitelist = ["T2_CH_CERN"]

#config.section_("Debug")
#config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
