from CRABClient.UserUtilities import config

config = config()

#label = "MuonGunGlobalCorGen_v29"
#label = "MuonGunGlobalCorRec_GBLNo2dStrip"
#label = "MuonGunGlobalCorRec_GBLAnalyticPropSimpleMF"
#label = "MuonGunGlobalCorRec_GBLAnalyticProp"
#label = "MuonGunDesignFwd_GeantInt4"
#label = "MuonGunDesignFwd150"
#label = "MuonGunDesignFwd150_0T"
#label = "MuonGunDesignFwd150_2T"
#label = "MuonGunDesignFwd150_0p2T"
#label = "MuonGunDesignFwd150_ZeroMaterial"
#label = "MuonGunGlobalCorRec_GBLchisq"
#label = "MuonGunGlobalCorGen_GBLchisq"
#label = "MuonGunGlobalCorRecRefit_v35"

#label = "MuonGunUL2016_v30_Gen"
#label = "MuonGunUL2016_v30_Rec"
#label = "MuonGUNUL2016Fwd_v3_Gen_noquality"
#label = "MuonGUNUL2016Fwd_v3_Gen_idealnoquality"
#label = "MuonGUNUL2016Fwd_v3_Rec_idealnoquality"
#label = "MuonGUNUL2016Fwd_v3_Rec_idealquality"
#label = "MuonGUNUL2016Fwd_v3_Rec_quality"
#label = "MuonGUNUL2016Fwd_v4_Rec_quality"
#label = "MuonGUNUL2016Fwd_v33_Gen_idealquality"
#label = "MuonGUNUL2016Fwd_v33_Rec_idealquality"
label = "MuonGUNUL2016Fwd_v33_Rec_idealnoquality_1"


#label = "MuonGunUL2016APV_v5_Gen"
#label = "MuonGunUL2016APV_v5_Rec"

config.General.requestName = label
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

#config.JobType.pluginName = 'PrivateMC'
#config.JobType.psetName = 'runGlobalCorGen.py'
#config.JobType.psetName = 'muongun_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_RECOSIM.py'
#config.JobType.psetName = 'muongun150.py'
#config.JobType.psetName = 'muongun150zerobfield.py'
#config.JobType.psetName = 'muongun150zeromaterial.py'
#config.JobType.psetName = 'muongun150_2T.py'
#config.JobType.psetName = 'muongun150_0p2T.py'
#config.JobType.psetName = 'muongunfevtsim.py'

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'RunGlobalCorGenIdeal.py'
config.JobType.psetName = 'RunGlobalCorRecIdeal.py'


ncores = 8
outfiles = []
for i in range(ncores):
    outfiles.append("globalcor_%i.root" % i)

print(outfiles)

#config.JobType.psetName = 'runGlobalCorRecRefit.py'
#config.JobType.outputFiles = ["trackTreeGrads.root"]
#config.JobType.eventsPerLumi = 1000

config.JobType.outputFiles = outfiles
config.JobType.maxMemoryMB = ncores*2000
config.JobType.numCores = ncores
config.JobType.maxJobRuntimeMin = 6*60

#config.Data.outputPrimaryDataset = label
#config.Data.splitting = 'EventBased'
##config.Data.unitsPerJob = 10000
##config.Data.unitsPerJob = 32000
#config.Data.unitsPerJob = 20000
#NJOBS = 250  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
##NJOBS = 20  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
#config.Data.publication = True
#config.Data.outputDatasetTag = label

#config.Data.inputDataset = '/DoubleMuonGun_Pt3To150/RunIISummer20UL16RECO-NoPU_106X_mcRun2_asymptotic_v13-v2/GEN-SIM-RECO'
#config.Data.inputDataset = '/DoubleMuonGun_Pt3To150/RunIISummer20UL16RECOAPV-NoPU_106X_mcRun2_asymptotic_preVFP_v8-v2/GEN-SIM-RECO'

config.Data.inputDataset = '/MuonGunDesignFwd_GeantInt/bendavid-MuonGunDesignFwd_GeantInt-7291c1dff7e6b6e708385462ffd0d31e/USER'
#config.Data.inputDataset = '/MuonGunFwd_GeantInt/bendavid-MuonGunFwd_GeantInt-4e53689e60e7ad899aed64c8d53db09e/USER'

config.Data.inputDBS = 'phys03'
#config.Data.inputDBS = 'global'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.unitsPerJob = 30
config.Data.publication = False
config.Data.outputDatasetTag = label


config.Site.storageSite = "T2_CH_CERN"
config.Data.outLFNDirBase = "/store/group/cmst3/group/wmass/bendavid/muoncal"
#config.Data.outLFNDirBase = "/store/group/phys_smp/bendavid"

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
