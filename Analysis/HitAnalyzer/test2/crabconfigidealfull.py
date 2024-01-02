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

#label = "MuonGunUL2016_v33_Gen_quality"
#label = "MuonGunUL2016_v33_Rec_quality"
#label = "MuonGunUL2016_v33_Rec_noquality"

#label = "MuonGunUL2016_v34_Gen_quality"
#label = "MuonGunUL2016_v34_Rec_quality"

#label = "MuonGunUL2016_v35_Gen_quality"
#label = "MuonGunUL2016_v35_Rec_quality"

##label = "MuonGunUL2016_v33_Gen_idealquality"
#label = "MuonGunUL2016_v33_Rec_idealquality"

#label = "MuonGunUL2016_v54_Gen_idealquality"
#label = "MuonGunUL2016_v54_Rec_idealquality"

#label = "MuonGunUL2016_v55_Rec_idealquality"
#label = "MuonGunUL2016_v56_Rec_idealquality"
#label = "MuonGunUL2016_v57_Rec_idealquality"
#label = "MuonGunUL2016_v121_Rec_idealquality_zeromaterial"

#label = "MuonGunUL2016_v123_RecJpsiPhotos_idealquality_zeromaterial"

#label = "MuonGunUL2016_v149_RecJpsiPhotosSingle_idealquality_biasm10"
#label = "MuonGunUL2016_v149_RecJpsiPhotos_idealquality_constraint_biasm10"
#label = "MuonGunUL2016_v149_Rec_idealquality_biasm10"

#label = "MuonGunUL2016_v179_RecJpsiPhotosSingle_idealquality"
#label = "MuonGunUL2016_v176_RecJpsiPhotos_idealquality_constraint"
#label = "MuonGunUL2016_v179_Rec_idealquality"

#label = "MuonGunUL2016_v202_Gen_idealquality_zeromaterial"
##label = "MuonGunUL2016_v202_GenJpsiPhotosSingle_idealquality_zeromaterial"

#label = "MuonGunUL2016_v259_Gen_idealquality"
#label = "MuonGunUL2016_v286_GenJpsiPhotosSingle_idealquality"
#label = "MuonGunUL2016_v321_GenSimJpsiPhotosSingle_idealquality"
#label = "MuonGunUL2016_v325a_GenJpsiPhotosSingle_idealquality"

# label = "MuonGunUL2016_v709_Rec_idealquality"
label = "MuonGunUL2016_v710_Gen_idealquality"

#label = "MuonGunUL2016_v567_RecJpsiPhotosSingle_idealquality"

#label = "MuonGunUL2016_v387_RecJpsiPhotos_idealquality_constraintnofsr"
#label = "MuonGunUL2016_v301_RecJpsiPhotos_idealquality_constraintfsr29"
#label = "MuonGunUL2016_v292_RecJpsiPhotos_idealquality_constraintfsr29_novtx_nograds"
#label = "MuonGunUL2016_v225_RecJpsiPhotos_idealquality_constraintnofsr_nograds"


#label = "MuonGunUL2016_v277_Rec_idealquality"
#label = "MuonGunUL2016_v277_RecJpsiPhotosSingle_idealquality"

#label = "MuonGunUL2016_v171_Gen_idealquality_zeromaterial_nograds"
#label = "MuonGunUL2016_v171_GenJpsiPhotosSingle_idealquality_zeromaterial_nograds"

#label = "MuonGunUL2016_v166_Gen_idealquality_nograds"
#label = "MuonGunUL2016_v166_GenJpsiPhotosSingle_idealquality_nograds"

#label = "MuonGunUL2016_v163_Gen_idealquality_zeromaterial_nograds"
#label = "MuonGunUL2016_v163_GenJpsiPhotosSingle_idealquality_zeromaterial_nograds"

#label = "MuonGunUL2016_v149_RecJpsiPhotosSingle_idealquality_zeromaterial"
#label = "MuonGunUL2016_v149_RecJpsiPhotos_idealquality_zeromaterial_constraint"
#label = "MuonGunUL2016_v149_Rec_idealquality_zeromaterial"

#label = "MuonGunUL2016_v123_RecJpsiPhotos_idealquality_zeromaterial_noconstraint"
#label = "MuonGunUL2016_v123a_RecJpsiPhotosSingle_idealquality_zeromaterial"
#label = "MuonGunUL2016_v121_Gen_idealquality_zeromaterial"
#label = "MuonGunUL2016_v131a_Gen_idealquality"


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
config.JobType.psetName = 'RunGlobalCorGenIdeal.py'
# config.JobType.psetName = 'RunGlobalCorRecIdeal.py'
#config.JobType.psetName = 'RunGlobalCorRecJpsiIdeal.py'


ncores = 1
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
#config.JobType.maxJobRuntimeMin = 6*60

# config.JobType.inputFiles = ["correctionResults_v560_idealgen_matm10.root"]


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

#config.Data.inputDataset = '/MuonGunDesign_GeantInt/bendavid-MuonGunDesign_GeantInt-69e9a86bc950310c8d202348f84bf938/USER'

#config.Data.inputDataset = '/MuonGunDesign_GeantIntZeroMaterial/bendavid-MuonGunDesign_GeantIntZeroMaterial-a58ca486694af8abd4f911ab5fa2b8fe/USER'
#config.Data.inputDataset = '/JPsiToMuMuGun_Pt0To30withTail-pythia8-photos_Ideal_ZeroMaterial/bendavid-JPsiToMuMuGun_Pt0To30withTail-pythia8-photos_Ideal_ZeroMaterial-b2ab5baffccc221e55480b3601a241ea/USER'


#config.Data.inputDataset = '/DoubleMuonGun_Pt3To150/RunIISummer20UL16RECO-NoPU_106X_mcRun2_asymptotic_v13-v2/GEN-SIM-RECO'
#config.Data.inputDataset = '/DoubleMuonGun_Pt3To150/RunIISummer20UL16RECOAPV-NoPU_106X_mcRun2_asymptotic_preVFP_v8-v2/GEN-SIM-RECO'

#config.Data.inputDataset = '/MuonGunDesign_GeantIntElossIoniOnly_v2/phys_smp-MuonGunDesign_GeantIntElossIoniOnly_v2-24fff29e0e65e7b0aebc5e3ec31e923f/USER'


#config.Data.inputDataset = '/JPsiToMuMuGun_Pt0To30withTail-pythia8-photos_Ideal/bendavid-JPsiToMuMuGun_Pt0To30withTail-pythia8-photos_Ideal-83067967246c6edbcda633d2e73221a8/USER'
config.Data.inputDataset = '/DoubleMuonGun_Pt3To150_Ideal/bendavid-DoubleMuonGun_Pt3To150_Ideal-69e9a86bc950310c8d202348f84bf938/USER'

#config.Data.inputDataset = "/DoubleMuonGun_Pt3To150_IdealNoScattering/bendavid-DoubleMuonGun_Pt3To150_IdealNoScattering-69e9a86bc950310c8d202348f84bf938/USER"
#config.Data.inputDataset = "/JPsiToMuMuGun_Pt0To30withTail-pythia8-photos_IdealNoScattering/bendavid-JPsiToMuMuGun_Pt0To30withTail-pythia8-photos_IdealNoScattering-83067967246c6edbcda633d2e73221a8/USER"


#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'


config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

#config.Data.splitting = 'EventAwareLumiBased'
#config.Data.unitsPerJob = 10000
#config.Data.totalUnits = 999*10000

#config.Data.unitsPerJob = 30
config.Data.publication = False
config.Data.outputDatasetTag = label


config.Site.storageSite = "T2_CH_CERN"
#config.Data.outLFNDirBase = "/store/group/cmst3/group/wmass/bendavid/muoncal"
config.Data.outLFNDirBase = "/store/group/phys_smp/ec/bendavid/muoncal"

config.Data.ignoreLocality = True
#config.Data.ignoreLocality = False
#config.Site.blacklist = ["T2_IT_Pisa"]
#config.Site.blacklist = ["T2_US_Florida"]
#config.Site.whitelist = ["T2_US_Nebraska", "T2_IT_Pisa"]
#config.Site.whitelist = ["T2_US_Nebraska"]
config.Site.whitelist = ["T2*"]
#config.Site.whitelist = ["T2_US_*"]
#config.Site.whitelist = ["T2_CH_CERN"]
#config.Site.whitelist = ["T2_US_MIT"]

#config.section_("Debug")
#config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
