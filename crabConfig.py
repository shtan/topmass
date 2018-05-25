from CRABClient.client_utilities import getBasicConfig
config = getBasicConfig()

config.General.requestName = 'test1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = none
outputFiles = fitresults.root
scriptExe = 'runcrab.sh'

config.Data.inputDataset = '/GenericTTbar/HC-CMSSW_5_3_1_START53_V5-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFN = '/store/user/testdircrab1' # or '/store/group/<subdir>'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

config.Site.storageSite = 'T3_US_FNAL'
