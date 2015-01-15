from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.workArea = 'MBPrimitivesProduction'
config.General.requestName = 'particleGunGenSim_MuMinus_v1' #also dest dir
config.General.transferLogs = True
#config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'particleGunGenSimCrab_MuMinus.py'
#config.JobType.allowNonProductionCMSSW = True

config.section_('Data')
config.Data.primaryDataset = 'SingleMuMinus'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 10000
NJOBS = 10
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'MBPrimitivesProduction_particleGunGenSim_MuMinus_v1'

config.section_('Site')
config.Site.storageSite = 'T3_IT_Bologna'
#config.Site.storageSite = 'T2_IT_Legnaro'

