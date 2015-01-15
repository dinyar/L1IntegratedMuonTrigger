from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.workArea = 'MBPrimitivesProduction'
config.General.requestName = 'Neutrino_Pt-2to20_gun-AVE20BX25_v1' #also dest dir
config.General.transferLogs = True
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_L1ITMuonBarrelLocalTrigger_72X_25ns.cfg.py'
#config.JobType.allowNonProductionCMSSW = True

config.section_('Data')
#config.Data.inputDataset = '/DYToMuMu_Tune4C_13TeV-pythia8/Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/GEN-SIM-RAW'
config.Data.inputDataset = '/Neutrino_Pt-2to20_gun/Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/GEN-SIM-RAW'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.inputDBS = 'global'
#config.Data.publishDBS = 'phys03'
#config.Data.publishDataName = 'MBPrimitivesProduction_Neutrino_Pt-2to20_gun_v1'

config.section_('Site')
config.Site.storageSite = 'T3_IT_Bologna'
#config.Site.storageSite = 'T2_IT_Legnaro'

