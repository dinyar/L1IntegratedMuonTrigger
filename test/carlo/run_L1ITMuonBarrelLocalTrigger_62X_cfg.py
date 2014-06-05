import FWCore.ParameterSet.Config as cms

process = cms.Process('L1ITMU')


process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")

from SLHCUpgradeSimulations.Configuration.postLS1Customs import*

process = customisePostLS1(process)

process.simCsctfTrackDigis.DTproducer = cms.untracked.InputTag("dttfDigis")
process.simRpcTriggerDigis.label      = cms.string("muonRPCDigis")

process.load('L1Trigger.L1IntegratedMuonTrigger.L1ITMuTriggerPrimitiveProducer_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1CSCTFTrackConverter_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1DTTFTrackConverter_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1RPCTFTrackConverter_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.MBLTProducer_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.MBTracksProducer_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1ITMuonBarrelPrimitiveProducer_cfi')

process.L1ITMuTriggerPrimitives.DT.src  = cms.InputTag('dttfDigis') # DT from unpacker
process.L1ITMuTriggerPrimitives.CSC.src = cms.InputTag('simCscTriggerPrimitiveDigis', 'MPCSORTED') # CSC are kept in 62X GEN-RAW
process.L1ITMuTriggerPrimitives.RPC.src = cms.InputTag('muonRPCDigis') # RPC from unpacker

process.L1CSCTFTrackConverter.CSCTrackSrc = cms.InputTag('simCsctfTrackDigis','') # FIXME CSCTF not there (but do we need it?)
process.L1DTTFTrackConverter.DTTrackSrc = cms.InputTag('dttfDigis','DATA')
#process.L1RPCbTFTrackConverter.RPCTrackSrc = cms.InputTag('gtDigis','RPCb')
#process.L1RPCbTFTrackConverter.RPCL1LinkSrc = cms.InputTag('gtDigis','RPCb')
#process.L1RPCfTFTrackConverter.RPCTrackSrc = cms.InputTag('gtDigis','RPCf')
#process.L1RPCfTFTrackConverter.RPCL1LinkSrc = cms.InputTag('gtDigis','RPCf')

process.MBTracksProducer.DTTrackSrc = cms.InputTag('dttfDigis','DATA')
process.MBTracksProducer.gmtDigis   = cms.InputTag('gtDigis','')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS162_V2::All', '')

#infile = ['file:SingleMuFlatPt_minusEta_1GeVto200GeV_GEN_SIM_DIGI_L1.root']
#infile.append('file:SingleMuFlatPt_plusEta_1GeVto200GeV_GEN_SIM_DIGI_L1.root')
#infile.append('file:SingleMuFlatPt_plusEta_1GeVto200GeV_GEN_SIM_DIGI_L1_2.root')
#infile.append('file:SingleMuFlatPt_minusEta_1GeVto200GeV_GEN_SIM_DIGI_L1_2.root')


process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.threshold = 'WARNING'
process.MessageLogger.cerr.threshold = 'ERROR'

infile = [
'file:/data2/battilan/L1Trigger/62X_RAW_RECO.root',
]

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

process.L1ITMUSequence = cms.Sequence( process.RawToDigi                +
                                       process.simCscTriggerPrimitiveDigis +
                                       process.simCsctfTrackDigis          +
                                       process.simRpcTriggerDigis          +
                                       process.L1ITMuTriggerPrimitives  +
                                       process.L1CSCTFTrackConverter    +
                                       process.L1DTTFTrackConverter     +
                                       process.L1RPCTFTrackConverters   +
                                       process.MBLTProducer             +
                                       process.L1ITMuonBarrelPrimitiveProducer +
                                       process.MBTracksProducer        )

process.L1ITMUPath = cms.Path(process.L1ITMUSequence)

outCommands = cms.untracked.vstring('drop *')
outCommands.append('keep *_muonDTDigis_*_*')
outCommands.append('keep *_simMuonDTDigis_*_*')
outCommands.append('keep *_genParticles_*_*')
outCommands.append('keep *_gtDigis_*_*')
outCommands.append('drop L1GlobalTriggerReadoutRecord_gtDigis_*_*')
outCommands.append('keep *_simMuonRPCDigis_RPCDigiSimLink_*')
outCommands.append('keep *_muonRPCDigis_*_*')
outCommands.append('keep *_dttfDigis_*_*')
outCommands.append('keep *_simRpcTriggerDigis_*_*')
outCommands.append('keep *_simCscTriggerPrimitiveDigis_*_*')
outCommands.append('keep *_L1ITMuTriggerPrimitives_*_*')
outCommands.append('keep *_MBLTProducer_*_*')
outCommands.append('keep *_MBTracksProducer_*_*')
outCommands.append('keep *_L1ITMuonBarrelPrimitiveProducer_*_*')
outCommands.append('keep *_*Converter_*_*')
outCommands.append('keep *_*Matcher_*_*')

process.FEVTDEBUGoutput = cms.OutputModule(
    "PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = outCommands,
    fileName = cms.untracked.string('L1ITMBLT.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    )
)

process.outPath = cms.EndPath(process.FEVTDEBUGoutput)

process.schedule = cms.Schedule(process.L1ITMUPath,
                                process.outPath)

f = file("debug.cfg", 'w')
f.write(process.dumpPython())
f.close()
