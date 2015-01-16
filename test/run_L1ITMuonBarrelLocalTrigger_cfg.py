import FWCore.ParameterSet.Config as cms

process = cms.Process('L1ITMU')


process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('L1Trigger.L1IntegratedMuonTrigger.L1ITMuTriggerPrimitiveProducer_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1CSCTFTrackConverter_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1DTTFTrackConverter_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1RPCTFTrackConverter_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.MBLTProducer_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.MBTracksProducer_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1ITMuonBarrelPrimitiveProducer_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START72_V1::All', '')

#infile = ['file:SingleMuFlatPt_minusEta_1GeVto200GeV_GEN_SIM_DIGI_L1.root']
#infile.append('file:SingleMuFlatPt_plusEta_1GeVto200GeV_GEN_SIM_DIGI_L1.root')
#infile.append('file:SingleMuFlatPt_plusEta_1GeVto200GeV_GEN_SIM_DIGI_L1_2.root')
#infile.append('file:SingleMuFlatPt_minusEta_1GeVto200GeV_GEN_SIM_DIGI_L1_2.root')


process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.threshold = 'WARNING'
process.MessageLogger.cerr.threshold = 'ERROR'

rawfiles = [
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/1A2F2E5F-D581-E411-83A0-002481E0DA4E.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/2CDF009F-D881-E411-A45B-002590AB3A70.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/56958F40-D781-E411-909F-002481E0D448.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/602AE1B3-F081-E411-96A9-0025904B0FE4.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/7078C897-D181-E411-A719-00266CFFA604.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/763EC846-D481-E411-86DD-002590D94F8E.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/7E7739A1-DB81-E411-81BB-003048D436B4.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/96DCF1F0-D381-E411-9000-002590DB9232.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/96E19D79-D681-E411-9370-003048D4DF80.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/ACEEA343-D981-E411-BBEA-0025904B1452.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/AE3EFB13-D381-E411-B203-002481E0DDE8.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/B27108BF-D481-E411-BAF9-0025901D42BC.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/BC92E6F4-D681-E411-B96B-00266CF33288.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/C03EA077-D681-E411-B45E-0025901D490C.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/CA957AFA-D681-E411-AEDA-00266CF9B420.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/DC778DBD-D581-E411-B695-002590AC4C74.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/DE0ED418-E181-E411-A72F-0025907DCA4A.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/E2590A79-D881-E411-97EC-0025901D490C.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/E4232FA1-D381-E411-AA15-00266CFFA7A8.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/E647BE4B-D381-E411-8B13-00266CF33288.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/E678D1D2-D081-E411-A9D9-002590AC4FEA.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/E8493540-D281-E411-8BC8-00266CF33288.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RAW/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/E851BA36-D681-E411-A301-0025901D4C74.root'
]
recofiles = [
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/044CA221-ED81-E411-B545-0025904B0F96.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/18ABBC4B-EA81-E411-8634-002590AC5074.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/1E748B94-EA81-E411-B990-002590DB9216.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/24CB1B26-E481-E411-BB28-00266CF32CD0.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/26918BE2-E581-E411-B639-00266CF32F18.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/2C50D1FD-0582-E411-A969-0025904B0F96.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/30FEEECB-EA81-E411-B037-0025904B12A4.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/34CAB61F-EB81-E411-A5AC-0025904B0FE4.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/42A68E8B-EA81-E411-939A-00266CF33118.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/4643C441-EA81-E411-981B-002590DB91A0.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/5A569D46-EA81-E411-A85F-002590AC5062.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/5AE6FE43-E781-E411-90DF-00266CF9B420.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/5C3216D1-EA81-E411-8FC3-00266CFFA124.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/64792538-E481-E411-9523-002481E0D480.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/661DC8E4-E481-E411-A149-00266CF32CD0.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/6EA05E77-EA81-E411-A44F-00266CF327C0.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/6ED410E1-E581-E411-BE33-00266CF32F18.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/789D3444-E681-E411-9132-00266CF2718C.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/78B23DE3-EC81-E411-A071-00266CF26450.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/8661B70C-E781-E411-8C5B-0025901D4C32.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/944801B5-E681-E411-9AF5-00266CF32CD0.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/98B87EFA-EB81-E411-8A31-008CFA104E64.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/9EAFD183-EA81-E411-9E14-00266CF32CD0.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/B2A46D70-EA81-E411-A644-0025901D490C.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/B665FA82-EA81-E411-9036-0025904B0F96.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/B69EBCC8-E081-E411-81AB-002590DB918C.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/CC172A8B-EA81-E411-BB49-00266CFFA780.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/CEDEF30A-EC81-E411-8FA3-0025904B1370.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/DA4A82F5-E181-E411-8ED5-002481E0D708.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/DAA26A13-E381-E411-9FC1-00266CF32BC4.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/DE81E382-EA81-E411-A978-0025904B0F96.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/EAA8B8C5-EA81-E411-9AF9-0025904B12A4.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/ECB6E7B2-E681-E411-9994-00266CF32920.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/F24EF5EC-F181-E411-8A24-0025904B12A8.root',
'/store/mc/Phys14DR/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/PU20bx25_tsg_141126_PHYS14_25_V1-v1/40000/FA70E0CF-EA81-E411-B50C-008CFA104E64.root'
]
process.source = cms.Source(
'PoolSource',
fileNames = cms.untracked.vstring(rawfiles),
secondaryFileNames = cms.untracked.vstring(recofiles)
)

process.L1ITMUSequence = cms.Sequence(  process.L1ITMuTriggerPrimitives +
                                        process.L1CSCTFTrackConverter   +
                                        process.L1DTTFTrackConverter    +
                                        process.L1RPCTFTrackConverters  +
                                        process.MBLTProducer            +
                                        process.L1ITMuonBarrelPrimitiveProducer +
                                        process.MBTracksProducer
                                    )

process.L1ITMUPath = cms.Path(process.L1ITMUSequence)

outCommands = cms.untracked.vstring('drop *')
outCommands.append('keep *_*uon*_*_RECO')
outCommands.append('keep *_simMuonDTDigis_*_*')
outCommands.append('keep *_genParticles_*_*')
outCommands.append('keep *_simCsctfDigis_*_*')
outCommands.append('keep *_simDttfDigis_*_*')
outCommands.append('keep *_simGmtDigis_*_*')
outCommands.append('keep *_simRpcTriggerDigis_*_*')
outCommands.append('keep *_simMuonRPCDigis_*_*')
outCommands.append('keep *_simDtTriggerPrimitiveDigis_*_*')
outCommands.append('keep *_simCscTriggerPrimitiveDigis_*_*')
outCommands.append('keep *_L1ITMuTriggerPrimitives_*_*')
outCommands.append('keep *_MBLTProducer_*_*')
outCommands.append('keep *_MBTracksProducer_*_*')
outCommands.append('keep *_L1ITMuonBarrelPrimitiveProducer_*_*')
outCommands.append('keep *_*Converter_*_*')
outCommands.append('keep *_*Matcher_*_*')

process.FEVTDEBUGoutput = cms.OutputModule(
    "PoolOutputModule",
    splitLevel                   = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands               = outCommands,
    fileName                     = cms.untracked.string('L1ITMBLT.root'),
    dataset                      = cms.untracked.PSet(
        filterName               = cms.untracked.string(''),
        dataTier                 = cms.untracked.string('')
        )
    )

process.outPath = cms.EndPath(process.FEVTDEBUGoutput)

process.schedule = cms.Schedule(process.L1ITMUPath,
                                process.outPath)

f = file("debug.cfg", 'w')
f.write(process.dumpPython())
f.close()
