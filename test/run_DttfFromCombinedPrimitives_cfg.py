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

process.load("RecoLocalMuon.Configuration.RecoLocalMuon_cff")
process.dt1DRecHits.dtDigiLabel = 'simMuonDTDigis'

#process.load('L1Trigger.DTTrackFinder.dttfDigis_cfi')
process.load('L1Trigger.GlobalMuonTrigger.gmtDigis_cfi')
from  L1Trigger.DTTrackFinder.dttfDigis_cfi import dttfDigis as tmpDttfDigis
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START72_V1::All', '')

process.phase1DttfDigis = process.tmpDttfDigis.clone(
#    DTDigi_Source = cms.InputTag("simDtTriggerPrimitiveDigis"),
    DTDigi_Source = cms.InputTag("L1ITMuonBarrelPrimitiveProducer"),
    DTDigi_Theta_Source = cms.InputTag("simDtTriggerPrimitiveDigis"),
    CSCStub_Source = cms.InputTag("simCsctfTrackDigis")
    )

process.phase1GmtDigis = process.gmtDigis.clone(
    DTCandidates = cms.InputTag("phase1DttfDigis","DT"),
    RPCbCandidates = cms.InputTag("simRpcTriggerDigis","RPCb"),
    #CSCCandidates = cms.InputTag("l1CscTfEmulDigis","CSC"),
    RPCfCandidates = cms.InputTag("simRpcTriggerDigis","RPCf")
    #MipIsoData = cms.InputTag("L1RCTRegionSumsEmCands"),
    )

process.phase1L1ITMuTriggerPrimitives = process.L1ITMuTriggerPrimitives.clone(
      DT = cms.PSet( collectorType = cms.string('DTCollector'),
                     src = cms.InputTag('L1ITMuonBarrelPrimitiveProducer'),
                     src2 = cms.InputTag('simDtTriggerPrimitiveDigis'),
                     BX_min = cms.int32(-10),
                     BX_max = cms.int32(10),
                     runBunchCrossingCleaner = cms.bool(False),
                     bxCleanerCfg = cms.PSet(bxWindowSize = cms.int32(1)) )
      )
process.phase1MBLTProducer = process.MBLTProducer.clone(
    TriggerPrimitiveSrc = cms.InputTag('phase1L1ITMuTriggerPrimitives')
    )
process.phase1MBTracksProducer = process.MBTracksProducer.clone(
    DTTrackSrc = cms.InputTag('phase1DttfDigis','DTTF'),
    MBLTCollectionSrc = cms.InputTag('phase1MBLTProducer'),
    gmtDigis = cms.InputTag('phase1GmtDigis')
    )

process.L1ITMuonBarrelPrimitiveProducer.qualityRemappingMode = cms.int32(2)
process.L1ITMuonBarrelPrimitiveProducer.is7QualityCodes = False

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'WARNING'
process.MessageLogger.cerr.threshold = 'ERROR'

process.dtLocalRecoSequence = cms.Sequence(process.dt1DRecHits + process.dt4DSegments)

process.L1ITMUSequence = cms.Sequence( process.L1ITMuTriggerPrimitives +  # L1IT primitives container
                                       process.L1CSCTFTrackConverter   +  
                                       process.L1DTTFTrackConverter    +
                                       process.L1RPCTFTrackConverters  +
                                       process.MBLTProducer            +  # chamber level primitive matching
                                       process.MBTracksProducer        +  # track-to-primitive matching
                                       process.L1ITMuonBarrelPrimitiveProducer + # new combined primitive prod.
                                       process.phase1L1ITMuTriggerPrimitives   + # make new TP from MBPrimitive
                                       process.phase1MBLTProducer              + # make new LT
                                       process.phase1DttfDigis                 + # make new DTTF
                                       process.phase1GmtDigis                  + # make new GMT
                                       process.phase1MBTracksProducer            # make new matching
                                       )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.L1ITMUPath = cms.Path(process.dtLocalRecoSequence + process.L1ITMUSequence)

outCommands = cms.untracked.vstring('keep *') #note keep *
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
outCommands.append('keep *_myDttfDigis_*_*')
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

infile = ['/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_1_1_96z.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_46_1_BrH.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_47_1_C1p.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_48_1_RAP.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_49_1_vkY.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_4_1_jSt.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_50_1_uho.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_51_1_m2R.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_52_1_dmi.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_53_1_XHm.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_54_1_lU5.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_55_1_zGf.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_56_1_Odz.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_56_1_abB.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_57_1_U8o.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_58_1_UIG.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_59_1_dp3.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_5_1_2eg.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_60_1_4TR.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_61_1_EJQ.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_62_1_aMr.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_63_1_pqL.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_64_1_5oc.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_65_1_XRP.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_66_1_XVH.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_67_1_wGt.root',
'/store/caf/user/battilan/L1Trigger/BarrelTF/Generation/v1/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_68_1_rh9.root']
notthisoneinfile=['root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_100_1_EST.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_101_1_Hue.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_101_1_TpF.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_102_1_TLM.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_102_1_uxv.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_103_1_Rcg.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_103_1_uVX.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_104_1_NtU.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_105_1_19l.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_107_1_nP4.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_108_1_GXk.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_109_1_Yvy.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_10_1_NJG.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_110_1_gR8.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_111_1_FfN.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_11_1_g6T.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_127_1_w39.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_128_1_dUi.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_129_1_4pV.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_12_1_hBI.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_130_1_uIB.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_131_1_rFx.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_13_1_S87.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_14_1_gKc.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_15_1_8V1.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_16_1_LAu.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_17_1_icI.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_18_1_Og9.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_1_1_P8z.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_20_1_6YG.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_22_1_ROg.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_23_1_3rd.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_24_1_PTb.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_25_1_bX4.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_26_1_zjN.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_27_1_TVF.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_28_1_AHz.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_29_1_CMF.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_2_1_9JS.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_30_1_pj6.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_31_1_6TR.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_32_1_0Qy.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_33_1_ZxI.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_35_1_FIU.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_36_1_Uun.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_37_1_rcm.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_38_1_Em8.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_38_1_aQm.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_39_1_bo6.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_39_1_i9D.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_3_1_Fcs.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_40_1_9dN.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_40_1_I5E.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_41_1_BQa.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_41_1_KNo.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_42_1_BM4.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_42_1_Vgv.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_43_1_YwF.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_44_1_POc.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_45_1_qxg.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_46_1_MsJ.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_47_1_vIZ.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_48_1_iY7.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_49_1_DOW.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_4_1_jdU.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_50_1_UmE.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_51_1_5t5.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_53_1_Uw6.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_54_1_n6C.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_55_1_bFj.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_56_1_DpH.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_57_1_KZy.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_57_1_y6j.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_58_1_ROo.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_58_1_Ucm.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_59_1_XTh.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_59_1_lI9.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_5_1_zTs.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_60_1_2ZQ.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_60_1_vlN.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_61_1_eQd.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_62_1_K7h.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_62_1_w5w.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_63_1_lG2.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_64_1_GIK.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_65_1_wu1.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_66_1_Z3o.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_67_1_Q8W.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_67_1_gcy.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_68_1_ZLP.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_68_1_fvo.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_69_1_WWj.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_69_1_aJz.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_6_1_L3o.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_70_1_91X.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_70_1_xO6.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_71_1_dVD.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_71_1_hLZ.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_72_1_4XF.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_72_1_4vb.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_73_1_Q0h.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_74_1_Ilx.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_75_1_cED.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_76_1_t8Q.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_77_1_WGb.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_78_1_eAX.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_79_1_Hfg.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_79_1_rDm.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_7_1_2XM.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_80_1_M9j.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_81_1_8Eh.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_81_1_txG.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_82_1_8TE.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_82_1_jzv.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_83_1_xul.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_83_1_z4D.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_84_1_OFY.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_84_1_QZg.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_85_1_tS3.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_86_1_l6X.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_87_1_iuC.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_88_1_tgM.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_89_1_ajg.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_8_1_mu1.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_90_1_X07.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_91_1_fCf.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_92_1_9Ah.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_92_1_jzZ.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_93_1_J3H.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_93_1_T1h.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_94_1_XTD.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_94_1_mia.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_95_1_LoQ.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_95_1_wWA.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_96_1_6ID.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_97_1_1vU.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_98_1_1th.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_99_1_2HA.root',
'root://eoscms.cern.ch//eos/cms/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_9_1_n7n.root']



process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    )


f = file("debug.cfg", 'w')
f.write(process.dumpPython())
f.close()
