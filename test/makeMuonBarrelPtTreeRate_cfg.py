import FWCore.ParameterSet.Config as cms

process = cms.Process('PLOTMAKER')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.TFileService = cms.Service(
    "TFileService",

   fileName=cms.string('L1ITMuonBarrelPtTreeRate.root') 
    )

process.L1ITMuPtPlotter = cms.EDAnalyzer(
    'L1ITMBPtTreeRate',
    MBTracksCollection = cms.InputTag('MBTracksProducer'),
    GenParticlesCollection = cms.InputTag('genParticles')
    )

infile1 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_100_1_EST.root'
infile2 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_101_1_Hue.root'
infile3 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_101_1_TpF.root'
infile4 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_102_1_TLM.root'
infile5 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_102_1_uxv.root'
infile6 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_103_1_Rcg.root'
infile7 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_103_1_uVX.root'
infile8 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_104_1_NtU.root'
infile9 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_105_1_19l.root'
infile10 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_107_1_nP4.root'
infile11 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_108_1_GXk.root'
infile12 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_109_1_Yvy.root'
infile13 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_10_1_NJG.root'
infile14 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_110_1_gR8.root'
infile15 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_111_1_FfN.root'
infile16 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_11_1_g6T.root'
infile17 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_127_1_w39.root'
infile18 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_128_1_dUi.root'
infile19 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_129_1_4pV.root'
infile20 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_12_1_hBI.root'
infile21 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_130_1_uIB.root'
infile22 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_131_1_rFx.root'
infile23 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_13_1_S87.root'
infile24 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_14_1_gKc.root'
infile25 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_15_1_8V1.root'
infile26 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_16_1_LAu.root'
infile27 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_17_1_icI.root'
infile28 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_18_1_Og9.root'
infile29 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_1_1_P8z.root'
infile30 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_20_1_6YG.root'
infile31 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_22_1_ROg.root'
infile32 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_23_1_3rd.root'
infile33 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_24_1_PTb.root'
infile34 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_25_1_bX4.root'
infile35 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_26_1_zjN.root'
infile36 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_27_1_TVF.root'
infile37 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_28_1_AHz.root'
infile38 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_29_1_CMF.root'
infile39 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_2_1_9JS.root'
infile40 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_30_1_pj6.root'
infile41 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_31_1_6TR.root'
infile42 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_9_1_n7n.root'
infile43 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_32_1_0Qy.root'
infile44 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_33_1_ZxI.root'
infile45 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_35_1_FIU.root'
infile46 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_36_1_Uun.root'
infile47 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_37_1_rcm.root'
infile48 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_38_1_Em8.root'
infile49 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_38_1_aQm.root'
infile50 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_39_1_bo6.root'
infile51 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_39_1_i9D.root'
infile52 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_3_1_Fcs.root'
infile53 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_40_1_9dN.root'
infile54 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_40_1_I5E.root'
infile55 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_41_1_BQa.root'
infile56 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_41_1_KNo.root'
infile57 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_42_1_BM4.root'
infile58 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_42_1_Vgv.root'
infile59 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_43_1_YwF.root'
infile60 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_44_1_POc.root'
infile61 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_45_1_qxg.root'
infile62 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_46_1_MsJ.root'
infile63 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_47_1_vIZ.root'
infile64 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_48_1_iY7.root'
infile65 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_49_1_DOW.root'
infile66 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_4_1_jdU.root'
infile67 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_50_1_UmE.root'
infile68 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_51_1_5t5.root'
infile69 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_53_1_Uw6.root'
infile70 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_54_1_n6C.root'
infile71 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_55_1_bFj.root'
infile72 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_56_1_DpH.root'
infile73 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_57_1_KZy.root'
infile74 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_57_1_y6j.root'
infile75 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_58_1_ROo.root'
infile76 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_58_1_Ucm.root'
infile77 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_59_1_XTh.root'
infile78 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_59_1_lI9.root'
infile79 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_5_1_zTs.root'
infile80 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_60_1_2ZQ.root'
infile81 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_60_1_vlN.root'
infile82 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_61_1_eQd.root'
infile83 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_62_1_K7h.root'
infile84 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_62_1_w5w.root'
infile85 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_63_1_lG2.root'
infile86 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_64_1_GIK.root'
infile87 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_65_1_wu1.root'
infile88 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_66_1_Z3o.root'
infile89 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_67_1_Q8W.root'
infile90 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_67_1_gcy.root'
infile91 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_68_1_ZLP.root'
infile92 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_68_1_fvo.root'
infile93 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_69_1_WWj.root'
infile94 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_69_1_aJz.root'
infile95 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_6_1_L3o.root'
infile96 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_70_1_91X.root'
infile97 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_70_1_xO6.root'
infile98 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_71_1_dVD.root'
infile99 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_71_1_hLZ.root'
infile100 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_72_1_4XF.root'
infile101 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_72_1_4vb.root'
infile102 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_73_1_Q0h.root'
infile103 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_74_1_Ilx.root'
infile104 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_75_1_cED.root'
infile105 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_76_1_t8Q.root'
infile106 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_77_1_WGb.root'
infile107 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_78_1_eAX.root'
infile108 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_79_1_Hfg.root'
infile109 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_79_1_rDm.root'
infile110 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_7_1_2XM.root'
infile111 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_80_1_M9j.root'
infile112 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_81_1_8Eh.root'
infile113 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_81_1_txG.root'
infile114 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_82_1_8TE.root'
infile115 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_82_1_jzv.root'
infile116 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_83_1_xul.root'
infile117 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_83_1_z4D.root'
infile118 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_84_1_OFY.root'
infile119 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_84_1_QZg.root'
infile120 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_85_1_tS3.root'
infile121 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_86_1_l6X.root'
infile122 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_87_1_iuC.root'
infile123 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_88_1_tgM.root'
infile124 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_99_1_2HA.root'
infile125 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_89_1_ajg.root'
infile126 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_8_1_mu1.root'
infile127 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_90_1_X07.root'
infile128 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_91_1_fCf.root'
infile129 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_92_1_9Ah.root'
infile130 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_92_1_jzZ.root'
infile131 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_93_1_J3H.root'
infile132 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_93_1_T1h.root'
infile133 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_94_1_XTD.root'
infile134 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_94_1_mia.root'
infile135 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_95_1_LoQ.root'
infile136 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_95_1_wWA.root'
infile137 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_96_1_6ID.root'
infile138 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_97_1_1vU.root'
infile139 = '/store/user/battilan/L1Trigger/NewBarrelTrigger/v1/ZeroBias/L1ITMBLT_98_1_1th.root'

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile1,infile2,infile3,infile4,infile5,infile6,infile7,infile8,infile9,infile10,infile11,infile12,infile13,infile14,infile15,infile16,infile17,infile18,infile19,infile20,infile21,infile22,infile23,infile24,infile25,infile26,infile27,infile28,infile29,infile30,infile31,infile32,infile33,infile34,infile35,infile36,infile37,infile38,infile39,infile40, infile41,infile42,infile43,infile44,infile45,infile46,infile47,infile48,infile49,infile50,infile51,infile52,infile53,infile54,infile55,infile56,infile57,infile58,infile59,infile60,infile61,infile62,infile63,infile64,infile65,infile66,infile67,infile68,infile69,infile70,infile71,infile72,infile73,infile74,infile75,infile76,infile77,infile78,infile79,infile80,infile81,infile82,infile83,infile84,infile85,infile86,infile87,infile88,infile89,infile90,infile91,infile92,infile93,infile94,infile95,infile96,infile97,infile98,infile99,infile100,infile101,infile102,infile103,infile104,infile105,infile106,infile107,infile108,infile109,infile110,infile111,infile112,infile113,infile114,infile115,infile116,infile117,infile118,infile119,infile120,infile121,infile122,infile123,infile124,infile125,infile126,infile127,infile128,infile129,infile130,infile131,infile132,infile133,infile134,infile135,infile136,infile137,infile138,infile139)
)    
process.L1ITMUSequence = cms.Sequence(process.L1ITMuPtPlotter)

process.L1ITMUPath = cms.Path(process.L1ITMUSequence)

process.schedule = cms.Schedule(process.L1ITMUPath)

