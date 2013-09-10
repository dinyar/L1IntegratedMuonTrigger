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
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V7A::All', '')

#infile = ['file:SingleMuFlatPt_minusEta_1GeVto200GeV_GEN_SIM_DIGI_L1.root']
#infile.append('file:SingleMuFlatPt_plusEta_1GeVto200GeV_GEN_SIM_DIGI_L1.root')
#infile.append('file:SingleMuFlatPt_plusEta_1GeVto200GeV_GEN_SIM_DIGI_L1_2.root')
#infile.append('file:SingleMuFlatPt_minusEta_1GeVto200GeV_GEN_SIM_DIGI_L1_2.root')


process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.threshold = 'WARNING'
process.MessageLogger.cerr.threshold = 'ERROR'

infile = [
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_100_1_0vY.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_10_1_0ut.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_101_1_5dL.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_102_1_fKp.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_103_1_6Xf.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_104_1_p0h.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_105_1_tTa.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_106_1_UA0.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_107_1_1SZ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_108_1_cXC.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_108_1_SoN.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_110_1_IGA.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_111_1_zkk.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_11_1_run.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_112_1_37U.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_113_1_Nd1.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_114_1_Mag.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_115_1_Ck0.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_116_1_GFa.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_117_1_ti2.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_118_1_6Hb.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_119_1_D5S.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_1_1_96z.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_120_1_VFh.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_121_1_ze6.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_12_1_vGN.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_122_1_TkJ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_123_1_PPB.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_124_1_HbO.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_125_1_gri.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_126_1_G05.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_127_1_9mR.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_128_1_y1a.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_129_1_arA.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_130_1_bVl.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_131_1_2jn.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_13_1_L18.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_132_1_V9h.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_133_1_RST.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_134_1_rC1.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_135_1_xnw.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_136_1_05Y.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_137_1_Igw.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_138_1_88F.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_139_1_4N2.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_141_1_Uua.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_14_1_gG8.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_142_1_O92.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_143_1_aS5.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_144_1_1Hr.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_145_1_2pD.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_146_1_DVi.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_148_1_RSJ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_149_1_u3y.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_150_1_3Mn.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_151_1_pYy.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_15_1_dXP.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_152_1_2k4.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_152_1_dgr.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_154_1_cgz.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_155_1_D4K.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_156_1_IOd.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_159_1_ODQ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_159_1_pj7.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_161_1_2NY.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_16_1_NoQ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_162_1_2NY.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_163_1_aiF.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_165_1_BWl.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_166_1_sp7.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_168_1_f3I.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_169_1_Tm5.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_170_1_vRH.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_171_1_o2g.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_17_1_3i0.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_172_1_tdW.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_174_1_vKw.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_175_1_731.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_176_1_l3W.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_177_1_xZX.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_178_1_4kq.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_179_1_QC7.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_180_1_D2q.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_181_1_kaQ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_18_1_yRS.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_182_1_aLn.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_183_1_X5t.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_184_1_U9n.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_185_1_9d2.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_186_1_PWl.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_187_1_8ZO.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_188_1_7S7.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_189_1_FsW.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_189_1_wTd.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_191_1_fbi.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_19_1_UNw.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_192_1_BMQ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_193_1_h9I.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_194_1_gzu.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_195_1_gG8.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_196_1_RfV.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_197_1_PXe.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_198_1_7ny.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_199_1_F3A.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_200_1_cor.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_201_1_tTX.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_20_1_9dK.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_203_1_JsA.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_203_1_V1Q.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_204_1_nBU.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_205_1_GQ0.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_206_1_Fai.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_207_1_GD2.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_208_1_pUG.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_209_1_isi.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_210_1_eUi.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_21_1_11Z.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_211_1_YJ7.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_212_1_zIy.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_214_1_5H9.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_215_1_yGl.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_216_1_ENb.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_217_1_msE.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_218_1_7W8.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_2_1_CX0.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_220_1_WT0.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_22_1_zsg.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_223_1_8Q8.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_224_1_XJU.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_225_1_QhX.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_227_1_A23.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_228_1_fBJ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_229_1_yvZ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_230_1_eYo.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_23_1_d5x.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_232_1_b6B.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_233_1_ARz.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_234_1_cuH.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_235_1_DVf.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_236_1_anY.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_237_1_KvH.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_238_1_h7J.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_239_1_0eR.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_240_1_ZcO.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_241_1_6kB.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_24_1_1XP.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_242_1_Rwe.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_243_1_3GC.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_244_1_wkS.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_245_1_ANE.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_246_1_QMt.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_247_1_CEh.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_248_1_eG2.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_249_1_VnB.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_250_1_Pfm.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_25_1_T7w.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_26_1_11X.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_26_1_HT2.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_27_1_nkX.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_28_1_6oH.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_29_1_Igt.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_30_1_r9k.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_31_1_xmX.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_3_1_WGA.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_32_1_h1b.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_33_1_1oQ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_34_1_vHY.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_35_1_BwQ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_36_1_GeJ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_37_1_WVC.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_38_1_9Y7.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_39_1_4qk.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_40_1_H4x.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_41_1_y0I.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_4_1_jSt.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_42_1_HXB.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_43_1_C4r.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_44_1_KIz.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_45_1_AJ0.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_46_1_BrH.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_47_1_C1p.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_48_1_RAP.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_49_1_vkY.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_50_1_uho.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_51_1_m2R.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_5_1_2eg.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_52_1_dmi.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_53_1_XHm.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_54_1_lU5.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_55_1_zGf.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_56_1_abB.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_56_1_Odz.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_57_1_U8o.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_58_1_UIG.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_59_1_dp3.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_60_1_4TR.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_61_1_EJQ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_6_1_v4T.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_62_1_aMr.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_63_1_pqL.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_64_1_5oc.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_65_1_XRP.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_66_1_XVH.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_67_1_wGt.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_68_1_rh9.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_69_1_A1t.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_70_1_Oaq.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_71_1_30x.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_7_1_8TG.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_72_1_CXY.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_73_1_qqm.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_74_1_IUH.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_75_1_hw3.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_76_1_TsB.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_77_1_H4z.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_78_1_km5.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_79_1_JzK.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_80_1_nph.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_81_1_Nrl.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_8_1_vOV.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_82_1_WGV.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_83_1_jQq.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_84_1_UZK.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_85_1_igQ.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_86_1_SGc.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_87_1_dOR.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_88_1_3fx.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_89_1_GYq.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_90_1_LZv.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_91_1_Gfo.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_9_1_tzs.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_92_1_n1j.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_93_1_EPK.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_94_1_VJP.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_95_1_y3l.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_96_1_QDB.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_97_1_yOI.root',
'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/Generation/SingleMuFlatPt_3GeVto140GeV_GEN_SIM_DIGI_L1_99_1_Sw0.root'
#'file:/afs/cern.ch/user/b/battilan/work/public/SingleMuFlatPt_3GeVto100GeV_GEN_SIM_DIGI_L1.root'
]

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

process.L1ITMUSequence = cms.Sequence( process.L1ITMuTriggerPrimitives +
                                       process.L1CSCTFTrackConverter   +
                                       process.L1DTTFTrackConverter    +
                                       process.L1RPCTFTrackConverters  +
                                       process.MBLTProducer            +
                                       process.L1ITMuonBarrelPrimitiveProducer +
                                       process.MBTracksProducer        )

process.L1ITMUPath = cms.Path(process.L1ITMUSequence)

outCommands = cms.untracked.vstring('drop *')
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
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = outCommands,
    fileName = cms.untracked.string('/data2/battilan/L1Trigger/BarrelDTTFUpgrade/L1BarrelMuon/L1ITMBLT_100K.root'),
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
