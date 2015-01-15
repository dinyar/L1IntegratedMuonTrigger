import FWCore.ParameterSet.Config as cms

process = cms.Process('PLOTMAKER')

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Geometry.DTGeometryBuilder.idealForDigiDtGeometry_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START72_V1::All', '')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("RecoLocalMuon.Configuration.RecoLocalMuon_cff")
process.dt1DRecHits.dtDigiLabel = 'simMuonDTDigis'

process.TFileService = cms.Service(
    "TFileService",
    #fileName=cms.string('L1ITMuonBarrelTimingPlots_befLuigi.root')
    #fileName=cms.string('L1ITMuonBarrelTimingPlots_afLuigi.root')
    fileName=cms.string('L1ITMuonBarrelTimingEffPlots_17aprilT0.root')
    )

process.L1ITMuTimingPlotter = cms.EDAnalyzer(
    'L1ITMBPrimitiveTimingEffPlots_T0',
    # labels of DDU/DCC data and 4D segments
    inputTagDCC = cms.untracked.InputTag("simDtTriggerPrimitiveDigis"),
    # replace the previous with the following line to run using "super primitives" as input 
 #   inputTagDCC = cms.untracked.InputTag("L1ITMuonBarrelPrimitiveProducer"),
    inputTagSEG = cms.untracked.InputTag("dt4DSegments"),
    # set outflows to boudaries
    rebinOutFlowsInGraph = cms.untracked.bool(True),
    # label of the geometry used to feed DTTrigGeomUtils
    geomLabel = cms.untracked.string("idealForDigi"),
    )

#infile = 'file:L1ITMBLT_befLuigi17feb2014.root'
infile = 'file:L1ITMBLT_afLuigi.root'

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    )

process.L1ITMUSequence = cms.Sequence(process.L1ITMuTimingPlotter)
process.dtLocalRecoSequence = cms.Sequence(process.dt1DRecHits + process.dt4DSegments)


process.L1ITMUPath = cms.Path(process.dtLocalRecoSequence + process.L1ITMUSequence)

process.schedule = cms.Schedule(process.L1ITMUPath)
