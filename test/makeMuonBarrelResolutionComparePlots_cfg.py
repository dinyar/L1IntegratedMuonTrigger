import FWCore.ParameterSet.Config as cms

process = cms.Process('PLOTMAKER')

process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Geometry.DTGeometryBuilder.idealForDigiDtGeometry_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START72_V1::All', '')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('RecoLocalMuon.Configuration.RecoLocalMuon_cff')
process.dt1DRecHits.dtDigiLabel = 'simMuonDTDigis'

process.TFileService = cms.Service(
	'TFileService',
	fileName=cms.string('L1ITMuonBarrelResolutionComparePlots.root')
)

process.L1ITMuResolutionPlotter = cms.EDAnalyzer(
	'L1ITMBPrimitiveResolutionComparePlots',
	# labels of DDU/DCC data and 4D segments
	inputTagDCC_old = cms.untracked.InputTag('simDtTriggerPrimitiveDigis'),
	inputTagDCC_new = cms.untracked.InputTag('L1ITMuonBarrelPrimitiveProducer'), # superprimitives

	inputTagSEG = cms.untracked.InputTag('dt4DSegments'),
	# set outflows to boudaries
	rebinOutFlowsInGraph = cms.untracked.bool(True),
	# label of the geometry used to feed DTTrigGeomUtils
	geomLabel = cms.untracked.string('idealForDigi'),
)

infile = 'L1ITMBLT.root'

process.source = cms.Source(
	'PoolSource',
	fileNames = cms.untracked.vstring(infile)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.L1ITMUSequence = cms.Sequence(process.L1ITMuResolutionPlotter)
process.dtLocalRecoSequence = cms.Sequence(process.dt1DRecHits + process.dt4DSegments)


process.L1ITMUPath = cms.Path(process.dtLocalRecoSequence + process.L1ITMUSequence)

process.schedule = cms.Schedule(process.L1ITMUPath)

