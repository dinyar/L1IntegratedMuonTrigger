import FWCore.ParameterSet.Config as cms

process = cms.Process('TEXTDUMP')

process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string('L1ITMuDtSuperStation.root')
    )

process.L1ITMuPlotter = cms.EDAnalyzer(
    'L1ITMuDtSuperStationPlots',
    L1ITMuDtSuperStation = cms.InputTag('L1ITMuDtSuperStationProducer')
)

infile = 'file:L1ITMU.root'

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    )

process.L1ITMUSequence = cms.Sequence(process.L1ITMuPlotter)

process.L1ITMUPath = cms.Path(process.L1ITMUSequence)

process.schedule = cms.Schedule(process.L1ITMUPath)
