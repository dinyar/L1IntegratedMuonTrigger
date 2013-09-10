import FWCore.ParameterSet.Config as cms

process = cms.Process('PLOTMAKER')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string('L1ITMuonBarrelPtPlots.root')
    )

process.L1ITMuPtPlotter = cms.EDAnalyzer(
    'L1ITMBPtLutPlots',
    MBTracksCollection = cms.InputTag('MBTracksProducer'),
    GenParticlesCollection = cms.InputTag('genParticles')
    )

infile = 'file:/data2/battilan/L1Trigger/BarrelDTTFUpgrade/L1BarrelMuon/L1ITMBLT.root'

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.L1ITMUSequence = cms.Sequence(process.L1ITMuPtPlotter)

process.L1ITMUPath = cms.Path(process.L1ITMUSequence)

process.schedule = cms.Schedule(process.L1ITMUPath)
