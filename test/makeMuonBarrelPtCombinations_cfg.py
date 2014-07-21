import FWCore.ParameterSet.Config as cms

process = cms.Process('PLOTMAKER')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.TFileService = cms.Service(
    "TFileService",
  
     fileName=cms.string('L1ITMuonBarrelPtStudies2.root')
    )

process.L1ITMuPtPlotter = cms.EDAnalyzer(
    'L1ITMBPtComb',
    MBTracksCollection = cms.InputTag('MBTracksProducer'),
    GenParticlesCollection = cms.InputTag('genParticles')
    )

infile = 'file:L1ITMBLT.root'
#infile = '/store/caf/user/battilan/L1Trigger/BarrelTF/L1ITBarrelMuon/v1/L1ITMBLT.root'

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    )

process.L1ITMUSequence = cms.Sequence(process.L1ITMuPtPlotter)

process.L1ITMUPath = cms.Path(process.L1ITMUSequence)

process.schedule = cms.Schedule(process.L1ITMUPath)

