import FWCore.ParameterSet.Config as cms

process = cms.Process('TEXTDUMP')

process.L1Tsubprime = cms.EDProducer(
    'SubsystemPrimitiveProducer',
    doGen = cms.untracked.bool(True),
    genSrc = cms.InputTag("genParticles"),
    trackSrcs = cms.VInputTag(
    cms.InputTag('L1ITMuSimpleDeltaEtaHitMatcher')
    )
)

infile = 'file:/afs/cern.ch/work/b/battilan/public/SingleMuFlatPt_3GeVto100GeV_GEN_SIM_DIGI_L1.root'

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    )

process.L1ITMUSequence = cms.Sequence(process.L1ITMuPlotter)

process.L1ITMUPath = cms.Path(process.L1ITMUSequence)

process.schedule = cms.Schedule(process.L1ITMUPath)
