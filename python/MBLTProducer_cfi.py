import FWCore.ParameterSet.Config as cms

MBLTProducer = cms.EDProducer(
    'MBLTProducer',
    TriggerPrimitiveSrc = cms.InputTag('L1ITMuTriggerPrimitives'),
    MaxDeltaPhi = cms.double( 0.05 )
    )
