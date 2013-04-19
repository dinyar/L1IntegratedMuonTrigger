import FWCore.ParameterSet.Config as cms

L1ITMuDtSuperStationProducer = cms.EDProducer(
    'L1ITMuDtSuperStationProducer',
    TriggerPrimitiveSrc = cms.InputTag('L1ITMuTriggerPrimitives'),
    MaxDeltaPhi = cms.double( 2 )
    )
