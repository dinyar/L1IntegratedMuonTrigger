import FWCore.ParameterSet.Config as cms

from L1Trigger.DTTrackFinder.dttfDigis_cfi import dttfDigis as tmpDttfDigis

L1DTTFTrackConverter = cms.EDProducer(
    'L1DTTFTrackConverter',
    DTTrackSrc = cms.InputTag('simDttfDigis','DTTF'),
    TriggerPrimitiveSrc = cms.InputTag('L1ITMuTriggerPrimitives',''),
    BX_min = cms.int32(tmpDttfDigis.BX_min.value()),
    BX_max = cms.int32(tmpDttfDigis.BX_max.value()) 
    )
