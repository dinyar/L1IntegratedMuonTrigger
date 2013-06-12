import FWCore.ParameterSet.Config as cms

from L1Trigger.DTTrackFinder.dttfDigis_cfi import dttfDigis

MBTracksProducer = cms.EDProducer(
    'MBTracksProducer',
    DTTrackSrc = cms.InputTag('simDttfDigis','DTTF'),
    gmtDigis = cms.InputTag('simGmtDigis',''),
    MBLTCollectionSrc = cms.InputTag('MBLTProducer'),
    BX_min = cms.int32(dttfDigis.BX_min.value()),
    BX_max = cms.int32(dttfDigis.BX_max.value()),
    MaxDeltaPhi = cms.double( 1e-3 ),
    MaxDeltaR   = cms.double( 0.3 )
    )
    