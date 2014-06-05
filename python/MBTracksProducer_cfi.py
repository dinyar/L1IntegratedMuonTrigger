import FWCore.ParameterSet.Config as cms

from L1Trigger.DTTrackFinder.dttfDigis_cfi import dttfDigis as tmpDttfDigis  

MBTracksProducer = cms.EDProducer(
    'MBTracksProducer',
    DTTrackSrc = cms.InputTag('simDttfDigis','DTTF'),
    gmtDigis = cms.InputTag('simGmtDigis',''),
    MBLTCollectionSrc = cms.InputTag('MBLTProducer'),
    BX_min = cms.int32(tmpDttfDigis.BX_min.value()),
    BX_max = cms.int32(tmpDttfDigis.BX_max.value()),
    MaxDeltaPhiGmtIn = cms.double( 1e-3 ),
    MaxDeltaPhiGmtOut = cms.double( 1e-3 )
    )
    
