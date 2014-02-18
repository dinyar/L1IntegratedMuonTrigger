import FWCore.ParameterSet.Config as cms

L1ITMuonBarrelPrimitiveProducer = cms.EDProducer(
    'L1ITMuonBarrelPrimitiveProducer',
    MBLTCollection = cms.InputTag('MBLTProducer'),
    xDtResol = cms.double( 1 ),
    xRpcResol = cms.double( 20 ),
    phibDtCorrResol = cms.double( 0.005 ),
    phibDtUnCorrResol = cms.double( 0.04 ),
    qualityRemappingMode = cms.int32( 2 )
    )
