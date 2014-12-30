import FWCore.ParameterSet.Config as cms

L1ITMuonBarrelPrimitiveProducer = cms.EDProducer(
    'L1ITMuonBarrelPrimitiveProducer',
    MBLTCollection = cms.InputTag('MBLTProducer'),
    xDtResol = cms.double( 0.1 ), # cm
    xRpcResol = cms.double( 2 ), # cm
    phibDtCorrResol = cms.double( 0.005 ), # rad
    phibDtUnCorrResol = cms.double( 0.04 ), # rad 
    qualityRemappingMode = cms.int32( 2 ),
    is7QualityCodes = cms.bool( True ),
    useRpcBxForDtBelowQuality = cms.int32(8)
    )
