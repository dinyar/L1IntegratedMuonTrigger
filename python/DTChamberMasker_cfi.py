import FWCore.ParameterSet.Config as cms

DTChamberMasker = cms.EDProducer('DTChamberMasker',
                                 digiTag = cms.InputTag('preDtDigis'),
                                 triggerPrimPhTag = cms.InputTag('preDttfDigis'),
                                 triggerPrimThTag = cms.InputTag('preDttfDigis'),
                                 doTriggerFromDDU = cms.bool(True),
                                 maskedChRegEx = cms.vstring("WH0_ST1_SEC1",
                                                             "WH0_ST1_SEC3",
                                                             "WH0_ST1_SEC5",
                                                             "WH0_ST1_SEC7",
                                                             "WH0_ST1_SEC9",
                                                             "WH0_ST1_SEC11")
)
