import FWCore.ParameterSet.Config as cms

from L1Trigger.L1IntegratedMuonTrigger.DTChamberMasker_cfi import DTChamberMasker

def appendChamberMaskerAtUnpacking(process, doDigis, doTrigger, chambRegEx):

    if doDigis and hasattr(process,'muonDTDigis') :

        print "[appendChamberMasker] : Found muonDTDigis, applying filter"

        process.preDtDigis = process.muonDTDigis.clone()
        process.muonDTDigis = DTChamberMasker.clone()
        if len(chambRegEx) > 0 :
            process.dtDigis.maskedChRegEx = chambRegEx
        process.muonDTDigis.triggerPrimPhTag = cms.InputTag('')
        process.muonDTDigis.triggerPrimThTag = cms.InputTag('')
        process.filteredDigiSequence = cms.Sequence(process.preDtDigis + process.muonDTDigis)
        process.RawToDigi.replace(process.muonDTDigis, process.filteredDigiSequence)

    if doDigis and hasattr(process,'dttfDigis') :

        print "[appendChamberMasker] : Found dttfDigis, applying filter"

        process.preDttfDigis = process.dttfDigis.clone()
        process.dttfDigis = DTChamberMasker.clone()
        if len(chambRegEx) > 0 :
            process.dtDigis.maskedChRegEx = chambRegEx
        process.dttfDigis.digiTag = cms.InputTag('')
        process.filteredTrigSequence = cms.Sequence(process.preDttfDigis + process.dttfDigis)
        process.RawToDigi.replace(process.dttfDigis, process.filteredTrigSequence)

        
