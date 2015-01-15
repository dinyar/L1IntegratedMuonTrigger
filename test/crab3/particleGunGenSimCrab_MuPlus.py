# Auto generated configuration file
# using: 
# Revision: 1.381.2.7 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: SingleMuFlatLogPt_100MeVto2TeV_cfi.py -s GEN,SIM,DIGI,L1 --conditions START72_V1::All --eventcontent FEVTDEBUG --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('L1')

## Configuration parameters

# The GlobalTag 
# must be defined according to the release and the configuration you need
# this is release/geometry depended, what is now here has been tested to make the
# workflow work in CMSSW_7_2_3_patch1 using 2012 geomery and was not used for 
# big sample production, only technical workflow was tested
globalTag = "START72_V1"

# The eta range for GEN muon production
# the present DTTF goes up to |eta|<1.04, the BarrelTF will go roughly up to
# |eta|<0.9, here putting 1.05 as limit to include scattering of muons 
# before reaching the muon chambers
minEta = -1.05
maxEta =  1.05

# The phi range for GEN muon production
# presently set to study the performance of one single sector plus neighbours
minPhi = - 3.14159265359/6.
maxPhi =   3.14159265359/6.

# The pT range for GEN muon production
# presently set using limits of pt for muons to reach the barrel spectrometer
# and the present DTTF pT scale overflow
minPt = 3
maxPt = 140

# The sign of the muon
# -1 for mu- and +1 for mu+
muonCharge = +1

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

if muonCharge > 0 :
	chargeTag='Plus'
else :
        chargeTag='Minus'

configTag = 'SingleMu' + chargeTag + '_FlatPt_' + str(minPt) + 'to' + str(maxPt) \
    + '_eta' + str(minEta) + 'to' + str(maxEta) + '_phi' + str(int(minPhi*100)/100.) \
    + 'to' + str(int(maxPhi*100)/100.)
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.7 $'),
    annotation = cms.untracked.string(configTag),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition
outname = configTag.replace( '.', 'p' ) + '_GEN_SIM_DIGI_L1.root'
process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string( outname ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag+'::All', '')

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
	PGunParameters = cms.PSet(
        MinPt = cms.double(minPt),
	MaxPt = cms.double(maxPt),
        PartID = cms.vint32(-13 * muonCharge),        
        MaxPhi = cms.double(maxPhi),
	MinPhi = cms.double(minPhi),
	MaxEta = cms.double(maxEta),
        MinEta = cms.double(minEta)        
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single mu pt ' + str(minPt) + 'to' + str(maxPt)),
    AddAntiParticle = cms.bool(False), #need *single* muons dammit
    firstRun = cms.untracked.uint32(1)
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

outCommands = cms.untracked.vstring('drop *')
outCommands.append('keep *_simMuonDTDigis_*_*')
outCommands.append('keep *_simMuonRPCDigis_*_*')
outCommands.append('keep *_simMuonCSCDigis_*_*')
outCommands.append('keep *_genParticles_*_*')
outCommands.append('keep *_simCsctfTrackDigis_*_*')
outCommands.append('keep *_simDttfDigis_*_*')
outCommands.append('keep *_simGmtDigis_*_*')
outCommands.append('keep *_simRpcTriggerDigis_*_*')
outCommands.append('keep *_simDtTriggerPrimitiveDigis_*_*')
outCommands.append('keep *_simCscTriggerPrimitiveDigis_*_*')
outCommands.append('keep *_L1ITMuTriggerPrimitives_*_*')
outCommands.append('keep *_MBLTProducer_*_*')
outCommands.append('keep *_MBTracksProducer_*_*')
outCommands.append('keep *_L1ITMuonBarrelPrimitiveProducer_*_*')
outCommands.append('keep *_*Converter_*_*')
outCommands.append('keep *_*Matcher_*_*')

process.FEVTDEBUGoutput.outputCommands = outCommands

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,
				process.genfiltersummary_step,
				process.simulation_step,
				process.digitisation_step,
				process.L1simulation_step,
				process.endjob_step,
				process.FEVTDEBUGoutput_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

