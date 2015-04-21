# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: miniAOD-prod -s PAT --eventcontent MINIAODSIM --runUnscheduled --fast --filein root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_AODSIM_100.root --conditions auto:run2_mc --no_exec -n 10
import FWCore.ParameterSet.Config as cms

process = cms.Process('PAT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
#process.load('FastSimulation.PileUpProducer.PileUpSimulator_NoPileUp_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator_E13TeV_AVE_20_inTimeOnly_cff')
process.load('FastSimulation.Configuration.Geometries_MC_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_101.root',
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_102.root',
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_103.root',
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_104.root',
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_105.root',
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_106.root',
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_107.root',
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_108.root',
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_109.root',
        'root://eoscms//eos/cms/store/user/taroni/723/ggHiggsETau/GGHetau_8TeVBeamSpot_20PU_AODSIM_110.root',
        )
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('miniAOD-prod nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    fileName = cms.untracked.string('GGHetau_8TeVBeamSpot_20PU_miniAOD.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

# Additional output definition

# Other statements
process.mix.playback= True
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# customisation of the process.

# Automatic addition of the customisation function from FastSimulation.Configuration.MixingModule_Full2Fast
from FastSimulation.Configuration.MixingModule_Full2Fast import setVertexGeneratorPileUpProducer 

#call to customisation function setVertexGeneratorPileUpProducer imported from FastSimulation.Configuration.MixingModule_Full2Fast
process = setVertexGeneratorPileUpProducer(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
process.load('Configuration.StandardSequences.PATMC_cff')

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff
from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff import miniAOD_customizeMETFiltersFastSim 

#call to customisation function miniAOD_customizeMETFiltersFastSim imported from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff
process = miniAOD_customizeMETFiltersFastSim(process)

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions
