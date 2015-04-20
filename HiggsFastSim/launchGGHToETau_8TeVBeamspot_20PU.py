#!/bin/env python


from os import popen

#numero iniziale
i=100

#inizio ciclo
while i<110:


    #inizio file cfg
    cfg = """
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator_E13TeV_AVE_20_inTimeOnly_cff')
process.load('FastSimulation.Configuration.Geometries_MC_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('HLTrigger.Configuration.HLT_GRun_Famos_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000)
)

# Input source
process.source = cms.Source("EmptySource",
               firstLuminosityBlock=cms.untracked.uint32("""+str(1+i)+"""),
               firstEvent=cms.untracked.uint32("""+str(1+1000*i)+""")

)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('GGToHtautau_13TeV_pythia8_cff nevts:50'),
    name = cms.untracked.string('Applications')
)

# Output definition


process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.AODSIMEventContent.outputCommands,
    fileName = cms.untracked.string('GGHetau_8TeVBeamSpot_20PU_AODSIM_"""+str(1+i)+""".root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

# Additional output definition

# Other statements


process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.HLTEndSequence = cms.Sequence( process.dummyModule )
process.simulation = cms.Sequence( process.dummyModule )
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring('pythia8CommonSettings',
            'pythia8CUEP8M1Settings',
            'processParameters'),
        processParameters = cms.vstring('HiggsSM:gg2H = on',
            'HiggsSM:ff2Hff(t:ZZ) = off',
            'HiggsSM:ff2Hff(t:WW) = off',
            'HiggsSM:ffbar2HZ = off',
            'HiggsSM:ffbar2HW = off',
            'TauDecays:mode = 2',
            'TauDecays:tauPolarization = 0',
            'TauDecays:tauMother = 25',
            '6:m0 = 172.04',
            '25:m0 = 125.0',
            '25:addChannel 1 0.1 100 15 -11',
            '25:addChannel 1 0.1 100 11 -15',
            '25:onMode = off',
            '25:onIfMatch 15 11'),
        pythia8CUEP8M1Settings = cms.vstring('Tune:pp 14',
            'Tune:ee 7',
            'MultipartonInteractions:pT0Ref=2.4024',
            'MultipartonInteractions:ecmPow=0.25208',
            'MultipartonInteractions:expPow=1.6'),
        pythia8CommonSettings = cms.vstring('Main:timesAllowErrors = 10000',
            'Check:epTolErr = 0.01',
            'Beams:setProductionScalesFromLHEF = on',
            'SLHA:keepSM = on',
            'SLHA:minMassSM = 1000.',
            'ParticleDecays:limitTau0 = on',
            'ParticleDecays:tau0Max = 10',
            'ParticleDecays:allowPhotonRadiation = on')
    ),
    comEnergy = cms.double(13000.0),
    filterEfficiency = cms.untracked.double(1),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(1)
)
process.RandomNumberGeneratorService.generator.initialSeed = """+str(123456789+1000000*i)+"""
process.ProductionFilterSequence = cms.Sequence(process.generator)

from FastSimulation.PileUpProducer.PileUpSimulator_E13TeV_AVE_20_inTimeOnly_cff import famosPileUp
process.famosPileUp.PileUpSimulator.averageNumber = 20.
process.famosPileUp.PileUpSimulator.usePoisson = True



# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.simulationWithFamos)
process.reconstruction_step = cms.Path(process.reconstructionWithFamos)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.AODoutput_step = cms.EndPath(process.AODoutput)


# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.reconstruction_step,process.eventinterpretaion_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.AODoutput_step ])
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from FastSimulation.Configuration.MixingModule_Full2Fast
from FastSimulation.Configuration.MixingModule_Full2Fast import setVertexGeneratorPileUpProducer 

#call to customisation function setVertexGeneratorPileUpProducer imported from FastSimulation.Configuration.MixingModule_Full2Fast
process = setVertexGeneratorPileUpProducer(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions


"""
    #scrittura del file di configurazione
    cfg_file = open("ggHiggsETau_8TeVBeamSpot_20PU-"+str(i)+"_cfg.py","w")
    cfg_file.write(cfg)
    cfg_file.close()


    #crea lo script per lanciare su cluster
    sh = """#!/bin/tcsh -f
            set W_DIR = \"/afs/cern.ch/user/t/taroni/scratch1/LFVH/CMSSW_7_2_3/src/FastSimulation/Event/test\"
            set CFG =   \"/afs/cern.ch/user/t/taroni/scratch1/LFVH/CMSSW_7_2_3/src/FastSimulation/Event/test/ggHiggsETau_8TeVBeamSpot_20PU-"""+str(i)+"""_cfg.py\"
            cd $W_DIR
            eval `scramv1 runtime -csh`
            cd -
            cmsRun $CFG
            cmsStage GGHetau_8TeVBeamSpot_20PU_AODSIM_"""+str(1+i)+""".root /store/user/taroni/723/ggHiggsETau/
            
            exit
    """

    #scrive script
    sh_file = open("submitggHiggsETau_8TeVBeamSpot_20PU-"+str(i)+".sh","w")
    sh_file.write(sh)
    sh_file.close()

    #sottomette script
    popen("chmod a+x submitggHiggsETau_8TeVBeamSpot_20PU-"+str(i)+".sh" )
    popen("bsub -q 1nw submitggHiggsETau_8TeVBeamSpot_20PU-"+str(i)+".sh")
    
    i+=1
