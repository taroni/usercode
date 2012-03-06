#!/bin/env python


from os import popen

#numero iniziale
i=0

#inizio ciclo
while i<100:


    #inizio file cfg
    cfg = """# Auto generated configuration file
import FWCore.ParameterSet.Config as cms

process = cms.Process('ALCARECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

## # Input source
## process.source = cms.Source(\"PoolSource\",
##     secondaryFileNames = cms.untracked.vstring(),
##     fileNames = cms.untracked.vstring(
## #'/store/data/Run2011B/MuOnia/RECO/19Nov2011-v1/0000/0A32E3F9-701A-E111-B5F7-003048678B0C.root',
##     'file:/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN2/HIP/taroni/validation/CMSSW_4_4_2_patch8/src/Alignment/CommonAlignmentProducer/test/muOnia_002AAB01-F81B-E111-853F-002618FDA259.root',
## ##     'file:/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN2/HIP/taroni/validation/CMSSW_4_4_2_patch8/src/Alignment/CommonAlignmentProducer/test/muOnia_003CC2E1-2C1A-E111-880E-003048FFD7D4.root',
## ##     'file:/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN2/HIP/taroni/validation/CMSSW_4_4_2_patch8/src/Alignment/CommonAlignmentProducer/test/muOnia_002C4A3A-751A-E111-B47D-00304867C0EA.root'
##     )
## )
process.source = cms.Source (\"PoolSource\",fileNames =  cms.untracked.vstring())
process.source.fileNames = []

process.options = cms.untracked.PSet(

)

# AlCaReco for track based alignment using Upsilon->MuMu events
import FWCore.ParameterSet.Config as cms

import  HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.ALCARECOTkAlUpsilonMuMuHLT = hlt.hltHighLevel.clone(
    andOr = True, ## choose logical OR between Triggerbits
    eventSetupPathsKey = 'TkAlUpsilonMuMu',
    throw = False # tolerate triggers stated above, but not available
    )

import DPGAnalysis.Skims.skim_detstatus_cfi as detstatus
process.ALCARECOTkAlUpsilonMuMuDCSFilter = detstatus.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX',
                               'DT0','DTp','DTm','CSCp','CSCm'),
    ApplyFilter  = cms.bool(True),
    AndOr        = cms.bool(True),
    DebugOn      = cms.untracked.bool(False)
)

## import Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi as TkAlMuonSelectors
## process.ALCARECOTkAlUpsilonMuMuGoodMuons = TkAlMuonSelectors.TkAlGoodIdMuonSelector.clone()
## process.ALCARECOTkAlUpsilonMuMuRelCombIsoMuons = TkAlMuonSelectors.TkAlRelCombIsoMuonSelector.clone(
##     src = 'ALCARECOTkAlUpsilonMuMuGoodMuons'
## )
process.ALCARECOTkAlUpsilonMuMuGoodMuons = cms.EDFilter(\"MuonSelector\",
    src = cms.InputTag('muons'),
    cut = cms.string('isGlobalMuon &'
                     'isTrackerMuon &'
                     'numberOfMatches > 1 &'
                     'globalTrack.hitPattern.numberOfValidMuonHits > 0 &'
                     'abs(eta) < 2.5 &'
                     'globalTrack.normalizedChi2 < 20.'),
    filter = cms.bool(True)
)

process.ALCARECOTkAlUpsilonMuMuRelCombIsoMuons = cms.EDFilter(\"MuonSelector\",
    src = cms.InputTag('ALCARECOTkAlUpsilonMuMuGoodMuons'),
    cut = cms.string('(isolationR03().sumPt + isolationR03().emEt + isolationR03().hadEt)/pt  < 0.3'),
    filter = cms.bool(True)
)

import Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi as TrackSelector
process.ALCARECOTkAlUpsilonMuMu = TrackSelector.AlignmentTrackSelector.clone()
process.ALCARECOTkAlUpsilonMuMu.filter = True ##do not store empty events

process.ALCARECOTkAlUpsilonMuMu.applyBasicCuts = True
process.ALCARECOTkAlUpsilonMuMu.ptMin = 0.8 ##GeV
process.ALCARECOTkAlUpsilonMuMu.etaMin = -3.5
process.ALCARECOTkAlUpsilonMuMu.etaMax = 3.5
process.ALCARECOTkAlUpsilonMuMu.nHitMin = 0

##process.ALCARECOTkAlUpsilonMuMu.GlobalSelector.muonSource = 'ALCARECOTkAlUpsilonMuMuRelCombIsoMuons'
process.ALCARECOTkAlUpsilonMuMu.GlobalSelector.muonSource = 'ALCARECOTkAlUpsilonMuMuGoodMuons'
# Isolation is shifted to the muon preselection, and then applied intrinsically if applyGlobalMuonFilter = True
process.ALCARECOTkAlUpsilonMuMu.GlobalSelector.applyIsolationtest = False
process.ALCARECOTkAlUpsilonMuMu.GlobalSelector.applyGlobalMuonFilter = True

process.ALCARECOTkAlUpsilonMuMu.TwoBodyDecaySelector.applyMassrangeFilter = True
process.ALCARECOTkAlUpsilonMuMu.TwoBodyDecaySelector.minXMass = 8.9 ##GeV
process.ALCARECOTkAlUpsilonMuMu.TwoBodyDecaySelector.maxXMass = 9.9 ##GeV
process.ALCARECOTkAlUpsilonMuMu.TwoBodyDecaySelector.daughterMass = 0.105 ##GeV (Muons)
process.ALCARECOTkAlUpsilonMuMu.TwoBodyDecaySelector.applyChargeFilter = True
process.ALCARECOTkAlUpsilonMuMu.TwoBodyDecaySelector.charge = 0
process.ALCARECOTkAlUpsilonMuMu.TwoBodyDecaySelector.applyAcoplanarityFilter = False
process.ALCARECOTkAlUpsilonMuMu.TwoBodyDecaySelector.acoplanarDistance = 1 ##radian



process.ALCARECOStreamTkAlUpsilonMuMu = cms.OutputModule(\"PoolOutputModule\",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlUpsilonMuMu')
    ),
    outputCommands = cms.untracked.vstring('drop *',
        'keep *_ALCARECOTkAlUpsilonMuMu_*_*',
        'keep L1AcceptBunchCrossings_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_TriggerResults_*_*',
        'keep DcsStatuss_scalersRawToDigi_*_*',
        'keep *_offlinePrimaryVertices_*_*'),
    fileName = cms.untracked.string('myTkAlUpsilonMuMu_"""+str(i)+""".root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('TkAlUpsilonMuMu'),
        dataTier = cms.untracked.string('ALCARECO')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880)
)

process.demo = cms.EDAnalyzer(\"NewTrackAnalyzer\",
                                    TkTag = cms.string ('ALCARECOTkAlUpsilonMuMu'),
                                    maxMass = cms.double(9.9),
                                    minMass = cms.double(8.9)
                                    )

process.TFileService = cms.Service(\"TFileService\",
#    fileName = cms.string('myUpsilonMuMu.root',')
#    fileName = cms.string('myUpsilonMuMu_all.root',')
    fileName = cms.string('myUpsilonMuMu_"""+str(i)+""".root')
)

process.seqALCARECOTkAlUpsilonMuMu = cms.Sequence(process.ALCARECOTkAlUpsilonMuMuHLT+process.ALCARECOTkAlUpsilonMuMuDCSFilter+process.ALCARECOTkAlUpsilonMuMuGoodMuons+process.ALCARECOTkAlUpsilonMuMu+ process.demo)

process.MessageLogger = cms.Service(\"MessageLogger\",
         destinations = cms.untracked.vstring(
                        \"cout\"
                        ),
         cout = cms.untracked.PSet(
         threshold = cms.untracked.string('DEBUG'),
         INFO = cms.untracked.PSet(
              reportEvery = cms.untracked.int32(1000))
         ),

)


# Other statements
process.GlobalTag.globaltag = 'GR_R_44_V12::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.dqmoffline_step = cms.Path(process.DQMOffline)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.ALCARECOoutput_step = cms.EndPath(process.ALCARECOoutput)
process.ALCARECOStreamTkAlUpsilonMuMuOutPath = cms.EndPath(process.ALCARECOStreamTkAlUpsilonMuMu)
process.pathALCARECOTkAlUpsilonMuMu = cms.Path(process.seqALCARECOTkAlUpsilonMuMu*process.ALCARECOTkAlUpsilonMuMuDQM)

# Schedule definition
##process.schedule = cms.Schedule(process.pathALCARECOTkAlMuonIsolated,process.endjob_step,process.ALCARECOoutput_step,process.ALCARECOStreamTkAlMuonIsolatedOutPath)

process.schedule = cms.Schedule(process.pathALCARECOTkAlUpsilonMuMu)
start=20*"""+str(i)+"""
end=20+20*"""+str(i)+"""
count=0
mf = open(\"Oniafiles.txt\")
for ii in mf:
    if (count >= start and count < end):

        ii = ii.strip(\"\\n\")

        process.source.fileNames.append(ii)

    count+=1
"""
    #scrittura del file di configurazione
    cfg_file = open("UpsilonNoIso-"+str(i)+"_cfg.py","w")
    cfg_file.write(cfg)
    cfg_file.close()


    #crea lo script per lanciare su cluster
    sh = """#!/bin/tcsh -f

            set W_DIR = \"/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN2/HIP/taroni/validation/CMSSW_4_4_2_patch10/src/StreamValidation/TrackAnalyzer/test/provaNoIso\"
            set CFG = \"/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN2/HIP/taroni/validation/CMSSW_4_4_2_patch10/src/StreamValidation/TrackAnalyzer/test/provaNoIso/UpsilonNoIso-"""+str(i)+"""_cfg.py\"


            set FILELIST = \"/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN2/HIP/taroni/validation/CMSSW_4_4_2_patch10/src/StreamValidation/TrackAnalyzer/test/provaNoIso/Oniafiles.txt\"

            cp $FILELIST .

            cd $W_DIR
            setenv SCRAM_ARCH slc5_ia32_gcc434
            eval `scramv1 runtime -csh`
            cd -
            cmsRun $CFG

            cmsStage myUpsilonMuMu_"""+str(i)+""".root /store/caf/user/taroni/UpsilonTest2/NoIso/
            exit


    """

    #scrive script
    sh_file = open("onia-"+str(i)+".sh","w")
    sh_file.write(sh)
    sh_file.close()

    #sottomette script
##     popen("chmod a+x QCD100-"+str(i)+".sh" )
##     popen("bsub -q cmscaf1nw QCD100-"+str(i)+".sh" )
    
    i+=1



