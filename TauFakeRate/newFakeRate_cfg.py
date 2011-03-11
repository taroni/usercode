import FWCore.ParameterSet.Config as cms

#isMC = True

process = cms.Process("Demo")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START38_V14::All'

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##-------------------- Disable the CondDB for the L1FastJet (until they are included in a new global tag) -------
#process.ak5PFL1Fastjet.useCondDB = False
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation -----------------------
#process.kt6PFJets.doRhoFastjet = True
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(4.5)


process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## process.source = cms.Source("PoolSource",
##     # replace 'myfile.root' with the source file you want to use
##     fileNames = cms.untracked.vstring(
## #        'file:myfile.root'
##     ''
##     )
## )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

process.demo = cms.EDAnalyzer('FakeRate',
   TauAgainstMuon= cms.InputTag("hpsPFTauDiscriminationAgainstMuon"),
   TauAgainstElectron   = cms.InputTag("hpsPFTauDiscriminationAgainstMuon"),
   TauLooseIso    = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
   TauMediumIso   = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation"),
   TauTightIso    = cms.InputTag("hpsPFTauDiscriminationByTightIsolation"),
   TauDecay       = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
   TauSource = cms.InputTag("hpsPFTauProducer"),
#   JetCorrectionService = cms.string('ak5PFL2L3')
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('new_fakeRate_ntuple_jetCorrected.root' )
                                   )
process.myPartons = cms.EDProducer("PartonSelector",
    withLeptons = cms.bool(True),
    src = cms.InputTag("genParticles")
)
process.flavourByRef = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5PFJets"),
##    jets = cms.InputTag("ak5PFJetsL2L3"),
    coneSizeToAssociate = cms.double(0.5),
    partons = cms.InputTag("myPartons")
)
process.load('RecoJets.Configuration.RecoPFJets_cff')

## process.p = cms.Path(process.ak5PFJets * process.myPartons*process.flavourByRef*process.demo)
##process.p = cms.Path(process.myPartons*process.flavourByRef*process.ak5PFJets *process.demo)
##process.p = cms.Path(process.ak5PFJetsL2L3*process.myPartons*process.flavourByRef*process.demo)
process.p = cms.Path(process.myPartons*process.flavourByRef*process.demo)

readFiles.extend( [
    'rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari2/F008894F-D9E1-DF11-BD4F-0002C90B7F2E.root'
]);


secFiles.extend( [
               ] )

