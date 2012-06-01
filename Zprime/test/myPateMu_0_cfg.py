import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patTemplate_cfg import *

process = cms.Process("PAT")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load('Configuration.StandardSequences.Services_cff')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.source.fileNames = cms.untracked.vstring(
'root://eoscms//eos/cms/store/caf/user/taroni/Zpm/243CE41A-53B6-E011-83D4-00215E21DD44.root',

)

process.genFilter = cms.EDFilter("ChannelFilter",
                                 ParticleID= cms.untracked.int32(15),
                                 DaughterIDs= cms.untracked.vint32(11,13)
                                 )


postfix = "PFlow"
runOnMC = True


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.load("Configuration.StandardSequences.MagneticField_cff")
if runOnMC:
    process.GlobalTag.globaltag = cms.string('START42_V11::All')

else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V14::All')

process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')


process.primaryVertexFilter = cms.EDFilter(
    "GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
    )

process.scraping = cms.EDFilter(
    "FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
    )

from PhysicsTools.PatAlgos.tools.coreTools import *

if not runOnMC:
    removeMCMatching(process,["All"])
    
removeSpecificPATObjects(process, ['Photons'],
                         outputInProcess=False)
removeCleaning(process,
               outputInProcess=False)

restrictInputToAOD(process, ['All'])

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, postfix)


from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute'])),
                 doType1MET   = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "ak5"
                 )


getattr(process,"patMuons").embedTrack = True
## getattr(process,"patMuons").embedHighLevel = True
## process.patMuons.usePV = True


getattr(process,"patElectrons").embedTrack = True
getattr(process,"patElectrons").embedGsfTrack = True
getattr(process,"patElectrons").embedSuperCluster = True

from Zprime.ZprimeMass.selectionPat import *
process.eleCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("selectedPatElectrons"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )


process.muCounter =  cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("selectedPatMuons"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )


process.atLeastOneGoodVertexSequence = cms.Sequence(
    process.primaryVertexFilter*process.scraping
    )

process.atLeastOneMuEleSequence = cms.Sequence(
    process.eleCounter*process.muCounter
    )

process.load("Zprime.ZprimeMass.diTausReconstruction_cff")
process.diTau = process.allElecMuPairs.clone()
process.diTau.srcLeg1 = cms.InputTag("selectedPatElectrons")
process.diTau.srcLeg2 = cms.InputTag("selectedPatMuons")
process.diTau.srcMET  = cms.InputTag("patMETsPFlow")
process.diTau.dRmin12  = cms.double(0.7)
process.diTau.doSVreco = cms.bool(True)
if not runOnMC:
    process.diTau.srcGenParticles = ""

process.selectedDiTau = cms.EDFilter(
    "PATElecMuPairSelector",
    src = cms.InputTag("diTau"),
    ##cut = cms.string("charge==0 && mt1MET<40")
##    cut = cms.string("charge ==0 || charge !=0")
    cut = cms.string("dR12 > 0.7 &&  cos(dPhi12) < -0.95 &&((leg1().pt > leg2().pt &&  cos(dPhi1MET)< -0.6)||(leg1().pt < leg2().pt  && cos(dPhi2MET)< -0.6) ) &&(pZeta - 1.5*pZetaVis) > -10. && met().pt > 20"),
    )


process.atLeast1selectedDiTau = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("selectedDiTau"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.diTauSequence = cms.Sequence(
    process.diTau*
    process.selectedDiTau*
    process.atLeast1selectedDiTau
    )


process.pat = cms.Sequence(
##    process.genFilter*
##   process.atLeastOneGoodVertexSequence*
    process.patDefaultSequence*
    process.atLeastOneMuEleSequence
##     process.diTauSequence

    )

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning

## process.demo = cms.EDAnalyzer("ZprimeRecoMass")
## process.TFileService = cms.Service("TFileService",
##                                    fileName = cms.string("myRecoHisto.root")
##                                    )
process.p = cms.Path(process.pat
#                     *process.demo
                     )
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands = cms.untracked.vstring('drop *', *patEventContentNoCleaning)
                               )

process.out.outputCommands.extend( cms.vstring(
    'keep *_TriggerResults_*_*',
    'keep *_hltTriggerSummaryAOD_*_*',
    'keep recoGenParticles_genParticles*_*_*',
    'keep *_patTriggerEvent_*_*',
    'keep *_patTrigger_*_*',
    'keep *_selectedPatJets_*_*',
    'keep *_ak5PFJets_*_*',
    'keep *_particleFlow__*',
    'keep *_offlinePrimaryVerticesDA_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep *_offlinePrimaryVerticesWithBS_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_patMETsPFlow_*_*',
    'keep *_tauGenJetsSelectorAllHadrons_*_*',
    'keep *_kt6PFJetsCentral_rho_*',
    'keep *_kt6PFJetsNeutral_rho_*',
    'keep *_muPtEtaID_*_*',
    'keep *_muPtEtaRelID_*_*',
    'keep *_addPileupInfo_*_*',
    'keep *_generalTracks_*_*',
    'keep *_electronGsfTracks_*_*',
    'keep recoTrackExtras_*_*_*',
    'keep recoGsfTrackExtras_*_*_*',
    'keep *_selectedDiTau_*_*',    
    'keep *_selectedPatElectrons_*_*',
    'keep *_selectedPatMuons_*_*',
    'keep *_selectedPatTaus_*_*',
    'keep *_selectedPatMuonsTriggerMatch_*_*',
    'keep *_selectedPatElectronsTriggerMatch_*_*',
    'keep *_selectedPatTausTriggerMatch_*_*',
    'keep *_selectedPatMuonsTriggerMatchUserEmbedded_*_*',
    'keep *_selectedPatTausTriggerMatchUserEmbedded_*_*',
    'keep *_PatElectrons_*_*',
    'keep *_patMuons_*_*',
    'keep *_patTaus_*_*',
    'keep *_patMuonsTriggerMatch_*_*',
    'keep *_patElectronsTriggerMatch_*_*',
    'keep *_patTausTriggerMatch_*_*',
##     'drop *_TriggerResults_*_HLT',
##     'drop *_TriggerResults_*_RECO',
##     'drop *_selectedPatElectrons_*_*',
##     'drop *_selectedPatMuons_*_*',
##     'drop *_selectedPatTaus_*_*',
##     'drop *_selectedPatMuonsTriggerMatch_*_*',
##     'drop *_selectedPatElectronsTriggerMatch_*_*',
##     'drop *_selectedPatTausTriggerMatch_*_*',
##     'drop *_selectedPatMuonsTriggerMatchUserEmbedded_*_*',
##     'drop *_selectedPatTausTriggerMatchUserEmbedded_*_*',
     )
                                  )
process.outpath = cms.EndPath(process.out)
