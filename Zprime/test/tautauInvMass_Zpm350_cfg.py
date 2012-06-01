import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patTemplate_cfg import *

process = cms.Process("daTau")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load('Configuration.StandardSequences.Services_cff')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from  Zprime.ZprimeMass.ZprimeSSMToTauTauM350PAT_cfi import *
process.source = cms.Source("PoolSource",
                 fileNames = readFiles,
        secondaryFileNames = secFiles
                                )
#
runOnMC = True



process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
if runOnMC:
    process.GlobalTag.globaltag = cms.string('START42_V11::All')

else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V14::All')


from Zprime.ZprimeMass.diTausReconstruction_cff import *
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
    cut = cms.string("dR12 > 0.7 &&  cos(dPhi12) < -0.95 &&((leg1().pt > leg2().pt &&  cos(dPhi1MET)< -0.6)||(leg1().pt < leg2().pt  && cos(dPhi2MET)< -0.6) ) &&(pZeta - 1.25*pZetaVis) > -10. && met().pt > 20"),
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


process.diTauSequence = cms.Sequence(
    process.diTau*
    process.selectedDiTau*
    process.atLeast1selectedDiTau
    )



process.pat = cms.Sequence(
    process.diTauSequence
    )

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning

process.demo = cms.EDAnalyzer("ZprimeRecoMass")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("myRecoHisto_Zpm350.root")
                                   )
process.p = cms.Path(process.pat*process.demo)
