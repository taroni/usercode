import FWCore.ParameterSet.Config as cms

process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_Candidate_forECALStudies', '')

# input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
outputFile = "electron_ntuple_2016H_MINIAOD.root"
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(
#'file:step3_2016H.root'
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/049041CA-029D-E611-9195-02163E012498.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/0A9AB335-019D-E611-BFCC-02163E0137B8.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/0CCBA92F-019D-E611-814E-02163E01445F.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/1A5185BA-019D-E611-94E6-02163E014379.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/2065012A-019D-E611-8F83-FA163EF2773E.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/206F6B06-019D-E611-B6DF-02163E0146F6.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/2870D72C-019D-E611-8EE0-FA163E212C65.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/3AA33AFA-039D-E611-BEC4-02163E01458F.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/5CED012E-019D-E611-838D-FA163EFE6BCB.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/5EE93F54-039D-E611-8E41-02163E014500.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/6426E32E-019D-E611-9AA3-02163E0133CF.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/6CB578F8-009D-E611-B7B9-02163E01436C.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/6E03F134-019D-E611-8BDE-02163E0119A2.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/6E1F8A2A-019D-E611-B4F2-FA163E2648F2.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/74C07396-049D-E611-9DC7-02163E0140F0.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/8089239E-049D-E611-AC20-02163E011FDA.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/84D481F8-019D-E611-9B73-02163E0144DE.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/94BFEBDD-029D-E611-AD09-02163E012AA0.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/9ADB1F33-029D-E611-9EA8-02163E0142AC.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/BC632366-039D-E611-862F-FA163E6198B8.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/CA15C2D5-029D-E611-99BA-02163E011A88.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/D0DB6828-019D-E611-89A6-FA163E35ED47.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/E40EAFF4-019D-E611-9326-02163E014672.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/E4C29A55-019D-E611-BB08-02163E0138EF.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/EC86217E-039D-E611-9033-02163E014559.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/F89B64E6-049D-E611-890F-02163E0144B5.root',
'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/877/00000/FEED3897-049D-E611-A9E1-FA163EFC9C91.root'
))
import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'goodList2016H.json').getVLuminosityBlockRange()



process.ntupler = cms.EDAnalyzer(
    'ElectronPlots',
    beamSpot = cms.InputTag('offlineBeamSpot'),
    genEventInfoProduct = cms.InputTag('generator'),
#    electrons    = cms.InputTag("gedGsfElectrons"),
    electrons    = cms.InputTag("slimmedElectrons"),
    genParticles = cms.InputTag("genParticles"),
#    vertices     = cms.InputTag("offlinePrimaryVertices"),
    vertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
    conversions  = cms.InputTag('allConversions'),
    isMC         = cms.bool(True)
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.load("DPGAnalysis/Skims/ZElectronSkim_cff") 
process.p = cms.Path(process.zdiElectronSequence*process.ntupler)
