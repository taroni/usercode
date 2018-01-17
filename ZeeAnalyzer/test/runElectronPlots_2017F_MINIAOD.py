import FWCore.ParameterSet.Config as cms

process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_Candidate_forECALStudies', '')

# input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
outputFile = "electron_ntuple_2017F_MINIAOD.root"
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(
'file:step3_2017F.root'
))
import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'goodList2017F.json').getVLuminosityBlockRange()



process.ntupler = cms.EDAnalyzer(
    'ElectronPlots',
    beamSpot = cms.InputTag('offlineBeamSpot'),
    genEventInfoProduct = cms.InputTag('generator'),
    electrons    = cms.InputTag("gedGsfElectrons"),
#    electrons    = cms.InputTag("slimmedElectrons"),
    genParticles = cms.InputTag("genParticles"),
    vertices     = cms.InputTag("offlinePrimaryVertices"),
#    vertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
    conversions  = cms.InputTag('allConversions'),
    isMC         = cms.bool(True)
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.load("DPGAnalysis/Skims/ZElectronSkim_cff") 
process.p = cms.Path(process.zdiElectronSequence*process.ntupler)
