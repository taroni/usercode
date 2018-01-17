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
inputFilesData = cms.untracked.vstring(
'/store/data/Run2017D/DoubleEG/MINIAOD/PromptReco-v1/000/302/030/00000/E69F63AA-EE8E-E711-8121-02163E019BAF.root',
'/store/data/Run2017D/DoubleEG/MINIAOD/PromptReco-v1/000/302/031/00000/008329E5-368F-E711-A1CD-02163E01A21D.root',
'/store/data/Run2017D/DoubleEG/MINIAOD/PromptReco-v1/000/302/031/00000/04BB9D82-398F-E711-B74B-02163E019BDF.root',
'/store/data/Run2017D/DoubleEG/MINIAOD/PromptReco-v1/000/302/031/00000/407638D4-4B8F-E711-AC24-02163E01437E.root',
'/store/data/Run2017D/DoubleEG/MINIAOD/PromptReco-v1/000/302/031/00000/44B91A0E-488F-E711-A372-02163E019CA5.root',
'/store/data/Run2017D/DoubleEG/MINIAOD/PromptReco-v1/000/302/031/00000/5479D9DF-3C8F-E711-BCF4-02163E01A5EB.root',
'/store/data/Run2017D/DoubleEG/MINIAOD/PromptReco-v1/000/302/031/00000/6496C386-518F-E711-B09E-02163E01341D.root'
)
inputFilesMC = cms.untracked.vstring(
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/00CDB4C7-5C93-E711-AF33-02163E0142CA.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/027E1441-3994-E711-BFBD-02163E01A6D8.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/02FD6F07-5D93-E711-85AC-02163E01A334.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/061B6C49-5793-E711-AF23-02163E011B7C.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/0A66322F-5793-E711-9184-02163E01A2BD.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/0EFBF8C4-5C93-E711-94C9-02163E012207.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/14FDD26B-7493-E711-8B21-001E67792532.root'
    )    

inputFiles = inputFilesMC
outputFile = "electron_ntuple_oldSkim.root"
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(
#'file:/eos/cms/store/user/taroni/ZeeMC17/MINIAOD_00F9D855-E293-E711-B625-02163E014200.root'
'file:/eos/cms/store/user/taroni/ZeeMC17/00CDB4C7-5C93-E711-AF33-02163E0142CA.root'
))

process.ntupler = cms.EDAnalyzer(
    'ElectronPlots',
    beamSpot = cms.InputTag('offlineBeamSpot'),
    genEventInfoProduct = cms.InputTag('generator'),
    electrons    = cms.InputTag("gedGsfElectrons"),
    genParticles = cms.InputTag("genParticles"),
    vertices     = cms.InputTag("offlinePrimaryVertices"),
    conversions  = cms.InputTag('allConversions'),
    isMC         = cms.bool(True)
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.load("DPGAnalysis/Skims/ZElectronSkim_cff") 
process.p = cms.Path(process.zdiElectronSequence*process.ntupler)
