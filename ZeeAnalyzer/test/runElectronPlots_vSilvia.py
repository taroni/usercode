import FWCore.ParameterSet.Config as cms

process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_Candidate_forECALStudies', '')

# input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

inputFilesMC = cms.untracked.vstring(
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/00CDB4C7-5C93-E711-AF33-02163E0142CA.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/027E1441-3994-E711-BFBD-02163E01A6D8.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/02FD6F07-5D93-E711-85AC-02163E01A334.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/061B6C49-5793-E711-AF23-02163E011B7C.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/0A66322F-5793-E711-9184-02163E01A2BD.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/0EFBF8C4-5C93-E711-94C9-02163E012207.root',
'/store/mc/RunIISummer17DRPremix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/14FDD26B-7493-E711-8B21-001E67792532.root'
    )    
inputFilesMINIAODSIM = cms.untracked.vstring(
'/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/00F9D855-E293-E711-B625-02163E014200.root',
'/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/06ED1148-1794-E711-9388-02163E013722.root',
'/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/0A8639B4-1D94-E711-9C68-02163E0135C6.root',
'/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/20C45D34-0C94-E711-9D21-02163E013722.root',
'/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/2E664004-8496-E711-A6B0-0025905C43EA.root',
'/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/38232F06-F694-E711-BCFB-02163E01185E.root',
'/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/3894CC9D-0094-E711-A0B7-02163E0140FE.root'
    )    

# xrdcp root://eoscms.cern.ch//eos/cms/store/mc/RunIISummer17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10_ext1-v2/10000/00F9D855-E293-E711-B625-02163E014200.root /tmp/taroni/
inputFiles = inputFilesMINIAODSIM
outputFile = "electron_ntuple.root"
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring('file:/tmp/taroni/00F9D855-E293-E711-B625-02163E014200.root') )                             

process.ntupler = cms.EDAnalyzer(
    'ElectronPlots',
    beamSpot = cms.InputTag('offlineBeamSpot'),
    genEventInfoProduct = cms.InputTag('generator'),
    electrons    = cms.InputTag("gedGsfElectrons"),
    patelectrons = cms.InputTag("slimmedElectronsWithUserData"),
    genParticles = cms.InputTag("genParticles"),
    vertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
    conversions  = cms.InputTag('allConversions'),
    isMC         = cms.bool(True)
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.load("SimpleAnalyzer/ZeeAnalyzer/zeeSkimVetoId2016_silvia") 
process.p = cms.Path(process.zdiElectronSequence*process.ntupler)

