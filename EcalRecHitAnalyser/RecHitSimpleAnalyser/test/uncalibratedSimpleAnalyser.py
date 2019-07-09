import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

readFiles=cms.untracked.vstring()
readFiles.extend([
'file:/afs/cern.ch/work/t/taroni/private/CMSSW_11_0_X_2019-07-07-2300/src/136.85_RunEGamma2018A+RunEGamma2018A+HLTDR2_2018+RECODR2_2018reHLT_skimEGamma_Offline_L1TEgDQM+HARVEST2018_L1TEgDQM_ON_UncalRH/step3_withUncalibRecHit.root'
#'file:/afs/cern.ch/work/t/taroni/private/CMSSW_11_0_X_2019-07-07-2300/src/tmpOFF/step3.root'
])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputFiles = readFiles
#outputFile = "zMuMu_prasanna_noTagNoRecov.root"
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             
process.simpleAnalyser = cms.EDAnalyzer(
    'RecHitSimpleAnalyser',
    EBuncalibRecHitCollection = cms.InputTag("ecalMultiFitUncalibRecHit","EcalUncalibRecHitsEB"),
    )

#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string( outputFile )
#                                   )

process.p = cms.Path(process.simpleAnalyser)

