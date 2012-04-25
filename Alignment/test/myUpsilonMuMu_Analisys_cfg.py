import FWCore.ParameterSet.Config as cms

process = cms.Process("AlcarecoAnalysis")

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### Standard Configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

### Conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_52_V7::All'

process.source = cms.Source ("PoolSource",fileNames =  cms.untracked.vstring())
process.source.fileNames = [
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/939/DC6147A9-0B8D-E111-9FEE-001D09F252E9.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/931/A2295F5E-F28C-E111-B0E9-001D09F2905B.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/930/90CAEA30-E28C-E111-AA7D-00215AEDFD98.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/929/0ABCF986-D58C-E111-B063-5404A63886EE.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/906/26F37585-D08C-E111-83A1-001D09F295FB.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/903/30B0718F-B88C-E111-90BC-00215AEDFD74.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/891/7C6193EE-B78C-E111-811C-001D09F2447F.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/884/4A8394EE-B78C-E111-A555-001D09F291D2.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/876/34CD7EF8-B78C-E111-8767-E0CB4E55365D.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/859/0C41D7B0-F18C-E111-BF51-001D09F25041.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/858/4E978EEE-F08C-E111-8F37-001D09F23174.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/857/7C150365-E68C-E111-82BE-001D09F252DA.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/856/7675476F-ED8C-E111-BEAE-001D09F242EF.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/849/4E2DED5D-E78C-E111-9671-001D09F25109.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/845/C4EBB532-E38C-E111-A8D6-003048D2C020.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/842/0CF71656-F28C-E111-80C4-0019B9F4A1D7.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/839/262F8680-E88C-E111-9443-001D09F2525D.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/837/AA619C08-F38C-E111-9C90-001D09F242EF.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/834/7AE0920C-FF8C-E111-B8AE-001D09F25041.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/833/2E60B986-D58C-E111-9196-BCAEC518FF65.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/830/F4BE3D4D-248D-E111-B9E5-001D09F29321.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/826/A02354EB-2D8C-E111-B984-BCAEC5364C93.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/811/E44A3D70-958C-E111-95F5-5404A63886C5.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/810/94C9522F-748C-E111-A33A-001D09F241B9.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/800/4EEAB509-7F8C-E111-A590-E0CB4E553651.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/798/86A10A2D-D78B-E111-9E99-002481E0DEC6.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/788/9431F41F-D08B-E111-9208-0015C5FDE067.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/783/30CF5E68-CA8B-E111-908F-001D09F295FB.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/749/02C6EF21-CB8B-E111-903D-001D09F244DE.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/726/A8ABD4B3-E68B-E111-BC95-5404A63886B2.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/721/3229EB2F-1A8C-E111-98FD-003048F110BE.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/720/26DC6184-208C-E111-B8DB-001D09F34488.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/718/C0CE0DCB-378C-E111-8402-001D09F291D2.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/716/94937597-708B-E111-BD02-001D09F2960F.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/701/BC6B2873-8C8B-E111-9938-003048D37666.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/700/C20EDF9E-0F8C-E111-9CA5-003048F1C424.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/697/C85C2AD0-218B-E111-9CA5-0025B32034EA.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/695/FEED9C68-278B-E111-83F3-002481E0D958.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/694/163A5BA8-928B-E111-A03E-0019B9F70468.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/693/B246E016-1A8B-E111-A26D-5404A63886AB.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/692/0279F61F-918B-E111-983D-E0CB4E4408C4.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/691/C87076A9-928B-E111-872A-0025901D62A0.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/689/A87A7498-168B-E111-ADA1-BCAEC5329709.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/682/34FF57C6-128B-E111-91B3-001D09F25041.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/679/D8EAC508-058B-E111-857C-001D09F24FBA.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/676/F4E35310-008B-E111-BF9F-001D09F2512C.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/612/38CC25F5-EA8A-E111-93F7-002481E0CC00.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/608/8E54A867-E58A-E111-867D-BCAEC518FF44.root',
'root://eoscms//eos/cms//store/data/Run2012A/MuOnia/ALCARECO/TkAlUpsilonMuMu-v1/000/191/597/2498612F-DB8A-E111-9A40-001D09F2305C.root'
    ]

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

process.myanalysis = cms.EDAnalyzer("TrackAnalyzerNewTwoBodyHisto",
                                    TkTag = cms.string ('ALCARECOTkAlUpsilonMuMu'),
                                    maxMass = cms.double(9.9),
                                    minMass = cms.double(8.9)
                                    )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('myUpsilonMuMu.root')
)


process.MessageLogger = cms.Service("MessageLogger",
         destinations = cms.untracked.vstring(
                        "cout"
                        ),
         cout = cms.untracked.PSet(
         threshold = cms.untracked.string('DEBUG'),
         INFO = cms.untracked.PSet(
              reportEvery = cms.untracked.int32(1000))
         ),
                                    
)

process.p1 = cms.Path(process.myanalysis)


 
