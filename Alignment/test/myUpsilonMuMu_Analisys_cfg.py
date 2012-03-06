import FWCore.ParameterSet.Config as cms

process = cms.Process("Refitting")

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### Standard Configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

### Conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_44_V12::All'

process.source = cms.Source ("PoolSource",fileNames =  cms.untracked.vstring())
process.source.fileNames = [
#    'file:myUpsilonFile/myTkAlUpsilonMuMu.root','
#    'file:myUpsilonFile/myTkAlUpsilonMuMu_all.root','
#    'file:myUpsilonFile/myTkAlUpsilonMuMu_modified.root','
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_379_1_pev.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_380_1_BgC.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_381_1_5qq.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_382_1_49G.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_383_1_xW7.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_384_1_Kxb.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_385_1_TwG.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_386_1_VAV.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_387_1_ik8.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_388_1_Adv.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_389_1_TIb.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_390_1_TOn.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_391_1_Q9L.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_392_1_Hx9.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_393_1_0CJ.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_394_1_V3z.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_395_1_vHg.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_396_1_ZmQ.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_397_1_gLb.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_398_1_KIk.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_399_1_YSq.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_400_1_Acm.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_401_1_rRX.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_402_1_ZRk.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_403_1_Ojr.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_404_1_ZaK.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_405_1_5E9.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_406_1_UFg.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_407_1_YHa.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_408_1_Z7H.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_409_1_UAZ.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_410_1_NEy.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_411_1_4Ed.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_412_1_M4R.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_413_1_nLV.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_414_1_gE0.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_415_1_7zc.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_416_1_rhE.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_417_1_260.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_418_1_hiP.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_419_1_QqT.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_420_1_aHM.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_421_1_XAK.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_422_1_O1H.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_423_1_Rzo.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_424_1_sjC.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_425_1_9Ey.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_426_1_QRm.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_427_1_DCC.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_428_1_I6u.root',
'root://eoscms//eos/cms/store/caf/user/taroni/Run2011A/myAlcareco/modified/myTkAlUpsilonMuMu_429_1_Kci.root'
    ]

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.myanalysis = cms.EDAnalyzer("TrackAnalyzerNewTwoBodyHisto",
                                    TkTag = cms.string ('ALCARECOTkAlUpsilonMuMu'),
                                    maxMass = cms.double(9.9),
                                    minMass = cms.double(8.9)
                                    )

process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('myUpsilonMuMu.root',')
#    fileName = cms.string('myUpsilonMuMu_all.root',')
    fileName = cms.string('myUpsilonMuMu_modified.root')
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


 
