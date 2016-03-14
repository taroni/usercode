import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eoscms//eos/cms/store/mc/RunIIFall15MiniAODv2/VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/08C704E1-FAB8-E511-BF98-003048C7A7C4.root',
        'root://eoscms//eos/cms/store/mc/RunIIFall15MiniAODv2/VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/2ADBE0CC-FCB8-E511-99B1-003048C7BA9C.root',
        'root://eoscms//eos/cms/store/mc/RunIIFall15MiniAODv2/VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/3268CE6A-FAB8-E511-8276-0025904E3FCC.root',
        'root://eoscms//eos/cms/store/mc/RunIIFall15MiniAODv2/VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/44638DD3-FCB8-E511-A297-0025904B1284.root',
        'root://eoscms//eos/cms/store/mc/RunIIFall15MiniAODv2/VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/6A3AEEA9-FAB8-E511-B8A2-003048C7B946.root',
        'root://eoscms//eos/cms/store/mc/RunIIFall15MiniAODv2/VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/7495718A-FFB8-E511-A41A-003048C7BA6C.root',
    )
)

process.demo = cms.EDAnalyzer('GenSystReco'
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TheoUnc_HiggsToMuTau_125_VBF.root")
)
process.p = cms.Path(process.demo)
