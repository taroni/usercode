import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupVIDSelection
from RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff import *
from RecoEgamma.ElectronIdentification.heepIdVarValueMapProducer_cfi import *
from math import ceil,log

# run on MIONAOD
RUN_ON_MINIAOD = True
print "ZEE SKIM. RUN ON MINIAOD = ",RUN_ON_MINIAOD

# cuts
ELECTRON_CUT=("pt > 10 && abs(eta)<2.5 && userInt('cutbasedID_veto')>0.5")
DIELECTRON_CUT=("mass > 70 && mass < 110 && daughter(0).pt>20 && daughter(1).pt()>10")

##https://cmssdt.cern.ch/lxr/source/PhysicsTools/NanoAOD/python/electrons_cff.py slimmedelectronwithuserdata
# single lepton selectors
egmGsfElectronIDs.physicsObjectIDs = cms.VPSet()
egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')

# define which IDs we want to produce
_electron_id_vid_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
_bitmapVIDForEle_WorkingPoints = cms.vstring(
    "egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto",
)
_bitmapVIDForEle_docstring = ''
for modname in _electron_id_vid_modules: 
    ids= __import__(modname, globals(), locals(), ['idName','cutFlow'])
    for name in dir(ids):
        _id = getattr(ids,name)
        if hasattr(_id,'idName') and hasattr(_id,'cutFlow'):
            setupVIDSelection(egmGsfElectronIDs,_id)
            if (len(_bitmapVIDForEle_WorkingPoints)>0 and _id.idName==_bitmapVIDForEle_WorkingPoints[0].split(':')[-1]):
                _bitmapVIDForEle_docstring = 'VID compressed bitmap (%s), %d bits per cut'%(','.join([cut.cutName.value() for cut in _id.cutFlow]),int(ceil(log(len(_bitmapVIDForEle_WorkingPoints)+1,2))))

if RUN_ON_MINIAOD:
    slimmedElectronsWithUserData = cms.EDProducer("PATElectronUserDataEmbedder",
                                                  src = cms.InputTag("slimmedElectrons"),
                                                  userIntFromBools = cms.PSet(
            cutbasedID_veto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
            ),
                                                  )
    goodZeeElectrons = cms.EDFilter("PATElectronRefSelector",
                                    src = cms.InputTag("slimmedElectronsWithUserData"),
                                    cut = cms.string(ELECTRON_CUT)
                                    )

else:
    goodZeeElectrons = cms.EDFilter("GsfElectronRefSelector",
                                    src = cms.InputTag("gedGsfElectrons"),
                                    cut = cms.string(ELECTRON_CUT)
                                    )


MinEleNumberFilter = cms.EDFilter("CandViewCountFilter",
                                   src =  cms.InputTag("goodZeeElectrons"),
                                   minNumber = cms.uint32(2)
                                   )
# dilepton selectors
diZeeElectrons = cms.EDProducer("CandViewShallowCloneCombiner",
                                decay       = cms.string("goodZeeElectrons goodZeeElectrons"),
                                checkCharge = cms.bool(False),
                                cut         = cms.string(DIELECTRON_CUT)
                                )
# dilepton counters
diZeeElectronsFilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("diZeeElectrons"),
                                    minNumber = cms.uint32(1)
                                    )

#sequences
zdiElectronSequence = cms.Sequence(  egmGsfElectronIDSequence*slimmedElectronsWithUserData* goodZeeElectrons * MinEleNumberFilter*diZeeElectrons* diZeeElectronsFilter )
