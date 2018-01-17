import FWCore.ParameterSet.Config as cms

# run on MIONAOD
RUN_ON_MINIAOD = False
#print "ZEE SKIM. RUN ON MINIAOD = ",RUN_ON_MINIAOD

# cuts
ELECTRON_CUT=("pt > 10 && abs(eta)<2.5")
DIELECTRON_CUT=("mass > 40 && mass < 140 && daughter(0).pt>20 && daughter(1).pt()>10")


# single lepton selectors
if RUN_ON_MINIAOD:
    goodZeeElectrons = cms.EDFilter("PATElectronRefSelector",
                                    src = cms.InputTag("slimmedElectrons"),
                                    cut = cms.string(ELECTRON_CUT)
                                    )
else:
    goodZeeElectrons = cms.EDFilter("GsfElectronRefSelector",
                                    src = cms.InputTag("gedGsfElectrons"),
                                    cut = cms.string(ELECTRON_CUT)
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
zdiElectronSequence = cms.Sequence( goodZeeElectrons  * diZeeElectrons * diZeeElectronsFilter )
