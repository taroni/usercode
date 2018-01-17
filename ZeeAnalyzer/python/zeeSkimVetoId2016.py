import FWCore.ParameterSet.Config as cms

# run on MIONAOD
RUN_ON_MINIAOD = True
print "ZEE SKIM. RUN ON MINIAOD = ",RUN_ON_MINIAOD

# cuts
ELECTRON_CUT=("pt > 10 && abs(eta)<2.5")
DIELECTRON_CUT=("mass > 70 && mass < 110 && daughter(0).pt>20 && daughter(1).pt()>10")


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


# electron ID (sync with the AlCaReco: https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_5_X/Calibration/EcalAlCaRecoProducers/python/WZElectronSkims_cff.py)
###RHO is wrong, I should use the effective area correction for the PU https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2#Rho_effective_area_corrections, not the SC rh. Rho to use fixedGridRhoFastjetAll
identifiedElectrons = goodZeeElectrons.clone(cut = cms.string(goodZeeElectrons.cut.value() +
                                                              " \
&& ( (isEB && (gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\')<=2) \
&& ( (abs(superCluster().position().eta()) <=1.)  \
&& ( (pfIsolationVariables().sumChargedHadronPt +  max(0.0,pfIsolationVariables().sumNeutralHadronEt +  pfIsolationVariables().sumPhotonEt -  0.1703 * superCluster().position().rho())) <0.175 )) \
|| ( ( 1.< abs(superCluster().position().eta()) <=1.479) \
&& ( (pfIsolationVariables().sumChargedHadronPt + max(0.0,pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt - 0.1715 *superCluster().position().rho())) <0.173  ))  \
&& (full5x5_sigmaIetaIeta<0.0115) \
&& ( - 0.228<deltaPhiSuperClusterTrackAtVtx< 0.228 ) \
&& ( -0.00749<deltaEtaSuperClusterTrackAtVtx<0.00749 ) \
&& (hadronicOverEm<0.346) ) \
|| (isEE  \
&& (gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\')<=3) \
&&  ((( 1.479< abs(superCluster().position().eta()) <=2.0) \
&& ( (pfIsolationVariables().sumChargedHadronPt + max(0.0,pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt - 0.1213 * superCluster().position().rho())) <0.159))) \
|| ( ( 2.0 < abs(superCluster().position().eta()) <=2.2) \
&& ( (pfIsolationVariables().sumChargedHadronPt + max(0.0,pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt - 0.1230 * superCluster().position().rho())) < 0.159 )) \
|| (( 2.2 < abs(superCluster().position().eta()) <=2.3) \
&& ( (pfIsolationVariables().sumChargedHadronPt + max(0.0,pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt - 0.1635 * superCluster().position().rho())) < 0.159 )) \
|| ( ( 2.3 < abs(superCluster().position().eta()) <=2.4) \
&& ( (pfIsolationVariables().sumChargedHadronPt + max(0.0,pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt - 0.1937 * superCluster().position().rho())) < 0.159 )) \
|| ( (2.4 < abs(superCluster().position().eta()) <=2.5) \
&& ( (pfIsolationVariables().sumChargedHadronPt + max(0.0,pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt - 0.2393 * superCluster().position().rho())) < 0.159 ))  \
&& (full5x5_sigmaIetaIeta<0.037)\
&& ( -0.213<deltaPhiSuperClusterTrackAtVtx<0.213 ) \
&& ( -0.00895<deltaEtaSuperClusterTrackAtVtx<0.00895 )\
&& (hadronicOverEm<0.211) \
))"
                                                              )
                                             )

# dilepton selectors
diZeeElectrons = cms.EDProducer("CandViewShallowCloneCombiner",
                                decay       = cms.string("identifiedElectrons identifiedElectrons"),
                                checkCharge = cms.bool(False),
                                cut         = cms.string(DIELECTRON_CUT)
                                )
# dilepton counters
diZeeElectronsFilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("diZeeElectrons"),
                                    minNumber = cms.uint32(1)
                                    )

#sequences
zdiElectronSequence = cms.Sequence( goodZeeElectrons * identifiedElectrons * diZeeElectrons * diZeeElectronsFilter )
