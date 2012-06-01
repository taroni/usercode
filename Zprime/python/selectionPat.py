import FWCore.ParameterSet.Config as cms
##&& abs(dB) < 0.02

from PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff import *
selectedPatMuons.cut = cms.string('isGlobalMuon && globalTrack().normalizedChi2 < 10. && globalTrack().hitPattern().numberOfValidTrackerHits>=10 '
                                  '&&pt>20 && abs(eta)<2.1 && isolationR05().sumPt < 3  && isolationR05().emEt <3.5 '
                                  )
selectedPatElectrons.cut = cms.string ('(ecalDrivenSeed ==1 && superCluster.eta() < 1.442&& hadronicOverEm < 0.05  && abs(deltaEtaSuperClusterTrackAtVtx) < 0.005 && abs(deltaPhiSuperClusterTrackAtVtx)< 0.09 && pt > 15 && abs(eta)<2.1 && dr03TkSumPt() < 3 && dr03EcalRecHitSumEt() <3.5 )|| '
                                               ##selection for endcap electrons
                                               '(ecalDrivenSeed ==1 && (superCluster.eta() > 1.56 && superCluster.eta() < 2.5 )&& hadronicOverEm < 0.05 && sigmaIetaIeta <0.03  && abs(deltaEtaSuperClusterTrackAtVtx) < 0.005 && abs(deltaPhiSuperClusterTrackAtVtx)< 0.09  && pt > 15 && abs(eta)<2.1 && dr03TkSumPt() < 3 && dr03EcalRecHitSumEt() <3.5 '
                                               ')')
##                                               ' &&gsfTrack().trackerExpectedHitsInner().numberOfHits ==0')
