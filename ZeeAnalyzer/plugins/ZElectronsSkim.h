#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
// system include files
#include <memory>
#include <algorithm>
#include <iostream>
#include <vector> 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
//
// class declaration
//

class ZElectronsSkim : public edm::global::EDFilter<> {
   public:
      explicit ZElectronsSkim(const edm::ParameterSet&);
      ~ZElectronsSkim();
      

   private:
      bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

      EffectiveAreas _effectiveAreas;
      edm::InputTag rhoTag;
      //int nSelectedEvents;

      std::vector<double> massRange;
      std::vector<double> ptThrs;
      edm::EDGetTokenT<double> theRhoToken;
      edm::EDGetTokenT<reco::GsfElectronCollection> theGsfEToken;
      // ----------member data ---------------------------
};

