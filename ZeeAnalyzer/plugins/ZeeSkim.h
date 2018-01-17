#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
// system include files
#include <memory>
#include <algorithm>
#include <iostream>
#include <vector> 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Utilities/interface/StreamID.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
//
// class declaration
//

class ZeeSkim : public edm::EDFilter {
   public:
      explicit ZeeSkim(const edm::ParameterSet&);
      ~ZeeSkim();
      
      //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);



   private:
      virtual void beginJob() override ;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override ;

      EffectiveAreas _effectiveAreas;
      edm::InputTag rhoTag;
      int nSelectedEvents;
      int nSelectedElectronsEB, nSelectedElectronsEE;
      std::vector<double> massRange;
      std::vector<double> ptThrs;
      edm::EDGetTokenT<double> theRhoToken;
      edm::EDGetTokenT<reco::GsfElectronCollection> theGsfEToken;
      // ----------member data ---------------------------
};

