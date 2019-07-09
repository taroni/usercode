// -*- C++ -*-
//
// Package:    Analyser/RecHitSimpleAnalyser
// Class:      RecHitSimpleAnalyser
//
/**\class RecHitSimpleAnalyser RecHitSimpleAnalyser.cc Analyser/RecHitSimpleAnalyser/plugins/RecHitSimpleAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Silvia Taroni
//         Created:  Mon, 08 Jul 2019 10:03:38 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"



#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class UncalibratedRecHitSimpleAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit UncalibratedRecHitSimpleAnalyser(const edm::ParameterSet&);
      ~UncalibratedRecHitSimpleAnalyser();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  std::vector<unsigned int> vRun, vLumi; 
  std::vector<unsigned long long> vEvt;
  std::vector<unsigned int> vRawId;

  edm::EDGetTokenT<EBUncalibratedRecHitCollection> ebUncalibRecHitToken_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
UncalibratedRecHitSimpleAnalyser::UncalibratedRecHitSimpleAnalyser(const edm::ParameterSet& iConfig)
{
  //ebRecHitToken_ = consumes<EBRechHitCollection>(ps.getParameter<std::string>("EBrechitCollection"));
  
  ebUncalibRecHitToken_ = consumes<EBUncalibratedRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBuncalibRecHitCollection"));
   //now do what ever initialization is needed

}


UncalibratedRecHitSimpleAnalyser::~UncalibratedRecHitSimpleAnalyser()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
UncalibratedRecHitSimpleAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   const EBUncalibratedRecHitCollection* ebUncalibRecHits = nullptr;

   Handle<EBUncalibratedRecHitCollection> pEBUncalibRecHits;
   iEvent.getByToken(ebUncalibRecHitToken_, pEBUncalibRecHits);
   ebUncalibRecHits = pEBUncalibRecHits.product();
 
   edm::ESHandle<EcalChannelStatus> pChannelStatus;
   iSetup.get<EcalChannelStatusRcd>().get(pChannelStatus);
   const EcalChannelStatus* ch_status = pChannelStatus.product();

   bool checkEvent=false;
 
   for (uint i=0; i<vRun.size(); i++){
     if (iEvent.eventAuxiliary().run()==vRun[i] && iEvent.eventAuxiliary().luminosityBlock()==vLumi[i] && iEvent.eventAuxiliary().event()==vEvt[i]) {
       checkEvent=true;
 
       break;
     }

   }
   if (checkEvent==false) return;
   for (EBUncalibratedRecHitCollection::const_iterator it = ebUncalibRecHits->begin(); it != ebUncalibRecHits->end();
     ++it) {
     bool checkRawId=false;

     unsigned int rawId_=0;
     for (uint i=0; i<vEvt.size(); i++){
       if (iEvent.eventAuxiliary().run()==vRun[i] && iEvent.eventAuxiliary().luminosityBlock()==vLumi[i] && iEvent.eventAuxiliary().event()==vEvt[i] && it->id().rawId()==vRawId[i]){
	 rawId_=it->id().rawId();
	 checkRawId=true;
	 break;
       }
     }
     if (checkRawId==false)continue;
     EcalChannelStatusMap::const_iterator chit;
     chit = ch_status->getMap().find(rawId_);
     if (chit != ch_status->getMap().end()) {
       EcalChannelStatusCode ch_code = (*chit);
       //11=kNonRespondingIsolated  from https://cmssdt.cern.ch/lxr/source/CondFormats/EcalObjects/interface/EcalChannelStatusCode.h?v=CMSSW_11_0_X_2019-07-04-2300#0018
       std::cout << iEvent.eventAuxiliary().run()<< " " << iEvent.eventAuxiliary().luminosityBlock()<< " " << iEvent.eventAuxiliary().event() << " "  << rawId_ << " EcalChannelStatus: " << std::setprecision(6) << ch_code.getStatusCode() << std::endl;
     } else {
       std::cout << "No channel status found for this xtal! something wrong with EcalChannelStatus in your DB? " << std::endl;
     }
   }

}


// ------------ method called once each job just before starting event loop  ------------
void
UncalibratedRecHitSimpleAnalyser::beginJob()
{
  std::fstream inFile;
  inFile.open("evt.csv");
  std::string evt_, rawId_, run_, lumi_; 
  std::string str;
  if(!inFile){
      std::cout << "Cannot open the File : evt.csv"<<std::endl;
  }
  while (std::getline(inFile, str)) {
    std::cout << str<< std::endl;
    int pos0=str.find(',');
    int pos1=str.find(',', pos0+1); 
    int pos2=str.find(',', pos1+1); 
    run_=str.substr(0, pos0);
    lumi_=str.substr(pos0+1, pos1-pos0-1);
    evt_=str.substr(pos1+1, pos2-pos1-1);
    rawId_=str.substr(pos2+1, str.length()-pos2);
    //std::getline(inFile, run_, ',');
    //std::getline(inFile, lumi_, ','); 
    //std::getline(inFile, evt_, ',');
    //std::getline(inFile, rawId_);
    vLumi.push_back((unsigned int)atoi(lumi_.c_str())) ;
    vEvt.push_back((unsigned long long)atoll(evt_.c_str())) ;
    vRawId.push_back((unsigned int)atoi(rawId_.c_str()));
    vRun.push_back((unsigned int)atoi(run_.c_str()));
    std::cout << run_ << " " <<  lumi_ << " " << evt_ << " " << rawId_ << std::endl;
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void
UncalibratedRecHitSimpleAnalyser::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
UncalibratedRecHitSimpleAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(UncalibratedRecHitSimpleAnalyser);
