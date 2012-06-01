// -*- C++ -*-
//
// Package:    ZprimeRecoMass
// Class:      ZprimeRecoMass
// 
/**\class ZprimeRecoMass ZprimeRecoMass.cc Zprime/ZprimeRecoMass/src/ZprimeRecoMass.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Silvia Taroni,21 1-007,+41227676459,
//         Created:  Tue May  8 12:25:26 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include <vector>
#include "TFile.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "TauAnalysis/CandidateTools/interface/CompositePtrCandidateT1T2MEtProducer.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>

#include "TH1D.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;
//
// class declaration
//

class ZprimeRecoMass : public edm::EDAnalyzer {
   public:
      explicit ZprimeRecoMass(const edm::ParameterSet&);
      ~ZprimeRecoMass();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      edm::Service<TFileService> fs;
      TH1F *hRecoMass_int ; 
      TH1F *hRecoMass_fit ; 
      TH1F *hVisMass_nsvfit; 
//       TH1F *hVisMass; 
  //  TH1F *hDeltaPhi;
//       TH1F *hDeltaEta;
//       TH1F *hDeltaR  ;
//   TH1F *hMuPt; 
//   TH1F *hElePt;
  TH1F *hZpMass;
//   TH1F *hDiTauMass;

       // ----------member data ---------------------------
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
ZprimeRecoMass::ZprimeRecoMass(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


ZprimeRecoMass::~ZprimeRecoMass()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZprimeRecoMass::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  

  edm::Handle<PATElecMuPairCollection> diTauHandle;
  iEvent.getByLabel("selectedDiTau",diTauHandle);
  if( !diTauHandle.isValid() )  
    edm::LogError("DataNotAvailable")
      << "No diTau label available \n";
  const PATElecMuPairCollection* diTaus = diTauHandle.product();

  edm::Handle<reco::GenParticleCollection> genParticles ;
  iEvent.getByLabel("genParticles", genParticles) ;


  for (PATElecMuPairCollection::const_iterator iditau = diTaus->begin(); iditau != diTaus->end() ; iditau ++){
//     const NSVfitResonanceHypothesisSummary* nSVfitSol(const std::string& algorithm, int* errorFlag = 0) const;
//     if (iditau->hasNSVFitSolutions() )    for ( std::vector<NSVfitResonanceHypothesisSummary>::const_iterator nSVfitSol= (iditau->nSVfitSolutions())->begin();
// 						sol != iditau->nSVfitSolutions()->end(); ++sol ) {
//       cout << __LINE__ << " " << sol->name() << endl;
//     }
    hVisMass_nsvfit->Fill(iditau->p4Vis().M());
    
    if (iditau->hasNSVFitSolutions() &&iditau->nSVfitSolution("psKine_MEt_logM_int",0)!=0){
//       cout << "Visible mass ==> " << (iditau->p4Vis()).M() << ",  Mass" <<iditau->nSVfitSolution("psKine_MEt_logM_int",0)->mass() << endl;  
      hRecoMass_int ->Fill (iditau->nSVfitSolution("psKine_MEt_logM_int",0)->mass() );
//       if (iditau->nSVfitSolution("psKine_MEt_logM_int",0)->mass() 
    } 
    if (iditau->hasNSVFitSolutions() &&iditau->nSVfitSolution("psKine_MEt_logM_fit",0)!=0){
//       cout << "Visible mass ==> " << (iditau->p4Vis()).M() << ",  Mass" <<iditau->nSVfitSolution("psKine_MEt_logM_fit",0)->mass() << endl;  
      hRecoMass_fit ->Fill (iditau->nSVfitSolution("psKine_MEt_logM_fit",0)->mass() );
    } 
    
  }//theDiTau->nSVfitSolution("psKine_MEt_logM_int",0)->mass()   
  const reco::Candidate *myZp = 0 ; 
  for (reco::GenParticleCollection::const_iterator mcIter=genParticles->begin() ; mcIter!=genParticles->end() ; mcIter++ ){
    if (abs(mcIter ->pdgId()) != 32) continue ;
    if (mcIter ->status() ==3) myZp=&*mcIter;
  }
  if (myZp !=0)  hZpMass -> Fill (myZp->mass()); 

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ZprimeRecoMass::beginJob()
{
  
  hRecoMass_int =  fs->make<TH1F>  ("hRecoMass_int", "hRecoMass_int", 1250,0,2500); 
  hRecoMass_fit =  fs->make<TH1F>  ("hRecoMass_fit", "hRecoMass_fit", 1250,0,2500); 
  hVisMass_nsvfit=  fs->make<TH1F>  ("hVisMass_nsvfit","hVisMass_nsvfit", 1250,0,2500); 

//   hDeltaPhi =  fs->make<TH1F> ("hDeltaPhi" , "hDeltaPhi", 320, -3.2, 3.2) ; 
//   hDeltaEta = fs->make<TH1F>  ("hDeltaEta","hDeltaEta",500,-5,5); 
//   hDeltaR = fs->make<TH1F> ("hDeltaR", "hDeltaR", 500,0,5); 

//   hDiTauMass = fs->make<TH1F> ("hDiTauMass","hDiTauMass",500, 0, 1000);
  hZpMass = fs->make<TH1F> ("hZpMass","hZpMass",500, 0, 1000);

//   hElePt  = fs->make<TH1F> ("hElePt","hElePt",500,0,500);
//   hMuPt  = fs->make<TH1F> ("hMuPt","hMuPt",500,0,500);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZprimeRecoMass::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ZprimeRecoMass::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZprimeRecoMass::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZprimeRecoMass::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZprimeRecoMass::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZprimeRecoMass::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZprimeRecoMass);
