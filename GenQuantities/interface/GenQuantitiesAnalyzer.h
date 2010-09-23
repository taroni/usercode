#ifndef GenQuantities_h
#define GenQuantities_h

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include<string>
#include <TFile.h>
#include <TH1.h>



//
// class declaration
//

class GenQuantities : public edm::EDAnalyzer {
   public:
      explicit GenQuantities(const edm::ParameterSet&);
      ~GenQuantities();

//   bool convertParticle(reco::GenParticle& cand, const HepMC::GenParticle * part);
//       int chargeTimesThree( int ) const;


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  TFile *file           ;
  std::string filename  ;

  //---Higgs
  TH1F * genHiggsPt;
  TH1F * genHiggsP;
  TH1F * genHiggsE;
  TH1F * genHiggsEta;
  TH1F * genHiggsPhi;

  //---W
  TH1F * genWPt;
  TH1F * genWP;
  TH1F * genWE;
  TH1F * genWEta;
  TH1F * genWPhi;
  
  //---Muons
  TH1F * genMuonPt;
  TH1F * genMuonP;
  TH1F * genMuonE;
  TH1F * genMuonEta;
  TH1F * genMuonPhi;
  
  int ievt;
  

};
#endif
