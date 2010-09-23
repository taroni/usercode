#ifndef GenQuantities_h
#define GenQuantities_h

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
/* #include "FWCore/Framework/interface/EDAnalyzer.h" */
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/Math/interface/Vector.h"


#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include <vector>
#include <string>
#include <sstream>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVector3.h>

typedef struct { 
/*   int run, event; */
  std::vector<float> muPt, muP,muE, muEta, muPhi, muPx, muPy, muPz;
  std::vector<int> muQ;
  std::vector<double> muConeChPt, muConeChP, muConeChEt, muConeChE, muConeChNPart, muConeNePt, muConeNeP, muConeNeEt, muConeNeE, muConeNeNPart;
} BRANCH;


//
// class declaration
//

class GenQuantities : public edm::EDProducer {
   public:
      explicit GenQuantities(const edm::ParameterSet&);
      ~GenQuantities();

//   bool convertParticle(reco::GenParticle& cand, const HepMC::GenParticle * part);
//       int chargeTimesThree( int ) const;


   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
/*       std::string decayChain_ ; */
/*       std::vector<std::string> valias; */
      BRANCH branch;
      TTree *ntuple;


      // ----------member data ---------------------------
  TFile *file           ;
  TFile *ntuplefile     ;
  std::string filename  ;
  std::string outfile  ;
/*   std::string decayChain; */

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
  TH1F * genW1E;
  TH1F * genW2E;
  TH1F * genWMass;
  TH1F * genWMassOn;
  TH1F * genWMassOff;
  TH1F * genWEta;
  TH1F * genWPhi;
  TH1F * genWhiE;
  TH1F * genWloE;
  TH1F * genWpE;
  TH1F * genWmE;
  TH1F * genWsameDirEnergy;
  
  //---Muons
  TH1F * genMuonPt;
  TH1F * genMuonP;
  TH1F * genMuonE;
  TH1F * genMuonEta;
  TH1F * genMuonPhi;
  
  //--- muonsIso
  TH1F * R2muH ;
  
  TH1F * mu1ConeChPtH     ;
  TH1F * mu1ConeChPH      ;
  TH1F * mu1ConeChEtH     ;
  TH1F * mu1ConeChEH      ;
  TH1F * mu1ConeChNPartH  ;
  TH1F * mu2ConeChPtH     ;
  TH1F * mu2ConeChPH      ;
  TH1F * mu2ConeChEtH     ;
  TH1F * mu2ConeChEH      ;
  TH1F * mu2ConeChNPartH  ;
  TH1F * mu1ConeNePtH     ;
  TH1F * mu1ConeNePH      ;
  TH1F * mu1ConeNeEtH     ;
  TH1F * mu1ConeNeEH      ;
  TH1F * mu1ConeNeNPartH  ;
  TH1F * mu2ConeNePtH     ;
  TH1F * mu2ConeNePH      ;
  TH1F * mu2ConeNeEtH     ;
  TH1F * mu2ConeNeEH      ;
  TH1F * mu2ConeNeNPartH  ;

  

};
#endif
