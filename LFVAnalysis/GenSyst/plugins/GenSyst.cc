// -*- C++ -*-
//
// Package:    LFVAnalysis/GenSyst
// Class:      GenSyst
// 
/**\class GenSyst GenSyst.cc LFVAnalysis/GenSyst/plugins/GenSyst.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Silvia Taroni
//         Created:  Sun, 28 Feb 2016 14:58:30 GMT
//
//

#define DEBUG 0
// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaPhi.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class GenSyst : public edm::one::EDAnalyzer<edm::one::WatchRuns>  {
   public:
      explicit GenSyst(const edm::ParameterSet&);
      ~GenSyst();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual void beginRun( edm::Run const&,  edm::EventSetup const&) override;
      virtual void endRun( edm::Run const&, edm::EventSetup const &) override;

      double dPhi(math::XYZTLorentzVector&, math::XYZTLorentzVector&);

      edm::EDGetTokenT<GenEventInfoProduct>    tok_gen_;
 
      edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken_;
      edm::EDGetTokenT<LHEEventProduct> tok_lhe_ ;
  //      edm::EDGetTokenT<reco::GenParticleCollection> tok_genPart_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> tok_genPart_; 
      edm::EDGetTokenT<reco::GenJetCollection> tok_genJet_;
      edm::EDGetTokenT<pat::JetCollection> tok_jet_;
      int nevent_run;  

      edm::Service<TFileService> fs;      
      TH1F * hMass_had; 
      std::vector<TH1F *> hMass_had_scale;
      std::vector<TH1F *> hMass_had_pdf;
      TH1F * hMass_had_0j; 
      std::vector<TH1F *> hMass_had_scale_0j;
      std::vector<TH1F *> hMass_had_pdf_0j;
      TH1F * hMass_had_1j; 
      std::vector<TH1F *> hMass_had_scale_1j;
      std::vector<TH1F *> hMass_had_pdf_1j;
      TH1F * hMass_had_2j; 
      std::vector<TH1F *> hMass_had_scale_2j;
      std::vector<TH1F *> hMass_had_pdf_2j;
      TH1F * hMass_lep; 
      std::vector<TH1F *> hMass_lep_scale;
      std::vector<TH1F *> hMass_lep_pdf;
      TH1F * hMass_lep_0j; 
      std::vector<TH1F *> hMass_lep_scale_0j;
      std::vector<TH1F *> hMass_lep_pdf_0j;
      TH1F * hMass_lep_1j; 
      std::vector<TH1F *> hMass_lep_scale_1j;
      std::vector<TH1F *> hMass_lep_pdf_1j;
      TH1F * hMass_lep_2j; 
      std::vector<TH1F *> hMass_lep_scale_2j;
      std::vector<TH1F *> hMass_lep_pdf_2j;


     std::vector<TH2F *>  hMuPt_vs_TauPt_had_0j; 
     std::vector<TH2F *>  hMuPt_vs_TauPt_had_1j; 
     std::vector<TH2F *>  hMuPt_vs_TauPt_had_2j; 
     std::vector<TH2F *>  hMuPt_vs_EPt_lep_0j; 
     std::vector<TH2F *>  hMuPt_vs_EPt_lep_1j; 
     std::vector<TH2F *>  hMuPt_vs_EPt_lep_2j;

     std::vector<TH1F *> hWeights; 
   
      TH1F * hMuPt_had; 
      TH1F * hMuMt_had;
      TH1F * hTauPt_had;
      TH1F * hTauMt_had;
      TH1F * hMuPt_had_0j; 
      TH1F * hMuMt_had_0j;
      TH1F * hTauPt_had_0j;
      TH1F * hTauMt_had_0j;
      TH1F * hMuPt_had_1j; 
      TH1F * hMuMt_had_1j;
      TH1F * hTauPt_had_1j;
      TH1F * hTauMt_had_1j;
      TH1F * hMuPt_had_2j; 
      TH1F * hMuMt_had_2j;
      TH1F * hTauPt_had_2j;
      TH1F * hTauMt_had_2j;
    
      TH1F * hMuPt_lep;
      TH1F * hMuMt_lep;
      TH1F * hEPt_lep;
      TH1F * hEMt_lep;
      TH1F * hMuPt_lep_0j;
      TH1F * hMuMt_lep_0j;
      TH1F * hEPt_lep_0j;
      TH1F * hEMt_lep_0j;
      TH1F * hMuPt_lep_1j;
      TH1F * hMuMt_lep_1j;
      TH1F * hEPt_lep_1j;
      TH1F * hEMt_lep_1j;
      TH1F * hMuPt_lep_2j;
      TH1F * hMuMt_lep_2j;
      TH1F * hEPt_lep_2j;
      TH1F * hEMt_lep_2j;

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
GenSyst::GenSyst(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   //usesResource("TFileService");
   tok_gen_       = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   tok_lhe_       = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer")) ;
   lheRunInfoToken_ = consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer",""));
   //tok_genPart_ =consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
   tok_genPart_ =consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"));
   tok_genJet_ =consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJets"));
   tok_jet_ =consumes<pat::JetCollection>(edm::InputTag("slimmedJets"));

}


GenSyst::~GenSyst()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenSyst::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  nevent_run++;
   using namespace edm;
   using namespace reco;
   using namespace std;
   
  if (DEBUG) std::cout << __LINE__ << std::endl;
   edm::Handle<GenEventInfoProduct> genEventInfo;
   iEvent.getByToken(tok_gen_, genEventInfo);


   //std::vector<double>& evtWeights = (std::vector<double>&) genEventInfo->weights();
   double theWeight = genEventInfo->weight();

   edm::Handle<LHEEventProduct> EvtHandle ;
   iEvent.getByToken( tok_lhe_ , EvtHandle ) ;

  if (DEBUG) std::cout << __LINE__ << std::endl;

  // std::string whichWeightId = "1002";
  // for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
  //   if (EvtHandle->weights()[i].id == whichWeightId) theWeight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
  //   if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
    
  // }

  edm::Handle<std::vector<reco::GenJet>> genjets;
  iEvent.getByToken(tok_genJet_,genjets);
  

  
  
  
  edm::Handle<pat::PackedGenParticleCollection> genParticles; 
  iEvent.getByToken(tok_genPart_,genParticles);
 
  const reco::Candidate* higgs  = 0;
  const reco::Candidate* mu     = 0;
  const reco::Candidate* tau     = 0;
  const reco::Candidate* etau   = 0;
  const reco::Candidate* nutau   = 0;
  const reco::Candidate* nue    = 0;
  vector < const reco::Candidate *> tauDau ; 

  for (pat::PackedGenParticleCollection::const_iterator p = genParticles->begin();p != genParticles->end(); ++p ) {
    bool isTau=false; 
    const reco::Candidate *mother = 0;
    const reco::Candidate *newtauDau=0;
    if( p->numberOfMothers()!=0) mother= p->mother(0);
    //if (DEBUG) cout<< __LINE__ << " " << p->pdgId()<< " " << p->status() << " ";
    for (unsigned int im =0; im<10 ; im++){
      if (mother->numberOfMothers()!=0 ) {
	mother=&*mother->mother(0);
	if (isTau==false && abs(mother->pdgId())==15){
	  isTau=true; 
	  newtauDau=&*p;
	}
	//cout << mother->pdgId() << " " ;
	if (mother->pdgId()==25){
	  if(fabs(p->pdgId())==13){
	    mu=&*p;
	  }
	  if(isTau==true) {
	    tau=&*p;
	    tauDau.push_back(&*newtauDau);
	  }
	  higgs = &*p;
	  break;
	}
      }else{
	continue; 
      }

    }
  }
  for (unsigned int ih = 0 ; ih < 222; ih++){
    int pdfset = 1001;
    if (ih<9) {
      pdfset=1001+ih;
    }else if (ih>=9  && ih< 111){
      pdfset=2000+ih-8; 
    }else if (ih>=111 && ih < 166){
      pdfset=3000+ih-110;
    } else if (ih>=166 && ih<222){
      pdfset=4000+ih-165;
    } else{
      cout << "ERROR: non existing pdfset" << endl;
    }
    std::string whichWeightId = std::to_string(pdfset);
    if (DEBUG) cout << __LINE__ << endl;
    double weight = theWeight;
    for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
      if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
    }
    if (DEBUG) cout << __LINE__ << endl;
    hWeights[ih] ->Fill (1, weight);
    if (DEBUG) cout << __LINE__ << endl;
    
  }
  
  int ngenjet=0;
  const reco::GenJet * genJet1=0; 
  const reco::GenJet * genJet2=0; 
  
  for (GenJetCollection::const_iterator jet = genjets->begin(); jet!=genjets->end(); jet++){
    if (jet->pt()< 30) continue; 
    if (abs(jet->eta())>4.7) continue ;
    if (DEBUG) cout<< __LINE__ << " " <<  jet->mother(0) << endl;
    if (tau!=0){
      if ( reco::deltaR(jet->eta(),jet->phi(),tau->eta(),tau->phi()) < 0.4) continue;
    }
    if(mu!=0) {
      if ( reco::deltaR(jet->eta(),jet->phi(),mu->eta(),mu->phi()) < 0.4) continue;
    }

    ngenjet++;
    if (ngenjet>2) continue;
    if (genJet1==0 ) {
      genJet1=&*jet; 

    }else {
      genJet2=&*jet;
    }
  }
  //cout << "number of gen jet " << ngenjet << endl;
 
  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(tok_jet_,jets);
  int njet=0;
  const pat::Jet * jet1=0; 
  const pat::Jet * jet2=0; 
  
  for (pat::JetCollection::const_iterator jet = jets->begin(); jet!=jets->end(); jet++){
    if (jet->pt()< 30) continue; 
    if (abs(jet->eta())>4.7) continue ;
    if (tau!=0){
      if ( reco::deltaR(jet->eta(),jet->phi(),tau->eta(),tau->phi()) < 0.4) continue;
    }
    if(mu!=0) {
      if ( reco::deltaR(jet->eta(),jet->phi(),mu->eta(),mu->phi()) < 0.4) continue;
    }
    if (DEBUG) cout<< __LINE__ << " " <<  jet->mother(0) << endl;
    njet++;
    if (njet>2) continue;
    if (jet1==0 ) {
      jet1=&*jet; 

    }else {
      jet2=&*jet;
    }
  }

  njet=ngenjet;
  
  if (DEBUG) cout << __LINE__ << endl;
  if (higgs == 0){
     return;
  }
  if (DEBUG) cout <<__LINE__ << " " << "tau daughter size " << tauDau.size() << endl;
  if (DEBUG) {
    for (unsigned int id = 0 ; id< tauDau.size(); id++){
      cout<< __LINE__ << " " << tauDau[id]->pdgId() << endl;
          }
  }
  if (tauDau.size() < 3) return;


  for (unsigned int it = 0 ; it< tauDau.size() ; it++) {
    if(fabs(tauDau[it]->pdgId())==11 && etau==0 ){
      etau=&*tauDau[it]; 
    }else if (fabs(tauDau[it]->pdgId())==11 && etau!=0 ){ 
      if(DEBUG)cout << "more than one electron" << endl;
      etau=0;
    }
    if (fabs(tauDau[it]->pdgId())==16 && nutau==0 ){
      nutau=&*tauDau[it]; 
    }else if (fabs(tauDau[it]->pdgId())==16 && nutau!=0){
      if(DEBUG)cout << "more than one tau neutrino" << endl;
      nutau=0;
    }
    if (fabs(tauDau[it]->pdgId())==12 && nue==0 ){
      nue=&*tauDau[it]; 
    }else if (fabs(tauDau[it]->pdgId())==12 && nue!=0 ) {
      if(DEBUG)cout << "more than one ele neutrino" << endl;
      nue=0; 
    }
  }
  if (DEBUG) cout << __LINE__ << endl;
  bool hadDecay=false; 
  if (etau==0 ) hadDecay=true;

  if (DEBUG && tau!=0 && mu!=0) cout<< __LINE__ << " " << hadDecay << " " << tau->pt() << " " << tau->eta() << " " << mu->pt() << " " << mu->phi() << endl; 
  

  if (DEBUG) cout << __LINE__ << endl;
  if (DEBUG) cout << __LINE__ << " " << njet << endl;
  if (mu==0 ) return;
  if ( abs(mu->eta())> 2.1) return;  
  if (DEBUG) cout << __LINE__ << endl;

  if (hadDecay) {
    if (DEBUG) cout << __LINE__ << endl;
    if (abs(tau->eta()) > 2.3) return; 
    if (DEBUG) cout << __LINE__ << endl;
  
  } else{
    if (DEBUG) cout << __LINE__ << endl;
    if (abs(etau->eta())> 2.1 ) return; 
    if (DEBUG) cout << __LINE__ << endl;
  }
  if (DEBUG) cout << __LINE__ << endl;
  
  if (njet > 2) return; 

  if (hadDecay && mu!=0 && tau!=0){

    if(njet==0) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
	int pdfset = 1001;
	if (ih<9) {
	  pdfset=1001+ih;
	}else if (ih>=9  && ih< 111){
	  pdfset=2000+ih-8; 
	}else if (ih>=111 && ih < 166){
	  pdfset=3000+ih-110;
	} else if (ih>=166 && ih<222){
	  pdfset=4000+ih-165;
	} else{
	  cout << "ERROR: non existing pdfset" << endl;
	}
	std::string whichWeightId = std::to_string(pdfset);
	if (DEBUG) cout << __LINE__ << endl;
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	hMuPt_vs_TauPt_had_0j[ih] ->Fill (tau->pt(), mu->pt(), weight);
	if (DEBUG) cout << __LINE__ << endl;

      }
    }

    if(njet==1) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
	int pdfset = 1001;
	if (ih<9) {
	  pdfset=1001+ih;
	}else if (ih>=9  && ih< 111){
	  pdfset=2000+ih-8; 
	}else if (ih>=111 && ih < 166){
	  pdfset=3000+ih-110;
	} else if (ih>=166 && ih<222){
	  pdfset=4000+ih-165;
	} else{
	  cout << "ERROR: non existing pdfset" << endl;
	}
	std::string whichWeightId = std::to_string(pdfset);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << " "<< ih<< " " << hMuPt_vs_TauPt_had_1j.size() << " " << EvtHandle->weights().size() << endl;
	hMuPt_vs_TauPt_had_1j[ih] ->Fill (tau->pt(), mu->pt(), weight);
	if (DEBUG) cout << __LINE__ << endl;
      }
      if (DEBUG) cout << __LINE__ << endl;

    }
    if (DEBUG) cout << __LINE__ << endl;
    if(njet==2) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
	int pdfset = 1001;
	if (ih<9) {
	  pdfset=1001+ih;
	}else if (ih>=9  && ih< 111){
	  pdfset=2000+ih-8; 
	}else if (ih>=111 && ih < 166){
	  pdfset=3000+ih-110;
	} else if (ih>=166 && ih<222){
	  pdfset=4000+ih-165;
	} else{
	  cout << "ERROR: non existing pdfset" << endl;
	}
	if (DEBUG) cout << __LINE__ << endl;

	std::string whichWeightId = std::to_string(pdfset);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	hMuPt_vs_TauPt_had_2j[ih] ->Fill (tau->pt(), mu->pt(), weight);
      }
    }

			 
    
    if (DEBUG) cout << __LINE__ << " " << njet << endl;
    if (mu->pt() < 30) return;
    if (tau->pt() < 35 ) return; 
    //if (tau->mt() > 50 ) return; 
    hMuPt_had->Fill(mu->pt());
    hMuMt_had->Fill(mu->mt());
    hTauPt_had->Fill(tau->pt());
    hTauMt_had->Fill(tau->mt());
    hMass_had->Fill((mu->p4()+tau->p4()).M());
    for (unsigned int ih = 0 ; ih < hMass_had_scale.size(); ih++){
      std::string whichWeightId = std::to_string(1001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	//if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
      }
      if (DEBUG) cout << __LINE__ << endl;
      hMass_had_scale[ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
    }
    for (unsigned int ih = 0 ; ih < 56; ih++){
      std::string whichWeightId = std::to_string(4002+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      hMass_had_pdf[ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
    }
    for (unsigned int ih = 0 ; ih < 55; ih++){
      std::string whichWeightId = std::to_string(3001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      hMass_had_pdf[56+ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
    }
    for (unsigned int ih = 0 ; ih < 102; ih++){
      std::string whichWeightId = std::to_string(2001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      hMass_had_pdf[111+ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
    }

    if (njet==0){
      if (mu->pt() < 45) return;
      if (tau->pt() < 35 ) return; 
      //if (tau->mt() > 50 ) return;
      //if (reco::deltaPhi(tau->phi() , mu->phi())<2.7)return;
      hMuPt_had_0j->Fill(mu->pt());
      hMuMt_had_0j->Fill(mu->mt());
      hTauPt_had_0j->Fill(tau->pt());
      hTauMt_had_0j->Fill(tau->mt());
      hMass_had_0j->Fill((mu->p4()+tau->p4()).M());
      for (unsigned int ih = 0 ; ih < hMass_had_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	hMass_had_scale_0j[ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_had_pdf_0j[ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_had_pdf_0j[56+ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_had_pdf_0j[111+ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
    }//njet==0
    if (njet==1){
      if (mu->pt() < 35) return;
      if (tau->pt() < 40 ) return; 
      //if (tau->mt() > 35 ) return;
      hMuPt_had_1j->Fill(mu->pt());
      hMuMt_had_1j->Fill(mu->mt());
      hTauPt_had_1j->Fill(tau->pt());
      hTauMt_had_1j->Fill(tau->mt());
      hMass_had_1j->Fill((mu->p4()+tau->p4()).M());
      for (unsigned int ih = 0 ; ih < hMass_had_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	hMass_had_scale_1j[ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_had_pdf_1j[ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_had_pdf_1j[56+ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_had_pdf_1j[111+ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
    }//njet==1
    if (njet==2){
      if (mu->pt() < 30) return;
      if (tau->pt() < 40 ) return; 
      //if (tau->mt() > 35 ) return;
      //if (abs(jet1->eta()-jet2->eta()) <2.5) return;
      //if ((jet1->p4()+jet2->p4()).mass() < 200 ) return; 
      hMuPt_had_2j->Fill(mu->pt());
      hMuMt_had_2j->Fill(mu->mt());
      hTauPt_had_2j->Fill(tau->pt());
      hTauMt_had_2j->Fill(tau->mt());
      hMass_had_2j->Fill((mu->p4()+tau->p4()).M());
      for (unsigned int ih = 0 ; ih < hMass_had_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	hMass_had_scale_2j[ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_had_pdf_2j[ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_had_pdf_2j[56+ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_had_pdf_2j[111+ih]->Fill((mu->p4()+tau->p4()).M(), weight); 
      }
    }//njet==2
    
  }//hadDecay
  if (DEBUG) cout << __LINE__ << endl;

  if (!hadDecay){
    if(njet==0) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
	int pdfset = 1001;
	if (ih<9) {
	  pdfset=1001+ih;
	}else if (ih>=9  && ih< 111){
	  pdfset=2000+ih-8; 
	}else if (ih>=111 && ih < 166){
	  pdfset=3000+ih-110;
	} else if (ih>=166 && ih<222){
	  pdfset=4000+ih-165;
	} else{
	  cout << "ERROR: non existing pdfset" << endl;
	}
	std::string whichWeightId = std::to_string(pdfset);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	hMuPt_vs_EPt_lep_0j[ih] ->Fill (etau->pt(), mu->pt(), weight);
      }
    }

    if(njet==1) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
	int pdfset = 1001;
	if (ih<9) {
	  pdfset=1001+ih;
	}else if (ih>=9  && ih< 111){
	  pdfset=2000+ih-8; 
	}else if (ih>=111 && ih < 166){
	  pdfset=3000+ih-110;
	} else if (ih>=166 && ih<222){
	  pdfset=4000+ih-165;
	} else{
	  cout << "ERROR: non existing pdfset" << endl;
	}
	std::string whichWeightId = std::to_string(pdfset);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	hMuPt_vs_EPt_lep_1j[ih] ->Fill (etau->pt(), mu->pt(), weight);
      }
    }
    if(njet==2) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
	int pdfset = 1001;
	if (ih<9) {
	  pdfset=1001+ih;
	}else if (ih>=9  && ih< 111){
	  pdfset=2000+ih-8; 
	}else if (ih>=111 && ih < 166){
	  pdfset=3000+ih-110;
	} else if (ih>=166 && ih<222){
	  pdfset=4000+ih-165;
	} else{
	  cout << "ERROR: non existing pdfset" << endl;
	}
	std::string whichWeightId = std::to_string(pdfset);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	hMuPt_vs_EPt_lep_2j[ih] ->Fill (etau->pt(), mu->pt(), weight);
      }
    }

    if (mu->pt() < 25) return; 
    if (etau->pt() < 15 ) return; 
    //if (etau->mt() > 40 ) return; 
    //if (mu->mt() < 15 ) return; 
    hMuPt_lep->Fill(mu->pt());
    hMuMt_lep->Fill(mu->mt());
    hEPt_lep->Fill(etau->pt());
    hEMt_lep->Fill(etau->mt());
    hMass_lep->Fill((mu->p4()+etau->p4()).M());
    for (unsigned int ih = 0 ; ih < hMass_lep_scale.size(); ih++){
      std::string whichWeightId = std::to_string(1001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	//if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
      hMass_lep_scale[ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
      //if (DEBUG) cout << __LINE__ << " HIGGS MASS " << higgs->mass() << endl;
      
    }
    for (unsigned int ih = 0 ; ih < 56; ih++){
      std::string whichWeightId = std::to_string(4001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      hMass_lep_pdf[ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
      if (DEBUG) cout << __LINE__ << endl;
    }
    for (unsigned int ih = 0 ; ih < 55; ih++){
      std::string whichWeightId = std::to_string(3001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      hMass_lep_pdf[56+ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
      if (DEBUG) cout << __LINE__ << endl;
    }
    for (unsigned int ih = 0 ; ih < 102; ih++){
      std::string whichWeightId = std::to_string(2001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      hMass_lep_pdf[111+ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
      if (DEBUG) cout << __LINE__ << endl;
    }
    if (njet==0){
      if (mu->pt() < 50) return; 
      if (etau->pt() < 15 ) return; 
      //if (etau->mt() > 65 ) return; 
      //if (mu->mt() < 50 ) return;
      
      hMuPt_lep_0j->Fill(mu->pt());
      hMuMt_lep_0j->Fill(mu->mt());
      hEPt_lep_0j->Fill(etau->pt());
      hEMt_lep_0j->Fill(etau->mt());
      hMass_lep_0j->Fill((mu->p4()+etau->p4()).M());
      for (unsigned int ih = 0 ; ih < hMass_lep_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	hMass_lep_scale_0j[ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	    //if (DEBUG) cout << __LINE__ << " HIGGS MASS " << higgs->mass() << endl;
	
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_lep_pdf_0j[ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	    if (DEBUG) cout << __LINE__ << endl;
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_lep_pdf_0j[56+ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	if (DEBUG) cout << __LINE__ << endl;
	  }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_lep_pdf_0j[111+ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	if (DEBUG) cout << __LINE__ << endl;
      }
    }//njets==0
    if (njet==1){
      if (mu->pt() < 45) return; 
      if (etau->pt() < 15 ) return; 
      //if (etau->mt() > 65 ) return; 
      //if (mu->mt() < 40 ) return;
      //if (reco::deltaPhi(mu->phi(), etau->phi())<1.) return;
      hMuPt_lep_1j->Fill(mu->pt());
      hMuMt_lep_1j->Fill(mu->mt());
      hEPt_lep_1j->Fill(etau->pt());
      hEMt_lep_1j->Fill(etau->mt());
      hMass_lep_1j->Fill((mu->p4()+etau->p4()).M());
      for (unsigned int ih = 0 ; ih < hMass_lep_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	hMass_lep_scale_1j[ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	//if (DEBUG) cout << __LINE__ << " HIGGS MASS " << higgs->mass() << endl;
	
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_lep_pdf_1j[ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	    if (DEBUG) cout << __LINE__ << endl;
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	    hMass_lep_pdf_1j[56+ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	    if (DEBUG) cout << __LINE__ << endl;
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_lep_pdf_1j[111+ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	if (DEBUG) cout << __LINE__ << endl;
      }
    }//njet==1
    if (njet==2){
      if (mu->pt() < 25) return; 
      if (etau->pt() < 15 ) return; 
      //if (etau->mt() > 40 ) return; 
      //if (mu->mt() < 15 ) return;
      //if (abs(jet1->eta()-jet2->eta()) <2.5) return;
      //if ((jet1->p4()+jet2->p4()).mass() < 200 ) return; 
      hMuPt_lep_2j->Fill(mu->pt());
      hMuMt_lep_2j->Fill(mu->mt());
      hEPt_lep_2j->Fill(etau->pt());
      hEMt_lep_2j->Fill(etau->mt());
      hMass_lep_2j->Fill((mu->p4()+etau->p4()).M());
      for (unsigned int ih = 0 ; ih < hMass_lep_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	hMass_lep_scale_2j[ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	//if (DEBUG) cout << __LINE__ << " HIGGS MASS " << higgs->mass() << endl;
	
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_lep_pdf_2j[ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	    if (DEBUG) cout << __LINE__ << endl;
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	    hMass_lep_pdf_2j[56+ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	    if (DEBUG) cout << __LINE__ << endl;
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	hMass_lep_pdf_2j[111+ih]->Fill((mu->p4()+etau->p4()).M(), weight); 
	if (DEBUG) cout << __LINE__ << endl;
      }
    }//njet==2
    
  }//!hadDecay
  
  



   
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenSyst::beginJob()
{

  std::stringstream name;
  
  
  hMass_had=fs->make<TH1F> ("HiggsMass_had",  "HiggsMass_had", 180, 0, 180); 
  hMass_had_0j=fs->make<TH1F> ("HiggsMass_had_0j",  "HiggsMass_had_0j", 180, 0, 180); 
  hMass_had_1j=fs->make<TH1F> ("HiggsMass_had_1j",  "HiggsMass_had_1j", 180, 0, 180); 
  hMass_had_2j=fs->make<TH1F> ("HiggsMass_had_2j",  "HiggsMass_had_2j", 180, 0, 180); 
  for (unsigned int ih = 0; ih < 9 ; ih++) {
    name.str("");
    name << "HiggsMass_had_scale_" << 1001+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_scale.push_back(histo);
    name.str("");
    name << "HiggsMass_had_scale_0j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_scale_0j.push_back(histo);
    name.str("");
    name << "HiggsMass_had_scale_1j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_scale_1j.push_back(histo);
    name.str("");
    name << "HiggsMass_had_scale_2j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_scale_2j.push_back(histo);

    name.str("");
    name<< "hWeights_"<<1001+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);
    
    name.str("");
    name << "hMuPt_vs_TauPt_had_0j" << 1001+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_0j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_TauPt_had_1j" << 1001+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_1j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_TauPt_had_2j" << 1001+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_2j.push_back(histo2D);
    
   
  }
  for (unsigned int ih = 1; ih < 57 ; ih++) {
    name.str("");
    name << "HiggsMass_had_pdf_" << 4000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf.push_back(histo);
    name.str("");
    name << "HiggsMass_had_pdf_0j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf_0j.push_back(histo);
    name.str("");
    name << "HiggsMass_had_pdf_1j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf_1j.push_back(histo);
    name.str("");
    name << "HiggsMass_had_pdf_2j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf_2j.push_back(histo);

    name.str("");
    name<< "hWeights_"<<4000+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);

    name.str("");
    name << "hMuPt_vs_TauPt_had_0j" << 4000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_0j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_TauPt_had_1j" << 4000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_1j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_TauPt_had_2j" << 4000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_2j.push_back(histo2D);
   
  }
  for (unsigned int ih = 1; ih < 56 ; ih++) {
    name.str("");
    name << "HiggsMass_had_pdf_" <<3000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf.push_back(histo);
    name.str("");
    name << "HiggsMass_had_pdf_0j" << 3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf_0j.push_back(histo);
    name.str("");
    name << "HiggsMass_had_pdf_1j" << 3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf_1j.push_back(histo);
    name.str("");
    name << "HiggsMass_had_pdf_2j" << 3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf_2j.push_back(histo);

    name.str("");
    name<< "hWeights_"<<3000+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);

    name.str("");
    name << "hMuPt_vs_TauPt_had_0j" << 3000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_0j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_TauPt_had_1j" << 3000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_1j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_TauPt_had_2j" << 3000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_2j.push_back(histo2D);
   
  }
  for (unsigned int ih = 1; ih < 103 ; ih++) {
    name.str("");
    name << "HiggsMass_had_pdf_" <<2000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf.push_back(histo);
    name.str("");
    name << "HiggsMass_had_pdf_0j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf_0j.push_back(histo);
    name.str("");
    name << "HiggsMass_had_pdf_1j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf_1j.push_back(histo);
    name.str("");
    name << "HiggsMass_had_pdf_2j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_had_pdf_2j.push_back(histo);

    name.str("");
    name<< "hWeights_"<<2000+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);

    name.str("");
    name << "hMuPt_vs_TauPt_had_0j" << 2000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_0j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_TauPt_had_1j" << 2000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_1j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_TauPt_had_2j" << 2000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_TauPt_had_2j.push_back(histo2D);
   
  }
  hMass_lep=fs->make<TH1F> ("HiggsVisMass_lep",  "HiggsVisMass_lep", 180, 0, 180); 
  hMass_lep_0j=fs->make<TH1F> ("HiggsVisMass_lep_0j",  "HiggsVisMass_lep_0j", 180, 0, 180); 
  hMass_lep_1j=fs->make<TH1F> ("HiggsVisMass_lep_1j",  "HiggsVisMass_lep_1j", 180, 0, 180); 
  hMass_lep_2j=fs->make<TH1F> ("HiggsVisMass_lep_2j",  "HiggsVisMass_lep_2j", 180, 0, 180); 
  for (unsigned int ih = 0; ih < 9 ; ih++) {
    name.str("");
    name << "HiggsVisMass_lep_scale_" << 1001+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_scale.push_back(histo);
    name.str("");
    name << "HiggsVisMass_lep_scale_0j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_scale_0j.push_back(histo);
    name.str("");
    name << "HiggsVisMass_lep_scale_1j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_scale_1j.push_back(histo);
    name.str("");
    name << "HiggsVisMass_lep_scale_2j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_scale_2j.push_back(histo);
    name.str("");
    name << "hMuPt_vs_EPt_lep_0j" << 1001+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_0j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_EPt_lep_1j" << 1001+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_1j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_EPt_lep_2j" << 1001+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_2j.push_back(histo2D);
   
  }
  for (unsigned int ih = 1; ih < 57 ; ih++) {
    name.str("");
    name << "HiggsMass_lep_pdf_" << 4000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf.push_back(histo);
    name.str("");
    name << "HiggsMass_lep_pdf_0j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf_0j.push_back(histo);
    name.str("");
    name << "HiggsMass_lep_pdf_1j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf_1j.push_back(histo);
    name.str("");
    name << "HiggsMass_lep_pdf_2j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf_2j.push_back(histo);
    name.str("");
    name << "hMuPt_vs_EPt_lep_0j" << 4000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_0j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_EPt_lep_1j" << 4000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_1j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_EPt_lep_2j" << 4000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_2j.push_back(histo2D);
  }
  for (unsigned int ih = 1; ih < 56 ; ih++) {
    name.str("");
    name << "HiggsMass_lep_pdf_" <<3000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf.push_back(histo);
    name.str("");
    name << "HiggsMass_lep_pdf_0j" <<3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf_0j.push_back(histo);
    name.str("");
    name << "HiggsMass_lep_pdf_1j" <<3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf_1j.push_back(histo);
    name.str("");
    name << "HiggsMass_lep_pdf_2j" <<3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf_2j.push_back(histo);
    name.str("");
    name << "hMuPt_vs_EPt_lep_0j" << 3000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_0j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_EPt_lep_1j" << 3000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_1j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_EPt_lep_2j" << 3000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_2j.push_back(histo2D);
  }
  for (unsigned int ih = 1; ih < 103 ; ih++) {
    name.str("");
    name << "HiggsMass_lep_pdf_" <<2000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf.push_back(histo);
    name.str("");
    name << "HiggsMass_lep_pdf_0j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf_0j.push_back(histo);
    name.str("");
    name << "HiggsMass_lep_pdf_1j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf_1j.push_back(histo);
    name.str("");
    name << "HiggsMass_lep_pdf_2j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    hMass_lep_pdf_2j.push_back(histo);
    name.str("");
    name << "hMuPt_vs_EPt_lep_0j" << 2000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_0j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_EPt_lep_1j" << 2000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_1j.push_back(histo2D);
    name.str("");
    name << "hMuPt_vs_EPt_lep_2j" << 2000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    hMuPt_vs_EPt_lep_2j.push_back(histo2D);
 }
  hMuPt_had = fs->make<TH1F> ("hMuPt_had", "hMuPt_had", 100, 0, 100 ) ; 
  hMuMt_had = fs->make<TH1F> ("hMuMt_had", "hMuMt_had", 100, 0, 100 ) ; 
  hTauPt_had= fs->make<TH1F> ("hTauPt_had", "hTauPt_had", 100, 0, 100 ) ; 
  hTauMt_had= fs->make<TH1F> ("hTauMt_had", "hTauMt_had", 100, 0, 100 ) ; 
  hMuPt_had_0j = fs->make<TH1F> ("hMuPt_had_0j", "hMuPt_had_0j", 100, 0, 100 ) ; 
  hMuMt_had_0j = fs->make<TH1F> ("hMuMt_had_0j", "hMuMt_had_0j", 100, 0, 100 ) ; 
  hTauPt_had_0j= fs->make<TH1F> ("hTauPt_had_0j", "hTauPt_had_0j", 100, 0, 100 ) ; 
  hTauMt_had_0j= fs->make<TH1F> ("hTauMt_had_0j", "hTauMt_had_0j", 100, 0, 100 ) ; 
  hMuPt_had_1j = fs->make<TH1F> ("hMuPt_had_1j", "hMuPt_had_1j", 100, 0, 100 ) ; 
  hMuMt_had_1j = fs->make<TH1F> ("hMuMt_had_1j", "hMuMt_had_1j", 100, 0, 100 ) ; 
  hTauPt_had_1j= fs->make<TH1F> ("hTauPt_had_1j", "hTauPt_had_1j", 100, 0, 100 ) ; 
  hTauMt_had_1j= fs->make<TH1F> ("hTauMt_had_1j", "hTauMt_had_1j", 100, 0, 100 ) ; 
  hMuPt_had_2j = fs->make<TH1F> ("hMuPt_had_2j", "hMuPt_had_2j", 100, 0, 100 ) ; 
  hMuMt_had_2j = fs->make<TH1F> ("hMuMt_had_2j", "hMuMt_had_2j", 100, 0, 100 ) ; 
  hTauPt_had_2j= fs->make<TH1F> ("hTauPt_had_2j", "hTauPt_had_2j", 100, 0, 100 ) ; 
  hTauMt_had_2j= fs->make<TH1F> ("hTauMt_had_2j", "hTauMt_had_2j", 100, 0, 100 ) ; 

  hMuPt_lep = fs->make<TH1F> ("hMuPt_lep", "hMuPt_lep", 100, 0, 100 ) ; 
  hMuMt_lep = fs->make<TH1F> ("hMuMt_lep", "hMuMt_lep", 100, 0, 100 ) ; 
  hEPt_lep= fs->make<TH1F> ("hEPt_lep", "hEPt_lep", 100, 0, 100 ) ; 
  hEMt_lep= fs->make<TH1F> ("hEMt_lep", "hEMt_lep", 100, 0, 100 ) ; 
  hMuPt_lep_0j = fs->make<TH1F> ("hMuPt_lep_0j", "hMuPt_lep_0j", 100, 0, 100 ) ; 
  hMuMt_lep_0j = fs->make<TH1F> ("hMuMt_lep_0j", "hMuMt_lep_0j", 100, 0, 100 ) ; 
  hEPt_lep_0j= fs->make<TH1F> ("hEPt_lep_0j", "hEPt_lep_0j", 100, 0, 100 ) ; 
  hEMt_lep_0j= fs->make<TH1F> ("hEMt_lep_0j", "hEMt_lep_0j", 100, 0, 100 ) ; 
  hMuPt_lep_1j = fs->make<TH1F> ("hMuPt_lep_1j", "hMuPt_lep_1j", 100, 0, 100 ) ; 
  hMuMt_lep_1j = fs->make<TH1F> ("hMuMt_lep_1j", "hMuMt_lep_1j", 100, 0, 100 ) ; 
  hEPt_lep_1j= fs->make<TH1F> ("hEPt_lep_1j", "hEPt_lep_1j", 100, 0, 100 ) ; 
  hEMt_lep_1j= fs->make<TH1F> ("hEMt_lep_1j", "hEMt_lep_1j", 100, 0, 100 ) ; 
  hMuPt_lep_2j = fs->make<TH1F> ("hMuPt_lep_2j", "hMuPt_lep_2j", 100, 0, 100 ) ; 
  hMuMt_lep_2j = fs->make<TH1F> ("hMuMt_lep_2j", "hMuMt_lep_2j", 100, 0, 100 ) ; 
  hEPt_lep_2j= fs->make<TH1F> ("hEPt_lep_2j", "hEPt_lep_2j", 100, 0, 100 ) ; 
  hEMt_lep_2j= fs->make<TH1F> ("hEMt_lep_2j", "hEMt_lep_2j", 100, 0, 100 ) ; 


}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenSyst::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenSyst::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
GenSyst::beginRun( edm::Run const & iRun,  edm::EventSetup const& iSetup) {
  nevent_run = 0;
}
void
GenSyst::endRun( const edm::Run& iRun, const edm::EventSetup& iSetup) {
  std::cout <<" Runnumber "<<iRun.run()<<" Nevents  "<<nevent_run;
  if (DEBUG) std::cout << __LINE__ << std::endl;
     edm::Handle<LHERunInfoProduct> run; 
     typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
     
  if (DEBUG) std::cout << __LINE__ << std::endl;
  iRun.getByToken(lheRunInfoToken_, run );
     LHERunInfoProduct myLHERunInfoProduct = *(run.product());
     
    if (DEBUG) std::cout << __LINE__ << std::endl;
    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
       std::cout << iter->tag() << std::endl;
       std::vector<std::string> lines = iter->lines();
       for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
  	 std::cout << lines.at(iLine);
       }
     }
}
 
double GenSyst::dPhi(math::XYZTLorentzVector& vec1, math::XYZTLorentzVector& vec2) {
  return reco::deltaPhi(vec1.phi(),vec2.phi());
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenSyst);
