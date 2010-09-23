// -*- C++ -*-
//
// Package:    GenQuantities
// Class:      GenQuantities
// 
/**\class GenQuantities GenQuantities.cc GluonFusion/GenQuantities/src/GenQuantities.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Tue Aug  3 10:31:24 CEST 2010
// $Id$
//
//


#include "GluonFusion/GenQuantities/interface/GenQuantities.h"

using namespace edm;
using namespace std;
using namespace reco;

//
GenQuantities::GenQuantities(const edm::ParameterSet& iConfig)
  :filename(iConfig.getUntrackedParameter<std::string>("filename")),
   outfile(iConfig.getUntrackedParameter<std::string>("outfile"))
  { 

  }


GenQuantities::~GenQuantities()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GenQuantities::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//   cout << "Processing event #: " << iEvent.id()<< endl;

//   branch.run=iEvent.run();
//   branch.event=iEvent.id().event();
  std::vector<float> muPt, muP, muE, muEta, muPhi, muPx, muPy, muPz;
  std::vector<int> muQ;
  std::vector<double> muConeChPt, muConeChP, muConeChEt, muConeChE, muConeChNPart, muConeNePt, muConeNeP, muConeNeEt, muConeNeE, muConeNeNPart;


  Handle<GenParticleCollection> genParticles; 
  iEvent.getByLabel("genParticles",genParticles);
 

  Handle<HepMCProduct> evt;
  iEvent.getByType(evt);

  const reco::Candidate* higgs  = 0;
  const reco::Candidate* w1     = 0;
  const reco::Candidate* w2     = 0;
  const reco::Candidate* muW1 = 0;
  const reco::Candidate* nuW1 = 0;
  const reco::Candidate* muW2 = 0;
  const reco::Candidate* nuW2 = 0;

  const HepMC::GenEvent * myGenEvent = evt->GetEvent();

  for (GenParticleCollection::const_iterator p = genParticles->begin();p != genParticles->end(); ++p ) {
    if (fabs(p->pdgId())==25 && p->status()==3) {
      higgs = &*p;
    }
  }
  if (higgs == 0){
    cout << __LINE__ << " No Higgs in the event" <<endl;
    return;
  }
  for (Candidate::const_iterator higgsD=higgs->begin(); higgsD!=higgs->end();++higgsD){
    if (higgsD->pdgId()==higgs->pdgId() && higgsD->status()==higgs->status()) 
      cout << __LINE__ << "daughter is the same particle of the mother" << endl;
    if (higgsD->pdgId()==higgs->pdgId()) continue;
    if (fabs(higgsD->pdgId())!=24) continue;
    if (w1 == 0){
      w1 = &* higgsD;
    } else {
      if (w2!=0) cout << __LINE__ <<" MORE THAN 2 W IN THE EVENT " << iEvent.id()<< endl;
      w2 = &*higgsD;
    }
  }
  if (w1 == 0 || w2==0) {
    cout << __LINE__ << " No W in the event " << endl;
    return ;
  } // else {
//     cout << __LINE__<< " pdgId="<<w1->pdgId() <<" mass="<< w1->mass() << ", " << " pdgId="<<w2->pdgId() <<" mass="<< w2->mass() << endl;
//   }
  for (Candidate::const_iterator w1D= w1->begin(); w1D!=w1->end(); ++w1D){
    if (fabs(w1D->pdgId())==24) continue;
//     if (muW1 == 0 && fabs(w1D->pdgId())==13) cout << __LINE__ << " " << w1D->pdgId() << " " << w1D->status()<< " "<< 
//       w1D->daughter(0)->pdgId() << " "<<  w1D->daughter(0)->status()<< endl;
    if (!w1D->daughter(0))continue;
    if (muW1 == 0 && fabs(w1D->pdgId())==13){
      for (unsigned int dau = 0; dau < w1D->numberOfDaughters(); dau++){
	if(w1D->daughter(dau)->pdgId() == w1D->pdgId() &&  w1D->daughter(dau)->status()==1) muW1 = &*(w1D->daughter(dau));
      }
    }
    if (nuW1 == 0 && fabs(w1D->pdgId())==14){
      for (unsigned int dau = 0; dau < w1D->numberOfDaughters(); dau++){
	if(w1D->daughter(dau)->pdgId() == w1D->pdgId() &&  w1D->daughter(dau)->status()==1) nuW1 = &*(w1D->daughter(dau));
      } 
    }
  }
//   if (w2 == 0){
//     cout << __LINE__ << " No W in the event " << endl;
//     return ;
//   }

  for (Candidate::const_iterator w2D= w2->begin(); w2D!=w2->end(); ++w2D){
    if (fabs(w2D->pdgId())==24) continue;
//     if (muW2 == 0 && fabs(w2D->pdgId())==13) cout << __LINE__ << " " << w2D->pdgId() << " " << w2D->status()<< " "<< 
//       w2D->daughter(0)->pdgId() << " "<<  w2D->daughter(0)->status()<< endl;
    if (muW2 == 0 && fabs(w2D->pdgId())==13) {
      for (unsigned int dau = 0; dau < w2D->numberOfDaughters(); dau++){
	if(w2D->daughter(dau)->pdgId() == w2D->pdgId() &&  w2D->daughter(dau)->status()==1) muW2 = &*(w2D->daughter(dau));
      }
    }
    if (nuW2 == 0 && fabs(w2D->pdgId())==14){
      for (unsigned int dau = 0; dau < w2D->numberOfDaughters(); dau++){
	if( w2D->daughter(dau)->pdgId() == w2D->pdgId() &&  w2D->daughter(dau)->status()==1) nuW2 = &*(w2D->daughter(dau));
      }
    }
  }
   //DeltaR muW1 muW2
  if (muW1 == 0 ||  muW2 ==0) {
    cout << __LINE__ << " No muons from Ws; " << muW1 << " mother of the mother: " <<w1->mother()->pdgId() << " " << muW2<< ", mother of the mother: " <<w2->mother()->pdgId()<< endl;
    return;
  }

  //filling ntuple
  muPt.push_back(muW1->pt());
  muPx.push_back(muW1->px());
  muPy.push_back(muW1->py());
  muPz.push_back(muW1->pz());  
  muP.push_back(muW1->p());
  muE.push_back(muW1->energy());
  muEta.push_back(muW1->eta());
  muPhi.push_back(muW1->phi());
  muQ.push_back(muW1->charge());

  muPt.push_back(muW2->pt());
  muPx.push_back(muW2->px());
  muPy.push_back(muW2->py());
  muPz.push_back(muW2->pz());  
  muP.push_back(muW2->p());
  muE.push_back(muW2->energy());
  muEta.push_back(muW2->eta());
  muPhi.push_back(muW2->phi());
  muQ.push_back(muW2->charge());

  branch.muPt=muPt;
  branch.muPx=muPx;
  branch.muPy=muPy;
  branch.muPz=muPz;
  branch.muP=muP;
  branch.muE=muE;
  branch.muEta=muEta;
  branch.muPhi=muPhi;
  branch.muQ=muQ;

  double R2mu= sqrt(pow(muW2->eta()-muW1->eta(),2)+pow(muW2->phi()-muW1->phi(),2));

  //isolation 
  //charged
  double mu1ConeChPt = 0;
  double mu1ConeChP  = 0;
  double mu1ConeChEt = 0;
  double mu1ConeChE  = 0;
  double mu1ConeChNPart = 0;
  double mu2ConeChPt = 0;
  double mu2ConeChP  = 0;
  double mu2ConeChEt = 0;
  double mu2ConeChE  = 0;
  double mu2ConeChNPart = 0;
  //neutral
  double mu1ConeNePt = 0;
  double mu1ConeNeP  = 0;
  double mu1ConeNeEt = 0;
  double mu1ConeNeE  = 0;
  double mu1ConeNeNPart = 0;
  double mu2ConeNePt = 0;
  double mu2ConeNeP  = 0;
  double mu2ConeNeEt = 0;
  double mu2ConeNeE  = 0;
  double mu2ConeNeNPart = 0;
  

  for (GenParticleCollection::const_iterator p = genParticles->begin();p != genParticles->end(); ++p ) {
    if (p->status()!=1) continue;
    if (fabs(p->pdgId())==14) continue;
    double mu1R = 10000;
    double mu2R = 10000;
    
    if (&*p!=muW1) mu1R = sqrt(pow(p->eta()-muW1->eta(),2)+pow(p->phi()-muW1->phi(),2));
    if (&*p!=muW2) mu2R = sqrt(pow(p->eta()-muW2->eta(),2)+pow(p->phi()-muW2->phi(),2));

    if (mu1R < 0.5 ) {
      if ( p->charge()!=0){
	mu1ConeChPt += p->pt();
	mu1ConeChP  += p->p();
	mu1ConeChEt += p->et();
	mu1ConeChE  += p->energy();
	mu1ConeChNPart++;
      }
      
      if (p->charge()==0 && fabs(p->pdgId())!=12 && fabs(p->pdgId())!=14 && fabs(p->pdgId())!=16){
	mu1ConeNePt += p->pt();
	mu1ConeNeP  += p->p();
	mu1ConeNeEt += p->et();
	mu1ConeNeE  += p->energy();
	mu1ConeNeNPart++;
      }
    }
    if (mu2R < 0.5 ) {
      if ( p->charge()!=0){
	mu2ConeChPt += p->pt();
	mu2ConeChP  += p->p();
	mu2ConeChEt += p->et();
	mu2ConeChE  += p->energy();
	mu2ConeChNPart++;
      }
	
      if (p->charge()==0&& fabs(p->pdgId())!=12 && fabs(p->pdgId())!=14 && fabs(p->pdgId())!=16){
	mu2ConeNePt += p->pt();
	mu2ConeNeP  += p->p();
	mu2ConeNeEt += p->et();
	mu2ConeNeE  += p->energy();
      	mu2ConeNeNPart++;
      }
    }
    
  }

  muConeChPt.push_back(mu1ConeChPt);
  muConeChP.push_back(mu1ConeChP);
  muConeChEt.push_back(mu1ConeChEt);
  muConeChE.push_back(mu1ConeChE);
  muConeChNPart.push_back(mu1ConeChNPart);
  muConeNePt.push_back(mu1ConeNePt);
  muConeNeP.push_back(mu1ConeNeP);
  muConeNeEt.push_back(mu1ConeNeEt);
  muConeNeE.push_back(mu1ConeNeE);
  muConeNeNPart.push_back(mu1ConeNeNPart);

  muConeChPt.push_back(mu2ConeChPt);
  muConeChP.push_back(mu2ConeChP);
  muConeChEt.push_back(mu2ConeChEt);
  muConeChE.push_back(mu2ConeChE);
  muConeChNPart.push_back(mu2ConeChNPart);
  muConeNePt.push_back(mu2ConeNePt);
  muConeNeP.push_back(mu2ConeNeP);
  muConeNeEt.push_back(mu2ConeNeEt);
  muConeNeE.push_back(mu2ConeNeE);
  muConeNeNPart.push_back(mu2ConeNeNPart);

  branch.muConeChPt=muConeChPt;
  branch.muConeChP =muConeChP;
  branch.muConeChEt=muConeChEt;
  branch.muConeChE=muConeChE;
  branch.muConeChNPart=muConeChNPart;
  branch.muConeNePt=muConeNePt;
  branch.muConeNeP=muConeNeP;
  branch.muConeNeEt=muConeNeEt;
  branch.muConeNeE=muConeNeE;
  branch.muConeNeNPart=muConeNeNPart;

  //-----Higgs
  genHiggsPt ->Fill(higgs->pt());
  genHiggsP  ->Fill(higgs->p());
  genHiggsE  ->Fill(higgs->energy());
  genHiggsEta->Fill(higgs->eta());
  genHiggsPhi->Fill(higgs->phi());

  //-----W 
  genWPt   ->Fill(w1->pt());
  genWP    ->Fill(w1->p());
  genWE    ->Fill(w1->energy());
  genWEta  ->Fill(w1->eta());
  genWPhi  ->Fill(w1->phi());
  
  genWPt   ->Fill(w2->pt());
  genWP    ->Fill(w2->p());
  genWE    ->Fill(w2->energy());
  genWEta  ->Fill(w2->eta());
  genWPhi  ->Fill(w2->phi());

  genWMass->Fill(w1->mass());
  genWMass->Fill(w2->mass());

  genW1E ->Fill (w1->energy());
  genW2E ->Fill (w2->energy());
  if (w1->mass()>w2->mass()){
    genWMassOn->Fill(w1->energy());
    genWMassOff->Fill(w2->energy());
  }else{
    genWMassOn->Fill(w2->energy());
    genWMassOff->Fill(w1->energy());
  }
  if (w1->energy()>w2->energy()){
    genWhiE->Fill(w1->energy());
    genWloE->Fill(w2->energy());
  }else{
    genWhiE->Fill(w2->energy());
    genWloE->Fill(w1->energy());
  }
  if (w1->charge() > 0){
    genWpE ->Fill(w2->energy());
    genWmE ->Fill(w1->energy());
  }
  if (w1->charge() < 0){
    genWmE ->Fill(w1->energy());
    genWpE ->Fill(w2->energy());
  }
  
  const TVector3 hmomentum = TVector3(higgs->momentum().x(), higgs->momentum().y(),higgs->momentum().z());
  const TVector3 w1momentum = TVector3 (w1->momentum().x(), w1->momentum().y(),w1->momentum().z());
  const TVector3 w2momentum = TVector3 (w2->momentum().x(), w2->momentum().y(),w2->momentum().z());
//    cout << __LINE__ << " W1 " <<  w1->momentum().x()<<" "<<  w1->momentum().y() << " " << w1->momentum().z() << " " << w1momentum.Angle(hmomentum)<< endl;
//    cout << __LINE__ << " W2 " <<  w2->momentum().x()<<" "<<  w2->momentum().y() << " " << w2->momentum().z() << " " << w2momentum.Angle(hmomentum)<<  endl;
  if (cos(w1momentum.Angle(hmomentum))>cos(w2momentum.Angle(hmomentum))) 
    genWsameDirEnergy ->Fill(w1->energy());
  if (cos(w2momentum.Angle(hmomentum))>cos(w1momentum.Angle(hmomentum))) 
    genWsameDirEnergy ->Fill(w2->energy());
  

  
  //---Muons
  genMuonPt  ->Fill(muW1->pt());
  genMuonP   ->Fill(muW1->p());
  genMuonE   ->Fill(muW1->energy());
  genMuonEta ->Fill(muW1->eta());
  genMuonPhi ->Fill(muW1->phi());
  
  genMuonPt  ->Fill(muW2->pt());
  genMuonP   ->Fill(muW2->p());
  genMuonE   ->Fill(muW2->energy());
  genMuonEta ->Fill(muW2->eta());
  genMuonPhi ->Fill(muW2->phi());
   
   //muonIso
  R2muH->Fill(R2mu);
  if (mu1ConeChNPart !=0) {
    mu1ConeChPtH    ->Fill(mu1ConeChPt);
    mu1ConeChPH     ->Fill(mu1ConeChP);
    mu1ConeChEtH    ->Fill(mu1ConeChEt);
    mu1ConeChEH     ->Fill(mu1ConeChE);
    mu1ConeChNPartH ->Fill(mu1ConeChNPart);
   }
  if (mu2ConeChNPart !=0) {
    mu2ConeChPtH    ->Fill(mu2ConeChPt);
    mu2ConeChPH     ->Fill(mu2ConeChP);
    mu2ConeChEtH    ->Fill(mu2ConeChEt);
    mu2ConeChEH     ->Fill(mu2ConeChE);
    mu2ConeChNPartH ->Fill(mu2ConeChNPart);
  }
  if (mu1ConeNeNPart !=0) {   
    mu1ConeNePtH    ->Fill(mu1ConeNePt);
    mu1ConeNePH     ->Fill(mu1ConeNeP);
    mu1ConeNeEtH    ->Fill(mu1ConeNeEt);
    mu1ConeNeEH     ->Fill(mu1ConeNeE);
    mu1ConeNeNPartH ->Fill(mu1ConeNeNPart);
  }
  if (mu2ConeNeNPart !=0) {   
    mu2ConeNePtH    ->Fill(mu2ConeNePt);
    mu2ConeNePH     ->Fill(mu2ConeNeP);
    mu2ConeNeEtH    ->Fill(mu2ConeNeEt);
    mu2ConeNeEH     ->Fill(mu2ConeNeE);
    mu2ConeNeNPartH ->Fill(mu2ConeNeNPart);
  }
  
  ntuple->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenQuantities::beginJob()
{
  file = new TFile(filename.c_str(),"recreate");
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);


  genHiggsPt  = new TH1F ("genHiggsPt" , "genHiggsPt" , 1000,0,1000);
  genHiggsP   = new TH1F ("genHiggsP"  , "genHiggsP " , 1000,0,1000);
  genHiggsE   = new TH1F ("genHiggsE"  , "genHiggsE " , 1000,0,1000);
  genHiggsEta = new TH1F ("genHiggsEta", "genHiggsEta", 100,-10,10);
  genHiggsPhi = new TH1F ("genHiggsPhi", "genHiggsPhi", 720,-3.2, 3.2);

  //---W
  genWPt  = new TH1F ("genWPt" , "genWPt" , 1000,0,1000);
  genWP   = new TH1F ("genWP"  , "genWP " , 1000,0,1000);
  genWE   = new TH1F ("genWE"  , "genWE " , 1000,0,1000);
  genWEta = new TH1F ("genWEta", "genWEta", 100,-10,10);
  genWPhi = new TH1F ("genWPhi", "genWPhi", 720,3.2, 3.2);

  genW1E   = new TH1F ("genW1E"  , "genW1E " , 1000,0,1000);
  genW2E   = new TH1F ("genW2E"  , "genW2E " , 1000,0,1000);

  genWMass   = new TH1F ("genWMass"  , "genWMass "   , 200,0,100);

  genWMassOn = new TH1F ("genWMassOnE" , "genWMassOnE" , 1000,0,1000);
  genWMassOff= new TH1F ("genWMassOffE", "genWMassOffE", 1000,0,1000);
  genWhiE  = new TH1F ("genWhiE"  , "genWhiE " , 1000,0,1000);
  genWloE  = new TH1F ("genWloE"  , "genWloE " , 1000,0,1000);
  genWpE   = new TH1F ("genWpE"  , "genWpE " , 1000,0,1000);
  genWmE   = new TH1F ("genWmE"  , "genWmE " , 1000,0,1000);
  genWsameDirEnergy = new TH1F ("genWsameDirEnergy","genWsameDirEnergy",1000,0,1000);

  //---Muons
  genMuonPt  = new TH1F ("genMuonPt","genMuonPt"   , 2000,0,1000);
  genMuonP   = new TH1F ("genMuonP"  , "genMuonP " , 1000,0,1000);
  genMuonE   = new TH1F ("genMuonE"  , "genMuonE " , 1000,0,1000);
  genMuonEta = new TH1F ("genMuonEta", "genMuonEta", 100, -10,10);
  genMuonPhi = new TH1F ("genMuonPhi", "genMuonPhi", 720,-3.2,3.2);
  
  R2muH = new TH1F ("R2mu","R2mu", 100, 0, 5);
  
  mu1ConeChPtH     = new TH1F ("mu1ConeChPt","mu1ConeChPt",1000,0, 500);
  mu1ConeChPH      = new TH1F ("mu1ConeChP ","mu1ConeChP ",1000,0, 500);
  mu1ConeChEtH     = new TH1F ("mu1ConeChEt","mu1ConeChEt",1000,0, 500);
  mu1ConeChEH      = new TH1F ("mu1ConeChE ","mu1ConeChE ",1000,0, 500);
  mu1ConeChNPartH  = new TH1F ("mu1ConeChNPart","mu1ConeChNPart",100,-0.5,99.5);
  mu2ConeChPtH     = new TH1F ("mu2ConeChPt","mu2ConeChPt",1000,0, 500);
  mu2ConeChPH      = new TH1F ("mu2ConeChP ","mu2ConeChP ",1000,0, 500);
  mu2ConeChEtH     = new TH1F ("mu2ConeChEt","mu2ConeChEt",1000,0, 500);
  mu2ConeChEH      = new TH1F ("mu2ConeChE ","mu2ConeChE ",1000,0, 500);
  mu2ConeChNPartH  = new TH1F ("mu2ConeChNPart","mu2ConeChNPart",100,-0.5,99.5);

  mu1ConeNePtH     = new TH1F ("mu1ConeNePt","mu1ConeNePt",1000,0, 500);
  mu1ConeNePH      = new TH1F ("mu1ConeNeP ","mu1ConeNeP ",1000,0, 500);
  mu1ConeNeEtH     = new TH1F ("mu1ConeNeEt","mu1ConeNeEt",1000,0, 500);
  mu1ConeNeEH      = new TH1F ("mu1ConeNeE ","mu1ConeNeE ",1000,0, 500);
  mu1ConeNeNPartH  = new TH1F ("mu1ConeNeNPart","mu1ConeNeNPart",100,-0.5,99.5);
  mu2ConeNePtH     = new TH1F ("mu2ConeNePt","mu2ConeNePt",1000,0, 500);
  mu2ConeNePH      = new TH1F ("mu2ConeNeP ","mu2ConeNeP ",1000,0, 500);
  mu2ConeNeEtH     = new TH1F ("mu2ConeNeEt","mu2ConeNeEt",1000,0, 500);
  mu2ConeNeEH      = new TH1F ("mu2ConeNeE ","mu2ConeNeE ",1000,0, 500);
  mu2ConeNeNPartH  = new TH1F ("mu2ConeNeNPart","mu2ConeNeNPart",100,0-0.5,99.5);

  
  TH1::AddDirectory(oldAddDir);

  ntuplefile = new TFile(outfile.c_str(),"recreate");
  
  const bool oldAddDir2 = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);

  ntuple = new TTree("ntuple","MuonsfromHtoWW");
//   ntuple->Branch("run",&(branch.run),"run/I");
//   ntuple->Branch("event",&(branch.event),"event/I");


  ntuple->Branch("muPt",&(branch.muPt));
  ntuple->Branch("muP",&(branch.muP));
  ntuple->Branch("muPx",&(branch.muPx));
  ntuple->Branch("muPy",&(branch.muPy));
  ntuple->Branch("muPz",&(branch.muPz));
  ntuple->Branch("muE",&(branch.muE));
  ntuple->Branch("muEta",&(branch.muEta));
  ntuple->Branch("muPhi",&(branch.muPhi));
  ntuple->Branch("muQ",&(branch.muQ));
  
  ntuple-> Branch("muConeChPt",&(branch.muConeChPt));
  ntuple-> Branch("muConeChP" ,&(branch.muConeChP));
  ntuple-> Branch("muConeChEt",&(branch.muConeChEt));
  ntuple-> Branch("muConeChE" ,&(branch.muConeChE));
  ntuple-> Branch("muConeChNPart",&(branch.muConeChNPart));
  ntuple-> Branch("muConeNePt",&(branch.muConeNePt));
  ntuple-> Branch("muConeNeP",&(branch.muConeNeP));
  ntuple-> Branch("muConeNeEt",&(branch.muConeNeEt));
  ntuple-> Branch("muConeNeE",&(branch.muConeNeE));
  ntuple-> Branch("muConeNeNPart",&(branch.muConeNeNPart));


  TH1::AddDirectory(oldAddDir2);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenQuantities::endJob() {
  file->Write();
  file->Close();
  ntuplefile->Write();
  ntuplefile->Close();

}

//define this as a plug-in
DEFINE_FWK_MODULE(GenQuantities);
