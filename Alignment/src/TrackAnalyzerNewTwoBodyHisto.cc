#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "IOMC/RandomEngine/src/RandomEngineStateProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"

#include "CondFormats/HLTObjects/interface/AlCaRecoTriggerBits.h"
#include "CondFormats/DataRecord/interface/AlCaRecoTriggerBitsRcd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
// #include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TBranch.h"
#include <map>

#define DEBUG 0

using namespace std;
using namespace edm;

class TrackAnalyzerNewTwoBodyHisto : public edm::EDAnalyzer {
 public:
  TrackAnalyzerNewTwoBodyHisto(const edm::ParameterSet& pset) {
    TkTag_ = pset.getParameter<string>("TkTag");
    theMinMass = pset.getParameter<double>( "minMass" );
    theMaxMass = pset.getParameter<double>( "maxMass" );
   }

  ~TrackAnalyzerNewTwoBodyHisto(){}

  edm::Service<TFileService> fs;
  
  TTree* AnalysisTree;
//   BRANCH branch;

  TH1D *hInvMassJpsi  ;
  TH1D *hInvMassUpsilon  ;
  TH1D *hInvMassZ ;
  vector <pair <string, bool> > triggerInfo;
 
  TH1D *hchi2  ;
  TH1D *hNtrk  ;
  TH1D *hP     ;
  TH1D *hPPlus ;
  TH1D *hPMinus;
  TH1D *hPt    ;
  TH1D *hMinPt ;
  TH1D *hPtPlus;
  TH1D *hPtMinus;
  TH1D *hHit;    
  TH1D *hHitPlus;
  TH1D *hHitMinus;


  TH1D *hEta;
  TH1D *hEtaPlus; 
  TH1D *hEtaMinus;
  TH1D *hPhi;
  TH1D * hPhiDT;
  TH1D * hPhiRPCPlus;
  TH1D * hPhiRPCMinus;
  TH1D * hPhiCSCPlus;
  TH1D * hPhiCSCMinus;
  TH1D *hPhiPlus;
  TH1D *hPhiMinus;

  TH1D *hDeltaPhi;
  TH1D *hDeltaEta;
  TH1D *hDeltaR  ;



  TH1D *hvx  ;
  TH1D *hvy  ;
  TH1D *hvz  ;
  TH1D *hd0  ;

  TH2D *hd0vsphi  ;
  TH2D *hd0vseta  ;
  TH2D *hd0vspt   ;

  TH1D *hnhpxb  ;
  TH1D *hnhpxe  ;
  TH1D *hnhTIB  ;
  TH1D *hnhTID  ;
  TH1D *hnhTOB  ;
  TH1D *hnhTEC  ;


  TH1D *hMultCand ;

  TH1D * hPx ;
  TH1D * hPy ;
  TH1D * hPz ;
  TH1D * hmuQ;
  TH1D * hdxy;
  TH1D * hdz ;

  TH1D * hEtaMother;
  TH1D * hPtMother ;
  TH1D * hPhiMother;


  std::string trigNames;
  int trigCount;
  int ievt;
  int mol;
  float InvMass;
//   double invMass;
  bool firstEvent_;
  const TrackerGeometry* trackerGeometry_;

  std::string TkTag_;
  double theMinMass, theMaxMass;

  virtual void analyze(const edm::Event& event, const edm::EventSetup& setup){
    double pigreco = 3.141592;
    
    if (DEBUG) cout << __LINE__ << endl;
    
    // extract tracker geometry
    //
//     edm::ESHandle<TrackerGeometry> theG;
//     setup.get<TrackerDigiGeometryRecord>().get(theG);

//     if(firstEvent_){
//      edm::ESHandle<TrackerGeometry> tmpTkGeometry;
//      setup.get<TrackerDigiGeometryRecord>().get(tmpTkGeometry);
//      trackerGeometry_=&(*tmpTkGeometry);
//      firstEvent_=false;
//     }

    using namespace std;
    vector <double> invMass;
    ievt++;
    if (DEBUG) cout << __LINE__ << endl;

    edm::Handle<reco::TrackCollection> trackCollection;
    event.getByLabel(TkTag_, trackCollection);
//     event.getByType(trackCollection);

//     edm::Handle<reco::TrackCollection> globalMuons;
//     event.getByLabel("globalMuons", globalMuons);

//     const reco::TrackCollection gMu= *(globalMuons.product());

    if (DEBUG) cout << __LINE__ << endl;
    
    const reco::TrackCollection tC = *(trackCollection.product());\
    
    int MultCand = 0 ;
    if (tC.size() > 2 ) MultCand =1 ; 
    if (DEBUG) cout << __LINE__ << endl;
    hNtrk->Fill(tC.size());
    if (MultCand==1) std::cout << "event#" << ievt << " Event ID = "<< event.id() << std::endl ;
    if (MultCand==1) std::cout << "Reconstructed "<< tC.size() << " tracks" << std::endl ;
    TLorentzVector mother(0.,0.,0.,0.);
    
    if (DEBUG) cout << __LINE__ << endl;
    int iCounter=0;

//     edm::Handle<edm::TriggerResults> hltresults;
//     InputTag tag("TriggerResults","","HLT");//InputTag tag("TriggerResults");
//     event.getByLabel(tag,hltresults);
//     edm::TriggerNames triggerNames_ = event.triggerNames(*hltresults);
//     int ntrigs=hltresults->size();
//     vector<string> triggernames = triggerNames_.triggerNames();

//     for (int itrig = 0; itrig != ntrigs; ++itrig){
//       string trigName=triggerNames_.triggerName(itrig);
//       bool accept = hltresults->accept(itrig);
//       if (ievt < 3 ) cout << trigName << " " << accept << endl;
//       if (accept== 1){
// 	triggerInfo.push_back(pair <string, bool > (trigName, accept));
//       }
//     }
    if (DEBUG) cout << __LINE__ << endl;

    invMass.clear();
    for (reco::TrackCollection::const_iterator track=tC.begin(); track!=tC.end(); track++){
      TLorentzVector track0(track->px(),track->py(),track->pz(),sqrt((track->p()*track->p())+(0.105*0.105))); //old 106
      InvMass=0; 
      if (DEBUG) cout << __LINE__ << endl;
      for (reco::TrackCollection::const_iterator track1=track+1; track1!=tC.end(); track1++){
	TLorentzVector track01(track1->px(),track1->py(),track1->pz(),sqrt((track1->p()*track1->p())+(0.105*0.105)));
      if (DEBUG) cout << __LINE__ << endl;
	mother=track0+track01;
      if (DEBUG) cout << __LINE__ << endl;
	InvMass=mother.M();
	if (MultCand==1) cout << __LINE__ << " " << InvMass <<" " << theMaxMass << " " << theMinMass  << ", pt " << mother.Pt() << endl;
      if (DEBUG) cout << __LINE__ << endl;
// 	if (InvMass < theMaxMass && InvMass > theMinMass ) {
	  if (DEBUG) cout << __LINE__ << endl;
	  invMass.push_back(InvMass);
	  if (DEBUG) cout << __LINE__ << endl;
	  if ( theMaxMass < 5. && theMinMass > 0. )hInvMassJpsi->Fill(InvMass);
	  if ( theMaxMass < 15. && theMinMass > 5. )hInvMassUpsilon->Fill(InvMass);
	  if ( theMaxMass < 150. && theMinMass > 50. )hInvMassZ->Fill(InvMass);
	  if (DEBUG) cout << __LINE__ << endl;
	  hEtaMother->Fill(mother.Eta());
	  if (DEBUG) cout << __LINE__ << endl;
	  hPtMother->Fill(mother.Pt());
	  if (DEBUG) cout << __LINE__ << endl;
	  hPhiMother->Fill(mother.Phi());
	  
	  double etaMu1 = track->eta();
	  double phiMu1 = track->phi();
 	  double ptMu1 = track->pt();
	  
	  double etaMu2 = track1->eta();
	  double phiMu2 = track1->phi();
 	  double ptMu2 = track1->pt();

	  double deltaEta = etaMu1 - etaMu2;
	  double deltaPhi = phiMu1 - phiMu2;
	  if (deltaPhi < - pigreco) deltaPhi = 2*pigreco + deltaPhi;
	  if (deltaPhi >   pigreco) deltaPhi =-2*pigreco + deltaPhi;
	  double R  = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
	  hDeltaPhi -> Fill (deltaPhi);
 	  hDeltaEta -> Fill (deltaEta);
	  hDeltaR   -> Fill (R) ;
// 	}//invmass
	if (DEBUG) cout << __LINE__ << endl;
	
      }//TRACK1
      if (invMass.size()==0) continue;
      if (DEBUG) cout << __LINE__ << endl;

      hHit->Fill(track->numberOfValidHits());
      hnhpxb->Fill( track->hitPattern().numberOfValidPixelBarrelHits());
      hnhpxe->Fill(track->hitPattern().numberOfValidPixelEndcapHits());
      hnhTIB->Fill(track->hitPattern().numberOfValidStripTIBHits());
      hnhTID->Fill(track->hitPattern().numberOfValidStripTIDHits());
      hnhTOB->Fill(track->hitPattern().numberOfValidStripTOBHits());
      hnhTEC->Fill(track->hitPattern().numberOfValidStripTECHits());
      hPt->Fill(track->pt());
      hP->Fill(track->p());
      hchi2->Fill(track->normalizedChi2());
      hEta->Fill(track->eta());
      hPhi->Fill(track->phi());
      if (fabs(track->eta() ) < 0.8 )  hPhiDT->Fill(track->phi());
      if (track->eta()>0.8 && track->eta()  < 1.4 )  hPhiRPCPlus->Fill(track->phi());
      if (track->eta()<-0.8 && track->eta() >- 1.4 )  hPhiRPCMinus->Fill(track->phi());
      if (track->eta()>1.4) hPhiCSCPlus->Fill(track->phi());
      if (track->eta()<-1.4) hPhiCSCMinus->Fill(track->phi());
      hd0->Fill(track->d0());
      hvx->Fill(track->vx());
      hvy->Fill(track->vy());
      hvz->Fill(track->vz());
      hd0vsphi ->Fill(track->phi(),track->d0());
      hd0vseta ->Fill(track->eta(), track->d0());
      hd0vspt  ->Fill(track->pt(), track->d0());
      if (DEBUG) cout << __LINE__ << endl;

      hPx->Fill(track->px());
      hPy->Fill(track->py());
      hPz->Fill(track->pz());
      hmuQ->Fill(track ->charge());     
      if (DEBUG) cout << __LINE__ << endl;

      hdxy->Fill(track->dxy());
      hdz->Fill(track->dz());

      if (DEBUG) cout << __LINE__ << endl;
      iCounter++;

    }

//     return; 


  }

  virtual void beginJob() {
    ievt=0;
    hInvMassJpsi = fs->make<TH1D>("hInvMassJpsi","hInvMassJpsi",400,2,4);
    hInvMassUpsilon = fs->make<TH1D>("hInvMassUpsilon","hInvMassUpsilon",600,8,11);
    hInvMassZ = fs->make<TH1D>("hInvMassZ","hInvMassZ",600,60,120);

    hchi2 = fs->make<TH1D>("hchi2","#chi^{2}/ndf",100,0,5.);
    hNtrk=fs->make<TH1D>("hNtrk","Number of Tracks",200,0.,200.);
    hP=fs->make<TH1D>("hP","Momentum",200,0.,200.);
    hPt=fs->make<TH1D>("hPt","Transverse Momentum",200,0.,200.);
    hHit=fs->make<TH1D>("hHit","Number of hit",30,0,30);
    hEta=fs->make<TH1D>("hEta","Eta",100,-4.,4.);
    hPhi=fs->make<TH1D>("hPhi","Phi",100,-4.,4.);
    hPhiDT=fs->make<TH1D>("hPhiDT","hPhiDT",100,-4.,4.);
    hPhiRPCPlus=fs->make<TH1D>("hPhiRPCPlus","hPhiRPCPlus",100,-4.,4.);
    hPhiRPCMinus=fs->make<TH1D>("hPhiRPCMinus","hPhiRPCMinus",100,-4.,4.);
    hPhiCSCPlus=fs->make<TH1D>("hPhiCSCPlus","hPhiCSCPlus",100,-4.,4.);
    hPhiCSCMinus=fs->make<TH1D>("hPhiCSCMinus","hPhiCSCMinus",100,-4.,4.);
    hvx  =fs->make<TH1D>("hvx"  ,"hvx"  ,100,0.,.5);
    hvy  =fs->make<TH1D>("hvy"  ,"hvy"  ,100,0.,.5);
    hvz  =fs->make<TH1D>("hvz"  ,"hvz"  ,100,-10.,10.);
    hd0  =fs->make<TH1D>("hd0"  ,"hd0"  ,100,-1.,1.);
    hd0vsphi  =fs->make<TH2D>("hd0vsphi"  ,"hd0vsphi"  ,160,-3.20,3.20,100,-1.,1.);
    hd0vseta  =fs->make<TH2D>("hd0vseta"  ,"hd0vseta"  ,160,-3.20,3.20,100,-1.,1.);
    hd0vspt   =fs->make<TH2D>("hd0vspt"   ,"hd0vspt"   ,50,0.,100.,100,-0.25,0.25);
    if (DEBUG) cout << __LINE__ << endl;

    hnhpxb  =fs->make<TH1D>("nhpxb"  ,"nhpxb"  ,10,0.,10.);
    hnhpxe  =fs->make<TH1D>("nhpxe"  ,"nhpxe"  ,10,0.,10.);
    hnhTIB  =fs->make<TH1D>("nhTIB"  ,"nhTIB"  ,20,0.,20.);
    hnhTID  =fs->make<TH1D>("nhTID"  ,"nhTID"  ,20,0.,20.);
    hnhTOB  =fs->make<TH1D>("nhTOB"  ,"nhTOB"  ,20,0.,20.);
    hnhTEC  =fs->make<TH1D>("nhTEC"  ,"nhTEC"  ,20,0.,20.);
    if (DEBUG) cout << __LINE__ << endl;
      

    hPx=fs->make<TH1D>("hPx","hPx",200,-100.,100.);
    hPy=fs->make<TH1D>("hPy","hPy",200,-100.,100.);
    hPz=fs->make<TH1D>("hPz","hPz",200,-100.,100.);
    
    hmuQ=fs->make<TH1D>("hmuQ","hmuQ",5,-2.5,2.5);
    
    hdxy=fs->make<TH1D>("hdxy","hdxy",200,-0.1,0.1);
    hdz =fs->make<TH1D>("hdz","hdz",200,-0.5,0.5);


    hEtaMother=fs->make<TH1D>("hEtaMother","EtaMother",100,-4.,4.);
    hPtMother=fs->make<TH1D>("hPtMother","Transverse Momentum Mother",200,0.,200.);
    hPhiMother=fs->make<TH1D>("hPhiMother","PhiMother",100,-4.,4.);

    hDeltaPhi= fs->make<TH1D>("hDeltaPhi","DeltaPhi",100,-10.,10.);
    hDeltaEta= fs->make<TH1D>("hDeltaEta","DeltaEta",100,-10.,10.);
    hDeltaR  = fs->make<TH1D>("R"  ,"R"  ,100,-1.,10.);

    
    firstEvent_=true;

  }//beginJob

  virtual void endJob()
  { 

  }

};


DEFINE_FWK_MODULE(TrackAnalyzerNewTwoBodyHisto);

