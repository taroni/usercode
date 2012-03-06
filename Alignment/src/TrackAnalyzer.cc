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
#include "TBranch.h"
#include <map>
using namespace std;
typedef struct { 
  vector <double>  Eta;
  vector <double>  Chi;
  vector <double>  P;
  vector <int>  Hit;
  vector <int>  nhpxb;
  vector <int>  nhpxe;
  vector <int>  nhTIB;
  vector <int>  nhTID;
  vector <int>  nhTOB;
  vector <int>  nhTEC;
  vector <double> Pt;	
  vector <double> E_Pt;
  vector <double>  Px;
  vector <double>  Py;
  vector <double>  Pz;
  vector <double> Phi;
  vector <int> muQ;
  vector <double> d0; // d0 = -dxy
  vector <double> dxy; 
  vector <double> dz;
  vector <double> d0Err;
  vector <double> dxyErr;
  vector <double> dzErr;
  vector <double> vx;
  vector <double> vy;
  vector <double> vz;
  vector <bool> globalMu;
//   vector <double> invMass;
//   double invMass;
  vector <pair <string, bool> > triggerInfo;
  vector <vector <double> > eHitFr;
} BRANCH;

using namespace edm;

class TrackAnalyzer : public edm::EDAnalyzer {
 public:
  TrackAnalyzer(const edm::ParameterSet& pset) {
    TkTag_ = pset.getParameter<string>("TkTag");
  }

  ~TrackAnalyzer(){}

  edm::Service<TFileService> fs;
  
  TTree* AnalysisTree;
  BRANCH branch;


  std::string trigNames;
  int trigCount;
  int ievt;
  int mol;
  float InvMass;
  double invMass;
  bool firstEvent_;
  const TrackerGeometry* trackerGeometry_;

  std::string TkTag_;

  virtual void analyze(const edm::Event& event, const edm::EventSetup& setup){
    
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
    vector <double>  Eta;
    vector <double>  Chi;
    vector <double>  P;
    vector <int>  Hit;
    vector <int>  nhpxb;
    vector <int>  nhpxe;
    vector <int>  nhTIB;
    vector <int>  nhTID;
    vector <int>  nhTOB;
    vector <int>  nhTEC;
    vector <double> Pt;	
    vector <double> E_Pt;
    vector <double>  Px;
    vector <double>  Py;
    vector <double>  Pz;
    vector <double> Phi;
    vector <int> muQ;
    vector <double> d0; // d0 = -dxy
    vector <double> dxy; 
    vector <double> dz;
    vector <double> d0Err;
    vector <double> dxyErr;
    vector <double> dzErr;
    vector <double> vx;
    vector <double> vy;
    vector <double> vz;
    vector <bool> globalMu;
    //    vector <double> invMass;
 //     double invMass;
    vector <pair <string, bool> > triggerInfo;
    vector <double> equalHitFraction;
    vector <vector <double> > eHitFr;
   ievt++;
    std::cout << "event#" << ievt << " Event ID = "<< event.id() << std::endl ;


    edm::Handle<reco::TrackCollection> trackCollection;
    event.getByLabel(TkTag_, trackCollection);
//     event.getByType(trackCollection);

//     edm::Handle<reco::TrackCollection> globalMuons;
//     event.getByLabel("globalMuons", globalMuons);

//     const reco::TrackCollection gMu= *(globalMuons.product());

    
    const reco::TrackCollection tC = *(trackCollection.product());

    std::cout << "Reconstructed "<< tC.size() << " tracks" << std::endl ;
    TLorentzVector mother(0.,0.,0.,0.);
    
    int iCounter=0;

    edm::Handle<edm::TriggerResults> hltresults;
    InputTag tag("TriggerResults","","HLT");//InputTag tag("TriggerResults");
    event.getByLabel(tag,hltresults);
    edm::TriggerNames triggerNames_ = event.triggerNames(*hltresults);
    int ntrigs=hltresults->size();
    vector<string> triggernames = triggerNames_.triggerNames();

    for (int itrig = 0; itrig != ntrigs; ++itrig){
      string trigName=triggerNames_.triggerName(itrig);
      bool accept = hltresults->accept(itrig);
//        cout << trigName << " " << accept << endl;
      if (accept== 1){
	triggerInfo.push_back(pair <string, int> (trigName, accept));
      }
    }


    for (reco::TrackCollection::const_iterator track=tC.begin(); track!=tC.end(); track++){
      TLorentzVector track0(track->px(),track->py(),track->pz(),sqrt((track->p()*track->p())+(0.106*0.106)));

        for (reco::TrackCollection::const_iterator track1=track; track1!=tC.end(); track1++){
	trackingRecHit_iterator trkHitsIt ;
	trackingRecHit_iterator trkHitsIt1 ;
	int nequalHits=0; 
	for ( trkHitsIt = track->recHitsBegin() ; trkHitsIt != track->recHitsEnd() ;trkHitsIt++){
	  int hitCheck=0;
	  for ( trkHitsIt1 = track1->recHitsBegin() ; trkHitsIt1 != track1->recHitsEnd() ;trkHitsIt1++){
	    if(!((**trkHitsIt).isValid())|| !((**trkHitsIt1).isValid())) {
	      continue ; 
	    }
	    if ( trkHitsIt1 == trkHitsIt1) hitCheck++; 
	  }
	  if (hitCheck==1) nequalHits++ ;
	}      
	equalHitFraction.push_back(nequalHits/track->numberOfValidHits());
      }//TRACK1
      eHitFr.push_back(equalHitFraction);
      Hit.push_back(track->numberOfValidHits());	
      nhpxb.push_back( track->hitPattern().numberOfValidPixelBarrelHits());
      nhpxe.push_back(track->hitPattern().numberOfValidPixelEndcapHits());
      nhTIB.push_back(track->hitPattern().numberOfValidStripTIBHits());
      nhTID.push_back(track->hitPattern().numberOfValidStripTIDHits());
      nhTOB.push_back(track->hitPattern().numberOfValidStripTOBHits());
      nhTEC.push_back(track->hitPattern().numberOfValidStripTECHits());
      Pt.push_back(track->pt());
      P.push_back(track->p());
      Chi.push_back(track->normalizedChi2());
      Eta.push_back(track->eta());
      E_Pt.push_back(track->ptError());
      Px.push_back(track->px());
      Py.push_back(track->py());
      Pz.push_back(track->pz());
      Phi.push_back(track->phi());
      muQ.push_back(track ->charge());     
      d0.push_back(track->d0());
      dxy.push_back(track->dxy());
      dz.push_back(track->dz());
      d0Err.push_back(track->d0Error());
      dxyErr.push_back(track->dxyError());
      dzErr.push_back(track->dzError());
      vx.push_back(track->vx());
      vy.push_back(track->vy());
      vz.push_back(track->vz());
      
      iCounter++;
      mother=track0+mother;
//       InvMass= -1.;
//       double iHalf, fractpart, intpart;
//       iHalf = (double)iCounter/2.;
//       fractpart = modf (iHalf , &intpart);
//       if ((double)(iCounter)/2. == intpart){
// 	InvMass=mother.M();
// // 	invMass.push_back(InvMass);
// 	invMass=(double)InvMass;
// 	cout << InvMass<<endl;
// 	mother.SetXYZM(0.,0.,0.,0.);
// 	cout << __LINE__ << " " << (double)(iCounter)/2. << " " << intpart << endl;
//       }

    }
    mol=iCounter;
//     if(iCounter==2) InvMass=mother.M();
//     cout << InvMass<<endl;


    branch.eHitFr=eHitFr;
    branch.Hit=Hit;
    branch.nhpxb=nhpxb;
    branch.nhpxe=nhpxe;
    branch.nhTIB=nhTIB;
    branch.nhTID=nhTID;
    branch.nhTOB=nhTOB;
    branch.nhTEC=nhTEC;
    branch.Pt=Pt;
    branch.P=P;
    branch.Chi=Chi;
    branch.Eta=Eta;
    branch.E_Pt=E_Pt;
    branch.Px=Px;
    branch.Py=Py;
    branch.Pz=Px;
    branch.Phi=Phi;
    branch.muQ=muQ;

    branch.d0=d0;
    branch.dxy=dxy;
    branch.dz=dz;
    branch.d0Err=d0Err;
    branch.dxyErr=dxyErr;
    branch.dzErr=dzErr;
    branch.vx=vx;
    branch.vy=vy;
    branch.vz=vz;

    branch.globalMu=globalMu;
    branch.triggerInfo=triggerInfo;

//     branch.invMass = invMass;


    AnalysisTree->Fill();
//     return; 


  }

  virtual void beginJob() {
    ievt=0;
    AnalysisTree = fs->make<TTree>("AnalysisTree","ntp1");
    
    AnalysisTree->Branch("Chi",&(branch.Chi));
    AnalysisTree->Branch("P",&(branch.P));
    AnalysisTree->Branch("eqHitFr",&(branch.eHitFr));
    AnalysisTree->Branch("Hit",&(branch.Hit));
    AnalysisTree->Branch("nhpxb",&(branch.nhpxb));
    AnalysisTree->Branch("nhpxe",&(branch.nhpxe));
    AnalysisTree->Branch("nhTIB",&(branch.nhTIB));
    AnalysisTree->Branch("nhTID",&(branch.nhTID));
    AnalysisTree->Branch("nhTOB",&(branch.nhTOB));
    AnalysisTree->Branch("nhTEC",&(branch.nhTEC));
    AnalysisTree->Branch("ErrorPt",&(branch.E_Pt));
    /*AnalysisTree->Branch("Pz",Pz,"Pz[10000]/D");*/
    AnalysisTree->Branch("Pt",&(branch.Pt));
    AnalysisTree->Branch("Px",&(branch.Px));
    AnalysisTree->Branch("Py",&(branch.Py));
    AnalysisTree->Branch("Pz",&(branch.Pz));
    AnalysisTree->Branch("Phi",&(branch.Phi));
    AnalysisTree->Branch("Eta",&(branch.Eta));
    AnalysisTree->Branch("mol",&mol,"mol/I");
    AnalysisTree->Branch("invMass",&invMass,"invMass/D");
    AnalysisTree->Branch("muQ",&(branch.muQ));
    AnalysisTree->Branch("d0",&(branch.d0));
    AnalysisTree->Branch("dxy",&(branch.dxy));
    AnalysisTree->Branch("dz",&(branch.dz));
    AnalysisTree->Branch("d0Err",&(branch.d0Err));
    AnalysisTree->Branch("dxyErr",&(branch.dxyErr));
    AnalysisTree->Branch("dzErr",&(branch.dzErr));
    AnalysisTree->Branch("vx",&(branch.vx));
    AnalysisTree->Branch("vy",&(branch.vy));
    AnalysisTree->Branch("vz",&(branch.vz));

    AnalysisTree->Branch("globalMu",&(branch.globalMu));
    AnalysisTree->Branch("triggerInfo",&(branch.triggerInfo));


    firstEvent_=true;

  }//beginJob

  virtual void endJob()
  { 

  }

};


DEFINE_FWK_MODULE(TrackAnalyzer);

