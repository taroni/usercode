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
using namespace std;

using namespace edm;

class TrackAnalyzerHisto : public edm::EDAnalyzer {
 public:
  TrackAnalyzerHisto(const edm::ParameterSet& pset) {
    TkTag_ = pset.getParameter<string>("TkTag");
  }

  ~TrackAnalyzerHisto(){}

  edm::Service<TFileService> fs;
  
  TTree* AnalysisTree;
//   BRANCH branch;

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

//     for (int itrig = 0; itrig != ntrigs; ++itrig){
//       string trigName=triggerNames_.triggerName(itrig);
//       bool accept = hltresults->accept(itrig);
// //        cout << trigName << " " << accept << endl;
//       if (accept== 1){
// 	triggerInfo.push_back(pair <string, int> (trigName, accept));
//       }
//     }


    for (reco::TrackCollection::const_iterator track=tC.begin(); track!=tC.end(); track++){

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
      hd0->Fill(track->d0());
      hvx->Fill(track->vx());
      hvy->Fill(track->vy());
      hvz->Fill(track->vz());
      hd0vsphi ->Fill(track->phi(),track->d0());
      hd0vseta ->Fill(track->eta(), track->d0());
      hd0vspt  ->Fill(track->pt(), track->d0());
 
//       if (track ->charge()> 0 ) {
// 	hHitMinus->Fill(track->numberOfValidHits());

//       } else if (track ->charge()< 0 ) {
// 	hHitPlus->Fill(track->numberOfValidHits());

//       }

    }
    hNtrk -> Fill (tC.size());


  }

  virtual void beginJob() {
    ievt=0;

    hchi2 = fs->make<TH1D>("hchi2","#chi^{2}/ndf",100,0,5.);
    hNtrk=fs->make<TH1D>("hNtrk","Number of Tracks",200,0.,200.);
    hP=fs->make<TH1D>("hP","Momentum",200,0.,200.);
//     hPPlus=fs->make<TH1D>("hPPlus","Momentum, #mu^{+}",200,0.,200.);
//     hPMinus=fs->make<TH1D>("hPMinus","Momentum, #mu^{-}",200,0.,200.);
    hPt=fs->make<TH1D>("hPt","Transverse Momentum",200,0.,200.);
//     hMinPt=fs->make<TH1D>("hMinPt","Transverse Momentum, lowest of the muon pair",200,0.,200.);
//     hPtPlus=fs->make<TH1D>("hPtPPlus","Transverse Momentum, #mu^{+}",200,0.,200.);
//     hPtMinus=fs->make<TH1D>("hPtMinus","Transverse Momentum, #mu^{-}",200,0.,200.);
    hHit=fs->make<TH1D>("hHit","Number of hit",30,0,30);
//     hHitPlus=fs->make<TH1D>("hHitPlus","Number of hit, #mu^{+}",30,0,30);
//     hHitMinus=fs->make<TH1D>("hHitMinus","Number of hit, #mu^{-}",30,0,30);

    
    hEta=fs->make<TH1D>("hEta","Eta",100,-4.,4.);
//     hEtaPlus=fs->make<TH1D>("hEtaPlus","#eta, #mu^{+}",100,-4.,4.);
//     hEtaMinus=fs->make<TH1D>("hEtaMinus","#eta, #mu^{-}",100,-4.,4.);
    hPhi=fs->make<TH1D>("hPhi","Phi",100,-4.,4.);
//     hPhiPlus=fs->make<TH1D>("hPhiPlus","#phi, #mu^{+}",100,-4.,4.);
//     hPhiMinus=fs->make<TH1D>("hPhiMinus","#phi, #mu^{-}",100,-4.,4.);
    
//     hDeltaPhi=fs->make<TH1D>("hDeltaPhi","DeltaPhi",100,-10.,10.);
//     hDeltaEta=fs->make<TH1D>("hDeltaEta","DeltaEta",100,-10.,10.);
//     hDeltaR  =fs->make<TH1D>("R"  ,"R"  ,100,-1.,10.);

    
    
    hvx  =fs->make<TH1D>("hvx"  ,"hvx"  ,100,0.,.5);
    hvy  =fs->make<TH1D>("hvy"  ,"hvy"  ,100,0.,.5);
    hvz  =fs->make<TH1D>("hvz"  ,"hvz"  ,100,-10.,10.);
    hd0  =fs->make<TH1D>("hd0"  ,"hd0"  ,100,-1.,1.);
    
    hd0vsphi  =fs->make<TH2D>("hd0vsphi"  ,"hd0vsphi"  ,160,-3.20,3.20,100,-1.,1.);
    hd0vseta  =fs->make<TH2D>("hd0vseta"  ,"hd0vseta"  ,160,-3.20,3.20,100,-1.,1.);
    hd0vspt   =fs->make<TH2D>("hd0vspt"   ,"hd0vspt"   ,50,0.,100.,100,-0.25,0.25);

    hnhpxb  =fs->make<TH1D>("nhpxb"  ,"nhpxb"  ,10,0.,10.);
    hnhpxe  =fs->make<TH1D>("nhpxe"  ,"nhpxe"  ,10,0.,10.);
    hnhTIB  =fs->make<TH1D>("nhTIB"  ,"nhTIB"  ,20,0.,20.);
    hnhTID  =fs->make<TH1D>("nhTID"  ,"nhTID"  ,20,0.,20.);
    hnhTOB  =fs->make<TH1D>("nhTOB"  ,"nhTOB"  ,20,0.,20.);
    hnhTEC  =fs->make<TH1D>("nhTEC"  ,"nhTEC"  ,20,0.,20.);
    
    
//     hMultCand =fs->make<TH1D>("hMultCand","MultipleCandidate",50,-.5,49.5);




    firstEvent_=true;

  }//beginJob

  virtual void endJob()
  { 

  }

};


DEFINE_FWK_MODULE(TrackAnalyzerHisto);

