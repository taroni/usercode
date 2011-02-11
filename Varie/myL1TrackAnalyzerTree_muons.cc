#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/myL1TrackAnalyzerTree.h"

//     CLASS IMPLEMENTATION     //

// CONSTRUCTOR
myL1TrackAnalyzerTree::myL1TrackAnalyzerTree(edm::ParameterSet const& conf) : 
  config(conf)
{
  /// Insert here what you need to initialize
}
// DESTRUCTOR
myL1TrackAnalyzerTree::~myL1TrackAnalyzerTree()
{
  /// Insert here what you need to delete
  /// when you close the class instance
}  

//////////
// ANALYZE
void myL1TrackAnalyzerTree::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  std::vector<bool> muRecoIsGlobal;
  std::vector<int> simTrkQ, L1TrkQ, recoTrkQ, simTrkPdgId, L1TrkQNoMatch;
  std::vector< pair <int, int > > genWQ, genMuQ, genWPdgId, genMuPdgId;
  std::vector<size_t> simTrkSize, L1TrkSize, recoTrkSize;
  std::vector<bool> recoTrkIsMuon, recoTrkIsHiggsMuon;
  std::vector<double> genHiggsPt, genHiggsPx,genHiggsPy, genHiggsPz,genHiggsEta,genHiggsPhi, genHiggsVtxX,genHiggsVtxY,genHiggsVtxZ,genHiggsMass,genHiggsE,
    simTrkPt, simTrkPx,simTrkPy, simTrkPz, simTrkEta, simTrkPhi, simTrkVtxX, simTrkVtxY, simTrkVtxZ, simTrkId, 
    L1TrkPt, L1TrkPx, L1TrkPy, L1TrkPz, L1TrkEta, L1TrkPhi, L1TrkVtxX, L1TrkVtxY, L1TrkVtxZ, L1TrkId, 
    L1TrkPtNoMatch, L1TrkPxNoMatch, L1TrkPyNoMatch, L1TrkPzNoMatch, L1TrkEtaNoMatch, L1TrkPhiNoMatch, L1TrkVtxXNoMatch, L1TrkVtxYNoMatch, L1TrkVtxZNoMatch, L1TrkIdNoMatch, 
    recoTrkPt, recoTrkPx, recoTrkPy, recoTrkPz, recoTrkEta, recoTrkPhi, recoTrkVtxX, recoTrkVtxY,recoTrkVtxZ, recoTrkId,
    muRecoTrkPt, muRecoTrkPx, muRecoTrkPy, muRecoTrkPz, muRecoTrkEta, muRecoTrkPhi, muRecoTrkVtxX, muRecoTrkVtxY, muRecoTrkVtxZ, muRecoTrkId;

  std::vector< pair <double, double > > genWPt, genWPx, genWPy, genWPz, genWEta, genWPhi, genWVtxX, genWVtxY, genWVtxZ, genWMass, genWE,
    genMuPt, genMuPx, genMuPy, genMuPz, genMuEta, genMuPhi, genMuVtxX, genMuVtxY, genMuVtxZ, genMuE,
    muHiggsRecoTrkPt, muHiggsRecoTrkPx, muHiggsRecoTrkPy, muHiggsRecoTrkPz, muHiggsRecoTrkEta, muHiggsRecoTrkPhi, muHiggsRecoTrkVtxX, muHiggsRecoTrkVtxY,
    muHiggsRecoTrkVtxZ, muHiggsRecoTrkId; 
  std::vector < pair <int, int > > muHiggsRecoIsGlobal;

  ///*** Fill Event NUmber Histogram
  h_EvtCnt->Fill(e.id().event()+0.2); /// The +0.2 is to be sure of being in the correct bin

  ////////////////////////
  // GET MAGNETIC FIELD //
  ////////////////////////
  edm::ESHandle<MagneticField> magnet;
  const MagneticField *magnetF;
  
  es.get<IdealMagneticFieldRecord>().get(magnet);
  magnetF = magnet.product();
  double mMagneticFieldStrength = magnetF->inTesla(GlobalPoint(0,0,0)).z();

  /////////////////////////////////////////
  // First of all, we need to look       //
  // at each vertex: this is to check    //
  // that we have that each single track //
  // event does not provide us with      //
  // (0,0,0) tracks but has really       //
  // a wide luminous region...           //
  // Moreover, here we put some counting //
  // of multiplicity and leading track   //
  /////////////////////////////////////////
  /// Get SimTracks and their Vertices
  /// NOTE: this is good also for later on,
  /// when checking efficiencies etc ...
  edm::Handle<edm::SimTrackContainer> theSimTracks;
  edm::Handle<edm::SimVertexContainer> theSimVtx;
  
  edm::InputTag simtrackHitsTag = config.getParameter<edm::InputTag>("simtrackHits");
  e.getByLabel( simtrackHitsTag, theSimTracks ); //e.getByLabel( "famosSimHits", theSimTracks );
  e.getByLabel( simtrackHitsTag, theSimVtx );    //e.getByLabel( "famosSimHits", theSimVtx );

  /// Get the L1 Tracks
  edm::Handle<cmsUpgrades::L1Track_PixelDigi_Collection> l1trackHandlePD;
  edm::InputTag l1tracksPixelDigisTag = config.getParameter<edm::InputTag>("l1tracksPixelDigis");
  e.getByLabel(l1tracksPixelDigisTag, l1trackHandlePD);
  //GetRecoTrack
  edm::Handle<reco::TrackCollection> recoTracks;
  e.getByLabel("generalTracks", recoTracks ); 
  edm::Handle<reco::MuonCollection> recoMuons;
  e.getByLabel("muons", recoMuons ); 
//   // cout << __PRETTY_FUNCTION__ << "  " << __LINE__ << endl;
//   edm::Handle<View<reco::Track> > theRecoTracks;
//   e.getByLabel("ctfWithMaterialTracks", theRecoTracks ); 
//   const View<Track> trackCol = *(theRecoTracks.product());
//   // cout << __PRETTY_FUNCTION__ << "  " << __LINE__ << endl;
  //GetGenParticle
  Handle<GenParticleCollection> genParticles; 
  e.getByLabel("genParticles",genParticles);

//   hrecoTrackNumber->Fill(recoTracks->size());

  const reco::Candidate* higgs  = 0;
  const reco::Candidate* w1     = 0;
  const reco::Candidate* w2     = 0;
  const reco::Candidate* muW1 = 0;
  const reco::Candidate* nuW1 = 0;
  const reco::Candidate* muW2 = 0;
  const reco::Candidate* nuW2 = 0;

  for (GenParticleCollection::const_iterator p = genParticles->begin();p != genParticles->end(); ++p ) {
    if (fabs(p->pdgId())==25 && p->status()==3) {
      higgs = &*p;
    }
  }

  if (higgs == 0){
    cout << __LINE__ << " No Higgs in the event" <<endl;
//     return;
  }
  
  // cout << __LINE__ << endl;
  
  if (higgs!=0){

    genHiggsPx.push_back(higgs->px());
    genHiggsPt.push_back(higgs->pt());
    genHiggsPy.push_back(higgs->py());
    genHiggsPz.push_back(higgs->pz());
    genHiggsEta.push_back(higgs->eta());
    genHiggsPhi.push_back(higgs->phi());
    genHiggsVtxX.push_back(higgs->vx());
    genHiggsVtxY.push_back(higgs->vy());
    genHiggsVtxZ.push_back(higgs->vz());
    genHiggsMass.push_back(higgs->mass());
    genHiggsE.push_back(higgs->energy());
    
    
    for (Candidate::const_iterator higgsD=higgs->begin(); higgsD!=higgs->end();++higgsD){
      if (higgsD->pdgId()==higgs->pdgId() && higgsD->status()==higgs->status()) 
	cout << __LINE__ << "daughter is the same particle of the mother" << endl;
      if (higgsD->pdgId()==higgs->pdgId()) continue;
      
      if (fabs(higgsD->pdgId())!=24) continue;
      if (w1 == 0){
	w1 = &* higgsD;
      } else {
	if (w2!=0) cout << __LINE__ <<" MORE THAN 2 W IN THE EVENT " << e.id()<< endl;
	w2 = &*higgsD;
      }
    }
    
    if (w1 == 0 || w2==0) {
      cout << __LINE__ << " No W in the event " << endl;
      return ;
    } 
    for (Candidate::const_iterator w1D= w1->begin(); w1D!=w1->end(); ++w1D){
      if (fabs(w1D->pdgId())==24) continue;
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
    
    for (Candidate::const_iterator w2D= w2->begin(); w2D!=w2->end(); ++w2D){
      if (fabs(w2D->pdgId())==24) continue;
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
    
    pair < double, double > wPt (w1->pt(),w2->pt()), wPx(w1->px(),w2->px()), wPy(w1->py(),w2->py()), wPz(w1->pz(),w2->pz()),
      wEta(w1->eta(),w2->eta()), wPhi(w1->phi(),w2->phi()), wVtxX (w1->vx(),w2->vx()), wVtxY(w1->vy(),w2->vy()), 
      wVtxZ(w1->vz(),w2->vz()), wMass(w1->mass(),w2->mass()), wE(w1->energy(),w2->energy()), wQ(w1->charge(),w2->charge()), 
      wPdgId(w1->pdgId(),w2->pdgId());

    genWPt.push_back(wPt);
    genWPx.push_back(wPx);
    genWPy.push_back(wPy);
    genWPz.push_back(wPz);
    genWEta.push_back(wEta);
    genWPhi.push_back(wPhi);
    genWVtxX.push_back(wVtxX);
    genWVtxY.push_back(wVtxY);
    genWVtxZ.push_back(wVtxZ);
    genWMass.push_back(wMass);
    genWE.push_back(wE);
    genWQ.push_back(wQ);
    genWPdgId.push_back(wPdgId);
    
    pair < double, double > muPt (muW1->pt(),muW2->pt()), muPx(muW1->px(),muW2->px()), muPy(muW1->py(),muW2->py()), muPz(muW1->pz(),muW2->pz()),
      muEta(muW1->eta(),muW2->eta()), muPhi(muW1->phi(),muW2->phi()), muVtxX (muW1->vx(),muW2->vx()), muVtxY(muW1->vy(),muW2->vy()), 
      muVtxZ(muW1->vz(),muW2->vz()), muE(muW1->energy(),muW2->energy()), muQ(muW1->charge(),muW2->charge()), muPdgId(muW1->pdgId(),muW2->pdgId());
    
    genMuPt.push_back(muPt);
    genMuPx.push_back(muPx);
    genMuPy.push_back(muPy);
    genMuPz.push_back(muPz);
    genMuEta.push_back(muEta);
    genMuPhi.push_back(muPhi);
    genMuVtxX.push_back(muVtxX);
    genMuVtxY.push_back(muVtxY);
    genMuVtxZ.push_back(muVtxZ);
    genMuE.push_back(muE);
    genMuQ.push_back(muQ);
    genMuPdgId.push_back(muPdgId);
  }/// if higgs!=0

  // cout << __LINE__ << endl;

  const reco::Track* mutrk1=0 ;
  const reco::Track* mutrk2=0;
  bool mu1Glb, mu2Glb; 
  double dRmin= 0.2;
  if (muW1 !=0 && muW2 !=0){
    for (reco::MuonCollection::const_iterator recoMu = recoMuons->begin(); recoMu != recoMuons->end(); recoMu++){
      const reco::Track *  t=0;
      if (!recoMu->track()) continue;
      t = &*(recoMu->track());
      if (t==0) continue;
      double dR  = deltaR( muW1->phi(), muW1->eta(), t->phi(), t->eta());
      if (dR<dRmin) {
	cout << __LINE__ << " " << dR << endl;
	dRmin=dR;
	mutrk1 = &*t;
	mu1Glb=recoMu->isGlobalMuon();
      }
    }
    double dRmin= 0.2;
    for (reco::MuonCollection::const_iterator recoMu = recoMuons->begin(); recoMu != recoMuons->end(); recoMu++){
      const reco::Track *  t=0;
      if (!recoMu->track()) continue;
      t = &*(recoMu->track());
      if (t==0) continue;
      if (&*t == mutrk1) continue;
      double dR  = deltaR( muW2->phi(), muW2->eta(), t->phi(), t->eta());
      if (dR<dRmin) {
	cout << __LINE__ << " " << dR << endl;
	dRmin=dR;
	mutrk2 = &*t;
	mu2Glb=recoMu->isGlobalMuon();
      }
    }
  }//muW1 !=0 && muW2 !=0

  // cout << __LINE__ << endl;

  if (mutrk1 !=0 && mutrk2 !=0){
    pair <double, double > muHiggsPt (mutrk1->pt(),mutrk2->pt()), muHiggsPx (mutrk1->px(),mutrk2->px()), muHiggsPy (mutrk1->py(),mutrk2->py()), 
      muHiggsPz (mutrk1->pz(),mutrk2->pz()), muHiggsEta (mutrk1->eta(),mutrk2->eta()),muHiggsPhi (mutrk1->phi(),mutrk2->phi()),muHiggsVx (mutrk1->vx(),mutrk2->vx()),
      muHiggsVy (mutrk1->vy(),mutrk2->vy()),muHiggsVz (mutrk1->vz(),mutrk2->vz());
    muHiggsRecoTrkPt.push_back(muHiggsPt);
    muHiggsRecoTrkPx.push_back(muHiggsPx);
    muHiggsRecoTrkPy.push_back(muHiggsPy);
    muHiggsRecoTrkPz.push_back(muHiggsPz);
    muHiggsRecoTrkEta.push_back(muHiggsEta);
    muHiggsRecoTrkPhi.push_back(muHiggsPhi);
    muHiggsRecoTrkVtxX.push_back(muHiggsVx);
    muHiggsRecoTrkVtxY.push_back(muHiggsVy);
    muHiggsRecoTrkVtxZ.push_back(muHiggsVz);
    pair < int, int > muHiggsIsGlobal(mu1Glb,mu2Glb);
    muHiggsRecoIsGlobal.push_back(muHiggsIsGlobal);

  }


  //Selecting gen muons
  vector < const reco::Candidate* > muons ;   
  for (GenParticleCollection::const_iterator p = genParticles->begin();p != genParticles->end(); ++p ) {
    if (fabs(p->pdgId())==13 && p->status()==1) {
      muons.push_back(&*p);
    }
  }
  vector < const reco::Track* > muonRecoTracks;

  // cout << __LINE__ << endl;

  for (reco::MuonCollection::const_iterator recoMu = recoMuons->begin(); recoMu != recoMuons->end(); recoMu++){
    const reco::Track *  trkTemp=0;
    if (!recoMu->track()) continue;
    trkTemp = &*(recoMu->track());

    
    if (trkTemp==0) continue;
    if (muonRecoTracks.size()==0) {
      muonRecoTracks.push_back(trkTemp);
      muRecoIsGlobal.push_back(recoMu->isGlobalMuon ());
      muRecoTrkPt.push_back(trkTemp->pt());
      muRecoTrkPx .push_back(trkTemp->px());
      muRecoTrkPy.push_back(trkTemp->py());
      muRecoTrkPz.push_back(trkTemp->pz());
      muRecoTrkEta.push_back(trkTemp->eta());
      muRecoTrkPhi.push_back(trkTemp->phi());
      muRecoTrkVtxX.push_back(trkTemp->vx());
      muRecoTrkVtxY.push_back(trkTemp->vy());
      muRecoTrkVtxZ.push_back(trkTemp->vz());
//       muRecoTrkId.push_back(trkTemp->);
    } else {
      int trackCheck=0; 
      for (unsigned int itrk =0; itrk <muonRecoTracks.size() ; itrk++){
	if (muonRecoTracks[itrk] == trkTemp) trackCheck++; 
      }
      if( trackCheck == 0){
	muonRecoTracks.push_back(trkTemp);
	muRecoIsGlobal.push_back(recoMu->isGlobalMuon ());
	muRecoTrkPt.push_back(trkTemp->pt());
	muRecoTrkPx .push_back(trkTemp->px());
	muRecoTrkPy.push_back(trkTemp->py());
	muRecoTrkPz.push_back(trkTemp->pz());
	muRecoTrkEta.push_back(trkTemp->eta());
	muRecoTrkPhi.push_back(trkTemp->phi());
	muRecoTrkVtxX.push_back(trkTemp->vx());
	muRecoTrkVtxY.push_back(trkTemp->vy());
	muRecoTrkVtxZ.push_back(trkTemp->vz());
// 	muRecoTrkId.push_back(trkTemp);
      }
    }
  }//recoMu


  // cout << __LINE__ << endl;
  ////////////////////////////////////////////////
  // Get L1Tracks from Pixel Digis and          //
  // start looking at relevant spectra and      //
  // correlations. Look also for Fake L1Tracks  //
  // which have a Fake Stub of Stubs coming     //
  // from different SimTracks. Check also       //
  // L1Tracks from the Luminous Region and try  //
  // to clean from those coming from outside    //
  ////////////////////////////////////////////////
  /// Go on only if there are L1Tracks from Pixel Digis

  vector <const cmsUpgrades::L1TrackFit *> L1Vector;
  if ( l1trackHandlePD->size() > 0 ) {
    L1TrkSize.push_back(l1trackHandlePD->size());
    hL1size->Fill( l1trackHandlePD->size());
    int validL1Tracks=0;
    int matchedL1Tracks=0;
    /// Loop over L1Tracks
//     cout << __LINE__ << endl;
    L1Track_PixelDigi_Collection::const_iterator iterL1Track;
    for ( iterL1Track = l1trackHandlePD->begin();  iterL1Track != l1trackHandlePD->end();  ++iterL1Track ) {
//     cout << __LINE__ << endl;
      /// Select only chosen window according to match probability
      if ( iterL1Track->probWindow() != probwindow ) continue;
//     cout << __LINE__ << endl;
      /// Select only beamspot-wise L1Tracks
//     cout << __LINE__ <<  " " << iterL1Track->isBeamSpot00() << " " << beamspot00 << endl;
      if ( iterL1Track->isBeamSpot00() != beamspot00 ) continue;
//     cout << __LINE__ << endl;
//     cout << __LINE__ << iterL1Track->whichSeed() << " " << seedsuperlayer << endl;
      /// Select only L1Tracks with chosen seed
      if ( iterL1Track->whichSeed() != seedsuperlayer ) continue;
           L1TrackFit iterL1TrackFit = iterL1Track->fitL1Track(true);
//       L1TrackFit iterL1TrackFit = iterL1Track->trickFitL1Track(true);


//      cout << __LINE__ << endl;
     /// Select only Tracks matching VTX
      if ( fabs(iterL1TrackFit.getVertex().z()) >= 20.0 ) continue;
//     cout << __LINE__ << endl;
//       /// Select only good L1Tracks
//       bool isFakeL1Track = false;
//       if ( iterL1Track->isFake() != 0 ) {
//         if (iterL1Track->numberStubs() == 2) {
//           /// Check if the Tracklet itself is fake
//           std::vector< GlobalStub<Ref_PixelDigi_> > theStubs = iterL1Track->getStubs();
//           if ( theStubs.at(0).isFake() || theStubs.at(1).isFake() ||
//                theStubs.at(0).trackID() != theStubs.at(1).trackID() ) {
//             isFakeL1Track = true;
//           }
//         }
//         else isFakeL1Track = true;
//       }

//       hL1Trk_Number_n->Fill(iterL1Track->numberStubs());

//       if (isFakeL1Track) {
//         hL1Trk_FakeRate_n->Fill(iterL1Track->numberStubs());
//         continue;
//       }


      /// Select only L1Tracks with chosen size
      if ( iterL1Track->numberStubs() != numberstubs ) continue;      
//     cout << __LINE__ << endl;

      /// Select only tight
      htightNumber->Fill(howmanytight);
      if (howmanytight==1) {
        bool tight = false;
        std::vector< GlobalStub<Ref_PixelDigi_> > theStubs = iterL1Track->getStubs();
        for (unsigned int p=0; p<theStubs.size(); p++ ) {
          if ( theStubs.at(p).localStub()->isTight() == true ) tight = true;
        }
        if (tight == false) continue;
      }
      else if (howmanytight==2) {
        bool tight = true;
        std::vector< GlobalStub<Ref_PixelDigi_> > theStubs = iterL1Track->getStubs();
        for (unsigned int p=0; p<theStubs.size(); p++ ) {
          if ( theStubs.at(p).localStub()->isTight() == false ) tight = false;
        }
        if (tight == false) continue;
      }
      validL1Tracks++;
      L1Vector.push_back(&iterL1TrackFit);

      // cout << __LINE__ << endl;

      L1TrkQNoMatch.push_back(iterL1TrackFit.getCharge());
      L1TrkPtNoMatch.push_back(iterL1TrackFit.getMomentum().perp());
      L1TrkPxNoMatch.push_back(cos(iterL1TrackFit.getMomentum().phi())*iterL1TrackFit.getMomentum().perp());
      L1TrkPyNoMatch.push_back(sin(iterL1TrackFit.getMomentum().phi())*iterL1TrackFit.getMomentum().perp());
      double tetha = 2.* atan(exp(-iterL1TrackFit.getMomentum().eta()));
      double pz = iterL1TrackFit.getMomentum().perp()/tan(tetha);
      L1TrkPzNoMatch.push_back(pz);
      L1TrkEtaNoMatch.push_back(iterL1TrackFit.getMomentum().eta());
      L1TrkPhiNoMatch.push_back(iterL1TrackFit.getMomentum().phi());
      L1TrkVtxXNoMatch.push_back(iterL1TrackFit.getVertex().x());
      L1TrkVtxYNoMatch.push_back(iterL1TrackFit.getVertex().y());
      L1TrkVtxZNoMatch.push_back(iterL1TrackFit.getVertex().z());
      L1TrkIdNoMatch.push_back(iterL1Track->simTrkId());

      // cout << __LINE__ << endl;

      /// Fit the Track already done in builder
      /// Compare with SimTracks (correlation)
      /// Go on only if there are SimTracks
      if ( theSimTracks->size() <= 0 )  continue;
      simTrkSize.push_back(theSimTracks->size());
      /// Loop over SimTracks
      

      // cout << __LINE__ << endl;
//       if ( iterL1Track->isFake() != 0 ) continue;
//     cout << __LINE__ << endl;

      SimTrackContainer::const_iterator iterSimTracks;
      for ( iterSimTracks = theSimTracks->begin();  iterSimTracks != theSimTracks->end();  ++iterSimTracks ) {
	int vertexIndex = iterSimTracks->vertIndex();
	const SimVertex& theSimVertex = (*theSimVtx)[vertexIndex];

      // cout << __LINE__ << endl;

	if (iterL1Track->isFake() == 0 && iterL1Track->simTrkId() == iterSimTracks->trackId() ) {
	  matchedL1Tracks++;
	  hSimTrackPt   -> Fill(iterSimTracks->momentum().pt());
	  hSimTrackPx   -> Fill(iterSimTracks->momentum().px());
	  hSimTrackPy   -> Fill(iterSimTracks->momentum().py());
	  hSimTrackPz   -> Fill(iterSimTracks->momentum().pz());
	  hSimTrackEta  -> Fill(iterSimTracks->momentum().eta());
	  hSimTrackPhi  -> Fill(iterSimTracks->momentum().phi());
	  hSimTrackvtxx -> Fill(theSimVertex.position().x());
	  hSimTrackvtxy -> Fill(theSimVertex.position().y());
	  hSimTrackvtxz -> Fill(theSimVertex.position().z());
	  hSimTrackq    -> Fill(iterSimTracks->charge());
	  
	  simTrkQ.push_back(iterSimTracks->charge());
 	  simTrkPdgId.push_back(iterSimTracks->type() );
	  simTrkPt.push_back(iterSimTracks->momentum().pt());
	  simTrkPx.push_back(iterSimTracks->momentum().px());
	  simTrkPy.push_back(iterSimTracks->momentum().py());
	  simTrkPz.push_back(iterSimTracks->momentum().pz());
	  simTrkEta.push_back(iterSimTracks->momentum().eta());
	  simTrkPhi.push_back(iterSimTracks->momentum().phi());
	  simTrkVtxX.push_back(theSimVertex.position().x());
	  simTrkVtxY.push_back(theSimVertex.position().y());
	  simTrkVtxZ.push_back(theSimVertex.position().z());
	  simTrkId.push_back(iterSimTracks->trackId());

      // cout << __LINE__ << endl;

	  hL1TrackPt   -> Fill(iterL1TrackFit.getMomentum().perp());
      // cout << __LINE__ << endl;
	  hL1TrackPx   -> Fill(cos(iterL1TrackFit.getMomentum().phi())*iterL1TrackFit.getMomentum().perp());
	  hL1TrackPy   -> Fill(sin(iterL1TrackFit.getMomentum().phi())*iterL1TrackFit.getMomentum().perp());
      // cout << __LINE__ << endl;
	  double tetha = 2.* atan(exp(-iterL1TrackFit.getMomentum().eta()));
	  double pz = iterL1TrackFit.getMomentum().perp()/tan(tetha);
	  //  	    double pz = - iterL1TrackFit->trkCharge() *iterL1TrackFit->Pt()/(iterL1TrackFit->fitRadius() * angcoeff);
	  hL1TrackPz   -> Fill(pz);
      // cout << __LINE__ << endl;
	  hL1TrackEta  -> Fill(iterL1TrackFit.getMomentum().eta());
	  hL1TrackPhi  -> Fill(iterL1TrackFit.getMomentum().phi());

      // cout << __LINE__ << endl;
// 	  hL1Trackvtxx -> Fill(iterL1TrackFit.getVertex().x());
// 	  hL1Trackvtxy -> Fill(iterL1TrackFit.getVertex().y());
 
      // cout << __LINE__ << endl;
	  hL1Trackvtxz -> Fill(iterL1TrackFit.getVertex().z());
	  hL1Trackq    -> Fill(iterL1TrackFit.getCharge());
	  hDeltaTrackPt   -> Fill(iterSimTracks->momentum().pt()-iterL1TrackFit.getMomentum().perp());
	  hDeltaTrackPx   -> Fill(iterSimTracks->momentum().px()-cos(iterL1TrackFit.getMomentum().phi())*iterL1TrackFit.getMomentum().perp());
	  hDeltaTrackPy   -> Fill(iterSimTracks->momentum().py()-sin(iterL1TrackFit.getMomentum().phi())*iterL1TrackFit.getMomentum().perp());
	  hDeltaTrackPz   -> Fill(iterSimTracks->momentum().pz()-pz);
	  hDeltaTrackEta  -> Fill(iterSimTracks->momentum().eta()-iterL1TrackFit.getMomentum().eta());
	  hDeltaTrackPhi  -> Fill(iterSimTracks->momentum().phi()-iterL1TrackFit.getMomentum().phi());
// 	  hDeltaTrackvtxx -> Fill(theSimVertex.position().x()-iterL1TrackFit.getVertex().x());
// 	  hDeltaTrackvtxy -> Fill(theSimVertex.position().y()-iterL1TrackFit.getVertex().y());
	  hDeltaTrackvtxz -> Fill(theSimVertex.position().z()-iterL1TrackFit.getVertex().z());
	  hDeltaTrackq    -> Fill(iterSimTracks->charge()-iterL1TrackFit.getCharge());
	    
	  hL1Trk_e_ST_e->Fill( iterSimTracks->momentum().eta(), iterL1TrackFit.getMomentum().eta() );
	  hL1Trk_p_ST_p->Fill( iterSimTracks->momentum().phi(), iterL1TrackFit.getMomentum().phi() );
	  hL1Trk_Pt_ST_Pt->Fill( iterSimTracks->momentum().pt(), iterL1TrackFit.getMomentum().perp() );
// 	  hL1Trk_xvtx_ST_xvtx->Fill( theSimVertex.position().x(), iterL1TrackFit.getVertex().x() );
// 	  hL1Trk_yvtx_ST_yvtx->Fill( theSimVertex.position().y(), iterL1TrackFit.getVertex().y() );
	  hL1Trk_zvtx_ST_zvtx->Fill( theSimVertex.position().z(), iterL1TrackFit.getVertex().z() );
	  hL1Trk_q_ST_q->Fill( iterSimTracks->charge(), iterL1TrackFit.getCharge() );
	  //             hL1Trk_chi2rphi_n->Fill( iterL1TrackFit->numberStubs(), iterL1TrackFit->fitChi2RPhi() );
	  //             hL1Trk_chi2rz_n->Fill( iterL1TrackFit->numberStubs(), iterL1TrackFit->fitChi2ZPhi() );
	  
	  L1TrkQ.push_back(iterL1TrackFit.getCharge());
	  L1TrkPt.push_back(iterL1TrackFit.getMomentum().perp());
	  L1TrkPx.push_back(cos(iterL1TrackFit.getMomentum().phi())*iterL1TrackFit.getMomentum().perp());
	  L1TrkPy.push_back(sin(iterL1TrackFit.getMomentum().phi())*iterL1TrackFit.getMomentum().perp());
	  L1TrkPz.push_back(pz);
	  L1TrkEta.push_back(iterL1TrackFit.getMomentum().eta());
	  L1TrkPhi.push_back(iterL1TrackFit.getMomentum().phi());
      // cout << __LINE__ << endl;
	  L1TrkVtxX.push_back(iterL1TrackFit.getVertex().x());
      // cout << __LINE__ << endl;
	  L1TrkVtxY.push_back(iterL1TrackFit.getVertex().y());
      // cout << __LINE__ << endl;
	  L1TrkVtxZ.push_back(iterL1TrackFit.getVertex().z());
      // cout << __LINE__ << endl;
	  L1TrkId.push_back(iterL1Track->simTrkId());

	  continue;
	} /// End of association to SimTrack
      } /// End of loop on SimTracks
    } /// End of loop over L1Tracks
    hValidL1TracksvsL1Number->Fill(l1trackHandlePD->size(), validL1Tracks);
    hMatchedL1TracksvsL1Number->Fill(l1trackHandlePD->size(), matchedL1Tracks);
  } /// End of l1trackHandlePD->size() > 0

      // cout << __LINE__ << endl;
  
  for (unsigned int imu=0 ; imu < muonRecoTracks.size() ; imu++ ) {
    hMuPt   -> Fill(muonRecoTracks[imu]->pt());
    hMuPx   -> Fill(muonRecoTracks[imu]->px());
    hMuPy   -> Fill(muonRecoTracks[imu]->py());
    hMuPz   -> Fill(muonRecoTracks[imu]->pz());
    hMuEta  -> Fill(muonRecoTracks[imu]->eta());
    hMuPhi  -> Fill(muonRecoTracks[imu]->phi());
    hMuvtxx -> Fill(muonRecoTracks[imu]->vx());
    hMuvtxy -> Fill(muonRecoTracks[imu]->vy());
    hMuvtxz -> Fill(muonRecoTracks[imu]->vz());
    hMuq    -> Fill(muonRecoTracks[imu]->charge());
  }
  if (muonRecoTracks.size() != 0) hMusize->Fill(muonRecoTracks.size());
  
  if (mutrk1!=0){
    hgenMuonPt -> Fill(muW1->pt());
    hMuFromHiggsPt   -> Fill(mutrk1->pt());
    hMuFromHiggsPx   -> Fill(mutrk1->px());
    hMuFromHiggsPy   -> Fill(mutrk1->py());
    hMuFromHiggsPz   -> Fill(mutrk1->pz());
    hMuFromHiggsEta  -> Fill(mutrk1->eta());
    hMuFromHiggsPhi  -> Fill(mutrk1->phi());
    hMuFromHiggsvtxx -> Fill(mutrk1->vx());
    hMuFromHiggsvtxy -> Fill(mutrk1->vy());
    hMuFromHiggsvtxz -> Fill(mutrk1->vz());
    hMuFromHiggsq    -> Fill(mutrk1->charge());
  }
  if (mutrk2!= 0) {
    hgenMuonPt -> Fill(muW1->pt());
    hMuFromHiggsPt   -> Fill(mutrk2->pt());
    hMuFromHiggsPx   -> Fill(mutrk2->px());
    hMuFromHiggsPy   -> Fill(mutrk2->py());
    hMuFromHiggsPz   -> Fill(mutrk2->pz());
    hMuFromHiggsEta  -> Fill(mutrk2->eta());
    hMuFromHiggsPhi  -> Fill(mutrk2->phi());
    hMuFromHiggsvtxx -> Fill(mutrk2->vx());
    hMuFromHiggsvtxy -> Fill(mutrk2->vy());
    hMuFromHiggsvtxz -> Fill(mutrk2->vz());
    hMuFromHiggsq    -> Fill(mutrk2->charge());
  }

  for (reco::TrackCollection::const_iterator t= recoTracks->begin();t != recoTracks->end(); t++ ) {
    if (&*t == mutrk2 || &*t == mutrk1 ) continue;
    if (mutrk1!=0) {
      double dR  = deltaR(mutrk1->phi(),mutrk1->eta(), t->phi(), t->eta());
      double isoPt;
      if (dR<0.3) {
	isoPt=+ t->pt() ;
      }
      hIsoPtRel->Fill(isoPt/mutrk1->pt());
      hIsoPt->Fill(isoPt);
    }
    if (mutrk2!=0) {
      double  dR  = deltaR(mutrk2->phi(),mutrk2->eta(), t->phi(), t->eta());
      double isoPt = 0; 
      if (dR<0.3) {
	isoPt=+ t->pt();
      }
      hIsoPtRel->Fill(isoPt/mutrk2->pt());
      hIsoPt->Fill(isoPt);
    }
  }

  dRmin = 0.2;   
//   for (unsigned int iL1Trk = 0; iL1Trk < L1Vector.size(); iL1Trk++){
//     for (unsigned int jL1Trk = 0; jL1Trk < L1Vector.size(); jL1Trk++){
//       if (iL1Trk == jL1Trk ) continue; 
//       double dR  = deltaR(L1Vector[iL1Trk]->getMomentum().phi(), L1Vector[iL1Trk]->getMomentum().eta(), L1Vector[jL1Trk]->getMomentum().phi(), L1Vector[jL1Trk]->getMomentum().eta());
//       double isoPt = 0; 
//       if (dR<0.3) {
// 	isoPt=+ L1Vector[jL1Trk]->getMomentum().perp();
//       }
//       hL1IsoPtRel->Fill(isoPt/L1Vector[iL1Trk]->getMomentum().perp());
//       hL1IsoPt->Fill(isoPt);
//     }
//     hL1TrackNoMatchPt   -> Fill(L1Vector[iL1Trk]->getMomentum().perp());
//     hL1TrackNoMatchPx   -> Fill(cos(L1Vector[iL1Trk]->getMomentum().phi())*L1Vector[iL1Trk]->getMomentum().perp());
//     hL1TrackNoMatchPy   -> Fill(sin(L1Vector[iL1Trk]->getMomentum().phi())*L1Vector[iL1Trk]->getMomentum().perp());
//     double tetha = 2.* atan(exp(-L1Vector[iL1Trk]->getMomentum().eta()));
//     double pz = L1Vector[iL1Trk]->getMomentum().perp()/tan(tetha);
//     //  	    double pz = - L1Vector[iL1Trk]->trkCharge() *L1Vector[iL1Trk]->Pt()/(L1Vector[iL1Trk]->fitRadius() * angcoeff);
//     hL1TrackNoMatchPz   -> Fill(pz);
//     hL1TrackNoMatchEta  -> Fill(L1Vector[iL1Trk]->getMomentum().eta());
//     hL1TrackNoMatchPhi  -> Fill(L1Vector[iL1Trk]->getMomentum().phi());

//     hL1TrackNoMatchvtxx -> Fill(L1Vector[iL1Trk]->gettVertex().x());
//     hL1TrackNoMatchvtxy -> Fill(L1Vector[iL1Trk]->fitVertex().y());
//     hL1TrackNoMatchvtxz -> Fill(L1Vector[iL1Trk]->fitVertex().z());
//     hL1TrackNoMatchq    -> Fill(L1Vector[iL1Trk]->getCharge());
    
//   }
  

  branch.simTrkQ= simTrkQ;
  branch.L1TrkQ= L1TrkQ;
  branch.L1TrkQNoMatch= L1TrkQNoMatch;
  branch.recoTrkQ= recoTrkQ;
  branch.simTrkPdgId= simTrkPdgId;
  branch.genWQ= genWQ; 
  branch.genMuQ= genMuQ; 
  branch.genWPdgId= genWPdgId; 
  branch.genMuPdgId=genMuPdgId;
  branch.simTrkSize= simTrkSize; 
  branch.L1TrkSize= L1TrkSize; 
  branch.recoTrkSize=recoTrkSize;
  branch.recoTrkIsMuon= recoTrkIsMuon; 
  branch.recoTrkIsHiggsMuon=recoTrkIsHiggsMuon;

  branch.genHiggsPt=genHiggsPt; 
  branch.genHiggsPx=genHiggsPx;
  branch.genHiggsPy=genHiggsPy; 
  branch.genHiggsPz=genHiggsPz;
  branch.genHiggsEta=genHiggsEta;
  branch.genHiggsPhi=genHiggsPhi;
  branch.genHiggsVtxX=genHiggsVtxX;
  branch.genHiggsVtxY=genHiggsVtxY;
  branch.genHiggsVtxZ=genHiggsVtxZ;
  branch.genHiggsMass=genHiggsMass;
  branch.genHiggsE=genHiggsE;
  branch.genWPt=genWPt; 
  branch.genWPx=genWPx;
  branch.genWPy=genWPy;
  branch.genWPz=genWPz; 
  branch.genWEta=genWEta; 
  branch.genWPhi=genWPhi; 
  branch.genWVtxX=genWVtxX;
  branch.genWVtxY=genWVtxY;
  branch.genWVtxZ=genWVtxZ; 
  branch.genWMass=genWMass; 
  branch.genWE=genWE;
  branch.genMuPt=genMuPt;
  branch.genMuPx=genMuPx;
  branch.genMuPy=genMuPy; 
  branch.genMuPz=genMuPz; 
  branch.genMuEta=genMuEta; 
  branch.genMuPhi=genMuPhi; 
  branch.genMuVtxX=genMuVtxX;
  branch.genMuVtxY=genMuVtxY;
  branch.genMuVtxZ=genMuVtxZ; 
  branch.genMuE=genMuE;
  branch.simTrkPt=simTrkPt;
  branch.simTrkPx=simTrkPx;
  branch.simTrkPy=simTrkPy;
  branch.simTrkPz=simTrkPz; 
  branch.simTrkEta=simTrkEta; 
  branch.simTrkPhi=simTrkPhi;
  branch.simTrkVtxX=simTrkVtxX;
  branch.simTrkVtxY=simTrkVtxY;
  branch.simTrkVtxZ=simTrkVtxZ; 
  branch.simTrkId=simTrkId; 

  branch.L1TrkPtNoMatch=L1TrkPtNoMatch; 
  branch.L1TrkPxNoMatch=L1TrkPxNoMatch;
  branch.L1TrkPyNoMatch=L1TrkPyNoMatch; 
  branch.L1TrkPzNoMatch=L1TrkPzNoMatch; 
  branch.L1TrkEtaNoMatch=L1TrkEtaNoMatch; 
  branch.L1TrkPhiNoMatch=L1TrkPhiNoMatch;
  branch.L1TrkVtxXNoMatch=L1TrkVtxXNoMatch;
  branch.L1TrkVtxYNoMatch=L1TrkVtxYNoMatch;
  branch.L1TrkVtxZNoMatch=L1TrkVtxZNoMatch;
  branch.L1TrkIdNoMatch=L1TrkIdNoMatch; 

  branch.L1TrkPt=L1TrkPt; 
  branch.L1TrkPx=L1TrkPx;
  branch.L1TrkPy=L1TrkPy; 
  branch.L1TrkPz=L1TrkPz; 
  branch.L1TrkEta=L1TrkEta; 
  branch.L1TrkPhi=L1TrkPhi;
  branch.L1TrkVtxX=L1TrkVtxX;
  branch.L1TrkVtxY=L1TrkVtxY;
  branch.L1TrkVtxZ=L1TrkVtxZ;
  branch.L1TrkId=L1TrkId; 

  branch.recoTrkPt=recoTrkPt; 
  branch.recoTrkPx=recoTrkPx;
  branch.recoTrkPy=recoTrkPy;
  branch.recoTrkPz=recoTrkPz; 
  branch.recoTrkEta=recoTrkEta; 
  branch.recoTrkPhi=recoTrkPhi; 
  branch.recoTrkVtxX=recoTrkVtxX;
  branch.recoTrkVtxY=recoTrkVtxY;
  branch.recoTrkVtxZ=recoTrkVtxZ;
  branch.recoTrkId=recoTrkId;

  branch.muRecoIsGlobal=muRecoIsGlobal;
  branch.muRecoTrkPt=muRecoTrkPt;
  branch.muRecoTrkPx=muRecoTrkPx;
  branch.muRecoTrkPy=muRecoTrkPy;
  branch.muRecoTrkPz=muRecoTrkPz;
  branch.muRecoTrkEta=muRecoTrkEta;
  branch.muRecoTrkPhi=muRecoTrkPhi;
  branch.muRecoTrkVtxX=muRecoTrkVtxX;
  branch.muRecoTrkVtxY=muRecoTrkVtxY;
  branch.muRecoTrkVtxZ=muRecoTrkVtxZ;

  branch.muHiggsRecoIsGlobal=muHiggsRecoIsGlobal;
  branch.muHiggsRecoTrkPx=muHiggsRecoTrkPx;
  branch.muHiggsRecoTrkPy=muHiggsRecoTrkPy;
  branch.muHiggsRecoTrkPz=muHiggsRecoTrkPz;
  branch.muHiggsRecoTrkEta=muHiggsRecoTrkEta;
  branch.muHiggsRecoTrkPhi=muHiggsRecoTrkPhi;
  branch.muHiggsRecoTrkVtxX=muHiggsRecoTrkVtxX;
  branch.muHiggsRecoTrkVtxY=muHiggsRecoTrkVtxY;
  branch.muHiggsRecoTrkVtxZ=muHiggsRecoTrkVtxZ;

  ntuple->Fill();
  

} /// End of analyze()

/////////////////////////////////
//                             //
// SOME OTHER CUSTOM FUNCTIONS //
//                             //
/////////////////////////////////


/////////////////////
// Calculate DeltaPhi
double myL1TrackAnalyzerTree::deltaPhi( double phi1, double phi2 ){
  double dp = phi1 - phi2;
  double pigreco = 4.0*atan(1.0);
  if (fabs(dp) < pigreco) return fabs(dp);
  else return (2*pigreco - fabs(dp));
}
double myL1TrackAnalyzerTree::deltaPhiNP( double phi1, double phi2 ){
  double deltaPhiNP = phi1 - phi2;
  double pigreco = 4.0*atan(1.0);
  if ( fabs(deltaPhiNP) >= pigreco) {
    if ( deltaPhiNP>0 ) deltaPhiNP = deltaPhiNP - 2*pigreco;
    else deltaPhiNP = 2*pigreco - fabs(deltaPhiNP);
  }
  return deltaPhiNP;
}

/////////////////////
// Calculate DeltaEta
double myL1TrackAnalyzerTree::deltaEta( double eta1, double eta2 ){
  return fabs(eta1 - eta2);
}

///////////////////
// Calculate DeltaR
double myL1TrackAnalyzerTree::deltaR( double phi1, double eta1, double phi2, double eta2 ){
  double dp2 = deltaPhi(phi1, phi2)*deltaPhi(phi1, phi2);
  double de2 = deltaEta(eta1, eta2)*deltaEta(eta1, eta2);
  return sqrt(dp2 + de2);
}


//////////////
// L1TRACKS //
//////////////

///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_ANOTHER_FWK_MODULE(myL1TrackAnalyzerTree);


// BEGIN JOB
void myL1TrackAnalyzerTree::beginJob(const edm::EventSetup& es)
{
  /// Initialize all slave variables
  /// mainly histogram ranges and resolution
  double trkMPt = 99999.9;
  int trkBn = 1;
  double trkltMPt = 99999.9;
  int trkltBn = 1;
  double rBnW = 0.087/2.0;
  int rBn = 40;
  PDGCode = 0;
  lumiZ = 10.0;
  beamspot00 = true;
  probwindow = -99.9;
  testedGeometry = false;
  seedsuperlayer = 0;
  numberstubs = -9;
  howmanytight = 0;
  /// Select the Superlayer you want to check
  if ( config.getParameter<int>("seedSuperLayer") >= 0  ) seedsuperlayer = config.getParameter<int>("seedSuperLayer");
  /// Select the number of Stubs of the L1Track
  if ( config.getParameter<int>("numberStubs") >= 0  ) numberstubs = config.getParameter<int>("numberStubs");
  /// Select the kind of Tracklets you want to handle
  if ( config.getParameter<string>("trackletVTX") == "offcenter" ) beamspot00 = false;
  /// Select the windows corresponding to chosen match probabilities
  if ( config.getParameter<double>("windowSize") >= 0 ) probwindow = config.getParameter<double>("windowSize");

  /// Select how many tight Stubs you want in the L1Track
  /// 0 means any
  /// 1 means at least 1
  /// 2 means all the Stubs in the L1Track must be tight
  if ( config.getParameter<int>("tightStubsL1Trk") >= 0 && config.getParameter<int>("tightStubsL1Trk") <= 2 ) howmanytight = config.getParameter<int>("tightStubsL1Trk");

  /// Things to be done before entering the event Loop
  /// Book histograms etc
  edm::Service<TFileService> fs;

  TFileDirectory histoDir = fs->mkdir( "histos" );
  TFileDirectory treeDir = fs->mkdir( "trees" );


  h_EvtCnt = histoDir.make<TH1D>( "h_EvtCnt", "MC Tau Decay", 1000, 0, 1000 );
  h_EvtCnt->Sumw2();

  hSimTrkBeforeCutEta =  histoDir.make<TH1D> ("hSimTrackBeforeCutEta","hSimTrackBeforeCutEta",300, -3, 3);

  hSimTrackPt  = histoDir.make<TH1D> ("hSimTrackPt","SimTrackPt",100, 0,200);
  hSimTrackPx  = histoDir.make<TH1D> ("hSimTrackPx","SimTrackPx",100, 0,200);
  hSimTrackPy  = histoDir.make<TH1D> ("hSimTrackPy","SimTrackPy",100, 0,200);
  hSimTrackPz  = histoDir.make<TH1D> ("hSimTrackPz","SimTrackPz",100, 0,200);
  hSimTrackEta = histoDir.make<TH1D> ("hSimTrackEta","hSimTrackEta",300, -3, 3);
  hSimTrackPhi = histoDir.make<TH1D> ("hSimTrackPhi","SimTrackPhi",320, -3.2, 3.2);
  hSimTrackvtxx= histoDir.make<TH1D> ("hSimTrackvtxx","SimTrackvtxx",100, -.1, .1);
  hSimTrackvtxy= histoDir.make<TH1D> ("hSimTrackvtxy","SimTrackvtxy",100, -.1, .1);
  hSimTrackvtxz= histoDir.make<TH1D> ("hSimTrackvtxz","SimTrackvtxz",100, -10, 10);
  hSimTrackq   = histoDir.make<TH1D> ("hSimTrackq","SimTrackq",10, -5, 5);
	      
  hL1TrackPt  = histoDir.make<TH1D> ("hL1TrackPt","L1TrackPt",100, 0,200);
  hL1TrackPx  = histoDir.make<TH1D> ("hL1TrackPx","L1TrackPx",100, 0,200);
  hL1TrackPy  = histoDir.make<TH1D> ("hL1TrackPy","L1TrackPy",100, 0,200);
  hL1TrackPz  = histoDir.make<TH1D> ("hL1TrackPz","L1TrackPz",100, 0,200);
  hL1TrackEta = histoDir.make<TH1D> ("hL1TrackEta","L1TrackEta",300, -3, 3);
  hL1TrackPhi = histoDir.make<TH1D> ("hL1TrackPhi","L1TrackPhi",320, -3.2, 3.2);
  hL1Trackvtxz= histoDir.make<TH1D> ("hL1Trackvtxx","L1Trackvtxx",100, -.1, .1);
  hL1Trackvtxz= histoDir.make<TH1D> ("hL1Trackvtxy","L1Trackvtxy",100, -.1, .1);
  hL1Trackvtxz= histoDir.make<TH1D> ("hL1Trackvtxz","L1Trackvtxz",100, -10, 10);
  hL1Trackq   = histoDir.make<TH1D> ("hL1Trackq","L1Trackq",10, -5, 5);

  hDeltaTrackPt  = histoDir.make<TH1D> ("hDeltaTrackPt" ,"#Delta Pt"  ,200, -1.,1.);
  hDeltaTrackPx  = histoDir.make<TH1D> ("hDeltaTrackPx" ,"#Delta Px"  ,200, -1.,1.);
  hDeltaTrackPy  = histoDir.make<TH1D> ("hDeltaTrackPy" ,"#Delta Py"  ,200, -1.,1.);
  hDeltaTrackPz  = histoDir.make<TH1D> ("hDeltaTrackPz" ,"#Delta Pz"  ,200, -1.,1.);
  hDeltaTrackEta = histoDir.make<TH1D> ("hDeltaTrackEta","#Delta #eta",200, -0.01, 0.01);
  hDeltaTrackPhi = histoDir.make<TH1D> ("hDeltaTrackPhi","#Delta #phi",200, -0.01, 0.01);
  hDeltaTrackvtxx= histoDir.make<TH1D> ("hDeltaTrackvtxx","#Delta vtxx",100, -0.2, 0.2);
  hDeltaTrackvtxy= histoDir.make<TH1D> ("hDeltaTrackvtxy","#Delta vtxy",100, -0.2, 0.2);
  hDeltaTrackvtxz= histoDir.make<TH1D> ("hDeltaTrackvtxz","#Delta vtxz",100, -0.2, 0.2);
  hDeltaTrackq   = histoDir.make<TH1D> ("hDeltaTrackq","#Delta q",7, -3.5, 3.5);

  hRecoTrackSize = histoDir.make<TH1D> ("recoTrackSize", "recoTrackSize", 100, 0, 100);
  hRecoTrackPt  = histoDir.make<TH1D> ("hRecoTrackPt" ,"reco Pt"  ,100, 0,200);
  hRecoTrackPx  = histoDir.make<TH1D> ("hRecoTrackPx" ,"reco Px"  ,100, 0,200);
  hRecoTrackPy  = histoDir.make<TH1D> ("hRecoTrackPy" ,"reco Py"  ,100, 0,200);
  hRecoTrackPz  = histoDir.make<TH1D> ("hRecoTrackPz" ,"reco Pz"  ,100, 0,200);
  hRecoTrackEta = histoDir.make<TH1D> ("hRecoTrackEta","reco #eta",300, -3, 3);
  hRecoTrackPhi = histoDir.make<TH1D> ("hRecoTrackPhi","reco #phi",320, -3.2, 3.2);
  hRecoTrackvtxx= histoDir.make<TH1D> ("hRecoTrackvtxx","reco vtxx",100, -.1, .1);
  hRecoTrackvtxy= histoDir.make<TH1D> ("hRecoTrackvtxy","reco vtxy",100, -.1, .1);
  hRecoTrackvtxz= histoDir.make<TH1D> ("hRecoTrackvtxz","reco vtxz",100, -10., 10.);
  hRecoTrackq   = histoDir.make<TH1D> ("hRecoTrackq","reco q",5, -2.5, 2.5);
  
  hMuPt    = histoDir.make<TH1D> ("hMuPt","MuPt",100, 0,200);
  hMuPx    = histoDir.make<TH1D> ("hMuPx","MuPx",100, 0,200);
  hMuPy    = histoDir.make<TH1D> ("hMuPy","MuPy",100, 0,200);
  hMuPz    = histoDir.make<TH1D> ("hMuPz","MuPz",100, 0,200);
  hMuEta   = histoDir.make<TH1D> ("hMuEta","MuEta",300, -3, 3);
  hMuPhi   = histoDir.make<TH1D> ("hMuPhi","MuPhi",320, -3.2, 3.2);
  hMuvtxx  = histoDir.make<TH1D> ("hMuvtxx","Muvtxx",100,-.1, .1);
  hMuvtxy  = histoDir.make<TH1D> ("hMuvtxy","Muvtxy",100,-.1, .1);
  hMuvtxz  = histoDir.make<TH1D> ("hMuvtxz","Muvtxz",100,-10., 10.);
  hMuq     = histoDir.make<TH1D> ("hMuq","Mu charge",5, -2.5, 2.5);
  hMusize  = histoDir.make<TH1D> ("hMusize","Mu size",100, 0, 100);

  hgenMuonPt = histoDir. make <TH1D> ("hgenMuonPt","genMuonPt",100,0,200);

  hMuFromHiggsPt    = histoDir.make<TH1D> ("hMuFromHiggsPt","MuFromHiggsPt",100, 0,200);
  hMuFromHiggsPx    = histoDir.make<TH1D> ("hMuFromHiggsPx","MuFromHiggsPx",100, 0,200);
  hMuFromHiggsPy    = histoDir.make<TH1D> ("hMuFromHiggsPy","MuFromHiggsPy",100, 0,200);
  hMuFromHiggsPz    = histoDir.make<TH1D> ("hMuFromHiggsPz","MuFromHiggsPz",100, 0,200);
  hMuFromHiggsEta   = histoDir.make<TH1D> ("hMuFromHiggsEta","MuFromHiggsEta",300, -3, 3);
  hMuFromHiggsPhi   = histoDir.make<TH1D> ("hMuFromHiggsPhi","MuFromHiggsPhi",320, -3.2, 3.2);
  hMuFromHiggsvtxx  = histoDir.make<TH1D> ("hMuFromHiggsvtxx","MuFromHiggsvtxx",100,-.1, .1);
  hMuFromHiggsvtxy  = histoDir.make<TH1D> ("hMuFromHiggsvtxy","MuFromHiggsvtxy",100,-.1, .1);
  hMuFromHiggsvtxz  = histoDir.make<TH1D> ("hMuFromHiggsvtxz","MuFromHiggsvtxz",100,-10., 10.);
  hMuFromHiggsq     = histoDir.make<TH1D> ("hMuFromHiggsq","MuFromHiggs charge",5, -2.5, 2.5);

  htightNumber = histoDir.make<TH1D> ("htightNumber","tightNumber",10,0,10);
  hValidL1TracksvsL1Number = histoDir.make<TH2D> ("hValidL1TracksvsL1Number","hValidL1TracksvsL1Number",20, 0, 20, 20, 0,20);
  hMatchedL1TracksvsL1Number = histoDir.make<TH2D> ("hMatchedL1TracksvsL1Number","hmatchedL1TracksvsL1Number",11, -0.5, 10.5, 101, -0.5,100);
  
  hrecoTrackNumber = histoDir.make<TH1D> ("hrecoTrackNumber","recoTrackNumber", 100, 0, 100);
  hIsoPtRel = histoDir.make<TH1D> ("hIsoPtRel","hIsoPtRel",200, 0, 10);
  hIsoPt    = histoDir.make<TH1D> ("hIsoPt"   ,"hIsoPt"   ,200, 0, 100);
  hL1IsoPtRel = histoDir.make <TH1D>("hL1IsoPtRel","hL1IsoPtRel",200,0,10); 
  hL1IsoPt    = histoDir.make <TH1D>("hL1IsoPt"   ,"hL1IsoPt"   ,200,0,100); 

  hL1TrackNoMatchPt  = histoDir.make<TH1D> ("hL1TrackNoMatchPt","L1TrackNoMatchPt",100, 0,200);
  hL1TrackNoMatchPx  = histoDir.make<TH1D> ("hL1TrackNoMatchPx","L1TrackNoMatchPx",100, 0,200);
  hL1TrackNoMatchPy  = histoDir.make<TH1D> ("hL1TrackNoMatchPy","L1TrackNoMatchPy",100, 0,200);
  hL1TrackNoMatchPz  = histoDir.make<TH1D> ("hL1TrackNoMatchPz","L1TrackNoMatchPz",100, 0,200);
  hL1TrackNoMatchEta = histoDir.make<TH1D> ("hL1TrackNoMatchEta","L1TrackNoMatchEta",300, -3, 3);
  hL1TrackNoMatchPhi = histoDir.make<TH1D> ("hL1TrackNoMatchPhi","L1TrackNoMatchPhi",320, -3.2, 3.2);
  hL1TrackNoMatchvtxz= histoDir.make<TH1D> ("hL1TrackNoMatchvtxx","L1TrackNoMatchvtxx",100, -.1, .1);
  hL1TrackNoMatchvtxz= histoDir.make<TH1D> ("hL1TrackNoMatchvtxy","L1TrackNoMatchvtxy",100, -.1, .1);
  hL1TrackNoMatchvtxz= histoDir.make<TH1D> ("hL1TrackNoMatchvtxz","L1TrackNoMatchvtxz",100, -10, 10);
  hL1TrackNoMatchq   = histoDir.make<TH1D> ("hL1TrackNoMatchq","L1TrackNoMatchq",10, -5, 5);
  hL1size= histoDir.make<TH1D> ("hL1size","L1size",100,0,100);

  hDET_LAYER_r = histoDir.make<TH2D>( "hDET_LAYER_r",  "Layer Number vs radius", 220, 0, 110, 12, 0, 12 );
  hDET_IPHI_p = histoDir.make<TH2D>( "hDET_IPHI_p",  "iPhi Number vs #phi", 660, -3.3, 3.3, 70, 0, 70 );
  hDET_IZ_z = histoDir.make<TH2D>( "hDET_IZ_z",  "iZ Number vs z", 580, -290, 290, 80, 0, 80 );
  hDET_LAYER_r->Sumw2();
  hDET_IPHI_p->Sumw2();
  hDET_IZ_z->Sumw2();

  hL1Trk_e_ST_e = histoDir.make<TH2D>( "hL1Trk_e_ST_e",  "L1Trk #eta vs SimTrack #eta", 90, -3, 3, 90, -3, 3 );
  hL1Trk_p_ST_p = histoDir.make<TH2D>( "hL1Trk_p_ST_p",  "L1Trk #phi vs SimTrack #phi", 99, -3.3, 3.3, 99, -3.3, 3.3 );
  hL1Trk_Pt_ST_Pt = histoDir.make<TH2D>( "hL1Trk_Pt_ST_Pt",  "L1Trk p_{T} vs SimTrack p_{T}", 440, 0, 110, 440, 0, 110 );
  hL1Trk_xvtx_ST_xvtx = histoDir.make<TH2D>( "hL1Trk_xvtx_ST_xvtx",  "L1Trk x_{vtx} vs SimTrack x_{vtx}", 99, -3.3, 3.3, 99, -3.3, 3.3 );
  hL1Trk_yvtx_ST_yvtx = histoDir.make<TH2D>( "hL1Trk_yvtx_ST_yvtx",  "L1Trk y_{vtx} vs SimTrack y_{vtx}", 99, -3.3, 3.3, 99, -3.3, 3.3 );
  hL1Trk_zvtx_ST_zvtx = histoDir.make<TH2D>( "hL1Trk_zvtx_ST_zvtx",  "L1Trk z_{vtx} vs SimTrack z_{vtx}", 300, -20, 20, 300, -20, 20 );
  hL1Trk_q_ST_q = histoDir.make<TH2D>( "hL1Trk_q_ST_q",  "L1Trk charge vs SimTrack charge", 9, -4.5, 4.5, 9, -4.5, 4.5 );
  hL1Trk_chi2rphi_n = histoDir.make<TH2D>( "hL1Trk_chi2rphi_n",  "L1Trk #chi^{2}_{r#phi} vs num of Stubs", 10, -0.5, 9.5, 150, 0, 0.03 );
  hL1Trk_chi2rz_n = histoDir.make<TH2D>( "hL1Trk_chi2rz_n",  "L1Trk #chi^{2}_{rz} vs num of Stubs", 10, -0.5, 9.5, 150, 0, 0.03 );
  hL1Trk_e_ST_e->Sumw2();
  hL1Trk_p_ST_p->Sumw2();
  hL1Trk_Pt_ST_Pt->Sumw2();
  hL1Trk_xvtx_ST_xvtx->Sumw2();
  hL1Trk_yvtx_ST_yvtx->Sumw2();
  hL1Trk_zvtx_ST_zvtx->Sumw2();
  hL1Trk_q_ST_q->Sumw2();
  hL1Trk_chi2rphi_n->Sumw2();
  hL1Trk_chi2rz_n->Sumw2();

  hL1Trk_FakeRate_n = histoDir.make<TH1D>( "hL1Trk_FakeRate_n", "Fake Rate vs num of Stubs", 10, -0.5, 9.5 );
  hL1Trk_Number_n = histoDir.make<TH1D>( "hL1Trk_Number_n", "Num of L1Tracks vs num of Stubs", 10, -0.5, 9.5 );
  hST_Pt_matched = histoDir.make<TH1D>( "hST_Pt_matched", "L1Tracks p_{T} matched", 440, 0, 110 );
  hST_Pt_matched_any = histoDir.make<TH1D>( "hST_Pt_matched_any", "L1Tracks p_{T} matched", 440, 0, 110 );
  hST_Pt = histoDir.make<TH1D>( "hST_Pt", "L1Tracks p_{T}", 440, 0, 110 );
  hL1Trk_FakeRate_n->Sumw2();
  hL1Trk_Number_n->Sumw2();
  hST_Pt_matched->Sumw2();
  hST_Pt_matched_any->Sumw2();
  hST_Pt->Sumw2();


  hST_e_Pt10_matched = histoDir.make<TH1D>( "hST_e_Pt10_matched", "L1Tracks #eta matched, p_{T}>10 GeV/c", 90, -3, 3 );
  hST_e_Pt10_matched_any = histoDir. make<TH1D>( "hST_e_Pt10_matched_any", "L1Tracks p_{T} #eta matched, p_{T}>10 GeV/c", 90, -3, 3 );
  hST_e_Pt10 = histoDir. make<TH1D>( "hST_e_Pt10", "L1Tracks #eta, p_{T}>10 GeV/c", 90, -3, 3 );
  hST_e_Pt10_matched->Sumw2();
  hST_e_Pt10_matched_any->Sumw2();
  hST_e_Pt10->Sumw2();
  
  /// Geometry setup
  /// Set pointers to Geometry
  es.get<TrackerDigiGeometryRecord>().get(geometryESH);
  theGeometry = &(*geometryESH);
  /// Set pointers to Stacked Modules
  es.get<StackedTrackerGeometryRecord>().get(stackedgeometryESH);
  theStackedGeometry = stackedgeometryESH.product(); /// Note this is different 
                                                     /// from the "global" geometry

  ntuple = treeDir.make<TTree>("ntuple","HtoWWvariable" ) ;

  ntuple ->Branch("simTrkQ", &(branch.simTrkQ));
  ntuple -> Branch("L1TrkQ", &(branch.L1TrkQ));
  ntuple -> Branch("L1TrkQNoMatch", &(branch.L1TrkQNoMatch));
  ntuple -> Branch("recoTrkQ", &(branch.recoTrkQ));
  ntuple -> Branch("simTrkPdgId", &(branch.simTrkPdgId));
  ntuple -> Branch("genWQ", &(branch.genWQ)); 
  ntuple -> Branch("genMuQ", &(branch.genMuQ)); 
  ntuple -> Branch("genWPdgId", &(branch.genWPdgId)); 
  ntuple -> Branch("genMuPdgId", &(branch.genMuPdgId));
  ntuple -> Branch("simTrkSize", &(branch.simTrkSize)); 
  ntuple -> Branch("L1TrkSize", &(branch.L1TrkSize)); 
  ntuple -> Branch("recoTrkSize", &(branch.recoTrkSize));
  ntuple -> Branch("recoTrkIsMuon", &(branch.recoTrkIsMuon)); 
  ntuple -> Branch("recoTrkIsHiggsMuon", &(branch.recoTrkIsHiggsMuon));
  ntuple -> Branch("genHiggsPt", &(branch.genHiggsPt)); 
  ntuple -> Branch("genHiggsPx", &(branch.genHiggsPx));
  ntuple -> Branch("genHiggsPy", &(branch.genHiggsPy)); 
  ntuple -> Branch("genHiggsPz", &(branch.genHiggsPz));
  ntuple -> Branch("genHiggsEta", &(branch.genHiggsEta));
  ntuple -> Branch("genHiggsPhi", &(branch.genHiggsPhi));
  ntuple -> Branch("genHiggsVtxX", &(branch.genHiggsVtxX));
  ntuple -> Branch("genHiggsVtxY", &(branch.genHiggsVtxY));
  ntuple -> Branch("genHiggsVtxZ", &(branch.genHiggsVtxZ));
  ntuple -> Branch("genHiggsMass", &(branch.genHiggsMass));
  ntuple -> Branch("genHiggsE", &(branch.genHiggsE));
  ntuple -> Branch("genWPt", &(branch.genWPt)); 
  ntuple -> Branch("genWPx", &(branch.genWPx));
  ntuple -> Branch("genWPy", &(branch.genWPy));
  ntuple -> Branch("genWPz", &(branch.genWPz)); 
  ntuple -> Branch("genWEta", &(branch.genWEta)); 
  ntuple -> Branch("genWPhi", &(branch.genWPhi)); 
  ntuple -> Branch("genWVtxX", &(branch.genWVtxX));
  ntuple -> Branch("genWVtxY", &(branch.genWVtxY));
  ntuple -> Branch("genWVtxZ", &(branch.genWVtxZ)); 
  ntuple -> Branch("genWMass", &(branch.genWMass)); 
  ntuple -> Branch("genWE", &(branch.genWE));
  ntuple -> Branch("genMuPt", &(branch.genMuPt));
  ntuple -> Branch("genMuPx", &(branch.genMuPx));
  ntuple -> Branch("genMuPy", &(branch.genMuPy)); 
  ntuple -> Branch("genMuPz", &(branch.genMuPz)); 
  ntuple -> Branch("genMuEta", &(branch.genMuEta)); 
  ntuple -> Branch("genMuPhi", &(branch.genMuPhi)); 
  ntuple -> Branch("genMuVtxX", &(branch.genMuVtxX));
  ntuple -> Branch("genMuVtxY", &(branch.genMuVtxY));
  ntuple -> Branch("genMuVtxZ", &(branch.genMuVtxZ)); 
  ntuple -> Branch("genMuE", &(branch.genMuE));
  ntuple -> Branch("simTrkPt", &(branch.simTrkPt));
  ntuple -> Branch("simTrkPx", &(branch.simTrkPx));
  ntuple -> Branch("simTrkPy", &(branch.simTrkPy));
  ntuple -> Branch("simTrkPz", &(branch.simTrkPz)); 
  ntuple -> Branch("simTrkEta", &(branch.simTrkEta)); 
  ntuple -> Branch("simTrkPhi", &(branch.simTrkPhi));
  ntuple -> Branch("simTrkVtxX", &(branch.simTrkVtxX));
  ntuple -> Branch("simTrkVtxY", &(branch.simTrkVtxY));
  ntuple -> Branch("simTrkVtxZ", &(branch.simTrkVtxZ)); 
  ntuple -> Branch("simTrkId", &(branch.simTrkId)); 

  ntuple -> Branch("L1TrkPt", &(branch.L1TrkPt)); 
  ntuple -> Branch("L1TrkPx", &(branch.L1TrkPx));
  ntuple -> Branch("L1TrkPy", &(branch.L1TrkPy)); 
  ntuple -> Branch("L1TrkPz", &(branch.L1TrkPz)); 
  ntuple -> Branch("L1TrkEta", &(branch.L1TrkEta)); 
  ntuple -> Branch("L1TrkPhi", &(branch.L1TrkPhi));
  ntuple -> Branch("L1TrkVtxX", &(branch.L1TrkVtxX));
  ntuple -> Branch("L1TrkVtxY", &(branch.L1TrkVtxY));
  ntuple -> Branch("L1TrkVtxZ", &(branch.L1TrkVtxZ));
  ntuple -> Branch("L1TrkPtNoMatch", &(branch.L1TrkPtNoMatch)); 
  ntuple -> Branch("L1TrkPxNoMatch", &(branch.L1TrkPxNoMatch));
  ntuple -> Branch("L1TrkPyNoMatch", &(branch.L1TrkPyNoMatch)); 
  ntuple -> Branch("L1TrkPzNoMatch", &(branch.L1TrkPzNoMatch)); 
  ntuple -> Branch("L1TrkEtaNoMatch", &(branch.L1TrkEtaNoMatch)); 
  ntuple -> Branch("L1TrkPhiNoMatch", &(branch.L1TrkPhiNoMatch));
  ntuple -> Branch("L1TrkVtxXNoMatch", &(branch.L1TrkVtxXNoMatch));
  ntuple -> Branch("L1TrkVtxYNoMatch", &(branch.L1TrkVtxYNoMatch));
  ntuple -> Branch("L1TrkVtxZNoMatch", &(branch.L1TrkVtxZNoMatch));
  ntuple -> Branch("L1TrkIdNoMatch", &(branch.L1TrkIdNoMatch)); 
  ntuple -> Branch("L1TrkIdNoMatch", &(branch.L1TrkIdNoMatch)); 

  ntuple -> Branch("recoTrkPt", &(branch.recoTrkPt)); 
  ntuple -> Branch("recoTrkPx", &(branch.recoTrkPx));
  ntuple -> Branch("recoTrkPy", &(branch.recoTrkPy));
  ntuple -> Branch("recoTrkPz", &(branch.recoTrkPz)); 
  ntuple -> Branch("recoTrkEta", &(branch.recoTrkEta)); 
  ntuple -> Branch("recoTrkPhi", &(branch.recoTrkPhi)); 
  ntuple -> Branch("recoTrkVtxX", &(branch.recoTrkVtxX));
  ntuple -> Branch("recoTrkVtxY", &(branch.recoTrkVtxY));
  ntuple -> Branch("recoTrkVtxZ", &(branch.recoTrkVtxZ));
  ntuple -> Branch("recoTrkId", &(branch.recoTrkId));

  ntuple -> Branch("muRecoIsGlobal",   &(branch.muRecoIsGlobal)); 
  ntuple -> Branch("muRecoTrkPt",   &(branch.muRecoTrkPt)); 
  ntuple -> Branch("muRecoTrkPx",   &(branch.muRecoTrkPx));
  ntuple -> Branch("muRecoTrkPy",   &(branch.muRecoTrkPy));
  ntuple -> Branch("muRecoTrkPz",   &(branch.muRecoTrkPz)); 
  ntuple -> Branch("muRecoTrkEta",  &(branch.muRecoTrkEta)); 
  ntuple -> Branch("muRecoTrkPhi",  &(branch.muRecoTrkPhi)); 
  ntuple -> Branch("muRecoTrkVtxX", &(branch.muRecoTrkVtxX));
  ntuple -> Branch("muRecoTrkVtxY", &(branch.muRecoTrkVtxY));
  ntuple -> Branch("muRecoTrkVtxZ", &(branch.muRecoTrkVtxZ));
//   ntuple -> Branch("muRecoTrkId",   &(branch.muRecoTrkId));

  ntuple -> Branch("muHiggsRecoIsGlobal",   &(branch.muHiggsRecoIsGlobal)); 
  ntuple -> Branch("muHiggsRecoTrkPt",   &(branch.muHiggsRecoTrkPt)); 
  ntuple -> Branch("muHiggsRecoTrkPx",   &(branch.muHiggsRecoTrkPx));
  ntuple -> Branch("muHiggsRecoTrkPy",   &(branch.muHiggsRecoTrkPy));
  ntuple -> Branch("muHiggsRecoTrkPz",   &(branch.muHiggsRecoTrkPz)); 
  ntuple -> Branch("muHiggsRecoTrkEta",  &(branch.muHiggsRecoTrkEta)); 
  ntuple -> Branch("muHiggsRecoTrkPhi",  &(branch.muHiggsRecoTrkPhi)); 
  ntuple -> Branch("muHiggsRecoTrkVtxX", &(branch.muHiggsRecoTrkVtxX));
  ntuple -> Branch("muHiggsRecoTrkVtxY", &(branch.muHiggsRecoTrkVtxY));
  ntuple -> Branch("muHiggsRecoTrkVtxZ", &(branch.muHiggsRecoTrkVtxZ));
//   ntuple -> Branch("muHiggsRecoTrkId",   &(branch.muHiggsRecoTrkId));


  /// End of things to be done before entering the event Loop

 
}

// END JOB
void myL1TrackAnalyzerTree::endJob() 
{
  /// Things to be done at the exit of the event Loop
  std::cout << " myL1TrackAnalyzerTree::endJob" << std::endl;
  /// End of things to be done at the exit from the event Loop
}

