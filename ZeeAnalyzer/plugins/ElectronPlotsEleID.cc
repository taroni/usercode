// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"


#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/RefToPtr.h"

#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"

class ElectronPlotsEleID : public edm::EDAnalyzer {
public:
  explicit ElectronPlotsEleID(const edm::ParameterSet&);
  ~ElectronPlotsEleID();
  const float getEffectiveArea(float eta) const;
  void printEffectiveAreas() const;

  enum ElectronMatchType {UNMATCHED = 0, 
			  TRUE_PROMPT_ELECTRON, 
			  TRUE_ELECTRON_FROM_TAU,
			  TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  int matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
		   const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);
  
  void findFirstNonElectronMother(const reco::Candidate *particle,
				  int &ancestorPID, int &ancestorStatus);
  
  // ----------member data ---------------------------
  
  bool isMC=false; 
  // Data members
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoProduct_;
  edm::EDGetToken electronsToken_;
  edm::EDGetToken patElectronsToken_;
 
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

  typedef edm::View<reco::GsfElectron>    ElectronView;
  edm::EDGetTokenT<double> theRhoToken;
 
  edm::Handle< double > _rhoHandle;

  std::vector<double> absEtaMin_; // low limit of the eta range
  std::vector<double> absEtaMax_; // upper limit of the eta range
  std::vector<double> effectiveAreaValues_; // effective area for this eta range

  edm::ParameterSet eleIDWP;

  std::vector<int> missHits ;
  std::vector<double> sigmaIEtaIEtaCut;
  std::vector<double> dEtaInSeedCut   ;
  std::vector<double> dPhiInCut       ;
  std::vector<double> hOverECut       ;
  std::vector<double> relCombIso      ;

  // Histos
  TH1F *h_eta, *h_phi, *h_isTrue, *h_eID, *h_Mass, *h_MassMCTruth;
  TH1F *h_EB_pt, *h_EE_pt;
  TH1F *h_EB_rawEne, *h_EE_rawEne;
  TH1F *h_EB_sigmaIeIe, *h_EE_sigmaIeIe;
  TH1F *h_EB_r9, *h_EE_r9;
  TH1F *h_EB_r9uz, *h_EE_r9uz;
  TH1F *h_EB_hoe, *h_EE_hoe;
  TH1F *h_EB_dz, *h_EE_dz;
  TH1F *h_EB_conv, *h_EE_conv;
};

ElectronPlotsEleID::ElectronPlotsEleID(const edm::ParameterSet& iConfig)
{
  absEtaMin_ = iConfig.getParameter<std::vector<double> >("absEtaMin");
  absEtaMax_ = iConfig.getParameter<std::vector<double> >("absEtaMax");
  effectiveAreaValues_= iConfig.getParameter<std::vector<double> >("effectiveAreaValues"); 
  printEffectiveAreas();
  eleIDWP = iConfig.getParameter<edm::ParameterSet>("eleID");

  missHits = eleIDWP.getParameter<std::vector<int> >("missingHitsCut");
  sigmaIEtaIEtaCut= eleIDWP.getParameter<std::vector<double> >("full5x5_sigmaIEtaIEtaCut");
  dEtaInSeedCut   = eleIDWP.getParameter<std::vector<double> >("dEtaInSeedCut");
  dPhiInCut       = eleIDWP.getParameter<std::vector<double> >("dPhiInCut");
  hOverECut       = eleIDWP.getParameter<std::vector<double> >("hOverECut");
  relCombIso      = eleIDWP.getParameter<std::vector<double> >("relCombIsolationWithEALowPtCut");

  std::cout << "missing hits: EB =  " << missHits[0] << ", EE "<< missHits[1] << std::endl;
  std::cout << "sigmaIEtaIEtaCut: EB =  " << sigmaIEtaIEtaCut[0] << ", EE "<< sigmaIEtaIEtaCut[1] << std::endl;
  std::cout << "dPhiInCut: EB =  " << dPhiInCut[0] << ", EE "<< dPhiInCut[1] << std::endl;
  std::cout << "hOverECut: EB =  " << hOverECut[0] << ", EE "<< hOverECut[1] << std::endl;
  std::cout << "relCombIso: EB = " << relCombIso[0] << ", EE "<< relCombIso[1]<< std::endl;
  isMC=iConfig.getParameter<bool>("isMC");

  theRhoToken= consumes <double> (iConfig.getParameter<edm::InputTag>("rho"));

  beamSpotToken_    = consumes<reco::BeamSpot> 
    (iConfig.getParameter <edm::InputTag>
     ("beamSpot"));

  genEventInfoProduct_ = consumes<GenEventInfoProduct> 
    (iConfig.getParameter <edm::InputTag>
     ("genEventInfoProduct"));

  electronsToken_    = mayConsume<ElectronView >
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));
  // electronsToken_    = mayConsume<reco::GsfElectronCollection >
  //   (iConfig.getParameter<edm::InputTag>
  //    ("electrons"));
  // patElectronsToken_    = mayConsume<pat::ElectronCollection >
  //   (iConfig.getParameter<edm::InputTag>
  //    ("patelectrons"));

  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));

  vtxToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("vertices"));

  conversionsToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversions"));


  edm::Service<TFileService> fs;
  h_Mass = fs->make<TH1F>("Mass", "Mass", 100,50.,150.);
  h_MassMCTruth = fs->make<TH1F>("MassMCTruth", "MassMCTruth", 100,50.,150.);
  h_eta = fs->make<TH1F>("eta", "eta", 50,-2.5,2.5);
  h_phi = fs->make<TH1F>("phi", "phi", 50,-2.5,2.5);
  h_isTrue = fs->make<TH1F>("isTrue", "isTrue", 2,-0.5,1.5);
  h_eID = fs->make<TH1F>("eIDSummer16VetoWP", "eIDSummer16VetoWP", 2, -0.5, 1.5);
  h_EB_pt = fs->make<TH1F>("EB_pt", "EB_pt", 50,0.,200.);
  h_EE_pt = fs->make<TH1F>("EE_pt", "EE_pt", 50,0.,200.);
  h_EB_rawEne = fs->make<TH1F>("EB_rawEne", "EB_rawEne", 50,0.,200.);
  h_EE_rawEne = fs->make<TH1F>("EE_rawEne", "EE_rawEne", 100,0.,400.);
  h_EB_sigmaIeIe = fs->make<TH1F>("EB_sigmaIeIe", "EB_sigmaIeIe", 50,0.,0.05);
  h_EE_sigmaIeIe = fs->make<TH1F>("EE_sigmaIeIe", "EE_sigmaIeIe", 50,0.,0.1);
  h_EB_r9 = fs->make<TH1F>("EB_r9", "EB_r9", 50,0.5,1.);
  h_EE_r9 = fs->make<TH1F>("EE_r9", "EE_r9", 50,0.5,1.);
  h_EB_r9uz = fs->make<TH1F>("EB_r9uz", "EB_r9uz", 100,0.,1.);
  h_EE_r9uz = fs->make<TH1F>("EE_r9uz", "EE_r9uz", 100,0.,1.);
  h_EB_hoe = fs->make<TH1F>("EB_hoe", "EB_hoe", 50,0.,0.1);
  h_EE_hoe = fs->make<TH1F>("EE_hoe", "EE_hoe", 50,0.,0.1);
  h_EB_dz = fs->make<TH1F>("EB_dz", "EB_dz", 50,0.,0.1);
  h_EE_dz = fs->make<TH1F>("EE_dz", "EE_dz", 50,0.,0.1);
  h_EB_conv = fs->make<TH1F>("EB_conv", "EB_conv", 2,-0.5,1.5);
  h_EE_conv = fs->make<TH1F>("EE_conv", "EE_conv", 2,-0.5,1.5); 

  typedef reco::Candidate::LorentzVector LorentzVector;
}


ElectronPlotsEleID::~ElectronPlotsEleID() { }

void ElectronPlotsEleID::printEffectiveAreas() const {
 
 
   printf("  eta_min   eta_max    effective area\n");
   uint nEtaBins = absEtaMin_.size();
   for(uint iEta = 0; iEta<nEtaBins; iEta++){
     printf("  %8.4f    %8.4f   %8.5f\n",
        absEtaMin_[iEta], absEtaMax_[iEta],
        effectiveAreaValues_[iEta]);
   }
 
}
const float ElectronPlotsEleID::getEffectiveArea(float eta) const{

  float effArea = 0;
  uint nEtaBins = absEtaMin_.size();
  for(uint iEta = 0; iEta<nEtaBins; iEta++){
    if( std::abs(eta) >= absEtaMin_[iEta]
    && std::abs(eta) < absEtaMax_[iEta] ){
      effArea = effectiveAreaValues_[iEta];
      break;
    }
  }

  return effArea;
}
void ElectronPlotsEleID::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace std;
  using namespace edm;
  using namespace reco;

  iEvent.getByToken(theRhoToken,_rhoHandle);

  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  
  
  // // Get electrons
  //   edm::Handle<reco::GsfElectronCollection>  electrons;
  //   iEvent.getByToken(electronsToken_, electrons);
  // // // Get electrons
  // //  edm::Handle<pat::ElectronCollection > electrons;
  // //  iEvent.getByToken(patElectronsToken_, electrons);
   
  edm::Handle<ElectronView> electrons;
  iEvent.getByToken(electronsToken_, electrons);

  // Get MC collection
  Handle<edm::View<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticlesToken_,genParticles);

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  

  // Find the first vertex in the collection that passes
  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin(); 
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    bool isFake = vtx->isFake();
    // Check the goodness
    if ( !isFake
	 &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }
  if ( firstGoodVertex==vertices->end() ) return; // skip event if there are no good PVs

  // // Get the conversions collection
  // edm::Handle<reco::ConversionCollection> conversions;
  // iEvent.getByToken(conversionsToken_, conversions);


  

  //  bool keepEle=false;
  int i=0; 
  for (ElectronView::const_iterator it = electrons->begin(); it != electrons->end(); it++){
  // for (pat::ElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); it++){
    //std::cout << "cut based id " << it->userInt("cutbasedID_veto")<< std::endl;
    auto etrack  = it->gsfTrack().get(); 
    if (it->isEB()){
      if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>missHits[0]) continue; 
      if (it->full5x5_sigmaIetaIeta() >  sigmaIEtaIEtaCut[0]) continue; 
      if (fabs(it->deltaPhiSuperClusterTrackAtVtx())> dPhiInCut[0]) continue; 
      if (fabs(it->deltaEtaSuperClusterTrackAtVtx())> dEtaInSeedCut[0]) continue; 
      if (it->hadronicOverEm()>hOverECut[0]) continue;
      float abseta = fabs((it->superCluster().get())->position().eta());
      if (abseta> 1.479) continue; // check if it is really needed
      const float  eA = getEffectiveArea( abseta );
      const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle) : 0; 
      if (( it->pfIsolationVariables().sumChargedHadronPt + 
	    std::max(float(0.0), it->pfIsolationVariables().sumNeutralHadronEt +  
		     it->pfIsolationVariables().sumPhotonEt -  eA*rho )
	    ) > relCombIso[0]*it->pt()) continue; 
      //keepEle=true;
    
    }else if (it->isEE()){
      if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>missHits[1]) continue; 
      if (it->full5x5_sigmaIetaIeta() >  sigmaIEtaIEtaCut[1]) continue; 
      if (fabs(it->deltaPhiSuperClusterTrackAtVtx())> dPhiInCut[1]) continue; 
      if (fabs(it->deltaEtaSuperClusterTrackAtVtx())> dEtaInSeedCut[1]) continue; 
      if (it->hadronicOverEm()>hOverECut[1]) continue;
      float abseta = fabs((it->superCluster())->position().eta());
      if (abseta< 1.479) continue; // check if it is really needed
      if (abseta>=2.5  ) continue; // check if it is really needed
      
      const float  eA = getEffectiveArea( abseta );
      const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle) : 0; 
      if ((it->pfIsolationVariables().sumChargedHadronPt + 
	   std::max(float(0.0), it->pfIsolationVariables().sumNeutralHadronEt + 
		    it->pfIsolationVariables().sumPhotonEt - eA*rho)
	   ) >relCombIso[1]*it->pt()) continue; 
      //keepEle=true;
    }
    
    
    int j=0; 
    for (ElectronView::const_iterator jt = it+1; jt != electrons->end(); jt++){
    // for (pat::ElectronCollection::const_iterator jt = it+1; jt != electrons->end(); jt++){
      auto etrack  = jt->gsfTrack().get(); 
      if (jt->isEB()){
	if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>missHits[0]) continue; 
	if (jt->full5x5_sigmaIetaIeta() >  sigmaIEtaIEtaCut[0]) continue; 
	if (fabs(jt->deltaPhiSuperClusterTrackAtVtx())> dPhiInCut[0]) continue; 
	if (fabs(jt->deltaEtaSuperClusterTrackAtVtx())> dEtaInSeedCut[0]) continue; 
	if (jt->hadronicOverEm()>hOverECut[0]) continue;
	float abseta = fabs((jt->superCluster().get())->position().eta());
	if (abseta> 1.479) continue; // check if jt is really needed
	const float  eA = getEffectiveArea( abseta );
	const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle) : 0; 
	if (( jt->pfIsolationVariables().sumChargedHadronPt + 
	      std::max(float(0.0), jt->pfIsolationVariables().sumNeutralHadronEt +  
		       jt->pfIsolationVariables().sumPhotonEt -  eA*rho )
		  ) > relCombIso[0]*jt->pt()) continue; 
	//keepEle=true;
	
      }else if (jt->isEE()){
	if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>missHits[1]) continue; 
	if (jt->full5x5_sigmaIetaIeta() >  sigmaIEtaIEtaCut[1]) continue; 
	if (fabs(jt->deltaPhiSuperClusterTrackAtVtx())> dPhiInCut[1]) continue; 
	if (fabs(jt->deltaEtaSuperClusterTrackAtVtx())> dEtaInSeedCut[1]) continue; 
	if (jt->hadronicOverEm()>hOverECut[1]) continue;
	float abseta = fabs((jt->superCluster())->position().eta());
	if (abseta< 1.479) continue; // check if jt is really needed
	if (abseta>=2.5  ) continue; // check if jt is really needed
	
	const float  eA = getEffectiveArea( abseta );
	const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle) : 0; 
	if ((jt->pfIsolationVariables().sumChargedHadronPt + 
	     std::max(float(0.0), jt->pfIsolationVariables().sumNeutralHadronEt + 
			  jt->pfIsolationVariables().sumPhotonEt - eA*rho)
	     ) >relCombIso[1]*jt->pt()) continue; 
	//keepEle=true;
      }
      //if (keepEle) std::cout << __LINE__ << " "<< keepEle << std::endl; 
      reco::CompositeCandidate diele;
      diele.addDaughter(*it);
      diele.addDaughter(*jt);
      reco::Candidate::LorentzVector zp4 = it->p4() + jt->p4();
      diele.setP4(zp4);
      if (abs( diele.p4().M() - 91)<20) {
	//std::cout << iEvent.id().event() << " electronsize "<< electrons->size()<<  ", mass  " << diele.p4().M() << ", electron "<< i<< ", pt "<< it->pt() << //" ID Veto " <<it->userInt("cutbasedID_veto") << ", electron "<< i+1 << jt->pt() //<< " ID Veto " <<jt->userInt("cutbasedID_veto") << std::endl;   
	h_Mass->Fill( diele.p4().M());
	//const reco::GsfElectron *el1= 0; el1= &*it; 
	//const reco::GsfElectron *el2= 0; el2= &*jt; 
	edm::Ptr<reco::GsfElectron> el1= (edm::Ptr<reco::GsfElectron>) electrons->ptrAt(i); 
	edm::Ptr<reco::GsfElectron> el2= (edm::Ptr<reco::GsfElectron>) electrons->ptrAt(j); ;
	if(bool(matchToTruth(el1, genParticles))==true && bool(matchToTruth( el2, genParticles))==true)
	  h_MassMCTruth->Fill(diele.p4().M());

      }
      j++;

    }
    i++;
    //h_eID ->Fill (it->userInt("cutbasedID_veto")); 
    
    //  }

    //  // Loop over electrons
    //  for (size_t i = 0; i < electrons->size(); ++i){
    //    const auto el = electrons->ptrAt(i);
    const auto el=it; 
    // acceptance
    if( el->pt() < 5 ) continue;
    if( fabs(el->eta()) > 2.5 ) continue;


    // MC truth match
    //    if (isMC==true) h_isTrue -> Fill(matchToTruth( el, genParticles));

    // Kinematics and shower shapes
    float scEta = el->superCluster()->eta();
    float scPhi = el->superCluster()->phi();
    h_eta -> Fill(scEta);
    h_phi -> Fill(scPhi);
    if (fabs(scEta)<1.5) {
      h_EB_pt -> Fill(el->pt());
      h_EB_rawEne -> Fill( el->superCluster()->rawEnergy() );
      h_EB_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ); 
      h_EB_r9 -> Fill( el->r9() );  
      h_EB_r9uz -> Fill( el->r9() );  
      h_EB_hoe -> Fill( el->full5x5_hcalOverEcal() );  
    } else {
      h_EE_pt -> Fill(el->pt());
      h_EE_rawEne -> Fill( el->superCluster()->rawEnergy() );
      h_EE_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ); 
      h_EE_r9 -> Fill( el->r9() );  
      h_EE_r9uz -> Fill( el->r9() );  
      h_EE_hoe -> Fill( el->full5x5_hcalOverEcal() );  
    }

    // Impact parameter
    reco::GsfTrackRef theTrack = el->gsfTrack();
    if (fabs(scEta)<1.5)
      h_EB_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) );
    else
      h_EE_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) );
    
    // // Conversion rejection
    // bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, conversions, theBeamSpot->position());
    // if (fabs(scEta)<1.5)   
    //   h_EB_conv -> Fill( (int) passConvVeto ); 
    // else
    //   h_EE_conv -> Fill( (int) passConvVeto ); 

  } // Loop over electrons
   

}


void ElectronPlotsEleID::beginJob() { }

void ElectronPlotsEleID::endJob() { }

int ElectronPlotsEleID::matchToTruth(const edm::Ptr<reco::GsfElectron> el, const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void ElectronPlotsEleID::findFirstNonElectronMother(const reco::Candidate *particle,
						 int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronPlotsEleID);
