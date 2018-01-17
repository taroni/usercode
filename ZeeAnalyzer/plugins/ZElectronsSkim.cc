// -*- C++ -*-
//
// Original Author:  Silvia Taroni
//         Created:  Wed, 29 Nov 2017 09:06:05 GMT
//
//



#include "SimpleAnalyzer/ZeeAnalyzer/plugins/ZElectronsSkim.h"
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZElectronsSkim::ZElectronsSkim(const edm::ParameterSet& iConfig):
  _effectiveAreas( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()),
  massRange ( iConfig.getParameter<std::vector<double>>("massRange")),
  ptThrs ( iConfig.getParameter<std::vector<double>>("ptcut")),
  theRhoToken(consumes <double> (iConfig.getParameter<edm::InputTag>("rho"))),
  theGsfEToken(consumes <reco::GsfElectronCollection>(iConfig.getParameter <edm::InputTag> ("electrons")))
{
  //produces<reco::GsfElectronCollection>();
}

ZElectronsSkim::~ZElectronsSkim()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZElectronsSkim::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{


  double massMin = massRange[0];
  double massMax = massRange[1];
  double pt_e1_thr = ptThrs[0];
  double pt_e2_thr = ptThrs[1];

  double ptthr = pt_e2_thr; 
  double etathr=2.5; 

  bool keepEvent   = false;
  int nSelectedElectronsEB =0 ; 
  int nSelectedElectronsEE =0 ; 
  using namespace edm;

  edm::Handle< double > _rhoHandle;
  iEvent.getByToken(theRhoToken,_rhoHandle);

  edm::Handle<reco::GsfElectronCollection> pTracks;
  iEvent.getByToken(theGsfEToken,pTracks);

  //edm::Handle<std::vector<double> > _rhoHandle;
  //iEvent.getByToken(rhoToken, _rhoHandle);

  std::unique_ptr< reco::GsfElectronCollection> IdentifiedElectrons (new reco::GsfElectronCollection());
  // std::vector<reco::GsfElectron> * IdentifiedElectrons=0; 
  std::cout << __LINE__ << " valid collection "<< pTracks.isValid() << std::endl;
  if ( bool(pTracks.isValid())==true ) {
    std::cout << __LINE__ << std::endl;
    const reco::GsfElectronCollection* eTracks = pTracks.product();
    std::cout << __LINE__ << std::endl;
   
    reco::GsfElectronCollection::const_iterator electrons;
    std::cout << __LINE__ << std::endl;

    for ( electrons = eTracks->begin(); electrons != eTracks->end(); ++electrons ) {
      float pt_e = electrons->pt();
      float eta_e = electrons->eta(); 
      std::cout << __LINE__ << std::endl;
    
      if (pt_e < ptthr ) continue; 
      std::cout << __LINE__ << std::endl;
      if (fabs(eta_e)>etathr) continue; 
      auto etrack  = electrons->gsfTrack().get(); 
      std::cout << __LINE__ << std::endl;

      if (electrons->isEB()){
	if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>2) continue; 
	if (electrons->full5x5_sigmaIetaIeta() >  0.0115) continue; 
	if (fabs(electrons->deltaPhiSuperClusterTrackAtVtx())> 0.228) continue; 
	if (fabs(electrons->deltaEtaSuperClusterTrackAtVtx())> 0.00749) continue; 
	if (electrons->hadronicOverEm()>0.346) continue;
	float abseta = fabs(electrons->superCluster()->position().eta());
	if (abseta> 1.479) continue; // check if it is really needed
	const float  eA = _effectiveAreas.getEffectiveArea( abseta );
	const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle) : 0; 
	std::cout << __LINE__ << " rhoHandle.isValid() " << _rhoHandle.isValid()  << std::endl; 
	std::cout << __LINE__ << " ept:" <<  pt_e << ", e eta: "<< eta_e << " " << eA<< " " << rho << std::endl; 
	if (( electrons->pfIsolationVariables().sumChargedHadronPt + 
	      std::max(float(0.0), electrons->pfIsolationVariables().sumNeutralHadronEt +  
		       electrons->pfIsolationVariables().sumPhotonEt -  eA*rho )
	      ) > 0.175*pt_e) continue; 
	nSelectedElectronsEB++; 
	IdentifiedElectrons->push_back(*electrons);

	
      }else if (electrons->isEE()){
	if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>3) continue; 
	if (electrons->full5x5_sigmaIetaIeta() > 0.037) continue; 
	if (fabs(electrons->deltaPhiSuperClusterTrackAtVtx())> 0.213) continue; 
	if (fabs(electrons->deltaEtaSuperClusterTrackAtVtx())> 0.00895) continue; 
	if (electrons->hadronicOverEm()>0.211) continue;
	float abseta = fabs(electrons->superCluster()->position().eta());
	if (abseta< 1.479) continue; // check if it is really needed
	if (abseta>=2.5  ) continue; // check if it is really needed

	const float  eA = _effectiveAreas.getEffectiveArea( abseta );
	const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle) : 0; 
	std::cout << __LINE__ << " rhoHandle.isValid() " << _rhoHandle.isValid()  << std::endl; 
	std::cout << __LINE__ << " ept:" <<  pt_e << ", e eta: "<< eta_e << " " << eA<< " " << rho << std::endl; 

	if ((electrons->pfIsolationVariables().sumChargedHadronPt + 
	     std::max(float(0.0), electrons->pfIsolationVariables().sumNeutralHadronEt + 
		 electrons->pfIsolationVariables().sumPhotonEt - eA*rho)
	     ) >0.159*pt_e) continue; 
	nSelectedElectronsEE++; 
	IdentifiedElectrons->push_back(*electrons);
      }//isEE
    }//electron loop
  }// valid collection

 

  if (IdentifiedElectrons->size() <2) 
    return false;
  std::cout << __LINE__<< std::endl; 

  for (reco::GsfElectronCollection::const_iterator ele1 = IdentifiedElectrons->begin() ; ele1!=IdentifiedElectrons->end(); ele1++){
    for (reco::GsfElectronCollection::const_iterator ele2 = ele1+1  ; ele2!=IdentifiedElectrons->end(); ele2++){
      if (ele1->pt()>ele2->pt()) {
	if (ele1->pt()< pt_e1_thr) continue;
	if (ele2->pt()< pt_e2_thr) continue;
      }else {
	if (ele1->pt()< pt_e2_thr) continue;
	if (ele2->pt()< pt_e1_thr) continue;
      }
      reco::Candidate::LorentzVector zp4= ele1->p4()+ele2->p4();
      double mass = zp4.M();
      std::cout <<__LINE__<< " Mass " << mass << std::endl;
      if (mass < massMin) continue;
      if (mass > massMax) continue; 
      keepEvent = true; 
    }
   }
  
  if (keepEvent==true){
    //nSelectedEvents++;
    iEvent.put(std::move(IdentifiedElectrons));
  }
  
  
  
    std::cout << __LINE__<< " filter result: "<< keepEvent  << " nSelectedElectrons "<< nSelectedElectronsEB << " + " << nSelectedElectronsEE << std::endl;
  return keepEvent;
}

//define this as a plug-in


DEFINE_FWK_MODULE( ZElectronsSkim );

