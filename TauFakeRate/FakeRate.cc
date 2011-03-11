// -*- C++ -*-
//
// Package:    FakeRate
// Class:      FakeRate
// 
/**\class FakeRate FakeRate.cc TauAnalysis/FakeRate/src/FakeRate.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Silvia Taroni,21 1-007,+41227676459,
//         Created:  Mon Feb 14 11:38:53 CET 2011
// $Id$
//
//



// user include files
#include "TauAnalysis/FakeRate/interface/FakeRate.h"
// system include files
#include <memory>
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
// constructors and destructor
//
FakeRate::FakeRate(const edm::ParameterSet& iConfig):
  discriminatorAgainstMuon_ ( iConfig.getParameter<edm::InputTag>("TauAgainstMuon")),
  discriminatorAgainstElectron_ ( iConfig.getParameter<edm::InputTag>("TauAgainstElectron")),
  discriminatorLooseIso_ ( iConfig.getParameter<edm::InputTag>("TauLooseIso")),
  discriminatorMediumIso_ ( iConfig.getParameter<edm::InputTag>("TauMediumIso")),
  discriminatorTightIso_ ( iConfig.getParameter<edm::InputTag>("TauTightIso")),
  discriminatorDecay_ ( iConfig.getParameter<edm::InputTag>("TauDecay")),
  tauSrc_ ( iConfig.getParameter<edm::InputTag>("TauSource"))
{

  //now do what ever initialization is needed

}


FakeRate::~FakeRate()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
namespace {
  reco::PFJetRef getJetRef(const reco::PFTau& tau) {
    if (tau.jetRef().isNonnull())
      return tau.jetRef();
    else if (tau.pfTauTagInfoRef()->pfjetRef().isNonnull())
      return tau.pfTauTagInfoRef()->pfjetRef();
    else throw cms::Exception("cant find jet ref");
  }
}

// ------------ method called to for each event  ------------
void
FakeRate::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  std::vector<double>  tauPt, tauEta, bTag, mcTruth, hT;
  std::vector<double>  looseIso, mediumIso, tightIso, againstEle, againstMu, byDecay;
  std::vector<int> selJetMulteplicity,jetMulteplicityNoTauSel,jetMulteplicity;
  std::vector<int> physId, st2Id, st3Id;

  edm::Handle<reco::TrackCollection> recoTracks;
  iEvent.getByLabel("generalTracks", recoTracks ); 
  edm::Handle<reco::PFJetCollection> pfjets;
//    iEvent.getByLabel("ak5PFJets", pfjets);
  iEvent.getByLabel("ak5PFJets", pfjets);
  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByLabel(tauSrc_, taus);
  edm::Handle<reco::PFTauDiscriminator> TauAgainstMuon;
  iEvent.getByLabel(discriminatorAgainstMuon_, TauAgainstMuon);
  edm::Handle<reco::PFTauDiscriminator> TauAgainstElectron;
  iEvent.getByLabel(discriminatorAgainstElectron_, TauAgainstElectron);
  edm::Handle<reco::PFTauDiscriminator> TauLooseIso;
  iEvent.getByLabel(discriminatorLooseIso_, TauLooseIso);
  edm::Handle<reco::PFTauDiscriminator> TauMediumIso;
  iEvent.getByLabel(discriminatorMediumIso_, TauMediumIso);
  edm::Handle<reco::PFTauDiscriminator> TauTightIso;
  iEvent.getByLabel(discriminatorTightIso_, TauTightIso);
  edm::Handle<reco::PFTauDiscriminator> TauDecay;
  iEvent.getByLabel(discriminatorDecay_, TauDecay);
  edm::Handle<reco::JetIDValueMap> hJetIDMap;
  iEvent.getByLabel( "ak5JetID", hJetIDMap );
  edm::Handle<reco::JetMatchedPartonsCollection> theTagByRef;
  iEvent.getByLabel ("flavourByRef", theTagByRef   );
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel ("ak5GenJets", genJets  );

  int theTauJetIndex = -1;
//   const reco::PFTau* theTauJet = getTheTauJet(*taus, tauJetEtaCut_, tauJetPtCut_, theTauJetIndex);

  edm::Handle<reco::PFJetCollection>  hJets; // uncorrected jets!
  iEvent.getByLabel("ak5PFJets", hJets );
  const JetCorrector* corrector = JetCorrector::getJetCorrector ("ak5PFL2L3", iSetup);   
  // jet ID selection
  PFJetIDSelectionFunctor jetIDFunctor( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE ); // this line can be done once, in a constructor or such
  pat::strbitset ret = jetIDFunctor.getBitTemplate(); // Only needed if you plan to use the detailed return values (see Selector docs)
  

  int jetIdPartPhy =0, jetIdPartSt2=0,jetIdPartSt3=0;
  int tauCounter = 0;
  h_nTracks->Fill(recoTracks->size()); 
  for ( unsigned iTau = 0; iTau < taus->size(); iTau++ ) {
     cout << __LINE__ << " PFJet size "<< pfjets->size() << ", tau size " << taus->size() << endl;
     jetIdPartPhy =0;
    jetIdPartSt2 =0;
    jetIdPartSt3 =0;
    reco::PFTauRef tauCandidate(taus, iTau);
    reco::PFJetRef jet = getJetRef(*tauCandidate);
    const Jet * myTauJet = (const Jet *)(&*jet);
    double dRmin=100.;
    const reco::PFJet * myJet=0;
    for (reco::PFJetCollection::const_iterator jpf = pfjets->begin(); jpf != pfjets->end(); jpf++){
      double deltaEta = tauCandidate->eta() - jpf->eta();
      double deltaPhi = tauCandidate->phi() - jpf->phi();
       double dR = sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);

     if (dR<dRmin){
	dRmin=dR;
	myJet=&*jpf;
      }
//        cout << deltaEta << " " << deltaPhi<< " " << dR<< " " << dRmin << endl;
    }
//     if (myGenJet !=0)cout  << myGenJet->print () << endl;;
//     const MatchedPartons * aMatch =0;
    const reco::Jet * myMatchedJet=0;
    double dRmatchmin = 1.;
    const MatchedPartons * aMatch =0;
    for ( JetMatchedPartonsCollection::const_iterator j  = theTagByRef->begin();
	  j != theTagByRef->end();
	  j ++ ) {
      const Jet *aJet       = (*j).first.get();
      double deltaEta = tauCandidate->eta() - aJet->eta();
      double deltaPhi = tauCandidate->phi() - aJet->phi();
      double dR = sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);

      if (dR<dRmatchmin){
	dRmatchmin=dR;
	aMatch = & j->second;
      }
    }
      // //       cout << "Line " << __LINE__ <<endl;
//       if (aJet != myTauJet) continue;
      //       cout << "Line " << __LINE__ <<endl;
    if (aMatch==0) continue;
//       const MatchedPartons * aMatch = & j->second;
    const GenParticleRef aPartPhy = aMatch->physicsDefinitionParton();
    const GenParticleRef aPartSt2 = aMatch->nearest_status2();
    const GenParticleRef aPartSt3 = aMatch->nearest_status3();
    
    if( aPartPhy.isNonnull()){
      jetIdPartPhy=aPartPhy.get()->pdgId();
      // 	cout << "Line " << __LINE__ << "phys ID "<< jetIdPartPhy <<  endl;
    }
    if(aPartSt2.isNonnull()){
      jetIdPartSt2=aPartSt2.get()->pdgId();
      // 	cout << "Line " << __LINE__ << ", jetIdPartSt2 "<<  jetIdPartSt2 << endl;
    }
    if(aPartSt3.isNonnull()) {
      jetIdPartSt3= aPartSt3.get()->pdgId();
      // 	cout << "Line " << __LINE__ << ", jetIdPartSt3 "<<  jetIdPartSt3 << endl
      ;
    }
//   }
//     cout << "Line " << __LINE__ << " phys ID "<< jetIdPartPhy << ", jetIdPartSt2 "<<  jetIdPartSt2 << ", jetIdPartSt3 "<<  jetIdPartSt3 << endl;
    
    bool partPhys = false;
    bool partSt2  = false;
    bool partSt3  = false; // Status 3 has to be used
    // from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate 
    // 3	 identifes the "hard part" of the interaction, i.e. the partons that are used in the matrix element calculation, including immediate decays of resonances. (documentation entry, defined separately from the event history. "This includes the two incoming colliding particles and partons produced in hard interaction." [ * ])


//     if (fabs(jetIdPartPhy)==1 || fabs(jetIdPartPhy)==2 ||
// 	fabs(jetIdPartPhy)==3 || fabs(jetIdPartPhy)==4 ||
// 	fabs(jetIdPartPhy)==5 || fabs(jetIdPartPhy)==6 ||
// 	fabs(jetIdPartPhy)==21 ) partPhys = true;

//     if (fabs(jetIdPartSt2)==1 || fabs(jetIdPartSt2)==2 ||
// 	fabs(jetIdPartSt2)==3 || fabs(jetIdPartSt2)==4 ||
// 	fabs(jetIdPartSt2)==5 || fabs(jetIdPartSt2)==6 ||
// 	fabs(jetIdPartSt2)==21) partSt2 = true;

//     if (fabs(jetIdPartSt3)==1 || fabs(jetIdPartSt3)==2 ||
// 	fabs(jetIdPartSt3)==3 || fabs(jetIdPartSt3)==4 ||
// 	fabs(jetIdPartSt3)==5 || fabs(jetIdPartSt3)==6 ||
// 	fabs(jetIdPartSt3)==21) partSt3 = true;

//     if (partPhys == false && partSt2 == false && partSt3==false ) continue;
//     if (partSt3 == false) continue;
    physId.push_back(jetIdPartPhy);
    st2Id.push_back(jetIdPartSt2);
    st3Id.push_back(jetIdPartSt3);	
    cout << "Line " << __LINE__ << " phys ID "<< jetIdPartPhy << ", jetIdPartSt2 "<<  jetIdPartSt2 << ", jetIdPartSt3 "<<  jetIdPartSt3 << endl;
    tauCounter++,
    h_ptvseta->Fill(tauCandidate->eta(), tauCandidate->pt());
    h_discriminator->Fill((*TauLooseIso)[tauCandidate]);
    h_byDecay->Fill((*TauDecay)[tauCandidate]);
    h_AgainstElectron->Fill((*TauAgainstElectron)[tauCandidate]);
    h_AgainstMuon->Fill((*TauAgainstMuon)[tauCandidate]);
    h_TightIsoSel->Fill((*TauTightIso)[tauCandidate]);
    h_MediumIsoSel->Fill((*TauMediumIso)[tauCandidate]);
    h_LooseIsoSel->Fill((*TauLooseIso)[tauCandidate]);
			     
    if( (*TauLooseIso)[tauCandidate] > 0.5 ){
      h_ptvseta_selTau->Fill(tauCandidate->eta(), tauCandidate->pt());
      h_discriminator_selTau->Fill((*TauLooseIso)[tauCandidate]);

      h_byDecay_LooseIsoSel->Fill((*TauDecay)[tauCandidate]);
      h_AgainstMuon_LooseIsoSel->Fill((*TauAgainstMuon)[tauCandidate]);
      h_AgainstElectron_LooseIsoSel->Fill((*TauAgainstElectron)[tauCandidate]);
    }//discriminator
    if( (*TauMediumIso)[tauCandidate] > 0.5 ){
      h_byDecay_MediumIsoSel->Fill((*TauDecay)[tauCandidate]);
      h_AgainstMuon_MediumIsoSel->Fill((*TauAgainstMuon)[tauCandidate]);
      h_AgainstElectron_MediumIsoSel->Fill((*TauAgainstElectron)[tauCandidate]);
    }//discriminator
    if( (*TauTightIso)[tauCandidate] > 0.5 ){
      h_byDecay_TightIsoSel->Fill((*TauDecay)[tauCandidate]);
      h_AgainstMuon_TightIsoSel->Fill((*TauAgainstMuon)[tauCandidate]);
      h_AgainstElectron_TightIsoSel->Fill((*TauAgainstElectron)[tauCandidate]);
    }//discriminator

    // cout << __LINE__ << endl;
    if ( (*TauAgainstMuon)[tauCandidate] > 0.5 ){
      h_byDecay_AgainstMuon->Fill((*TauDecay)[tauCandidate]);
      h_AgainstElectron_AgainstMuon->Fill((*TauAgainstElectron)[tauCandidate]);
      h_TightIsoSel_AgainstMuon->Fill((*TauTightIso)[tauCandidate]);
      h_MediumIsoSel_AgainstMuon->Fill((*TauMediumIso)[tauCandidate]);
      h_LooseIsoSel_AgainstMuon->Fill((*TauLooseIso)[tauCandidate]);
    }
    if ( (*TauAgainstElectron)[tauCandidate] > 0.5 ){
      h_byDecay_AgainstElectron->Fill((*TauDecay)[tauCandidate]);
      h_AgainstMuon_AgainstElectron->Fill((*TauAgainstMuon)[tauCandidate]);
      h_TightIsoSel_AgainstElectron->Fill((*TauTightIso)[tauCandidate]);
      h_MediumIsoSel_AgainstElectron->Fill((*TauMediumIso)[tauCandidate]);
      h_LooseIsoSel_AgainstElectron->Fill((*TauLooseIso)[tauCandidate]);
    }

    if ( (*TauDecay)[tauCandidate] > 0.5 ){
      h_AgainstElectron_byDecay->Fill((*TauAgainstElectron)[tauCandidate]);
      h_AgainstMuon_byDecay->Fill((*TauAgainstMuon)[tauCandidate]);
      h_TightIsoSel_byDecay->Fill((*TauTightIso)[tauCandidate]);
      h_MediumIsoSel_byDecay->Fill((*TauMediumIso)[tauCandidate]);
      h_LooseIsoSel_byDecay->Fill((*TauLooseIso)[tauCandidate]);
    }
    
    double htVal=0;
//     pat::strbitset ret = jetIDLoose.getBitTemplate();
//     unsigned int idx;
//     for (edm::View<reco::PFJet>::const_iterator jet = pfjets->begin(); jet!=pfjets->end(); jet++){
    int jetCountNoTauSel =0;
    int jetCount = 0 ;
    for (reco::PFJetCollection::const_iterator jj = pfjets->begin(); jj!=pfjets->end(); jj++){
      const reco::PFJet * jjet= &*jj;
      ret.set(false);
      bool passed = jetIDFunctor( *jjet,  ret );
      if (!passed) continue;
      jetCountNoTauSel++;
      if ( jjet == myJet) continue;
      reco::PFJet correctedJet = *jjet;
      if ( corrector ) {
	double scale = corrector->correction(jjet->p4());
	correctedJet.scaleEnergy(scale);
      }
      htVal+=correctedJet.pt();
      jetCount++;

    // cout << __LINE__ << endl;
    }
 //    cout << __LINE__ << " taujetEta " << jet.get()->eta() << " taujetpt " << jet.get()->pt()  << endl;
    //     // cout <<"================" << endl;    
    // cout << __LINE__ << endl;
    tauPt.push_back(tauCandidate->pt());
    tauEta.push_back(tauCandidate->eta());
//     mcTruth.push_back();
    hT.push_back(htVal);
    // cout << __LINE__ << endl;
    looseIso.push_back((*TauLooseIso)[tauCandidate]);
    mediumIso.push_back((*TauMediumIso)[tauCandidate]);
    tightIso.push_back((*TauTightIso)[tauCandidate]);
    againstEle.push_back((*TauAgainstElectron)[tauCandidate]);
    // cout << __LINE__ << endl;
    againstMu.push_back((*TauAgainstMuon)[tauCandidate]);
    byDecay.push_back((*TauDecay)[tauCandidate]);
    selJetMulteplicity.push_back(jetCount);
    jetMulteplicityNoTauSel.push_back(jetCountNoTauSel);
    jetMulteplicity.push_back(pfjets->size());
    cout << __LINE__ << " PFJet size "<< pfjets->size() << ", jetNoTauSel " << jetCountNoTauSel << ", selected jets " << jetCount <<", tau size " << taus->size() << ", selected taus "<< tauCounter << endl;
    cout << " ------------------- " << endl ;
  }//iTau

  // cout << __LINE__ << endl;
  branch.physId = physId;
  branch.st2Id = st2Id;
  branch.st3Id = st3Id;
  branch.tauPt  = tauPt;
  branch.tauEta = tauEta;
    // cout << __LINE__ << endl;
//   branch.mcTruth= mcTruth;
//   branch.bTag = bTag;
  branch.hT = hT;
    // cout << __LINE__ << endl;
  branch.looseIso = looseIso;
  branch.mediumIso= mediumIso;
    // cout << __LINE__ << endl;
  branch.tightIso = tightIso;
  branch.againstEle = againstEle;
    // cout << __LINE__ << endl;
  branch.againstMu = againstMu;
  branch.byDecay = byDecay;
  branch.selJetMulteplicity= selJetMulteplicity;
  branch.jetMulteplicityNoTauSel= jetMulteplicityNoTauSel;
  branch.jetMulteplicity=jetMulteplicity;
  ntuple->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
FakeRate::beginJob()
{
  edm::Service<TFileService> fs;
  TFileDirectory histoDir = fs->mkdir( "histos" );
  TFileDirectory treeDir = fs->mkdir( "trees" );

  h_nTracks = histoDir.make<TH1F>( "h_nTracks", "nTracks", 1000, 0, 1000 );

  h_discriminator = histoDir.make<TH1D>( "h_LooseIsoSel", "LooseIsoSel",100, -2,2);
  h_discriminator_selTau = histoDir.make<TH1D>( "h_LooseIsoSel_selTau", "LooseIsoSel_selTau",100, -2,2);
  h_byDecay = histoDir.make<TH1D>( "h_byDecay", "byDecay",100, -2,2);

  h_byDecay_LooseIsoSel = histoDir.make<TH1D>( "h_byDecay_LooseIsoSel", "byDecay_LooseIsoSel",100, -2,2);
  h_AgainstMuon_LooseIsoSel= histoDir.make<TH1D>( "h_AgainstMuon_LooseIsoSel","AgainstMuon if LooseIsoSel",100,-2,2);
  h_AgainstElectron_LooseIsoSel= histoDir.make<TH1D>( "h_AgainstElectron_LooseIsoSel","AgainstElectron if LooseIsoSel",100,-2,2);

  h_byDecay_MediumIsoSel = histoDir.make<TH1D>( "h_byDecay_MediumIsoSel", "byDecay_MediumIsoSel",100, -2,2);
  h_AgainstMuon_MediumIsoSel= histoDir.make<TH1D>( "h_AgainstMuon_MediumIsoSel","AgainstMuon if MediumIsoSel",100,-2,2);
  h_AgainstElectron_MediumIsoSel= histoDir.make<TH1D>( "h_AgainstElectron_MediumIsoSel","AgainstElectron if MediumIsoSel",100,-2,2);

  h_byDecay_TightIsoSel = histoDir.make<TH1D>( "h_byDecay_TightIsoSel", "byDecay_TightIsoSel",100, -2,2);
  h_AgainstMuon_TightIsoSel= histoDir.make<TH1D>( "h_AgainstMuon_TightIsoSel","AgainstMuon if TightIsoSel",100,-2,2);
  h_AgainstElectron_TightIsoSel= histoDir.make<TH1D>( "h_AgainstElectron_TightIsoSel","AgainstElectron if TightIsoSel",100,-2,2);


  h_AgainstElectron_byDecay = histoDir.make<TH1D>("h_AgainstElectron_byDecay", "AgainstElectron if byDecay", 100, -2,2);
  h_AgainstMuon_byDecay     = histoDir.make<TH1D>("h_AgainstMuon_byDecay"    , "AgainstMuon if byDecay", 100, -2,2);
  h_TightIsoSel_byDecay     = histoDir.make<TH1D>("h_TightIsoSel_byDecay"    , "TightIsoSel if byDecay", 100, -2,2);
  h_MediumIsoSel_byDecay    = histoDir.make<TH1D>("h_MediumIsoSel_byDecay"   , "MediumIsoSel if byDecay", 100, -2,2);
  h_LooseIsoSel_byDecay     = histoDir.make<TH1D>("h_LooseIsoSel_byDecay"    , "LooseIsoSel if byDecay", 100, -2,2);

  h_byDecay_AgainstElectron = histoDir.make<TH1D>("h_byDecay_AgainstElectron", "byDecay if AgainstElectron", 100, -2,2);
  h_AgainstMuon_AgainstElectron = histoDir.make<TH1D>("h_AgainstMuon_AgainstElectron"    , "AgainstMuon if AgainstElectron", 100, -2,2);
  h_TightIsoSel_AgainstElectron = histoDir.make<TH1D>("h_TightIsoSel_AgainstElectron"    , "TightIsoSel if AgainstElectron", 100, -2,2);
  h_MediumIsoSel_AgainstElectron= histoDir.make<TH1D>("h_MediumIsoSel_AgainstElectro"   , "MediumIsoSel if AgainstElectron", 100, -2,2);
  h_LooseIsoSel_AgainstElectron = histoDir.make<TH1D>("h_LooseIsoSel_AgainstElectron"    , "LooseIsoSel if AgainstElectron", 100, -2,2);

  h_byDecay_AgainstMuon = histoDir.make<TH1D>("h_byDecay_AgainstMuon", "byDecay if AgainstMuon", 100, -2,2);
  h_AgainstElectron_AgainstMuon = histoDir.make<TH1D>("h_AgainstMuon_AgainstMuon"    , "AgainstMuon if AgainstMuon", 100, -2,2);
  h_TightIsoSel_AgainstMuon = histoDir.make<TH1D>("h_TightIsoSel_AgainstMuon"    , "TightIsoSel if AgainstMuon", 100, -2,2);
  h_MediumIsoSel_AgainstMuon= histoDir.make<TH1D>("h_MediumIsoSel_AgainstElectro"   , "MediumIsoSel if AgainstMuon", 100, -2,2);
  h_LooseIsoSel_AgainstMuon = histoDir.make<TH1D>("h_LooseIsoSel_AgainstMuon"    , "LooseIsoSel if AgainstMuon", 100, -2,2);

  h_byDecay_AgainstMuon = histoDir.make<TH1D>("h_byDecay_AgainstMuon", "byDecay if AgainstMuon", 100, -2,2);
  h_AgainstElectron_AgainstMuon = histoDir.make<TH1D>("h_AgainstMuon_AgainstMuon"    , "AgainstMuon if AgainstMuon", 100, -2,2);
  h_TightIsoSel_AgainstMuon = histoDir.make<TH1D>("h_TightIsoSel_AgainstMuon"    , "TightIsoSel if AgainstMuon", 100, -2,2);
  h_MediumIsoSel_AgainstMuon= histoDir.make<TH1D>("h_MediumIsoSel_AgainstElectro"   , "MediumIsoSel if AgainstMuon", 100, -2,2);
  h_LooseIsoSel_AgainstMuon = histoDir.make<TH1D>("h_LooseIsoSel_AgainstMuon"    , "LooseIsoSel if AgainstMuon", 100, -2,2);

  h_AgainstElectron= histoDir.make<TH1D>("h_AgainstElectron","againstElectron",100,-2,2);
  h_AgainstMuon= histoDir.make<TH1D>("h_AgainstMuon","againstMuon",100,-2,2);
  h_TightIsoSel= histoDir.make<TH1D>("h_TightIsoSel","tightIsoSel",100,-2,2);
  h_MediumIsoSel= histoDir.make<TH1D>("h_MediumIsoSel","mediumIsoSel",100,-2,2);
  h_LooseIsoSel= histoDir.make<TH1D>("h_LooseIsoSel","looseIsoSel",100,-2,2);

  h_ptvseta = histoDir.make<TH2D>( "h_ptvseta", "ptvseta", 250,-2.5,2.5,100,0,200);
  h_ptvseta_selTau = histoDir.make<TH2D>( "h_ptvseta_selTau", "ptvseta_selTau", 250,-2.5,2.5,100,0,200);



  ntuple = treeDir.make<TTree>("ntuple","tauVariables" ) ;
  ntuple ->Branch("partonPhys", &(branch.physId));
  ntuple ->Branch("partonSt2" , &(branch.st2Id));
  ntuple ->Branch("partonSt3" , &(branch.st3Id));
  ntuple ->Branch("tauPt"    , &(branch.tauPt));
  ntuple ->Branch("tauEta"   , &(branch.tauEta));
  ntuple ->Branch("mcTruth"  , &(branch.mcTruth));
  ntuple ->Branch("bTag"     , &(branch.bTag));
  ntuple ->Branch("hT"       , &(branch.hT));
  ntuple ->Branch("looseIso" , &(branch.looseIso));
  ntuple ->Branch("mediumIso", &(branch.mediumIso));
  ntuple ->Branch("tightIso" , &(branch.tightIso));
  ntuple ->Branch("againstEle", &(branch.againstEle));
  ntuple ->Branch("againstMu" , &(branch.againstMu));
  ntuple ->Branch("byDecay"   , &(branch.byDecay));
  ntuple ->Branch("selJetMulteplicity"   , &(branch.selJetMulteplicity));
  ntuple ->Branch("jetMulteplicityNoTauSel" , &(branch.jetMulteplicityNoTauSel));
  ntuple ->Branch("jetMulteplicity"      , &(branch.jetMulteplicity));


}

// ------------ method called once each job just after ending the event loop  ------------
void 
FakeRate::endJob() {

}


//define this as a plug-in
DEFINE_FWK_MODULE(FakeRate);
