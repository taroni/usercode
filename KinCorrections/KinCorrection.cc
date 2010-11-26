#include "KinCorrection.hh"
#include "CommonOps.hh"
#include "EventData.hh"
#include "KinSuite.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "Types.hh"
#include "GenObject.hh"
#include "TRandom3.h"
#include "TString.h"
#include "Math/GenVector/DisplacementVector3D.h"

#include "ChAsymTemplateHistos.hh"


using namespace std;
using namespace Event;
using namespace ROOT::Math;


// -----------------------------------------------------------------------------
//
KinCorrection::KinCorrection(const std::string & filename, const Utils::ParameterSet& pset ) :
  // Misc
//   dirName_(   filename)
// {}
dirName_(filename)
{
	zMassMin_ = pset.Get<double>("zMassMin");
	zMassMax_ = pset.Get<double>("zMassMax");
	BarrelEtaMax_ = pset.Get<double>("BarrelEtaMax");
	EndCapEtaMin_ = pset.Get<double>("EndCapEtaMin");
	EndCapEtaMax_ = pset.Get<double>("EndCapEtaMax");
	ElecPtMin_  = pset.Get<double>("MinElecRescaledPt");
 	ReconPtMin_ = pset.Get<double>("c_ErNu_pt");
        wID_ = pset.Get<double>("wID");
        mW_ = pset.Get<double>("mW");
        mZ_ = pset.Get<double>("mZ");
//////////////////////////////////////////////////////////////////////////////
	chChk      = pset.Get<bool>("ChCheck");
	corrEE     = pset.Get<bool>("CorrEEMisalig");
	lepPtVeto  = pset.Get<double>("ElePtVeto");
	convChk    = pset.Get<bool>("ConvCheck");
	lepVeto    = pset.Get<bool>("EleVeto");	
	wp         = pset.Get<int>("WorkingPoint");
	dataset    = pset.Get<bool>("DataSet");//True= W, False= Z
	datatype   = pset.Get<bool>("DataType");//  True= Data, False= Monte Carlo
	ReWeightHistos_ 	= pset.Get<bool>("ReWeightHistos");	
	ReWeightCorrectionFile_ = pset.Get<std::string>("ReWeightCorrectionFile");
    
	//(0=95%,1=90%,2=85%,3=80%,4=70%,5=60%)
	ieff=-1;
	if (wp==95) ieff=0;
	if (wp==90) ieff=1;
	if (wp==85) ieff=2;
	if (wp==80) ieff=3;
	if (wp==70) ieff=4;
	if (wp==60) ieff=5;
	if (ieff<0){
	  
	  cout<<"Working point "<<wp<<" unknown"<<endl
	      <<"Working points accepted are 60,70,80,85,90,95"<<endl
	      <<"Job is going to stop"<<endl;
	  assert(0);
	} 
        mMissingHits=1;
        mDCot=0.02;
        mDist=0.02;
	cout << " ========> Analyzing " ;
	if ( dataset == 1 ) {
	  cout << "W to e nu ";
	} else {
	  cout << "Z to e e ";
	}
	if (datatype ==1 ) {
	  cout << "from LHC collisions" << endl;
	} else {
	  cout << "from MC simulation" << endl;
	}
}

// -----------------------------------------------------------------------------
//
KinCorrection::~KinCorrection() {}

// -----------------------------------------------------------------------------
//
void KinCorrection::Start( Event::Data& ev ) {
  initDir( ev.OutputFile(), dirName_.c_str() );
  BookHistos();
  if (ReWeightHistos_){
    WeightsFile_ = TFile::Open(ReWeightCorrectionFile_.c_str(), "READ");
    if (!WeightsFile_){
      cout << " no weight file " << endl;
      return;
    }
    if (wID_=24)  PtYWeight_   = (TH3F*) WeightsFile_->Get("KinCorrWPlus");
    if (wID_=-24) PtYWeight_   = (TH3F*) WeightsFile_->Get("KinCorrWMinus");
    
    if (!PtYWeight_){
      cout << " no weight histo " << endl;
      return;
    }
  }
 
 
  ievt=0;
}

// -----------------------------------------------------------------------------
//
void KinCorrection::BookHistos() {

  TString ib[10]={"_0","_1","_2","_3","_4","_5","_6","_7","_8","_9"};

// 	h_tagMatch_dR_   = new TH1F ("tagMatch_dR", "#Delta R of Match of Tag to MC Electron;#Delta R;Number of Events", 100, 0., 0.2);
// 	h_probeMatch_dR_ = new TH1F("probeMatch_dR", "#Delta R of Match of Supercluster Probe to MC Electron;#Delta R;Number of Events", 100, 0., 0.2);

// 	h_mc_Et_ = new TH1F("mc_Et", "MC Electron #Et;MC Electron E_{T}^{MC};", 100, 0, 100);
// 	h_mc_Eta_ = new TH1F("mc_Eta", "MC Electron #eta;MC Electron #eta^{MC};", 100, -2.5, 2.5);
// 	h_mc_fabsEta_ = new TH1F("mc_fabsEta", "MC Electron #eta;MC Electron #eta^{MC};", 50, 0., 2.5);
// 	h_Et_ = new TH1F("Et", "Electron #Et; Electron E_{T}^{MC};", 100, 0, 100);
// 	h_Eta_ = new TH1F("Eta", "Electron #eta; Electron #eta^{MC};", 100, -2.5, 2.5);
// 	h_fabsEta_ = new TH1F("fabsEta", "Electron #eta; Electron #eta^{MC};", 50, 0., 2.5);


    hist_Total 				= new TH1I("Total Number of Events", "Total # of events",10,0,10);
    hist_Total_Up 			= new TH1I("Total Number of Up Type Events", "Total # of Up Type events",10,0,10);
    hist_Total_Down 			= new TH1I("Total Number of Down Type Events", "Total # of Down Type events",10,0,10);
    hist_Total_gq 			= new TH1I("Total Number of Gluon-Quark Type Events", "Total # of gluon-quark Type events",10,0,10);
    hist_Total_Other 			= new TH1I("Total Number Other Events", "Total # of Other events",10,0,10);
    hist_Mother_PdgId 			= new TH1I("Mothers Pdg Id", "Mothers Pdg ID",46,-23,23);

    hist_PTvsETA_MassCut 		= new TH2F("PTvsETA_MassCut","Log PTvs ETA_MassCut",100,0.,100.,100,-8,8);
    hist_PTvsY_MassCut 			= new TH2F("PTvsY_MassCut","Log PT vs absY_MassCut",50,-2.,4.,50,0.,5.);
    hist_YvsETA_MassCut	 		= new TH2F("YvsETA_MassCut","YvsETA_MassCut",100,-8.,8.,100,-8,8);
    hist_PTvsYvsEta 			= new TH3F("PTvsYvsEta","PTvsYvsEta;Log_{10}(Pt);|Y|;Electron |#eta|", 100,-1.,3.,100,0.,6.,100,0.,6.);

    hist_Beta 				= new TH1F("Beta","Beta",100,0.,1.);	
    hist_Mass 				= new TH1F("Mass","Mass;Mass GeV",100,0.,100.);
    hist_Eta 				= new TH1F("Eta","Eta;#eta",100,-8.,8.);
    hist_Pt 				= new TH1F("Pt","PT;Pt",100,0.,100.);
    hist_P 				= new TH1F("P","P;P",100,0.,100.);
    hist_Rapidity 			= new TH1F("Rapidity","Rapidity;y",100,-8.,8.);
    
    hist_Electron_Pt 			= new TH1F("Electron_Pt","E PT;Pt",100,0.,100.);
    hist_Neutrino_Pt 			= new TH1F("Neutrino_Pt","Nu PT;Pt",100,0.,100.);
    hist_Electron_Eta 			= new TH1F("Electron_Eta","E Eta;#eta",100,-8.,8.);
    hist_Neutrino_Eta 			= new TH1F("Neutrino_Eta","Nu Eta;#eta",100,-8.,8.);
    hist_Electron_Rapidity 		= new TH1F("Electron_Rapidity","E y;y",100,-8.,8.);
    hist_Neutrino_Rapidity 		= new TH1F("Neutrino_Rapidity","Nu y;y",100,-8.,8.);

    hist_Beta_MassCut 			= new TH1F("Beta_MassCut","Beta_MassCut",100,0.,1.);
    hist_Mass_MassCut 			= new TH1F("Mass_MassCut","Mass_MassCut;Mass GeV",100,0.,100.);
    hist_Eta_MassCut 			= new TH1F("Eta_MassCut","Eta_MassCut;#eta",100,-8.,8.);
    hist_Pt_MassCut 			= new TH1F("Pt_MassCut","PT_MassCut;Pt",100,0.,100.);
    hist_P_MassCut 			= new TH1F("P_MassCut","P_MassCut;P",100,0.,100.);
    hist_Rapidity_MassCut 		= new TH1F("Rapidity_MassCut","Rapidity_MassCut;y",100,-8.,8.);
    
    hist_Electron_Pt_MassCut 		= new TH1F("Electron_Pt_MassCut","E PT MassCut;Pt",100,0.,100.);
    hist_Neutrino_Pt_MassCut 		= new TH1F("Neutrino_Pt_MassCut","Nu PT MassCut;Pt",100,0.,100.);
    hist_Electron_Eta_MassCut 		= new TH1F("Electron_Eta_MassCut","E Eta MassCut;#eta",100,-8.,8.);
    hist_Neutrino_Eta_MassCut 		= new TH1F("Neutrino_Eta_MassCut","Nu Eta MassCut;#eta",100,-8.,8.);
    hist_Electron_Rapidity_MassCut 	= new TH1F("Electron_Rapidity_MassCut","E y MassCut;y",100,-8.,8.);
    hist_Neutrino_Rapidity_MassCut 	= new TH1F("Neutrino_Rapidity_MassCut","Nu y MassCut;y",100,-8.,8.);
    

    hist_Beta_ElectronID 		= new TH1F("Beta_ElectronID","Beta_ElectronID",100,0.,1.);
    hist_Mass_ElectronID 		= new TH1F("Mass_ElectronID","Mass_ElectronID;Mass GeV",100,0.,100.);
    hist_Eta_ElectronID 		= new TH1F("Eta_ElectronID","Eta_ElectronID;#eta",100,-8.,8.);
    hist_Pt_ElectronID 			= new TH1F("Pt_ElectronID","PT_ElectronID;Pt",100,0.,100.);
    hist_P_ElectronID 			= new TH1F("P_ElectronID","P_ElectronID;P",100,0.,100.);
    hist_Rapidity_ElectronID 		= new TH1F("Rapidity_ElectronID","Rapidity_ElectronID;y",100,-8.,8.);
    
    hist_Electron_Pt_ElectronID 	= new TH1F("Electron_Pt_ElectronID","E PT ElectronID;Pt",100,0.,100.);
    hist_Neutrino_Pt_ElectronID 	= new TH1F("Neutrino_Pt_ElectronID","Nu PT ElectronID;Pt",100,0.,100.);
    hist_Electron_Eta_ElectronID 	= new TH1F("Electron_Eta_ElectronID","E Eta ElectronID;#eta",100,-8.,8.);
    hist_Neutrino_Eta_ElectronID 	= new TH1F("Neutrino_Eta_ElectronID","Nu Eta ElectronID;#eta",100,-8.,8.);
    hist_Electron_Rapidity_ElectronID 	= new TH1F("Electron_Rapidity_ElectronID","E y ElectronID;y",100,-8.,8.);
    hist_Neutrino_Rapidity_ElectronID 	= new TH1F("Neutrino_Rapidity_ElectronID","Nu y ElectronID;y",100,-8.,8.);
    
    hist_Beta_Selected_NuInAcc 		= new TH1F("Beta_Selected_NuInAcc","Beta_Selected_NuInAcc",100,0.,1.);
    hist_Mass_Selected_NuInAcc 		= new TH1F("Mass_Selected_NuInAcc","Mass_Selected_NuInAcc;Mass GeV",100,0.,100.);
    hist_Eta_Selected_NuInAcc 		= new TH1F("Eta_Selected_NuInAcc","Eta_Selected_NuInAcc;#eta",100,-8.,8.);
    hist_Pt_Selected_NuInAcc 		= new TH1F("Pt_Selected_NuInAcc","PT_Selected_NuInAcc;Pt",100,0.,100.);
    hist_P_Selected_NuInAcc 		= new TH1F("P_Selected_NuInAcc","P_Selected_NuInAcc;P",100,0.,100.);
    hist_Rapidity_Selected_NuInAcc 	= new TH1F("Rapidity_Selected_NuInAcc","Rapidity_Selected_NuInAcc;y",100,-8.,8.);
    
    hist_Electron_Pt_Selected_NuInAcc 	= new TH1F("Electron_Pt_Selected_NuInAcc","E PT Selected_NuInAcc;Pt",100,0.,100.);
    hist_Neutrino_Pt_Selected_NuInAcc 	= new TH1F("Neutrino_Pt_Selected_NuInAcc","Nu PT Selected_NuInAcc;Pt",100,0.,100.);
    hist_Electron_Eta_Selected_NuInAcc 	= new TH1F("Electron_Eta_Selected_NuInAcc","E Eta Selected_NuInAcc;#eta",100,-8.,8.);
    hist_Neutrino_Eta_Selected_NuInAcc 	= new TH1F("Neutrino_Eta_Selected_NuInAcc","Nu Eta Selected_NuInAcc;#eta",100,-8.,8.);
    hist_Electron_Rapidity_Selected_NuInAcc = new TH1F("Electron_Rapidity_Selected_NuInAcc","E y Selected_NuInAcc;y",100,-8.,8.);
    hist_Neutrino_Rapidity_Selected_NuInAcc = new TH1F("Neutrino_Rapidity_Selected_NuInAcc","Nu y Selected_NuInAcc;y",100,-8.,8.);

    hist_Beta_Selected 			= new TH1F("Beta_Selected","Beta_Selected",100,0.,1.);
    hist_Mass_Selected 			= new TH1F("Mass_Selected","Mass_Selected;Mass GeV",100,0.,100.);
    hist_Eta_Selected 			= new TH1F("Eta_Selected","Eta_Selected;#eta",100,-8.,8.);
    hist_Pt_Selected 			= new TH1F("Pt_Selected","PT_Selected;Pt",100,0.,100.);
    hist_P_Selected 			= new TH1F("P_Selected","P_Selected;P",100,0.,100.);
    hist_Rapidity_Selected 		= new TH1F("Rapidity_Selected","Rapidity_Selected;y",100,-8.,8.);
    
    hist_Electron_Pt_Selected 		= new TH1F("Electron_Pt_Selected","E PT Selected;Pt",100,0.,100.);
    hist_Neutrino_Pt_Selected 		= new TH1F("Neutrino_Pt_Selected","Nu PT Selected;Pt",100,0.,100.);
    hist_Electron_Eta_Selected 		= new TH1F("Electron_Eta_Selected","E Eta Selected;#eta",100,-8.,8.);
    hist_Neutrino_Eta_Selected 		= new TH1F("Neutrino_Eta_Selected","Nu Eta Selected;#eta",100,-8.,8.);
    hist_Electron_Rapidity_Selected 	= new TH1F("Electron_Rapidity_Selected","E y Selected;y",100,-8.,8.);
    hist_Neutrino_Rapidity_Selected 	= new TH1F("Neutrino_Rapidity_Selected","Nu y Selected;y",100,-8.,8.);

// 	h_PfSel_mcEta_ = new TH1F("PfSel_mcEta", "Selected PF Probe: #eta^{MC};#eta^{MC};Number of Events", 100, -2.5, 2.5);
// 	h_PfSel_fabsMCEta_ = new TH1F("PfSel_fabsMCEta", "Selected PF Probe: #eta^{MC};#eta^{MC};Number of Events", 50, 0., 2.5);
// 	h_PfSel_fabsEta_ = new TH1F("PfSel_fabsEta", "Selected PF Probe: #eta^{};#eta^{};Number of Events", 50, 0., 2.5);
// 	h_PfSel_mcEt_ = new TH1F("PfSel_mcEt", "Selected PF Probe: p_{T}^{MC};p_{T}^{MC} (GeV);Number of Events", 100, 0., 100.);

// 	h_PfSel_Eta_ = new TH1F("PfSel_Eta", "Selected PF Probe: #eta;#eta;Number of Events", 100, -2.5, 2.5);
// 	h_PfSel_Et_ = new TH1F("PfSel_Et", "Selected PF Probe: p_{T};p_{T} (GeV);Number of Events", 100, 0., 100.);

}

// -----------------------------------------------------------------------------
// PROCESS
bool KinCorrection::Process( Event::Data& ev ) {
  //  float w=ev.GetEventWeight();
  std::vector<Event::GenObject> Mj_Electrons;  //define Electron Container
  std::vector<Event::GenObject> Mj_Neutrinos;  //define Netrino Container
  std::vector<Event::GenObject> Mj_W;          //define W Container
  std::vector<Event::GenObject> Mj_Z;          //define Z Container
  

  int ngoodEle=0;
  
  std::vector<Lepton const *>::const_iterator goodLep;
  std::vector< pair <Lepton const *, int > > myEle;
  std::vector< pair <Lepton const *, int > > myNu;
 

  double zMassMin= zMassMin_;
  double zMassMax= zMassMax_;

  MeeMin_ = zMassMin;
  
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > Daughter1, Daughter2;

  double pigreco = 3.14159265358979323846;
  ievt++;
  cout << "event #  "<< ievt<< endl;

  //----------
  //Fill Particle Vectors
  //----------
  for (std::vector<Event::GenObject>::const_iterator j = ev.GenParticles().begin();
       j != ev.GenParticles().end();
       ++j) {      //Loop all Gen particles
    if (j->GetStatus()==3){                           //All stable particles
       int ID = abs(j->GetID());
      //      int ID = abs(j->GetID());
      switch (ID){
      case 11: //Electron
	Mj_Electrons.push_back(*j);
	break;
      case 12: //Neutrino 
	Mj_Neutrinos.push_back(*j); 
	break;
      case 24: //W Boson
	Mj_W.push_back(*j);
	break;
      case 23:
	Mj_Z.push_back(*j);
      }//End Switch
      
    }//End If
  }//End Gen Loop
  
  // Selection Bools ----------
  bool ptselect 		= false;
  bool neuptselect 		= false;
  bool ElectronInaccPreErsatz 	= false;
  bool NeutrinoInaccPreErsatz 	= false;
  bool ElectronInacc 		= false;
  bool NeutrinoInacc 		= false;
  bool massselect 		= false;
  bool reconstructed 		= false;
  // --------------------------
  
  if( Mj_Z.size()==1&&
      Mj_Electrons.size()==2){
    
    hist_Total -> Fill(5);
    
    int elec_daughter=-5000, pos_daughter=-5000;
    for(unsigned int iEle = 0; iEle < Mj_Electrons.size();iEle++){
      if (Mj_Electrons[iEle].GetID() ==  11) elec_daughter = iEle;
      if (Mj_Electrons[iEle].GetID() == -11) pos_daughter = iEle;
    }
    
    TRandom3 *random = new TRandom3();
    if(random-> Uniform()<0.5){
      Daughter1.SetXYZT( Mj_Electrons[elec_daughter].Px(), Mj_Electrons[elec_daughter].Py(), Mj_Electrons[elec_daughter].Pz(), Mj_Electrons[elec_daughter].E());
      Daughter2.SetXYZT( Mj_Electrons[pos_daughter].Px(), Mj_Electrons[pos_daughter].Py(), Mj_Electrons[pos_daughter].Pz(), Mj_Electrons[pos_daughter].E());
    }
    else{
      Daughter1.SetXYZT( Mj_Electrons[pos_daughter].Px(), Mj_Electrons[pos_daughter].Py(), Mj_Electrons[pos_daughter].Pz(), Mj_Electrons[pos_daughter].E());
      Daughter2.SetXYZT( Mj_Electrons[elec_daughter].Px(), Mj_Electrons[elec_daughter].Py(), Mj_Electrons[elec_daughter].Pz(), Mj_Electrons[elec_daughter].E());
    }
    

    if (fabs(Daughter1.Eta()) < EndCapEtaMax_){
      if(fabs(Daughter1.Eta()) <= BarrelEtaMax_ || fabs(Daughter1.Eta()) >= EndCapEtaMin_)
	{ElectronInaccPreErsatz = true;}
    }
    
    if (fabs(Daughter2.Eta()) < EndCapEtaMax_){
      if(fabs(Daughter2.Eta())<= BarrelEtaMax_ || fabs(Daughter2.Eta()) >= EndCapEtaMin_)
	{NeutrinoInaccPreErsatz = true;}
    }
    
    if(Mj_Z[0].M() > MeeMin_) massselect = true;


  
    // Finally find out if it will be reconstructed
    if (Daughter1.Pt() > ReconPtMin_ && Daughter2.Pt() > ReconPtMin_)  reconstructed = true;
    
    // Now I'm doing the Ersatz
    RescaleZee(&Mj_Z[0],&Daughter1,&Daughter2);
    
    //find out if Pt cut satisfied (after rescale selection) (Daughter 1 is my Ersatz wElectron)
    if(Daughter1.Pt() > ElecPtMin_) ptselect = true;
    if(Daughter2.Pt() > ReconPtMin_) neuptselect = true;
    
    if (fabs(Daughter1.Eta()) < EndCapEtaMax_){
      if(fabs(Daughter1.Eta())<= BarrelEtaMax_ || fabs(Daughter1.Eta()) >= EndCapEtaMin_)
	{ElectronInacc = true;}
    }
    
    if (fabs(Daughter2.Eta()) < EndCapEtaMax_){
      if(fabs(Daughter2.Eta())<= BarrelEtaMax_ || fabs(Daughter1.Eta()) >= EndCapEtaMin_)
	{NeutrinoInacc = true;}
    }
  
    // And Fill some histograms
    
  // Plots before Selection
    hist_Beta->Fill(Mj_Z[0].Beta());
    hist_Mass->Fill(Mj_Z[0].mass());
    hist_Eta->Fill(Mj_Z[0].Eta());
    hist_Pt->Fill(Mj_Z[0].Pt());
    hist_P->Fill(Mj_Z[0].P());
    hist_Rapidity->Fill(Mj_Z[0].Rapidity());
    hist_Electron_Pt -> Fill(Daughter1.Pt());
    hist_Neutrino_Pt -> Fill(Daughter2.Pt());
    hist_Electron_Eta -> Fill(Daughter1.Eta());
    hist_Neutrino_Eta-> Fill(Daughter2.Eta());
    hist_Electron_Rapidity-> Fill(Daughter1.Rapidity());
    hist_Neutrino_Rapidity-> Fill(Daughter2.Rapidity());

    double weight;
    
    if (ReWeightHistos_){
      cout << "reweighting..." << endl;
      weight = GetWeight3D(Mj_Z[0],Daughter1.Eta());
       cout << "reweighted" << endl;
   }
    else  {
      weight = 1.0;
    }
    
    // Plots After MassCut
    
    if (massselect){
    
      //Fill 3D histo here
      
      hist_PTvsYvsEta->Fill(TMath::Log10(Mj_Z[0].Pt()),fabs(Mj_Z[0].Rapidity()),fabs(Daughter1.Eta()));
      hist_PTvsETA_MassCut->Fill(Mj_Z[0].Pt(),Mj_Z[0].Eta());
      hist_PTvsY_MassCut->Fill(TMath::Log10(Mj_Z[0].Pt()),fabs(Mj_Z[0].Rapidity()));
      hist_YvsETA_MassCut->Fill(Mj_Z[0].Rapidity(),Mj_Z[0].Eta());
      
      hist_Beta_MassCut->Fill(Mj_Z[0].Beta(),weight);
      hist_Mass_MassCut->Fill(Mj_Z[0].mass(),weight);
      hist_Eta_MassCut->Fill(Mj_Z[0].Eta(),weight);
      hist_Pt_MassCut->Fill(Mj_Z[0].Pt(),weight);
      hist_P_MassCut->Fill(Mj_Z[0].P(),weight);
      hist_Rapidity_MassCut->Fill(Mj_Z[0].Rapidity(),weight);
      hist_Electron_Pt_MassCut -> Fill(Daughter1.Pt(),weight);
      hist_Neutrino_Pt_MassCut -> Fill(Daughter2.Pt(),weight);
      hist_Electron_Eta_MassCut -> Fill(Daughter1.Eta(),weight);
      hist_Neutrino_Eta_MassCut-> Fill(Daughter2.Eta(),weight);
      hist_Electron_Rapidity_MassCut -> Fill(Daughter1.Rapidity(),weight);
      hist_Neutrino_Rapidity_MassCut-> Fill(Daughter2.Rapidity(),weight);
    }
  
  // Plots after W Selection (dont worry if electrons were reconstructed or not_
  
    if (ptselect && massselect && ElectronInacc){
      hist_Beta_ElectronID->Fill(Mj_Z[0].Beta(),weight);
      hist_Mass_ElectronID->Fill(Mj_Z[0].mass(),weight);
      hist_Eta_ElectronID->Fill(Mj_Z[0].Eta(),weight);
      hist_Pt_ElectronID->Fill(Mj_Z[0].Pt(),weight);
      hist_P_ElectronID->Fill(Mj_Z[0].P(),weight);
      hist_Rapidity_ElectronID->Fill(Mj_Z[0].Rapidity(),weight);
      hist_Electron_Pt_ElectronID -> Fill(Daughter1.Pt(),weight);
      hist_Neutrino_Pt_ElectronID -> Fill(Daughter2.Pt(),weight);
      hist_Electron_Eta_ElectronID -> Fill(Daughter1.Eta(),weight);
      hist_Neutrino_Eta_ElectronID-> Fill(Daughter2.Eta(),weight);
      hist_Electron_Rapidity_ElectronID -> Fill(Daughter1.Rapidity(),weight);
      hist_Neutrino_Rapidity_ElectronID-> Fill(Daughter2.Rapidity(),weight);
    
    }
    
    // Plots After Electron In Fiducial
    
    if (ptselect && massselect && NeutrinoInacc && ElectronInacc && neuptselect){
      hist_Beta_Selected_NuInAcc->Fill(Mj_Z[0].Beta(),weight);
      hist_Mass_Selected_NuInAcc->Fill(Mj_Z[0].mass(),weight);
      hist_Eta_Selected_NuInAcc->Fill(Mj_Z[0].Eta(),weight);
      hist_Pt_Selected_NuInAcc->Fill(Mj_Z[0].Pt(),weight);
      hist_P_Selected_NuInAcc->Fill(Mj_Z[0].P(),weight);
      hist_Rapidity_Selected_NuInAcc->Fill(Mj_Z[0].Rapidity(),weight);
      hist_Electron_Pt_Selected_NuInAcc -> Fill(Daughter1.Pt(),weight);
      hist_Neutrino_Pt_Selected_NuInAcc -> Fill(Daughter2.Pt(),weight);
      hist_Electron_Eta_Selected_NuInAcc -> Fill(Daughter1.Eta(),weight);
      hist_Neutrino_Eta_Selected_NuInAcc-> Fill(Daughter2.Eta(),weight);
      hist_Electron_Rapidity_Selected_NuInAcc -> Fill(Daughter1.Rapidity(),weight);
      hist_Neutrino_Rapidity_Selected_NuInAcc-> Fill(Daughter2.Rapidity(),weight);
      
    }

		  // Plots after
    if (ptselect && massselect && NeutrinoInacc && ElectronInacc &&
	NeutrinoInaccPreErsatz &&  ElectronInaccPreErsatz && 
	reconstructed && neuptselect){    
      hist_Beta_Selected->Fill(Mj_Z[0].Beta(),weight);
      hist_Mass_Selected->Fill(Mj_Z[0].mass(),weight);
      hist_Eta_Selected->Fill(Mj_Z[0].Eta(),weight);
      hist_Pt_Selected->Fill(Mj_Z[0].Pt(),weight);
      hist_P_Selected->Fill(Mj_Z[0].P(),weight);
      hist_Rapidity_Selected->Fill(Mj_Z[0].Rapidity(),weight);
      hist_Electron_Pt_Selected -> Fill(Daughter1.Pt(),weight);
      hist_Neutrino_Pt_Selected -> Fill(Daughter2.Pt(),weight);
      hist_Electron_Eta_Selected -> Fill(Daughter1.Eta(),weight);
      hist_Neutrino_Eta_Selected-> Fill(Daughter2.Eta(),weight);
      hist_Electron_Rapidity_Selected -> Fill(Daughter1.Rapidity(),weight);
      hist_Neutrino_Rapidity_Selected-> Fill(Daughter2.Rapidity(),weight);
    }
    
  }//Mj_Z==1

  if( Mj_W.size()==1&&
      Mj_Electrons.size()==1 && 
      Mj_Neutrinos.size()==1){

    hist_Total -> Fill(5);
    		  
    double eta = fabs(Mj_Electrons[0].Eta());
    
    
		  //std::cout << EleCosThCS << std::endl;  
    
    if(eta < 2.5 && (eta <=1.4442 || eta >=1.566)) ElectronInacc = true;
    if(Mj_Electrons[0].Pt() > ReconPtMin_) reconstructed = true;
    if(Mj_Electrons[0].Pt() > ElecPtMin_) ptselect = true;
    
    double nueta = Mj_Neutrinos[0].Eta();
    if(nueta < 2.5 && (nueta <=1.4442 || nueta >=1.566)) NeutrinoInacc = true; 
    if(Mj_Neutrinos[0].Pt() > ReconPtMin_) neuptselect = true;


    if(Mj_W[0].GetID()==wID_) {
    
      hist_Beta->Fill(Mj_W[0].Beta());
      hist_Mass->Fill(Mj_W[0].mass());
      hist_Eta->Fill(Mj_W[0].Eta());
      hist_Pt->Fill(Mj_W[0].Pt());
      hist_P->Fill(Mj_W[0].P());
      hist_Rapidity->Fill(Mj_W[0].Rapidity());
      hist_Electron_Pt ->Fill(Mj_Electrons[0].pt());
      hist_Neutrino_Pt ->Fill(Mj_Neutrinos[0].pt());
      hist_Electron_Eta ->Fill(Mj_Electrons[0].eta());
      hist_Neutrino_Eta ->Fill(Mj_Neutrinos[0].eta());
      hist_Electron_Rapidity ->Fill(Mj_Electrons[0].Rapidity());
      hist_Neutrino_Rapidity ->Fill(Mj_Neutrinos[0].Rapidity());
      
    
      hist_PTvsYvsEta->Fill(TMath::Log10(Mj_W[0].Pt()),fabs(Mj_W[0].Rapidity()),
			    fabs(Mj_Electrons[0].eta()));
      hist_Beta_MassCut->Fill(Mj_W[0].Beta());
      hist_Mass_MassCut->Fill(Mj_W[0].mass());
      hist_Eta_MassCut->Fill(Mj_W[0].Eta());
      hist_Pt_MassCut->Fill(Mj_W[0].Pt());
      hist_P_MassCut->Fill(Mj_W[0].P());
      hist_Rapidity_MassCut->Fill(Mj_W[0].Rapidity());
      hist_Electron_Pt_MassCut ->Fill(Mj_Electrons[0].pt());
      hist_Neutrino_Pt_MassCut ->Fill(Mj_Neutrinos[0].pt());
      hist_Electron_Eta_MassCut ->Fill(Mj_Electrons[0].eta());
      hist_Neutrino_Eta_MassCut ->Fill(Mj_Neutrinos[0].eta());
      hist_Electron_Rapidity_MassCut ->Fill(Mj_Electrons[0].Rapidity());
      hist_Neutrino_Rapidity_MassCut ->Fill(Mj_Neutrinos[0].Rapidity());
      
      
      hist_PTvsETA_MassCut->Fill(Mj_W[0].Pt(),Mj_W[0].Eta());
      hist_PTvsY_MassCut->Fill(TMath::Log10(Mj_W[0].Pt()),fabs(Mj_W[0].Rapidity()));
      hist_YvsETA_MassCut->Fill(Mj_W[0].Rapidity(),Mj_W[0].Eta());
      

    
    
      if (ptselect && ElectronInacc){
	
	hist_Beta_ElectronID->Fill(Mj_W[0].Beta());
	hist_Mass_ElectronID->Fill(Mj_W[0].mass());
	hist_Eta_ElectronID->Fill(Mj_W[0].Eta());
	hist_Pt_ElectronID->Fill(Mj_W[0].Pt());
	hist_P_ElectronID->Fill(Mj_W[0].P());
	hist_Rapidity_ElectronID->Fill(Mj_W[0].Rapidity());
	hist_Electron_Pt_ElectronID ->Fill(Mj_Electrons[0].pt());
	hist_Neutrino_Pt_ElectronID ->Fill(Mj_Neutrinos[0].pt());
	hist_Electron_Eta_ElectronID ->Fill(Mj_Electrons[0].eta());
	hist_Neutrino_Eta_ElectronID ->Fill(Mj_Neutrinos[0].eta());
	hist_Electron_Rapidity_ElectronID ->Fill(Mj_Electrons[0].Rapidity());
	hist_Neutrino_Rapidity_ElectronID ->Fill(Mj_Neutrinos[0].Rapidity());
	
      }
      
      
      if(ElectronInacc && ptselect && NeutrinoInacc && neuptselect){
	
	hist_Beta_Selected_NuInAcc->Fill(Mj_W[0].Beta());
	hist_Mass_Selected_NuInAcc->Fill(Mj_W[0].mass());
	hist_Eta_Selected_NuInAcc->Fill(Mj_W[0].Eta());
	hist_Pt_Selected_NuInAcc->Fill(Mj_W[0].Pt());
	hist_P_Selected_NuInAcc->Fill(Mj_W[0].P());
	hist_Rapidity_Selected_NuInAcc->Fill(Mj_W[0].Rapidity());
	hist_Electron_Pt_Selected_NuInAcc ->Fill(Mj_Electrons[0].pt());
	hist_Neutrino_Pt_Selected_NuInAcc ->Fill(Mj_Neutrinos[0].pt());
	hist_Electron_Eta_Selected_NuInAcc ->Fill(Mj_Electrons[0].eta());
	hist_Neutrino_Eta_Selected_NuInAcc ->Fill(Mj_Neutrinos[0].eta());
	hist_Electron_Rapidity_Selected_NuInAcc ->Fill(Mj_Electrons[0].Rapidity());
	hist_Neutrino_Rapidity_Selected_NuInAcc ->Fill(Mj_Neutrinos[0].Rapidity());
	
	
	hist_Beta_Selected->Fill(Mj_W[0].Beta());
	hist_Mass_Selected->Fill(Mj_W[0].mass());
	hist_Eta_Selected->Fill(Mj_W[0].Eta());
	hist_Pt_Selected->Fill(Mj_W[0].Pt());
	hist_P_Selected->Fill(Mj_W[0].P());
	hist_Rapidity_Selected->Fill(Mj_W[0].Rapidity());
	hist_Electron_Pt_Selected ->Fill(Mj_Electrons[0].pt());
	hist_Neutrino_Pt_Selected ->Fill(Mj_Neutrinos[0].pt());
	hist_Electron_Eta_Selected ->Fill(Mj_Electrons[0].eta());
	hist_Neutrino_Eta_Selected ->Fill(Mj_Neutrinos[0].eta());
	hist_Electron_Rapidity_Selected ->Fill(Mj_Electrons[0].Rapidity());
	hist_Neutrino_Rapidity_Selected ->Fill(Mj_Neutrinos[0].Rapidity());
	
      
      }
      
    }//wCharge
  


  }// Mj_W==1

//   cout << " Mj_Z.size() " <<  Mj_Z.size() << ", Mj_W.size() " << Mj_W.size() << ", Mj_Electrons.size() " <<  Mj_Electrons.size() 
//        << ", Mj_Neutrinos.size() " << Mj_Neutrinos.size() << endl;
//   if( (Mj_Z.size()==1&&
//        Mj_Electrons.size()==2) || 
//       ( Mj_W.size()==1&&
// 	Mj_Electrons.size()==1 && 
// 	Mj_Neutrinos.size()==1)  ){

//     // Match Gen - Reco 
//     double dRmin = 0.2, dR = 0.5;
//     LorentzVector<ROOT::Math::PxPyPzE4D<float> > mcTag, mcProbe;
//     for (unsigned int iGen = 0; iGen < Mj_Electrons.size(); iGen++){	
//       double mcEt = Mj_Electrons[iGen].Et();
// //       double mcEta = Mj_Electrons[iGen].eta();
//       ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> >  momentum;
//       momentum.SetXYZ(Mj_Electrons[iGen].Px(), Mj_Electrons[iGen].Py(),Mj_Electrons[iGen].Pz());
//       double mcEta = ecalEta(momentum,ev.GetvertexPosition(0));
//       cout << "ecalEta " << mcEta << " eta " << Mj_Electrons[iGen].eta()<< endl;
//       double fabsMcEta = fabs(mcEta);
//       h_mc_Et_->Fill(mcEt);
//       h_mc_Eta_->Fill(mcEta);
//       h_mc_fabsEta_->Fill(fabsMcEta);
      

// //       for (unsigned int iGen = 0; iGen < Mj_Electrons.size(); iGen++){	
// // 	double mcEt = Mj_Electrons[iGen].Et();
// // 	double mcEta = Mj_Electrons[iGen].eta();
// // 	double fabsMcEta = fabs(mcEta);
// // 	h_mc_Et_->Fill(mcEt);
// // 	h_mc_Eta_->Fill(mcEta);
// // 	h_mc_fabsEta_->Fill(fabsMcEta);
// //       }

//     }
//     for (std::vector<Lepton const *>::const_iterator lep=ev.LD_CommonElectrons().accepted.begin();
// 	 lep != ev.LD_CommonElectrons().accepted.end();
// 	 ++lep){
//       int i = (*lep)->GetIndex();	
//       bool chargeOk=(chChk) ? (((*lep)->GetCharge()==ev.GetElectronSCCharge(i)) && ((*lep)->GetCharge()==ev.GetElectronKFCharge(i))) : true;
//       bool convOk= (convChk) ? ((ev.GetElectronGsfTrackTrackerExpectedHitsInner(i) <= mMissingHits) &&
// 				(fabs(ev.GetElectronDCot(i)) > mDCot || fabs(ev.GetElectronDist(i)) > mDist)) : true;
      
      
//       if (!fid((*lep)->Eta())) continue;
      
      
//       bool iso = passIsolation((*lep)->GetTrkIsolation()/(*lep)->Pt(),
// 			       (*lep)->GetEcalIsolation()/(*lep)->Pt(),
// 			       (*lep)->GetHcalIsolation()/(*lep)->Pt(),(*lep)->Eta(),ieff);
//       bool id = passID(ev.GetElectronSigmaIetaIeta(i),
// 		       ev.GetElectronDeltaPhiAtVtx(i),
// 		       ev.GetElectronDeltaEtaAtVtx(i), 
// 		       ev.GetElectronHoE(i),(*lep)->Eta(),ieff);
//       int index = 0, tagIndex = -1, probeIndex = -1;
//       if (iso && id && chargeOk && convOk){
// 	if ((*lep)->Et()>elecET){
// 	  for (unsigned int iGen = 0; iGen < Mj_Electrons.size(); iGen++){	
// 	    double mcEt = Mj_Electrons[iGen].Et();
// 	    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> >  momentum;
// 	    momentum.SetXYZ(Mj_Electrons[iGen].Px(), Mj_Electrons[iGen].Py(),Mj_Electrons[iGen].Pz());
// 	    double mcEta = ecalEta(momentum,ev.GetvertexPosition(0));
//       //	    double mcEta = Mj_Electrons[iGen].eta();
// 	    double fabsMcEta = fabs(mcEta);
// 	    h_mc_Et_->Fill(mcEt);
// 	    h_mc_Eta_->Fill(mcEta);
// 	    h_mc_fabsEta_->Fill(fabsMcEta);
// 	    double deltaEta = mcEta - (*lep)->Eta();
// 	    double deltaPhi = Mj_Electrons[iGen].phi() - (*lep)->Phi();
// 	    while (deltaPhi > pigreco) deltaPhi -= 2*pigreco;
// 	    while (deltaPhi <= -pigreco) deltaPhi += 2*pigreco;
// 	    double deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
// 	    double deltaR = sqrt(deltaR2);
// 	    dR=deltaR;
// 	    if (deltaR < dRmin){
// 	      mcTag = Mj_Electrons[iGen];
// 	      dRmin = dR;
// 	      tagIndex = index;
// 	    }
// 	    index++;
// 	  }//iGen
// 	  if (tagIndex !=-1) {
// 	    h_tagMatch_dR_->Fill(dR);
// // 	    probeIndex = (tagIndex+1)%2;
// // 	    mcProbe = Mj_Electrons[probeIndex];
// 	    cout << " tag index " << tagIndex << ", index " << index <<  endl;
// 	  }
	  
// 	}//Et
//       }//iso, id, charge, convOk 
//       cout << " .-------------------- " << endl;
      
//       if (tagIndex ==-1) continue;
// //       double mcEt = mcProbe.Et();
// //       double mcEta = mcProbe.eta();
// //       double fabsMcEta = fabs(mcEta);

//       h_Et_->Fill((*lep)->Et());
//       h_Eta_->Fill((*lep)->Eta());
//       h_fabsEta_->Fill(fabs((*lep)->Eta()));
//     }//lep
  




//   }// end of if( Mj_Z.size()==1&&Mj_Electrons.size()==2)









  return true; //No WENU in this event
}//End Method

// -----------------------------------------------------------------------------
std::ostream& KinCorrection::Description( std::ostream& ostrm ) {
  ostrm << "KinCorrection";
  return ostrm;
}
//------------------------------------------------------------------------------
bool KinCorrection::isInBarrel(double eta)
{
	return (fabs(eta) < BarrelEtaMax_);
}

bool KinCorrection::isInEndCap(double eta)
{
	return (fabs(eta) < EndCapEtaMax_ && fabs(eta) > EndCapEtaMin_);
}

bool KinCorrection::isInFiducial(double eta)
{
	return isInBarrel(eta) || isInEndCap(eta);
}

//------------------
bool KinCorrection::passIsolation (double track, double ecal, double hcal, double eta, int ieff)
{
  
  return CheckCuts(track, ecal, hcal, 0., 0., 0., 0.,  eta, ieff );
  
}

bool KinCorrection::passID (double sihih, double dfi, double dhi,double hoe, double eta, int ieff)
{
  
  return CheckCuts(0., 0., 0., sihih, dfi, dhi, hoe,  eta, ieff );

}


bool KinCorrection::fid(double eta) {
  return  (fabs(eta)<2.4 && (  fabs(eta) < 1.4442  || fabs(eta) > 1.56 ));
}

bool KinCorrection::CheckCuts(double v_trk, double v_ecal, double v_hcal, 
				     double v_sihih, double v_dfi, double v_dhi, double v_hoe,
				     double eta, int ieff)
{
  
  
  if (fabs(eta)< 1.479) {	  
    if ( v_trk  <  Trk[ieff]    && 
	 v_ecal <  Ecal[ieff]   &&
	 v_hcal <  Hcal[ieff]   &&
	 v_sihih < sihih[ieff]  &&
	 fabs(v_dfi) < Dphi[ieff]   &&
	 fabs(v_dhi) < Deta[ieff]   &&
	 fabs(v_hoe)< HoE[ieff]    
	 ) return true;
  }
  else {

    if ( v_trk <  Trk_ee[ieff]    && 
	 v_ecal < Ecal_ee[ieff]   &&
	 v_hcal < Hcal_ee[ieff]   &&
	 v_sihih <sihih_ee[ieff]  &&
	 fabs(v_dfi) <Dphi_ee[ieff]   &&
	 //MICHELE DA SCOMMENTARE   
	 //	fabs(v_dhi) <Deta_ee[ieff]     &&
	 fabs(v_hoe)< HoE_ee[ieff] 
	 ) return true;

  }
  return false;
}


// -----------------------------------------------------------------------------
//
double KinCorrection::ecalEta(ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> momentum, ICF_XYZPoint vertex)
{
static const float R_ECAL           = 136.5;
static const float Z_Endcap         = 328.0;
static const float etaBarrelEndcap  = 1.479;
static const float pi = acos(-1.);
        float etaParticle = momentum.eta();
        float vZ = vertex.z();
        float vRho = vertex.Rho();

        float theta = pi/2.;
        float zEcal = vZ;

        if(etaParticle != 0.0) zEcal = (R_ECAL-vRho)*sinh(etaParticle)+vZ;

        if(zEcal != 0.0) theta = atan(R_ECAL/zEcal);
        if(theta < 0.0) theta = theta + pi ;

        float ETA = 0.0;
        if(theta != 0.0) ETA = - log(tan(0.5*theta));

        if( fabs(ETA) > etaBarrelEndcap )
        {
                float Zend = Z_Endcap ;
                if(etaParticle < 0.0 )  Zend = - Zend;
                float Zlen = Zend - vZ ;
                float RR = Zlen/sinh(etaParticle);
                theta = atan((RR+vRho)/Zend);
                if(theta < 0.0) theta = theta+pi;
                ETA = -log(tan(0.5*theta));
        }
        return ETA;
}

// -----------------------------------------------------------------------------
//

void
KinCorrection::RescaleZee(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *Z, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *e1, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *e2)
{
// Method to rescale a Z event to "look like" a W event"

  ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > Booster = (ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> >)Z->BoostToCM();
  TVector3 B(Booster.X(),Booster.Y(),Booster.Z());
  
  TLorentzVector ZP(Z->Px(),Z->Py(),Z->Pz(),Z->E());
  TLorentzVector E1P(e1->Px(),e1->Py(),e1->Pz(),e1->E());
  TLorentzVector E2P(e2->Px(),e2->Py(),e2->Pz(),e2->E());
  
  // Boost all 3 into rest frame of Z
  ZP.Boost(B);
  E1P.Boost(B);
  E2P.Boost(B);
  
  // Rescale the mass of the mother
  // This seems like a long way of doing it but this way any COM Energy rescale
  // can be performed.
  
  double E = ZP.E();
  //E = (E-MZ_)*(GW_/GZ_) + MW_;
  E = E*mW_/mZ_;
  ZP.SetE(E);
  double pe1 = E1P.P();
  double pe2 = E2P.P();
  E1P.SetPxPyPzE(E1P.Px()*E/(2*pe1),E1P.Py()*E/(2*pe1),E1P.Pz()*E/(2*pe1), E/2 );
  E2P.SetPxPyPzE(E2P.Px()*E/(2*pe2),E2P.Py()*E/(2*pe2),E2P.Pz()*E/(2*pe2), E/2 );
  
  TVector3 reboostB;
  TLorentzVector W(Z->Px(),Z->Py(), Z->Pz(),
		   TMath::Sqrt(E*E+Z->Px()*Z->Px()+ Z->Py()*Z->Py() + Z->Pz()*Z->Pz()));
  reboostB = W.BoostVector();
  
  ZP.Boost(reboostB);
  E1P.Boost(reboostB);
  E2P.Boost(reboostB);
  
  Z->SetPxPyPzE(ZP.Px(),ZP.Py(),ZP.Pz(),ZP.E());
  e1->SetPxPyPzE(E1P.Px(),E1P.Py(),E1P.Pz(),E1P.E());
  e2->SetPxPyPzE(E2P.Px(),E2P.Py(),E2P.Pz(),E2P.E());
  
}

// -----------------------------------------------------------------------------
double
KinCorrection::GetWeight3D(Event::GenObject Z,double eta){

  double pt = TMath::Log10(Z.Pt());
  cout << " pt " << pt << endl;
  double y  = fabs(Z.Rapidity());
  cout << " Rapidity " << y << endl;
  double weight;
  if (pt > 3. || pt < -1. || y > 6. || fabs(eta) > 6.){
    weight = 1.0;
  }  else{
    cout << PtYWeight_ ->GetName()<< endl;
    weight = PtYWeight_->GetBinContent(PtYWeight_->FindBin(pt,y,fabs(eta)));
    cout << " weight " << weight << endl;
  }
  return weight;

}
//


