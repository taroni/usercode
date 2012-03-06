//Plot the validation histograms
// last update 9th Sep 2010
//----------------------------------------------
// to run  in ROOT 
// root [1] .L newhisto.C+
// root [2] newhisto("decay", saveFig)
// where decay is ZMuMu, JpsiMuMu, UpsilonMuMU
// and saveFig = 1 or 0, if 1 the canvases will 
// be saved in png and eps format
//-----------------------------------------------
#include <memory>

#include "Riostream.h" 
#include <math.h>
#include <sstream>
#include <vector>

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

using namespace std;
void myMass(string decay, bool saveFig= 0 , bool onlyFirstCandidate= 0){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gROOT->ForceStyle();
  double pigreco = 3.141592;

  stringstream file;
  file.str("");
  file << "variables"<<decay<<".root";

  cout << __LINE__ << " Opening file: " << file.str().c_str() << endl;
  
  TFile * f = new TFile();
  
  f = TFile::Open(file.str().c_str());
  if (!f){
    cout << "FILE NOT FOUND" << endl;
    return;
  }
  cout << __LINE__ << " file " << file.str().c_str() <<  " opened" << endl;
  stringstream outstring;
  outstring.str("");
  outstring << "outfileMyMass" << decay << ".root";
  
  TFile* outfile = new TFile(outstring.str().c_str(),"RECREATE");
  
  outfile -> cd();

  double etaPlus= 0;
  double etaMinus =0;
  double phiPlus = 0;
  double phiMinus = 0;

  double etaMu1= 0;
  double etaMu2 =0;
  double phiMu1 = 0;
  double phiMu2 = 0;
  double ptMu1 = 0;
  double ptMu2 = 0;
   
  double nbins = 0;
  double  firstbin = 0;
  double  lastbin  = 0;

  double maxMass =0, minMass=0;

  string str3;
  size_t pos;
  pos = decay.find("MuMu");    // position of "live" in str
  str3 = decay.substr (0, pos);   // get from "live" to the end
  cout << str3 << endl;

//   if (decay == "ZMuMu" || decay == "ZMuMuMC" ) {
  if (str3 == "Z" ) {
    nbins = 140;
    firstbin = 50;
    lastbin  = 120;
    maxMass = 115;
    minMass = 65;
  }  else if (str3 == "Jpsi") {
//     nbins = 80;
    nbins = 300;
//     firstbin = 2;
//     lastbin  = 4;
    firstbin = 2.;
    lastbin  = 4.1;
    maxMass = 3.4;
    minMass = 2.7;
  }  else if (str3 == "Upsilon") {
    nbins = 80;
    firstbin = 8.5;
    lastbin  = 10.5;
    maxMass = 9.9;
    minMass = 8.9;
  }


  
  
  TH1D *hchi2 = new TH1D("hchi2","#chi^{2}/ndf",100,0,5.);
  TH1D *hNtrk= new TH1D("hNtrk","Number of Tracks",10,0.,10.);
  TH1D *hMass= new TH1D("hMass","Invariant Mass",nbins,firstbin,lastbin);
  TH1D * hMassFirstCand  = new TH1D("hMassFirstCand","hMassFirstCand",nbins,firstbin,lastbin);
  TH1D * hMassSecondCand = new TH1D("hMassSecondCand","hMassSecondCand",nbins,firstbin,lastbin);
  TH1D * hMassThirdCand  = new TH1D("hMassThirdCand","hMassThirdCand",nbins,firstbin,lastbin);
  TH1D * hMassForthCand  = new TH1D("hMassForthCand","hMassForthCand",nbins,firstbin,lastbin);
  TH1D * hMassFifthCand  = new TH1D("hMassFifthCand","hMassFifthCand",nbins,firstbin,lastbin);
  TH1D * hMassOneCand = new TH1D("hMassOneCand","hMassOneCand",nbins,firstbin,lastbin);
  TH1D * hMassTwoCand = new TH1D("hMassTwoCand","hMassTwoCand",nbins,firstbin,lastbin);
  TH1D * hMassThreeCand = new TH1D("hMassThreeCand","hMassThreeCand",nbins,firstbin,lastbin);
  TH1D * hMassFourCand = new TH1D("hMassFourCand","hMassFourCand",nbins,firstbin,lastbin);
  TH1D * hMassFiveCand = new TH1D("hMassFiveCand","hMassFiveCand",nbins,firstbin,lastbin);

  TH1D *hP= new TH1D("hP","Momentum",200,0.,200.);
  TH1D *hPPlus= new TH1D("hPPlus","Momentum, #mu^{+}",200,0.,200.);
  TH1D *hPMinus= new TH1D("hPMinus","Momentum, #mu^{-}",200,0.,200.);
  TH1D *hPt= new TH1D("hPt","Transverse Momentum",200,0.,200.);
  TH1D *hMinPt= new TH1D("hMinPt","Transverse Momentum, lowest of the muon pair",200,0.,200.);
  TH1D *hPtPlus= new TH1D("hPtPPlus","Transverse Momentum, #mu^{+}",200,0.,200.);
  TH1D *hPtMinus= new TH1D("hPtMinus","Transverse Momentum, #mu^{-}",200,0.,200.);
  TH1I *hHit= new TH1I("hHit","Number of hit",30,0,30);
  TH1I *hHitPlus= new TH1I("hHitPlus","Number of hit, #mu^{+}",30,0,30);
  TH1I *hHitMinus= new TH1I("hHitMinus","Number of hit, #mu^{-}",30,0,30);


  TH1D *hEta= new TH1D("hEta","Eta",100,-4.,4.);
  TH1D *hEtaPlus= new TH1D("hEtaPlus","#eta, #mu^{+}",100,-4.,4.);
  TH1D *hEtaMinus= new TH1D("hEtaMinus","#eta, #mu^{-}",100,-4.,4.);
  TH1D *hPhi= new TH1D("hPhi","Phi",100,-4.,4.);
  TH1D *hPhiPlus= new TH1D("hPhiPlus","#phi, #mu^{+}",100,-4.,4.);
  TH1D *hPhiMinus= new TH1D("hPhiMinus","#phi, #mu^{-}",100,-4.,4.);

  TH1D *hDeltaPhi= new TH1D("hDeltaPhi","DeltaPhi",100,-10.,10.);
  TH1D *hDeltaEta= new TH1D("hDeltaEta","DeltaEta",100,-10.,10.);
  TH1D *hDeltaR  = new TH1D("R"  ,"R"  ,100,-1.,10.);

  TH2D *hd0DeltaPhivsPtMu1  = new TH2D("hd0DeltaPhivsPtMu1"  ,"hd0DeltaPhivsPtMu1"  ,100,0.,100.,100,-10.,10.);
  TH2D *hd0DeltaEtavsPtMu1  = new TH2D("hd0DeltaEtavsPtMu1"  ,"hd0DeltaEtavsPtMu1"  ,100,0.,100.,100,-5.,5.);
  TH2D *hd0DeltaPhivsPtMu2  = new TH2D("hd0DeltaPhivsPtMu2"  ,"hd0DeltaPhivsPtMu2"  ,100,0.,100.,100,-10.,10.);
  TH2D *hd0DeltaEtavsPtMu2  = new TH2D("hd0DeltaEtavsPtMu2"  ,"hd0DeltaEtavsPtMu2"  ,100,0.,100.,100,-5.,5.);


  TH1D *hvx  = new TH1D("hvx"  ,"hvx"  ,100,0.,.5);
  TH1D *hvy  = new TH1D("hvy"  ,"hvy"  ,100,0.,.5);
  TH1D *hvz  = new TH1D("hvz"  ,"hvz"  ,100,-10.,10.);
  TH1D *hd0  = new TH1D("hd0"  ,"hd0"  ,100,-1.,1.);


  TH2D *eta1vseta2 = new TH2D("eta1vseta2", "eta1vseta2", 250, -2.5, 2.5, 250,-2.5,2.5);

  TH2D *hd0vsphi  = new TH2D("hd0vsphi"  ,"hd0vsphi"  ,160,-3.20,3.20,100,-1.,1.);
  TH2D *hd0vseta  = new TH2D("hd0vseta"  ,"hd0vseta"  ,160,-3.20,3.20,100,-1.,1.);
  TH2D *hd0vspt   = new TH2D("hd0vspt"   ,"hd0vspt"   ,50,0.,100.,100,-0.25,0.25);

  TH1D *hnhpxb  = new TH1D("nhpxb"  ,"nhpxb"  ,10,0.,10.);
  TH1D *hnhpxe  = new TH1D("nhpxe"  ,"nhpxe"  ,10,0.,10.);
  TH1D *hnhTIB  = new TH1D("nhTIB"  ,"nhTIB"  ,20,0.,20.);
  TH1D *hnhTID  = new TH1D("nhTID"  ,"nhTID"  ,20,0.,20.);
  TH1D *hnhTOB  = new TH1D("nhTOB"  ,"nhTOB"  ,20,0.,20.);
  TH1D *hnhTEC  = new TH1D("nhTEC"  ,"nhTEC"  ,20,0.,20.);

  TH1D *hnhpxbe  = new TH1D("nhpxbe"  ,"nhpxbe"  ,20,0.,20.);
  TH1D *hnhTIBTOB  = new TH1D("nhTIBTOB"  ,"nhTIBTOB"  ,30,0.,30.);
  TH1D *hnhTIBTOBTECTID= new TH1D("hnhTIBTOBTECTID","hnhTIBTOBTECTID",40,0,40);

  TH2D *hnhpxbvspxe = new TH2D ("hnhpxbvspxe","hnhpxbvspxe",10,0,10,10,0,10);
  TH2D *hnhTOBvsTIB = new TH2D ("hnhTOBvsTIB","hnhTOBvsTIB",20,0,20,20,0,20);
  TH2D *hnhTIDvsTIB = new TH2D ("hnhTIDvsTIB","hnhTIDvsTIB",20,0,20,20,0,20);
  TH2D *hnhTECvsTIB = new TH2D ("hnhTECvsTIB","hnhTECvsTIB",20,0,20,20,0,20);
  TH2D *hnhTIDvsTOB = new TH2D ("hnhTIDvsTOB","hnhTIDvsTOB",20,0,20,20,0,20);
  TH2D *hnhTECvsTOB = new TH2D ("hnhTECvsTOB","hnhTECvsTOB",20,0,20,20,0,20);
  TH2D *hnhTECvsTID = new TH2D ("hnhTECvsTID","hnhTECvsTID",20,0,20,20,0,20);

  TH1D *hnhpxbMuSameCharge  = new TH1D("nhpxbMuSameCharge"  ,"nhpxbMuSameCharge"  ,10,0.,10.);
  TH1D *hnhpxeMuSameCharge  = new TH1D("nhpxeMuSameCharge"  ,"nhpxeMuSameCharge"  ,10,0.,10.);
  TH1D *hnhTIBMuSameCharge  = new TH1D("nhTIBMuSameCharge"  ,"nhTIBMuSameCharge"  ,20,0.,20.);
  TH1D *hnhTIDMuSameCharge  = new TH1D("nhTIDMuSameCharge"  ,"nhTIDMuSameCharge"  ,20,0.,20.);
  TH1D *hnhTOBMuSameCharge  = new TH1D("nhTOBMuSameCharge"  ,"nhTOBMuSameCharge"  ,20,0.,20.);
  TH1D *hnhTECMuSameCharge  = new TH1D("nhTECMuSameCharge"  ,"nhTECMuSameCharge"  ,20,0.,20.);
  TH1D *hEtaMuSameCharge = new TH1D("hEtaMuSameCharge","EtaMuSameCharge",100,-4.,4.);
  TH1D *hPhiMuSameCharge = new TH1D("hPhiMuSameCharge","PhiMuSameCharge",100,-4.,4.);
  TH1D *hPtMuSameCharge = new TH1D("hPtMuSameCharge","PtMuSameCharge",100,0.,100.);
  TH1D *hEtaBothMuP = new TH1D("hEtaBothMuP","EtaBothMuP",100,-4.,4.);
  TH1D *hPhiBothMuP = new TH1D("hPhiBothMuP","PhiBothMuP",100,-4.,4.);
  TH1D *hEtaBothMuM = new TH1D("hEtaBothMuM","EtaBothMuM",100,-4.,4.);
  TH1D *hPhiBothMuM = new TH1D("hPhiBothMuM","PhiBothMuM",100,-4.,4.);

  TH1I *hMultCand = new TH1I("hMultCand","MultipleCandidate",50,-.5,49.5);

  f->cd();

  TTree *t2 = (TTree*)f->Get("myanalysis/AnalysisTree");
  
  TBranch *chi2 = t2->GetBranch("Chi");
  TBranch *Ntk  = t2->GetBranch("mol");
  TBranch *M= t2->GetBranch("invMass");
  TBranch *P= t2->GetBranch("P");
  TBranch *Pt= t2->GetBranch("Pt");
  TBranch *Px= t2->GetBranch("Px");
  TBranch *Py= t2->GetBranch("Py");
  TBranch *Pz= t2->GetBranch("Pz");
  TBranch *Hit= t2->GetBranch("Hit");
  TBranch *muQ= t2->GetBranch("muQ");
  TBranch *Eta= t2->GetBranch("Eta");
  TBranch *Phi= t2->GetBranch("Phi");
  
  TBranch *Vx= t2->GetBranch("vx");
  TBranch *Vy= t2->GetBranch("vy");
  TBranch *Vz= t2->GetBranch("vz");
  TBranch *D0= t2->GetBranch("d0");

  TBranch *nHitPixB= t2->GetBranch("nhpxb");
  TBranch *nHitPixE= t2->GetBranch("nhpxe");
  TBranch *nHitTIB= t2->GetBranch("nhTIB");
  TBranch *nHitTID= t2->GetBranch("nhTID");
  TBranch *nHitTOB= t2->GetBranch("nhTOB");
  TBranch *nHitTEC= t2->GetBranch("nhTEC");

  vector <double> *chi2pt=0;
  int Ntrk;
  vector <int> *Hitpt=0, *muCharge=0, *nhpxb=0, *nhpxe=0, *nhTIB=0, *nhTID=0, *nhTOB=0, *nhTEC=0;
  //    float mass;
  vector <double> *Ppt=0, *Ptpt=0, *Pxpt=0, *Pypt=0, *Pzpt=0, *Etapt=0, *Phipt=0, *vx=0, *vy=0,*vz=0,*d0=0, *mass=0;

  Ntk->SetAddress(&Ntrk);
  M->SetAddress(&mass);

  chi2->SetAddress(&chi2pt);
  P->SetAddress(&Ppt);
  Pt->SetAddress(&Ptpt);
  Px->SetAddress(&Pxpt);
  Py->SetAddress(&Pypt);
  Pz->SetAddress(&Pzpt);
  Hit->SetAddress(&Hitpt);
  muQ->SetAddress(&muCharge);
  Eta->SetAddress(&Etapt);
  Phi->SetAddress(&Phipt);

  Vx ->SetAddress(&vx);
  Vy ->SetAddress(&vy);
  Vz ->SetAddress(&vz);
  D0 ->SetAddress(&d0);
  
  nHitPixB->SetAddress(&nhpxb);
  nHitPixE->SetAddress(&nhpxe);
  nHitTIB->SetAddress(&nhTIB);
  nHitTID->SetAddress(&nhTID);
  nHitTOB->SetAddress(&nhTOB);
  nHitTEC->SetAddress(&nhTEC);



  Int_t nentries = (Int_t)t2->GetEntries();
  
  cout << __LINE__ << " Entries " << nentries << endl;
  
  for ( int i=0;i<nentries;i++) {
    t2->GetEntry(i);      
//       if (mass < 3. ||  mass > 3.2 )continue;

//       if (chi2pt->size()!=2) {
// 	hMultCand->Fill( mass->size());
// 	cout << __LINE__ << " " <<i << " MORE THAN 2 MUONS " << endl;
// 	continue;
//       }
    hNtrk->Fill(Ntrk);
    unsigned int MassSize = mass->size();
    unsigned int nCand = 0;
    if (onlyFirstCandidate == 1 ) MassSize=1;
    hMultCand->Fill(MassSize);
    for (unsigned int ivec=0; ivec< MassSize;ivec++){
      if (mass->size()>0){
// 	if (mass->at(ivec) < maxMass+0.02 && mass->at(ivec) > minMass-0.02) {
	  hMass->Fill(mass->at(ivec));
	  if (mass->size()>1){
	    if (ivec==0) hMassFirstCand->Fill(mass->at(ivec));
	    if (ivec==1) hMassSecondCand->Fill(mass->at(ivec));
	    if (ivec==2) hMassThirdCand->Fill(mass->at(ivec));
	    if (ivec==3) hMassForthCand->Fill(mass->at(ivec));
	    if (ivec==4) hMassFifthCand->Fill(mass->at(ivec));
	  }
	  if (mass->size()==1)hMassOneCand->Fill(mass->at(ivec));
	  if (mass->size()==2)hMassTwoCand->Fill(mass->at(ivec));
	  if (mass->size()==3)hMassThreeCand->Fill(mass->at(ivec));
	  if (mass->size()==4)hMassFourCand->Fill(mass->at(ivec));
	  if (mass->size()==5)hMassFiveCand->Fill(mass->at(ivec));
// 	}
       }
    }
    if (chi2pt->size() >2 ) cout << __LINE__ << "================================" << endl;

    unsigned int muonSize = chi2pt->size();
    if (onlyFirstCandidate == 1) muonSize=2;
      

    for (unsigned int ivec=0; ivec< muonSize;ivec++){
      hchi2->Fill(chi2pt->at(ivec));
      hP->Fill(Ppt->at(ivec));
      hPt->Fill(Ptpt->at(ivec));
      
      hnhpxbvspxe ->Fill(nhpxe->at(ivec),nhpxb->at(ivec));
      
      hnhTOBvsTIB ->Fill(nhTIB->at(ivec),nhTOB->at(ivec));
      hnhTIDvsTIB ->Fill(nhTIB->at(ivec),nhTID->at(ivec));
      hnhTECvsTIB ->Fill(nhTIB->at(ivec),nhTEC->at(ivec));
      hnhTIDvsTOB ->Fill(nhTID->at(ivec),nhTOB->at(ivec));
      hnhTECvsTOB ->Fill(nhTOB->at(ivec),nhTEC->at(ivec));
      hnhTECvsTID ->Fill(nhTID->at(ivec),nhTEC->at(ivec));

      for (unsigned int jvec=ivec+1; jvec< muonSize;jvec++){
	
	etaMu1 = Etapt->at(ivec);
	phiMu1 = Phipt->at(ivec);
	ptMu1 = Ptpt->at(ivec);
	
	etaMu2 = Etapt->at(jvec);
	phiMu2 = Phipt->at(jvec);
	ptMu2 = Ptpt->at(jvec);
	
	double deltaEta = etaMu1 - etaMu2;
	double deltaPhi = phiMu1 - phiMu2;
	if (deltaPhi < - pigreco) deltaPhi = 2*pigreco + deltaPhi;
	if (deltaPhi >   pigreco) deltaPhi =-2*pigreco + deltaPhi;
	double R  = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
	double invMass =  sqrt(2*pow(0.105,2)+2*(sqrt((pow(0.105,2)+pow(Ppt->at(ivec),2))*(pow(0.105,2)+pow(Ppt->at(jvec),2)))
						 -(Pxpt->at(ivec)*Pxpt->at(jvec)+Pypt->at(ivec)*Pypt->at(jvec)+Pzpt->at(ivec)*Pzpt->at(jvec))));
	if ( muonSize >2 ) cout << __LINE__ << " mass " << invMass 
				<< ", mu1pt " << Ptpt->at(ivec) << ", ivec " << ivec
				<< ", mu2pt " << Ptpt->at(jvec) << ", jvec " << jvec
				<< endl;
	if (invMass > maxMass || invMass < minMass) continue;
// 	hMass-> Fill (invMass);
	if (Ptpt->at(ivec)<=Ptpt->at(jvec)){
	  hMinPt->Fill(Ptpt->at(ivec));
	  }else {
	  hMinPt->Fill(Ptpt->at(jvec));
	}
	eta1vseta2->Fill(etaMu2,etaMu1);
	hd0DeltaPhivsPtMu1->Fill(ptMu1,deltaPhi);
	hd0DeltaEtavsPtMu1->Fill(ptMu1,deltaEta);
	hd0DeltaPhivsPtMu2->Fill(ptMu2,deltaPhi);
	hd0DeltaEtavsPtMu2->Fill(ptMu2,deltaEta);

	hDeltaPhi -> Fill (deltaPhi);
	hDeltaEta -> Fill (fabs(fabs(deltaEta) - 2.5));
	hDeltaR   -> Fill (R) ;

	if ( muCharge->at(ivec)*muCharge->at(jvec)==+1){
	  hEtaMuSameCharge ->Fill(Etapt->at(ivec));
	  hPhiMuSameCharge ->Fill(Phipt->at(ivec));
	  hEtaMuSameCharge ->Fill(Etapt->at(jvec));
	  hPhiMuSameCharge ->Fill(Phipt->at(jvec));
	  if (nhpxb->at(ivec)!=0)hnhpxbMuSameCharge ->Fill(nhpxb->at(ivec));
	  if (nhpxe->at(ivec)!=0)hnhpxeMuSameCharge ->Fill(nhpxe->at(ivec));
	  if (nhTIB->at(ivec)!=0)hnhTIBMuSameCharge ->Fill(nhTIB->at(ivec));
	  if (nhTID->at(ivec)!=0)hnhTIDMuSameCharge ->Fill(nhTID->at(ivec));
	  if (nhTOB->at(ivec)!=0)hnhTOBMuSameCharge ->Fill(nhTOB->at(ivec));
	  if (nhTEC->at(ivec)!=0)hnhTECMuSameCharge ->Fill(nhTEC->at(ivec));
	  if (nhpxb->at(jvec)!=0)hnhpxbMuSameCharge ->Fill(nhpxb->at(jvec));
	  if (nhpxe->at(jvec)!=0)hnhpxeMuSameCharge ->Fill(nhpxe->at(jvec));
	  if (nhTIB->at(jvec)!=0)hnhTIBMuSameCharge ->Fill(nhTIB->at(jvec));
	  if (nhTID->at(jvec)!=0)hnhTIDMuSameCharge ->Fill(nhTID->at(jvec));
	  if (nhTOB->at(jvec)!=0)hnhTOBMuSameCharge ->Fill(nhTOB->at(jvec));
	  if (nhTEC->at(jvec)!=0)hnhTECMuSameCharge ->Fill(nhTEC->at(jvec));
	  
	  hPtMuSameCharge ->Fill(Ptpt->at(ivec));
	  hPtMuSameCharge ->Fill(Ptpt->at(ivec));
	  hPtMuSameCharge ->Fill(Ptpt->at(jvec));
	  hPtMuSameCharge ->Fill(Ptpt->at(jvec));
	  if ( muCharge->at(ivec)+muCharge->at(jvec)>0){
	    hEtaBothMuP ->Fill(Etapt->at(ivec));
	    hPhiBothMuP ->Fill(Phipt->at(ivec));
	    hEtaBothMuP ->Fill(Etapt->at(jvec));
	    hPhiBothMuP ->Fill(Phipt->at(jvec));
	  } else if ( muCharge->at(ivec)+muCharge->at(jvec)<0){
	    hEtaBothMuM ->Fill(Etapt->at(ivec));
	    hPhiBothMuM ->Fill(Phipt->at(ivec));
	    hEtaBothMuM ->Fill(Etapt->at(jvec));
	    hPhiBothMuM ->Fill(Phipt->at(jvec));
	  }
	}
      }//jvec
      
      hHit->Fill(Hitpt->at(ivec));
      hEta->Fill(Etapt->at(ivec));
      hPhi->Fill(Phipt->at(ivec));
	
      hnhpxbe->Fill(nhpxb->at(ivec)+nhpxe->at(ivec));
      hnhTIBTOB->Fill(nhTIB->at(ivec)+nhTOB->at(ivec));
      hnhTIBTOBTECTID->Fill(nhTIB->at(ivec)+nhTOB->at(ivec)+nhTEC->at(ivec)+nhTEC->at(ivec));

      hnhpxb ->Fill(nhpxb->at(ivec));
      hnhpxe ->Fill(nhpxe->at(ivec));
      hnhTIB ->Fill(nhTIB->at(ivec));
      hnhTID ->Fill(nhTID->at(ivec));
      hnhTOB ->Fill(nhTOB->at(ivec));
      hnhTEC ->Fill(nhTEC->at(ivec));
     
      hvx -> Fill(vx->at(ivec));
      hvy -> Fill(vy->at(ivec));
      hvz -> Fill(vz->at(ivec));
      hd0  -> Fill(d0->at(ivec));

      hd0vsphi ->Fill(Phipt->at(ivec),d0->at(ivec));
      hd0vseta ->Fill(Etapt->at(ivec),d0->at(ivec));
      hd0vspt  ->Fill(Ptpt->at(ivec) ,d0->at(ivec));
      
      if (muCharge->at(ivec)>0){
	hPPlus->Fill(Ppt->at(ivec));
	hPtPlus->Fill(Ptpt->at(ivec));
	hHitPlus->Fill(Hitpt->at(ivec));
	hEtaPlus->Fill(Etapt->at(ivec));
	hPhiPlus->Fill(Phipt->at(ivec));
	etaPlus = Etapt->at(ivec);
	phiPlus = Phipt->at(ivec);
      }
      if (muCharge->at(ivec)<0){
	hPMinus->Fill(Ppt->at(ivec));
	hPtMinus->Fill(Ptpt->at(ivec));
	hHitMinus->Fill(Hitpt->at(ivec));
	hEtaMinus->Fill(Etapt->at(ivec));
	hPhiMinus->Fill(Phipt->at(ivec));
	etaMinus = Etapt->at(ivec);
	phiMinus = Phipt->at(ivec);
      }

    }

      //   }//if
    
  }//entries
      
  
  f->Close();
  cout << __LINE__ << " file " << file.str().c_str() <<  " closed" << endl;
  outfile->cd();
  
 /////////////////////////////////////////////////////////////////////////////
   TCanvas *canv1= new TCanvas("canv1","firstAnalysis",1000,700);
   TCanvas *canv2= new TCanvas("canv2","firstPhiEtaAnalysis",1000,700);
   TCanvas *canv3= new TCanvas("canv3","DeltaPhiDeltaEtaR",1000,700);
   TCanvas *canv1Minus= new TCanvas("canv1Minus","MomentumMu+Mu-",1000,700);
   TCanvas *canv1Plus= new TCanvas("canv1Plus","nHitMu+Mu-",500,700);
   TCanvas *canv2Plus= new TCanvas("canv2Plus","PhiEtaMu+Mu-",1000,700);
   TCanvas *canv4=  new TCanvas("canv4","d0",1000,700);
   TCanvas *canv5=  new TCanvas("canv5","nhitDetector",1200,600);
   TCanvas *canv6=  new TCanvas("canv6","nhitDetectorSameCharge",1200,600);
   TCanvas *canv7=  new TCanvas("canv7","minPtofthePair",600,400);

   vector <TCanvas*> cVector;
   cVector.push_back(canv1);
   cVector.push_back(canv2);
   cVector.push_back(canv3);
   cVector.push_back(canv1Plus);
   cVector.push_back(canv1Minus);
   cVector.push_back(canv2Plus);
   cVector.push_back(canv4);
   cVector.push_back(canv5);
   cVector.push_back(canv6);
   cVector.push_back(canv7);

   canv1->Draw();
   canv1->Divide(2,2);
      
   canv1->cd(1);
   
   hchi2->SetFillStyle(3001);
   hchi2->SetFillColor(62);	
   hchi2->GetXaxis()->SetTitle("#chi^{2}/ndof");
   hchi2->Draw();

   canv1->cd(2);
   hP->SetFillStyle(3010);
   hP->SetFillColor(80);
   hP->GetXaxis()->SetTitle("P (GeV/c)");
   int maxBin = 0;
   for (int ibin = 1 ; ibin <= hP->GetNbinsX() ; ibin++){
     if (hP->GetBinContent(ibin) != 0) {
       maxBin = ibin;
     }
   }
   double maxX = (int )(hP->GetBinCenter(maxBin)/10)*10  + 11  ;
   hP->GetXaxis()->SetRangeUser(0,maxX);
   hP->Draw();

   canv1->cd(3);
   hPt->SetFillStyle(3033);
   hPt->SetFillColor(99);
   hPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
   maxBin = 0;
   for (int ibin = 1 ; ibin <= hPt->GetNbinsX() ; ibin++){
     if (hPt->GetBinContent(ibin) != 0) {
       maxBin = ibin;
     }
   }
   maxX = (int )(hPt->GetBinCenter(maxBin)/10)*10  + 11  ;
   hPt->GetXaxis()->SetRangeUser(0,maxX);
   hPt->Draw();
   
   canv1->cd(4);  
   hHit->SetFillStyle(3015);
   hHit->SetFillColor(51);
   hHit->GetXaxis()->SetTitle("Nhit");
   hHit->Draw();

   
   
   canv2->Divide(2,2);
   
   canv2->cd(1);
   hNtrk->SetFillStyle(3009);
   hNtrk->SetFillColor(1);
   hNtrk->GetXaxis()->SetTitle("Ntrk");
   hNtrk->Draw();
   
   canv2->cd(2);
   hMass->SetFillStyle(3020);
   hMass->SetFillColor(3);
   hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
   hMass->Draw();
   
   canv2->cd(3);
   hPhi->SetFillStyle(3009);
   hPhi->SetFillColor(2);
   hPhi->GetXaxis()->SetTitle("#phi (rad)");
   hPhi->Draw();
   
   canv2->cd(4);
   hEta->SetFillStyle(3020);
   hEta->SetFillColor(4);
   hEta->GetXaxis()->SetTitle("#eta");
   hEta->Draw();
 
   canv3->Divide(2,2);
   
   canv3->cd(1);
   hDeltaPhi->SetFillStyle(3009);
   hDeltaPhi->SetFillColor(1);
   hDeltaPhi->GetXaxis()->SetTitle("#phi");
   hDeltaPhi->Draw();
   
   canv3->cd(2);
   hDeltaEta->SetFillStyle(3020);
   hDeltaEta->SetFillColor(3);
   hDeltaEta->GetXaxis()->SetTitle("#eta");
   hDeltaEta->Draw();
   
   canv3->cd(3);
   hDeltaR->SetFillStyle(3009);
   hDeltaR->SetFillColor(2);
   hDeltaR->GetXaxis()->SetTitle("R");
   hDeltaR->Draw();

   canv1Minus->Divide(2,2);
   
   canv1Minus->cd(1);
   hPMinus->SetFillStyle(3010);
   hPMinus->SetFillColor(80);
   hPMinus->GetXaxis()->SetTitle("P (GeV/c)");
   maxBin = 0;
   for (int ibin = 1 ; ibin <= hPMinus->GetNbinsX() ; ibin++){
     if (hPMinus->GetBinContent(ibin) != 0) {
       maxBin = ibin;
     }
   }
   maxX = (int )(hPMinus->GetBinCenter(maxBin)/10)*10  + 11  ;
   hPMinus->GetXaxis()->SetRangeUser(0,maxX);
   hPMinus->Draw();

   canv1Minus->cd(2);
   hPtMinus->SetFillStyle(3033);
   hPtMinus->SetFillColor(99);
   hPtMinus->GetXaxis()->SetTitle("P_{T} (GeV/c)");
   maxBin = 0;
   for (int ibin = 1 ; ibin <= hPtMinus->GetNbinsX() ; ibin++){
     if (hPtMinus->GetBinContent(ibin) != 0) {
       maxBin = ibin;
     }
   }
   maxX = (int )(hPtMinus->GetBinCenter(maxBin)/10)*10  + 11  ;
   hPtMinus->GetXaxis()->SetRangeUser(0,maxX);
   hPtMinus->Draw();

   canv1Minus->cd(3);
   hPPlus->SetFillStyle(3010);
   hPPlus->SetFillColor(80);
   hPPlus->GetXaxis()->SetTitle("P (GeV/c)");
   maxBin = 0;
   for (int ibin = 1 ; ibin <= hPPlus->GetNbinsX() ; ibin++){
     if (hPPlus->GetBinContent(ibin) != 0) {
       maxBin = ibin;
     }
   }
   maxX = (int )(hPPlus->GetBinCenter(maxBin)/10)*10  + 11  ;
   hPPlus->GetXaxis()->SetRangeUser(0,maxX);
   hPPlus->Draw();

   canv1Minus->cd(4);
   hPtPlus->SetFillStyle(3033);
   hPtPlus->SetFillColor(99);
   hPtPlus->GetXaxis()->SetTitle("P_{T} (GeV/c)");
   maxBin = 0;
   for (int ibin = 1 ; ibin <= hPtPlus->GetNbinsX() ; ibin++){
     if (hPtPlus->GetBinContent(ibin) != 0) {
       maxBin = ibin;
     }
   }
   maxX = (int )(hPtPlus->GetBinCenter(maxBin)/10)*10  + 11  ;
   hPtPlus->GetXaxis()->SetRangeUser(0,maxX);
   hPtPlus->Draw();

   canv1Plus->Divide(1,2);
   canv1Plus->cd(1);  
   hHitMinus->SetFillStyle(3015);
   hHitMinus->SetFillColor(51);
   hHitMinus->GetXaxis()->SetTitle("Nhit");
   hHitMinus->Draw();
   canv1Plus->cd(2);  
   hHitPlus->SetFillStyle(3015);
   hHitPlus->SetFillColor(51);
   hHitPlus->GetXaxis()->SetTitle("Nhit");
   hHitPlus->Draw();
     
   canv2Plus->Divide(2,2);

   canv2Plus->cd(1);
   hPhiMinus->SetFillStyle(3009);
   hPhiMinus->SetFillColor(2);
   hPhiMinus->GetXaxis()->SetTitle("#phi (rad)");
   hPhiMinus->Draw();
   
   canv2Plus->cd(2);
   hEtaMinus->SetFillStyle(3020);
   hEtaMinus->SetFillColor(4);
   hEtaMinus->GetXaxis()->SetTitle("#eta");
   hEtaMinus->Draw();

   canv2Plus->cd(3);
   hPhiPlus->SetFillStyle(3009);
   hPhiPlus->SetFillColor(2);
   hPhiPlus->GetXaxis()->SetTitle("#phi (rad)");
   hPhiPlus->Draw();
   
   canv2Plus->cd(4);
   hEtaPlus->SetFillStyle(3020);
   hEtaPlus->SetFillColor(4);
   hEtaPlus->GetXaxis()->SetTitle("#eta");
   hEtaPlus->Draw();
   
   canv4->Draw();
   canv4->Divide(2,2);
   canv4->cd(1);
   hd0->Draw();
   canv4->cd(2);
   hd0vspt->Draw();
   canv4->cd(3);
   hd0vsphi->Draw();
   canv4->cd(4);
   hd0vseta->Draw();
   
   canv5->Draw();
   canv5->Divide(3,2);
   canv5->cd(1);
   hnhpxb->Draw();
   canv5->cd(2);
   hnhpxe->Draw();
   canv5->cd(3);
   hnhTIB->Draw();
   canv5->cd(4);
   hnhTID->Draw();
   canv5->cd(5);
   hnhTOB->Draw();
   canv5->cd(6);
   hnhTEC->Draw();

   canv6->Draw();
   canv6->Divide(3,2);
   canv6->cd(1);
   hnhpxbMuSameCharge->Draw();
   canv6->cd(2);
   hnhpxeMuSameCharge->Draw();
   canv6->cd(3);
   hnhTIBMuSameCharge->Draw();
   canv6->cd(4);
   hnhTIDMuSameCharge->Draw();
   canv6->cd(5);
   hnhTOBMuSameCharge->Draw();
   canv6->cd(6);
   hnhTECMuSameCharge->Draw();
   
   canv7->Draw();
   canv7->cd();
   hMinPt->Draw();
   hMinPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
   hMinPt->GetXaxis()->SetRangeUser(0.,30.);
   canv7->Update();
   
   if (saveFig == 1 ) {
     for (unsigned icanv = 0 ; icanv < cVector.size(); icanv++){
       stringstream imageName;
     // imageName.str("");
     // imageName <<decay<< cVector[icanv]->GetTitle()<<".jpg";
     // cVector[icanv]->SaveAs(imageName.str().c_str());
     // imageName.str("");
     // imageName <<decay<< cVector[icanv]->GetTitle()<<".gif";
     // cVector[icanv]->SaveAs(imageName.str().c_str());
       imageName.str("");
       imageName << "plots/"<<decay<< cVector[icanv]->GetTitle()<<".png";
       cVector[icanv]->SaveAs(imageName.str().c_str());
       imageName.str("");
       imageName <<"plots/"<<decay<< cVector[icanv]->GetTitle()<<".eps";
       cVector[icanv]->SaveAs(imageName.str().c_str());
     }
   }


   for (unsigned icanv = 0 ; icanv < cVector.size(); icanv++){
     cVector[icanv]->Write();
   }

   hchi2 ->Write();
   hNtrk->Write();
   hMass->Write();

   hMassFirstCand->Write();
   hMassSecondCand ->Write();
   hMassThirdCand  ->Write();
   hMassForthCand  ->Write();
   hMassFifthCand  ->Write();
   hMassOneCand ->Write();
   hMassTwoCand ->Write();
   hMassThreeCand ->Write();
   hMassFourCand ->Write();
   hMassFiveCand ->Write();

   hP->Write();
   hPPlus->Write();
   hPMinus->Write();
   hPt->Write();
   hMinPt->Write();
   hPtPlus->Write();
   hPtMinus->Write();
   hHit->Write();
   hHitPlus->Write();
   hHitMinus->Write();
   hEta->Write();
   hEtaPlus->Write();
   hEtaMinus->Write();
   hPhi->Write();
   hPhiPlus->Write();
   hPhiMinus->Write();

   hDeltaPhi->Write();
   hDeltaEta->Write();
   hDeltaR  ->Write();

   hvx -> Write();
   hvy -> Write();
   hvz -> Write();
   hd0  ->Write();

   hPtMuSameCharge ->Write();
   hEtaMuSameCharge ->Write();
   hPhiMuSameCharge ->Write();
   hEtaBothMuP ->Write();
   hPhiBothMuP ->Write();
   hEtaBothMuM ->Write();
   hPhiBothMuM ->Write();

   hnhpxbMuSameCharge ->Write();
   hnhpxeMuSameCharge ->Write();
   hnhTIBMuSameCharge ->Write();
   hnhTIDMuSameCharge ->Write();
   hnhTOBMuSameCharge ->Write();
   hnhTECMuSameCharge ->Write();

   hnhpxb ->Write();
   hnhpxe ->Write();
   hnhTIB ->Write();
   hnhTID ->Write();
   hnhTOB ->Write();
   hnhTEC ->Write();
   hnhTOBvsTIB ->Write();
   hnhTIDvsTIB ->Write();
   hnhTECvsTIB ->Write();
   hnhTIDvsTOB ->Write();
   hnhTECvsTOB ->Write();
   hnhTECvsTID ->Write();
   hnhpxbvspxe->Write();
   hnhpxbe ->Write();
   hnhTIBTOB ->Write();
   hnhTIBTOBTECTID ->Write();

   hd0vsphi ->Write();
   hd0vseta ->Write();
   hd0vspt  ->Write();

   hMultCand->Write();
   eta1vseta2->Write();
   cout << __LINE__ << endl;


   for (unsigned icanv = 0 ; icanv < cVector.size(); icanv++){
     cVector[icanv]->~TCanvas();
   }


   hchi2 ->~TH1D();
   hNtrk->~TH1D();
   hMass->~TH1D();

   hMassFirstCand->~TH1D();
   hMassSecondCand ->~TH1D();
   hMassThirdCand  ->~TH1D();
   hMassForthCand  ->~TH1D();
   hMassFifthCand  ->~TH1D();
   hMassOneCand ->~TH1D();
   hMassTwoCand ->~TH1D();
   hMassThreeCand ->~TH1D();
   hMassFourCand ->~TH1D();
   hMassFiveCand ->~TH1D();

   hP->~TH1D();
   hPPlus->~TH1D();
   hPMinus->~TH1D();
   hPt->~TH1D();
   hMinPt->~TH1D();
   hPtPlus->~TH1D();
   hPtMinus->~TH1D();
   hHit->~TH1I();
   hHitPlus->~TH1I();
   hHitMinus->~TH1I();
   hEta->~TH1D();
   hEtaPlus->~TH1D();
   hEtaMinus->~TH1D();
   hPhi->~TH1D();
   hPhiPlus->~TH1D();
   hPhiMinus->~TH1D();

   hDeltaPhi->~TH1D();
   hDeltaEta->~TH1D();
   hDeltaR  ->~TH1D();

   hPtMuSameCharge ->~TH1D();
   hEtaMuSameCharge ->~TH1D();
   hPhiMuSameCharge ->~TH1D();
   hEtaBothMuP ->~TH1D();
   hPhiBothMuP ->~TH1D();
   hEtaBothMuM ->~TH1D();
   hPhiBothMuM ->~TH1D();

   hnhpxbMuSameCharge ->~TH1D();
   hnhpxeMuSameCharge ->~TH1D();
   hnhTIBMuSameCharge ->~TH1D();
   hnhTIDMuSameCharge ->~TH1D();
   hnhTOBMuSameCharge ->~TH1D();
   hnhTECMuSameCharge ->~TH1D();


   hnhpxb ->~TH1D();
   hnhpxe ->~TH1D();
   hnhTIB ->~TH1D();
   hnhTID ->~TH1D();
   hnhTOB ->~TH1D();
   hnhTEC ->~TH1D();

   hd0vsphi ->~TH2D();
   hd0vseta ->~TH2D();
   hd0vspt  ->~TH2D();

   hvx -> ~TH1D();
   hvy -> ~TH1D();
   hvz -> ~TH1D();
   hd0  -> ~TH1D();

   hMultCand->~TH1I();

   outfile->Close();


}
