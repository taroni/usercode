//g++ -o histoTrbyTr histoTrbyTr.cc `root-config --cflags --libs`

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
// provare a fare il new dell'histo direttamente nella map
#define DEBUG 1

#include <memory>
#include <stdio.h>
#include <stdlib.h>

#include "Riostream.h" 
#include <math.h>
#include <sstream>
#include <map>
#include <vector>

#include <TApplication.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

using namespace std;

 
void histoTrbyTr(string decay, int saveFig= 0 , int onlyFirstCandidate= 0){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gStyle->SetHistLineWidth(2);
  //  gStyle->SetOptStat(111110);
  //gStyle->SetOptFit(111);

  gROOT->ForceStyle();
 //  double pigreco = 3.141592;

   stringstream file;
  file.str("");
  file << "rfdir:///castor/cern.ch/cms/store/caf/user/taroni/TkAlData2011-20110721/variables"<<decay<<".root";


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
   outstring << "triggerFile" << decay << ".root";

   TFile* outfile = new TFile(outstring.str().c_str(),"RECREATE");

   outfile -> cd();

//    double etaMu1;
//    double phiMu1;
//    double ptMu1;
//    double etaMu2;
//    double phiMu2;
//    double ptMu2;
   
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
    nbins = 80;
    firstbin = 2;
    lastbin  = 4;
    maxMass = 3.4;
    minMass = 2.7;
  }  else if (str3 == "Upsilon") {
    nbins = 80;
    firstbin = 8.5;
    lastbin  = 10.5;
    maxMass = 9.9;
    minMass = 8.9;
  }

  if (DEBUG)    cout<< __LINE__ << endl;

  TH1D * hTrigger = new TH1D("trigger","trigger",100,0.,100.); 
  
 //   TH1D *hMass= new TH1D("hMass","Invariant Mass",nbins,firstbin,lastbin);
//   TH1D *hPtMother= new TH1D("hPtMother","MotherTransverse Momentum",200,0.,200.);
//   TH1D *hEtaMother= new TH1D("hEtaMother","Mother Eta",100,-5.,5.);
//   TH1D *hPhiMother= new TH1D("hPhiMother","Mother Phi",100,-4.,4.);


  TH1D *hPt= new TH1D("hPt","Transverse Momentum",200,0.,200.);
  TH1D *hEta= new TH1D("hEta","Eta",100,-4.,4.);
  TH1D *hPhi= new TH1D("hPhi","Phi",100,-4.,4.);

//   TH1D *hDeltaPhi= new TH1D("hDeltaPhi","DeltaPhi",100,-10.,10.);
//   TH1D *hDeltaEta= new TH1D("hDeltaEta","DeltaEta",100,-10.,10.);
//   TH1D *hDeltaR  = new TH1D("R"  ,"R"  ,100,-1.,10.);

  TH1I *hMultCand = new TH1I("hMultCand","MultipleCandidate",50,-.5,49.5);

//   TH2D * hEta1vsDeltaPhi = new TH2D ("hEta1vsDeltaPhi", "EtaMu1vsDeltaPhi", 100, -2., 2. , 100, -4., 4.) ;
//   TH2D * hEta2vsDeltaPhi = new TH2D ("hEta2vsDeltaPhi", "EtaMu2vsDeltaPhi", 100, -2., 2. , 100, -4., 4.) ;

  if (DEBUG)  cout<< __LINE__ << endl;


  map <string,  vector <TH1D *> > hTrMap;
  f->cd();

    TTree *t2 = (TTree*)f->Get("myanalysis/AnalysisTree");
    if (DEBUG) cout<< __LINE__ << endl;

//     TBranch *EtaMother= t2->GetBranch("EtaMother");
//     TBranch *PhiMother= t2->GetBranch("PhiMother");
//     TBranch *PtMother= t2->GetBranch("PtMother");
   
    TBranch *chi2 = t2->GetBranch("Chi");
//     TBranch *Ntk  = t2->GetBranch("mol");
    TBranch *M= t2->GetBranch("invMass");
    TBranch *P= t2->GetBranch("P");
    TBranch *Pt= t2->GetBranch("Pt");
    TBranch *Px= t2->GetBranch("Px");
    TBranch *Py= t2->GetBranch("Py");
    TBranch *Pz= t2->GetBranch("Pz");
    TBranch *Eta= t2->GetBranch("Eta");
    TBranch *Phi= t2->GetBranch("Phi");

    TBranch *trigger= t2->GetBranch("triggerInfo");

      if (DEBUG)cout<< __LINE__ << endl;



    vector <double> *chi2pt=0;

    vector <double> *Ppt=0, *Ptpt=0, *Pxpt=0, *Pypt=0, *Pzpt=0, *Etapt=0, *Phipt=0;
//     vector<double> *mass=0, *PtMotherpt=0,*EtaMotherpt=0,*PhiMotherpt=0;

    vector <pair <string, bool> > *triggerInfo =0;

//     M->SetAddress(&mass);
//     PtMother->SetAddress(&PtMotherpt);
//     EtaMother->SetAddress(&EtaMotherpt);
//     PhiMother->SetAddress(&PhiMotherpt);

    P->SetAddress(&Ppt);
    Pt->SetAddress(&Ptpt);
    Px->SetAddress(&Pxpt);
    Py->SetAddress(&Pypt);
    Pz->SetAddress(&Pzpt);
    Eta->SetAddress(&Etapt);
    Phi->SetAddress(&Phipt);

    trigger->SetAddress(&triggerInfo);

    Int_t nentries = (Int_t)t2->GetEntries();
//     nentries = 100;
    cout << __LINE__ << " Entries " << nentries << endl;
    vector < string > triggerNames;

//     nentries = 100 ;
    outfile->cd();
   
   for ( int i=0;i<nentries;i++) {
      t2->GetEntry(i);
      hTrigger->Fill(triggerInfo->size());

      for (unsigned int itr = 0; itr < triggerInfo->size(); itr++){
	if (triggerInfo->at(itr).second== 1) {
	  string triggerName= (triggerInfo->at(itr)).first ;
	  bool checkName = false;
	  for( unsigned int iname =0; iname<triggerNames.size(); iname++){

	    if (triggerNames.at(iname).compare(triggerName)==0){
 	      checkName = true;
	    }
	  }
	  if (checkName == false ) {
	    triggerNames.push_back(triggerName);

	  }
	}
      }
  
   }
    cout << "triggersize" << " " << triggerNames.size()<< endl;

   for (unsigned int itr=0; itr<triggerNames.size(); itr++){
      
      vector < TH1D *> hTrVector; 
      stringstream histoEtaName, histoPhiName, histoPtName, histoEtaMotherName, histoPhiMotherName;
      histoEtaName.str("");
      histoPhiName.str("");
      histoPtName.str("");
      histoEtaName << triggerNames.at(itr) << "TrackEta";
      histoPhiName << triggerNames.at(itr) << "TrackPhi";
      histoPtName  << triggerNames.at(itr) << "TrackPt";
      if (DEBUG)cout<< __LINE__ << endl;
      if (DEBUG) cout<< __LINE__ << endl;
      hTrMap[triggerNames.at(itr)].push_back(new TH1D (histoEtaName.str().c_str(),histoEtaName.str().c_str(),500,-5,5));
      hTrMap[triggerNames.at(itr)].push_back( new TH1D (histoPhiName.str().c_str(),histoPhiName.str().c_str(),640,-3.2,3.2));
      hTrMap[triggerNames.at(itr)].push_back( new TH1D (histoPtName.str().c_str(),histoPtName.str().c_str(),200,0.,0.));
      
      hTrVector.clear();
   }

      if (DEBUG)cout<< __LINE__ << endl;
    f->cd();
    t2->ResetBranchAddresses();

//     M->SetAddress(&mass);

//     PtMother->SetAddress(&PtMotherpt);
//     EtaMother->SetAddress(&EtaMotherpt);
//     PhiMother->SetAddress(&PhiMotherpt);

    P->SetAddress(&Ppt);
    Pt->SetAddress(&Ptpt);
    Px->SetAddress(&Pxpt);
    Py->SetAddress(&Pypt);
    Pz->SetAddress(&Pzpt);
    Eta->SetAddress(&Etapt);
    Phi->SetAddress(&Phipt);

    trigger->SetAddress(&triggerInfo);
    outfile->cd();
    for ( int i=0;i<nentries;i++) {

      t2->GetEntry(i); 

      unsigned int muonSize = Ptpt->size();
      if (DEBUG)cout<< __LINE__ << endl;
      for (unsigned int ivec=0; ivec< muonSize;ivec++){
	hPt->Fill(Ptpt->at(ivec));
	if (DEBUG) cout<< __LINE__ << endl;
 
	hEta->Fill(Etapt->at(ivec));
	hPhi->Fill(Phipt->at(ivec));

 	for (unsigned int itr = 0; itr < triggerInfo->size(); itr++){
//	for (unsigned int itr = 0; itr < 100; itr++){   
	string name = triggerInfo->at(itr).first;
// 	  cout << __LINE__ << " trigger "  << triggerInfo->at(itr).first << " " << triggerInfo->at(itr).second << endl;
	  if (triggerInfo->at(itr).second== 1) {
// 	     if (DEBUG) cout<< __LINE__ << endl;
	    (hTrMap[name][0])->Fill(Etapt->at(ivec));
// 	     if (DEBUG) cout<< __LINE__ << endl;
	    hTrMap[name][1]->Fill(Phipt->at(ivec));
// 	     if (DEBUG) cout<< __LINE__ << endl;
	    hTrMap[name][2]->Fill(Ptpt->at(ivec));
	  }
	}
      
        // cout<< __LINE__ << endl;
    
    
//     for (unsigned int ivec=0; ivec< MassSize;ivec++){
// 	for (unsigned int itr = 0; itr < triggerInfo->size(); itr++){
// 	  string name = triggerInfo->at(itr).first;
// 	  hTrMap[name][2]->Fill(EtaMotherpt->at(ivec));
// 	  hTrMap[name][3]->Fill(PhiMotherpt->at(ivec));
// 	}
	//     }
      }

	if (DEBUG) cout<< __LINE__ << endl;
      //   }//if

    }//entries
    if (DEBUG) cout<< __LINE__ << endl;
      
   outfile->cd();

 /////////////////////////////////////////////////////////////////////////////
   TCanvas *canv1= new TCanvas("canv1","firstAnalysis",1000,700);
//    TCanvas *canv3= new TCanvas("canv3","DeltaPhiDeltaEtaR",1000,700);
//    TCanvas *canv4= new TCanvas("canv4","EtaMuVsDeltaPhi",600,800);
   TCanvas *canv7= new TCanvas("canv7","trigger",600,400);
//    TCanvas *canv2= new TCanvas("canv2","MotherAnalysis",1000,700);

   if (DEBUG) cout<< __LINE__ << endl;
   map <string , TCanvas *> canvasMap;
   for (map < string , vector <TH1D*> >::const_iterator it=hTrMap.begin(); it!=hTrMap.end(); it++){
   if (DEBUG) cout<< __LINE__ << endl;
     if (it->second.size()<1) continue ;
   if (DEBUG) cout<< __LINE__ << endl;
     if (it->second.at(1)->GetEntries() < 2*nentries/10) continue;
   if (DEBUG) cout<< __LINE__ << endl;
     string name = it->first;
     if (DEBUG) cout << __LINE__ << ", entries  "<<   it->second.at(1)->GetEntries() << endl;
   if (DEBUG) cout<< __LINE__ << endl;
     stringstream canvasName; 
     canvasName << "c_"<< decay<< "_" << name;
     canvasMap[name] = new TCanvas (canvasName.str().c_str(),canvasName.str().c_str(),1000,700);
   if (DEBUG) cout<< __LINE__ << endl;
     canvasMap[name]->Divide(2,2);
     canvasMap[name]->cd(1); 
     it->second[0]->Draw(); 
     canvasMap[name]->cd(2); 
     it->second[1]->Draw(); 
   if (DEBUG) cout<< __LINE__ << endl;
     canvasMap[name]->cd(3); 
     canvasMap[name]->GetPad(3)->SetLogy();
     it->second[2]->Draw(); 
     it->second[2]->GetXaxis()->SetTitle("GeV/c"); 
//      canvasMap[name]->cd(4); 
//      it->second[3]->Draw(); 
   if (DEBUG) cout<< __LINE__ << endl;
     stringstream canvasFileName;
     canvasFileName << "plots/"<< canvasName.str().c_str()<< ".png" ; 
   if (DEBUG) cout<< __LINE__ << endl;
     if (saveFig ==1){
       if (it->second.at(1)->GetEntries() < 2*nentries/10) continue;
       if (DEBUG) cout<< __LINE__ << endl;
       canvasMap[name]->SaveAs(canvasFileName.str().c_str());
     }
   }
   if (DEBUG) cout<< __LINE__ << endl;
   for (map < string , vector <TH1D*> >::const_iterator it=hTrMap.begin(); it!=hTrMap.end(); it++){
     if (it->second.size()<1) continue ;
     if (it->second.at(1)->GetEntries() < 2*nentries/10) continue;
     string name = it->first;
     canvasMap[name]->Write();
   }
   if (DEBUG) cout<< __LINE__ << endl;

   vector <TCanvas*> cVector;
   cVector.push_back(canv1);
//    cVector.push_back(canv2);
//    cVector.push_back(canv3);
//    cVector.push_back(canv4);
   cVector.push_back(canv7);

   canv1->Draw();
   canv1->Divide(2,2);
      
   canv1->cd(1);
   if (DEBUG) cout<< __LINE__ << endl;
	
   hPt->SetFillStyle(3033);
   hPt->SetFillColor(99);
   hPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
   int maxBin = 0;
   for (int ibin = 1 ; ibin <= hPt->GetNbinsX() ; ibin++){
     if (hPt->GetBinContent(ibin) != 0) {
       maxBin = ibin;
     }
   }
   double maxX = (int )(hPt->GetBinCenter(maxBin)/10)*10  + 11  ;
   hPt->GetXaxis()->SetRangeUser(0,maxX);
   canv1->GetPad(1)->SetLogy();
   hPt->Draw();
   if (DEBUG) cout<< __LINE__ << endl;
	
   canv1->cd(2);
//    hMass->SetFillStyle(3020);
//    hMass->SetFillColor(3);
//    hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
//    hMass->Draw();
   
   canv1->cd(3);
   hPhi->SetFillStyle(3009);
   hPhi->SetFillColor(2);
   hPhi->GetXaxis()->SetTitle("#phi (rad)");
   hPhi->Draw();
   
   canv1->cd(4);
   hEta->SetFillStyle(3020);
   hEta->SetFillColor(4);
   hEta->GetXaxis()->SetTitle("#eta");
   hEta->Draw();
   canv1->Update();
   if (DEBUG) cout<< __LINE__ << endl;

     
//    canv2->Draw();
//    canv2->Divide(2,2);
//    canv2->cd(1);
//    hMass->Draw();
//    canv2->cd(2);
//    hPhiMother->Draw();
//    hPhiMother->GetXaxis()->SetTitle("#phi (rad)");
//    canv2->cd(3);
//    hEtaMother->Draw();
//    hEtaMother->GetXaxis()->SetTitle("#eta");
//    canv2->cd(4);
//    hPtMother->Draw();
//    hPtMother->GetXaxis()->SetTitle("p_{T} (GeV/c)");
//    canv2->Update();
 
//    canv3->Divide(2,2);
   
//    canv3->cd(1);
//    hDeltaPhi->SetFillStyle(3009);
//    hDeltaPhi->SetFillColor(1);
//    hDeltaPhi->GetXaxis()->SetTitle("#phi");
//    hDeltaPhi->Draw();
   
//    canv3->cd(2);
//    hDeltaEta->SetFillStyle(3020);
//    hDeltaEta->SetFillColor(3);
//    hDeltaEta->GetXaxis()->SetTitle("#eta");
//    hDeltaEta->Draw();
   
//    canv3->cd(3);
//    hDeltaR->SetFillStyle(3009);
//    hDeltaR->SetFillColor(2);
//    hDeltaR->GetXaxis()->SetTitle("R");
//    hDeltaR->Draw();

//    canv3->Update();
   
   canv7->Draw();
   canv7->cd();
   hTrigger->Draw();
   canv7->Update();
   
//    canv4->Draw();
//    canv4->Divide(1,2);
//    canv4->cd(1);
//    hEta1vsDeltaPhi->Draw();
//    canv4->cd(2);
//    hEta2vsDeltaPhi->Draw();
//    canv4->Update();
   if (DEBUG) cout<< __LINE__ << endl;



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

//    hTriggCount -> Write();
   if (DEBUG) cout<< __LINE__ << endl;

//    hEta1vsDeltaPhi->Write();
//    hEta2vsDeltaPhi->Write();
//    if (DEBUG) cout<< __LINE__ << endl;
	
//    hMass->Write();
//    hPtMother->Write();
//    hEtaMother->Write();
//    hPhiMother->Write();
   hPt->Write();
   hEta->Write();
   hPhi->Write();

//    hDeltaPhi->Write();
//    hDeltaEta->Write();
//    hDeltaR  ->Write();

   hMultCand->Write();
   hTrigger ->Write();
   // cout<< __LINE__ << endl;

   for (map < string , vector <TH1D*> >::const_iterator it=hTrMap.begin(); it!=hTrMap.end(); it++){
     if (it->second.size()<1) continue ;
     if ((it->second).at(1)->GetEntries() < 2*nentries/10) continue;
     (it->second).at(0)->Write();
     (it->second).at(1)->Write();
     (it->second).at(2)->Write();
//      (it->second).at(3)->Write();
     
   }
   if (DEBUG) cout<< __LINE__ << endl;

   // cout<< __LINE__ << endl;


   for (unsigned icanv = 0 ; icanv < cVector.size(); icanv++){
     cVector[icanv]->~TCanvas();
   }


   if (DEBUG) cout<< __LINE__ << endl;

//    hMass->~TH1D();
   hPt->~TH1D();
   hEta->~TH1D();
   hPhi->~TH1D();

//    hDeltaPhi->~TH1D();
//    hDeltaEta->~TH1D();
//    hDeltaR  ->~TH1D();

   hTrigger ->~TH1D();

   hMultCand->~TH1I();
//    for (map<string, vector <TH1D*> >::iterator it=hTrMap.begin(); it!=hTrMap.end(); it++){
//      (it->second).at(0)->~TH1D();
//      (it->second).at(1)->~TH1D();
//      (it->second).at(2)->~TH1D();
//      (it->second).at(3)->~TH1D();
//    }

   outfile->Close();

   f->Close();
   cout << __LINE__ << " file " << file.str().c_str() <<  " closed" << endl;


}
