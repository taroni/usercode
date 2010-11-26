#include "Riostream.h" 
#include <memory>
#include <math.h>
#include <sstream>

#include <TBranch.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

void computeKinCorrection (){

  gROOT->SetStyle("Plain");

  TFile * wFile = new TFile();  
  TFile * zFile = new TFile();  
  wFile = TFile::Open("KinCorrection_wjets_madgraph_vols.root");
  if (!wFile){
    cout << "FILE NOT FOUND" << endl;
    return;
  }
  zFile = TFile::Open("KinCorrection_zjets_madgraph_vols.root");
  if (!zFile){
    cout << "FILE NOT FOUND" << endl;
    return;
  }

  TDirectory * dirPlus = (TDirectory *) wFile->FindObjectAny("KinCorrectionPlus_WP80");
  if (!dirPlus) {
    cout << "directory plus not found" << endl;
    return;
  }
  TDirectory * dirMinus = (TDirectory *) wFile->FindObjectAny("KinCorrectionMinus_WP80");
  if (!dirMinus) {
    cout << "directory minus not found" << endl;
    return;
  }

  TH3F * wh3Plus =  (TH3F *) dirPlus->FindObjectAny("PTvsYvsEta");
  if(!wh3Plus){
    cout << "histo not found" << endl;
    return;
  }
  TH3F * wh3Minus =  (TH3F *) dirMinus->FindObjectAny("PTvsYvsEta");
  if(!wh3Minus){
    cout << "histo not found" << endl;
    return;
  }
  TDirectory * dirZPlus = (TDirectory *) zFile->FindObjectAny("KinCorrectionPlus_WP80");
  if (!dirZPlus) {
    cout << "directory plus not found" << endl;
    return;
  }
  TDirectory * dirZMinus = (TDirectory *) zFile->FindObjectAny("KinCorrectionMinus_WP80");
  if (!dirZMinus) {
    cout << "directory minus not found" << endl;
    return;
  }
  
  TH3F * zh3Plus =  (TH3F *) dirZPlus->FindObjectAny("PTvsYvsEta");


  if(!zh3Plus){
    cout << "histo not found" << endl;
    return;
  }
  TH3F * zh3Minus =  (TH3F *) dirZMinus->FindObjectAny("PTvsYvsEta");


  if(!zh3Minus){
    cout << "histo not found" << endl;
    return;
  }


  TH3F * numPlus =(TH3F *) wh3Plus->Clone();
  TH3F * numMinus =(TH3F *) wh3Minus->Clone();
  TH3F * denPlus =(TH3F *) zh3Plus->Clone();
  TH3F * denMinus=(TH3F *) zh3Minus->Clone();
  
  cout << "wPlus integral " <<  numPlus->GetEntries() << "wMinus integral " <<  numMinus->GetEntries() << " z Integral (+) " << denPlus->GetEntries() <<" z Integral (-) " << denMinus->GetEntries() << endl;

  numPlus ->Scale(1./numPlus ->GetEntries());
  numMinus->Scale(1./numMinus->GetEntries());
  denPlus->Scale(1./denPlus->GetEntries());
  denMinus->Scale(1./denMinus->GetEntries());

  numPlus->Divide(denPlus);
  numPlus->SetName("KinCorrWPlus");
  numMinus->Divide(denMinus);
  numMinus->SetName("KinCorrWMinus");

  TFile *OutFile  = new TFile("KinematicReWeight.root","RECREATE");
  OutFile -> cd();
  
  numPlus->Write();
  numMinus->Write();


  wFile->Close();
  zFile->Close();
  OutFile->Close();

}
