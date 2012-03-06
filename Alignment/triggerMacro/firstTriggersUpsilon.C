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
#define DEBUG 0

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

void firstTriggersUpsilon(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gStyle->SetHistLineWidth(2);
  //  gStyle->SetOptStat(111110);
  //gStyle->SetOptFit(111);

  gROOT->ForceStyle();
 //  double pigreco = 3.141592;

  
  TFile * f = new TFile();
  
  f = TFile::Open("triggerFileUpsilonMuMu.root");
  if (!f){
    cout << "FILE NOT FOUND" << endl;
    return;
  }
//   cout << __LINE__ << " file " << file.str().c_str() <<  " opened" << endl;
//    stringstream outstring;
//    outstring.str("");
//    outstring << "triggerFile" << decay << ".root";

//    TFile* outfile = new TFile(outstring.str().c_str(),"RECREATE");

//   outfile -> cd();


  TH1D * h1_phi = (TH1D *) f->Get("HLT_Dimuon0_Barrel_Upsilon_v1TrackPhi");
  TH1D * h1_eta = (TH1D *) f->Get("HLT_Dimuon0_Barrel_Upsilon_v1TrackEta");
  TH1D * h1_pt  = (TH1D *) f->Get("HLT_Dimuon0_Barrel_Upsilon_v1TrackPt");

  TH1D * h2_phi = (TH1D *) f->Get("AlCa_RPCMuonNormalisation_v3TrackPhi");
  TH1D * h2_eta = (TH1D *) f->Get("AlCa_RPCMuonNormalisation_v3TrackEta");
  TH1D * h2_pt  = (TH1D *) f->Get("AlCa_RPCMuonNormalisation_v3TrackPt");

  TH1D * h3_phi = (TH1D *) f->Get("HLT_Physics_NanoDST_v1TrackPhi");
  TH1D * h3_eta = (TH1D *) f->Get("HLT_Physics_NanoDST_v1TrackEta");
  TH1D * h3_pt  = (TH1D *) f->Get("HLT_Physics_NanoDST_v1TrackPt");

  TH1D * h4_phi = (TH1D *) f->Get("HLT_DoubleMu3_Quarkonium_v2TrackPhi");
  TH1D * h4_eta = (TH1D *) f->Get("HLT_DoubleMu3_Quarkonium_v2TrackEta");
  TH1D * h4_pt  = (TH1D *) f->Get("HLT_DoubleMu3_Quarkonium_v2TrackPt");

  TH1D * h5_phi = (TH1D *) f->Get("HLT_DoubleMu3_Upsilon_v1TrackPhi");
  TH1D * h5_eta = (TH1D *) f->Get("HLT_DoubleMu3_Upsilon_v1TrackEta");
  TH1D * h5_pt  = (TH1D *) f->Get("HLT_DoubleMu3_Upsilon_v1TrackPt");


  TCanvas *c = new TCanvas ("c", "c", 1200, 800);
  if (DEBUG) cout << __LINE__ << endl;
  c->Divide(2,2);
  c->cd(1);
//   c->GetPad(1)->SetLogy();
  h1_phi ->Draw();
  h2_phi ->SetLineColor(2);
  h2_phi ->Draw("SAME");
  h1_phi ->GetXaxis()->SetTitle("#phi");
  if (DEBUG) cout << __LINE__ << endl;
  if (DEBUG) cout << __LINE__ << endl;
  h3_phi ->SetLineColor(3);
  h3_phi ->Draw("SAME");
  h4_phi ->SetLineColor(4);
  h4_phi ->Draw("SAME");
  if (DEBUG) cout << __LINE__ << endl;
  h5_phi ->SetLineColor(6);
  h5_phi ->Draw("SAME");
  if (DEBUG) cout << __LINE__ << endl;

  c->cd(2);
//   c->GetPad(2)->SetLogy();
  h1_eta ->Draw();
  h2_eta ->SetLineColor(2);
  if (DEBUG) cout << __LINE__ << endl;
  h2_eta ->Draw("SAME");
  h1_eta ->GetXaxis()->SetTitle("#eta");
  h3_eta ->SetLineColor(3);
  h3_eta ->Draw("SAME");
  if (DEBUG) cout << __LINE__ << endl;
  h4_eta ->SetLineColor(4);
  h4_eta ->Draw("SAME");
  h5_eta ->SetLineColor(6);
  h5_eta ->Draw("SAME");
  if (DEBUG) cout << __LINE__ << endl;

  c->cd(3);
  c->GetPad(3)->SetLogy();
  h1_pt ->Draw();
  h2_pt ->SetLineColor(2);
  h2_pt ->Draw("SAME");
  h1_pt ->GetXaxis()->SetTitle("p_{T} GeV/c");
  if (DEBUG) cout << __LINE__ << endl;
  h3_pt ->SetLineColor(3);
  h3_pt ->Draw("SAME");
  h4_pt ->SetLineColor(4);
  h4_pt ->Draw("SAME");
  if (DEBUG) cout << __LINE__ << endl;
  h5_pt ->SetLineColor(6);
  h5_pt ->Draw("SAME");

  if (DEBUG) cout << __LINE__ << endl;
  c->Update();
  c->SaveAs("firstTriggerUpsilon.png");
  c->SaveAs("firstTriggerUpsilon.root");

//   outfile->Close();

  f->Close();
//   cout << __LINE__ << " file " << file.str().c_str() <<  " closed" << endl;


}
