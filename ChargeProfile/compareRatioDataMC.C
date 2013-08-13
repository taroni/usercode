#include <memory>

#include "Riostream.h" 
#include <math.h>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <stdio.h>
#include <cstdio>

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TKey.h>

#define DEBUG 1

using namespace std;
void compareRatioDataMC(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  gStyle->SetPadGridY(1); 
  gStyle->SetPadGridX(1); 

  stringstream name; 
  int bias = 150; 
//   TFile * file1 = new TFile("chargeRatio_lorentzangleALCARECO.root", "READ"); 
//   TFile * file2 = new TFile("chargeRatio_lorentzangleALCARECOMC.root", "READ");
   TFile * file1 = new TFile("chargeRatio_driftdepthRun208391.root", "READ"); 
  TFile * file2 = new TFile("chargeRatio_driftdepth_MC.root", "READ"); 
 
  TCanvas * c = new TCanvas() ; 
  c->Draw() ; 
  if (DEBUG) cout <<__LINE__ << endl;
  TLegend * leg = new TLegend(0.35, 0.15, 0.65, 0.25);
  leg->SetFillColor(0); 
  for (int ilayer = 1; ilayer < 4; ilayer++){
    if (DEBUG) cout <<__LINE__ << endl;
    name.str(""); 
    name << "h_ratiovsmodule_layer"<<ilayer << "_bias" << bias; 
    if (DEBUG) cout <<__LINE__ <<" "<< name.str() << endl;  
    TH1F * h1 = (TH1F*) file1->Get(name.str().c_str());
    TH1F * h2 = (TH1F*) file2->Get(name.str().c_str());
    if (DEBUG) cout <<__LINE__ << endl;  
    c->cd() ; 
    h1->Draw("E") ;
    if (DEBUG) cout <<__LINE__ << endl;  
    h1->SetMarkerStyle(20); 
    if (DEBUG) cout <<__LINE__ << endl;  
    h2->Draw("ESAME"); 
    h2->SetMarkerColor(2); 
    h2->SetMarkerStyle(20); 
    if (DEBUG) cout <<__LINE__ << endl;
    if (ilayer ==1) {
      leg->AddEntry(h1, "Data, Run 208391", "lpf");
      leg->AddEntry(h2, "Simulation", "lpf");
    }
    leg->Draw(); 
    TPaveText *pt = new TPaveText(.66,.53,1.6,.8);
    pt->SetFillColor(0); 
    pt->AddText("1=ring8/ring1");
    pt->AddText("2=ring7/ring2");
    pt->AddText("3=ring6/ring3");
    pt->AddText("4=ring5/ring4");
      
    pt->Draw("SAME");

    name << ".png" ;
    c->SaveAs(name.str().c_str()); 
     
      
        

  }//ilayer
 
  file1->Close();
  file2->Close(); 
}
