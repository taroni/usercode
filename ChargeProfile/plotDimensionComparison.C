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
void plotDimensionComparison(string file1, string file2){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  gStyle->SetPadGridY(1); 
  gStyle->SetPadGridX(1); 

  stringstream name; 
  int bias = 150; 


  if (DEBUG) cout << __LINE__ << endl;
  
  TFile *_file0 = TFile::Open(file1.c_str());
  TFile *_file1 = TFile::Open(file2.c_str());
  TCanvas * c = new TCanvas(); 
  c->Draw(); 
  c->SetLogy(1);
  if (DEBUG) cout << __LINE__ << endl;
  for (int ilayer=1; ilayer < 4; ilayer++){
    for (int imodule =1; imodule<9; imodule++){
      if (DEBUG) cout << __LINE__ << endl;
      name.str("");
      name << "clusterSize_layer" << ilayer<<"_module"<< imodule << "_bias150"; 
      cout << __LINE__ << " " << name.str().c_str()<< endl;
      TH1F * h = (TH1F *) _file0->Get(name.str().c_str());
      TH1F * h1 = (TH1F *) _file1->Get(name.str().c_str());
      if (DEBUG) cout << __LINE__ << endl;
      c->cd(); 
      h->Sumw2();
      h1->Sumw2();
      h->Scale(1./h->Integral());
      h1->Scale(1./h1->Integral());
      h->Draw("E");
      h->GetXaxis()->SetRangeUser(0,50);
      h->Draw("EP");
      h->SetMarkerStyle(20);
      h->Draw("EP");
      h1->SetMarkerStyle(20);
      h->Draw("EPSAME");
      h1->Draw("EPSAME");
      h1->SetMarkerColor(2);
      h1->Draw("EPSAME");

      name << "comparison.png";
      c->SaveAs(name.str().c_str());
      h->~TH1F();
      h1->~TH1F(); 

      name.str("");
      name << "clusterSizeX_layer" << ilayer<<"_module"<< imodule << "_bias150"; 
      cout << __LINE__ << " " << name.str().c_str()<< endl;
      TH1F * hx = (TH1F *) _file0->Get(name.str().c_str());
      TH1F * hx1 = (TH1F *) _file1->Get(name.str().c_str());
      if (DEBUG) cout << __LINE__ << endl;
      c->cd(); 
      hx->Sumw2();
      hx1->Sumw2();
      hx->Scale(1./hx->Integral());
      hx1->Scale(1./hx1->Integral());
      hx->Draw("E");
      hx->GetXaxis()->SetRangeUser(0,50);
      hx->Draw("EP");
      hx->SetMarkerStyle(20);
      hx->Draw("EP");
      hx1->SetMarkerStyle(20);
      hx->Draw("EPSAME");
      hx1->Draw("EPSAME");
      hx1->SetMarkerColor(2);
      hx1->Draw("EPSAME");

      name << "comparison.png";
      c->SaveAs(name.str().c_str());
      hx->~TH1F();
      hx1->~TH1F(); 

      name.str("");
      name << "clusterSizeY_layer" << ilayer<<"_module"<< imodule << "_bias150"; 
      cout << __LINE__ << " " << name.str().c_str()<< endl;
      TH1F * hy = (TH1F *) _file0->Get(name.str().c_str());
      TH1F * hy1 = (TH1F *) _file1->Get(name.str().c_str());
      if (DEBUG) cout << __LINE__ << endl;
      c->cd(); 
      hy->Sumw2();
      hy1->Sumw2();
      hy->Scale(1./hy->Integral());
      hy1->Scale(1./hy1->Integral());
      hy->Draw("E");
      hy->GetXaxis()->SetRangeUser(0,50);
      hy->Draw("EP");
      hy->SetMarkerStyle(20);
      hy->Draw("EP");
      hy1->SetMarkerStyle(20);
      hy->Draw("EPSAME");
      hy1->Draw("EPSAME");
      hy1->SetMarkerColor(2);
      hy1->Draw("EPSAME");

      name << "comparison.png";
      c->SaveAs(name.str().c_str());
      hy->~TH1F();
      hy1->~TH1F(); 

    }//imodule
  }//ilayer

}
