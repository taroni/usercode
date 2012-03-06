//g++ -o  gaussYieldBin gaussYieldBin.cc `root-config --cflags --libs` -lRooFit
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <TBranch.h>
#include <TBasket.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFrame.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TEventList.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include "TPostScript.h"

void doFit(void);
double fitf(Double_t *x, Double_t *par);
TH1F *histo;
TH1F *histo2;
// TH1F * lifetime;
TCanvas * canvas;
TCanvas * canvas2;
TCanvas * canvas3;

// using namespace std; 
// using namespace RooFit; 

int main(int argc, char**argv)
{
  TApplication app("App",&argc, argv);
	
  doFit() ;
	
  app.Run (); 
 return 0;
}
//===========================//
void doFit(void)
{
  gSystem->Load("libRooFit");
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1111); 
  gStyle->SetOptStat(11); 
//   TFile * file = new TFile("fullAnalysisSignalDouble_0_10-ls3-3.root","READ") ;
  ofstream outFiletxt;
  outFiletxt.open("data.dat", ofstream::out);

  TFile * file = new TFile("fullAnalysisSignal_Mass08_0_10-ls1.root","READ") ;

  //  TFile * file = new TFile("fullAnalysisSignal_200pbDouble_0.root","READ") ;
//   TFile * file = new TFile("fullAnalysisSignal_0_10-ls3-1fb.root","READ") ;
//  TFile * file = new TFile("fullAnalysisSignalDouble_0_10-ls3-1fb.root","READ") ;
//   TFile * file = new TFile("fullAnalysisSignalDouble_0_10-ls3.root","READ");
  file->cd();
  TFile * outfile = new TFile("gaussYieldBin.root","RECREATE") ;
  TH1F * h = new TH1F ("BcInvariantMass","B_{c} Invariant Mass",20,6.1,6.5);
  
  TTree * tree = (TTree*) file->Get("MassBc2D");
//   tree->SetBranchStatus("*",1);
  tree ->Print();
  TBranch * Mass = tree -> GetBranch("Mass2D");
  double massAdd=0;
  Mass->SetAddress(&massAdd);
  
  for (int i=0; i<tree -> GetEntries(); i++){
      Mass->GetEntry(i);
      cout << massAdd << endl;
      h->Fill(massAdd);
      outFiletxt<< massAdd<< endl;
 }
  //  h->Fill();
//   tree->Draw("Mass2D>>hMass(15,6.08,6.52)");
  
//   TH1F * h1 = (TH1F*)gDirectory->Get("hMass");
//   TH1F *h2 = (TH1F*)h1->Clone();
 
// //    TH1F * h = (TH1F*)file->Get("MassBcCosLxyCovCut");
// // //  TH1F * h = (TH1F*)file->Get("MassB_c");
// //   TH1F * h2 = (TH1F*)file->Get("properTime2DCovH");
  
  histo =(TH1F*) h->Clone();
//   histo2=(TH1F*) h2->Clone();
  outfile->cd();
  file->Close();
  TH1F * hB = new TH1F ("hB", "hB", 200,3.99,7.99);
  int  iB=0;
  double iMeV =0;
  double mGeV =0;
  while (iB<582){
//     int seed = time (NULL);
//     int seed = 250051982;
    int seed = 20051982;
    if (iB!=0)  seed = (iMeV*mGeV+1000*iB);
    srand ( seed );
    
    iMeV = rand() % 600 + 6000;
//     cout<< iMeV << endl;
    mGeV = iMeV/1000.;
    if (mGeV>6.0 && mGeV <6.6){
      iB++;
      histo->Fill(mGeV);
      outFiletxt<< mGeV << endl;
//       h->Fill(mGeV);
      hB->Fill(mGeV);
//       cout << iB << " " << mGeV<< endl;
    }
  }
//   h2->Draw("Same");
//   h2->SetLineColor(2);

  canvas  = new TCanvas("c1","c1",700,500);
  histo->Draw();
//   histo->GetXaxis()->SetRangeUser(6.,6.6);
  canvas->Update();

  TF1 *fittu = new TF1("fit",fitf,6.1,6.5,4);
  fittu->SetParName(0,"Yieldgau");
  fittu->SetParName(1,"Meangau");
  fittu->SetParName(2,"Sigmagau");
  fittu->SetParName(3,"p0");

  fittu->SetParameter(0,320.);
  fittu->SetParameter(1,6.286);
  fittu->SetParameter(2,0.03);
  fittu->SetParameter(3,10.);

  
  fittu->SetParLimits(2,0.002,0.2);
  fittu->SetParLimits(3,0.001,50.);
  fittu->FixParameter(1,6.286);
   histo->Fit("fit","RL");
  fittu->ReleaseParameter(1);
  histo->Fit("fit","RL");
  histo->SetLineColor(1);
  fittu->SetLineColor(2);


  histo->GetXaxis()->SetRangeUser(6.1,6.5);
  histo->SetTitle("B_{c} invariant mass");
  histo->SetName("B_{c} invariant mass");
  histo->GetXaxis()->SetTitle("GeV");
  canvas->Update();
  canvas->SaveAs("Mass.eps");
  outFiletxt.close();

// // //-------------------------------------------------------
// // //endProgram

//   canvas->Write();
//   histo->Write();
//   hB->Write();
//   outfile->Close();

}
Double_t fitf(Double_t *x, Double_t *par)
   { 
//    Gaussian plus a polynomial
//     
//    Sqrt(2*pi)
      Double_t Rt2pi = 2.506628274922;
//            
      Double_t Binwid = histo->GetBinWidth(0);
      Double_t arg = 0;
      arg = -0.5*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]);
      if( arg < -70. ) arg = -70.;
      Double_t gaus = Binwid*par[0] / (par[2]*Rt2pi)*TMath::Exp(arg);
      Double_t poli = par[3] + par[4]*(x[0]-par[1]) + par[5]*
      (x[0]-par[1])*(x[0]-par[1]); 
      poli += par[6]*(x[0]-par[1])*(x[0]-par[1])*(x[0]-par[1]) +
       par[7]*(x[0]-par[1])*(x[0]-par[1])*(x[0]-par[1])*(x[0]-par[1]);

      Double_t fitval = gaus + poli;
      return fitval;
   }

