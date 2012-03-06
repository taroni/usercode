//g++ -o  gaussYield gaussYield.cc `root-config --cflags --libs` -lRooFit
#ifndef __CINT__
//#include "RooGlobalFunc.h"
#endif
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFrame.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
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
//using namespace RooFit; 

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
  gStyle->SetOptFit(1111111); 
 
//  TFile * file = new TFile("outfileJpsiMuMuNewMass2.root","READ") ;
  TFile * file = new TFile("outfileJpsiMuMuNewMass_noGlbMuSel2.root","READ") ;

  file->cd();
  TFile * outfile = new TFile("gaussYieldJpsi.root","RECREATE") ;
  

   TH1F * h = (TH1F*)file->Get("hMass");

  histo =(TH1F*) h->Clone();

  outfile->cd();

  canvas  = new TCanvas("c","c",700,500);
  histo->Draw();
  histo->GetXaxis()->SetRangeUser(2.,4.2);
  canvas->Update();

  TF1 *fittu = new TF1("fit",fitf,2.5,4,5);
  fittu->SetParName(0,"Yieldgau");
  fittu->SetParName(1,"Meangau");
  fittu->SetParName(2,"Sigmagau");
  fittu->SetParName(3,"p0");
  fittu->SetParName(4,"p1");
  fittu->SetParName(5,"p2");
  fittu->SetParameter(0,16000.);
  fittu->SetParameter(1,3.1);
  fittu->SetParameter(2,0.06);
  fittu->SetParameter(3,50000.);
  fittu->SetParameter(4,-15000.);
  
  fittu->SetParLimits(2,0.002,0.2);
  fittu->FixParameter(1,3.1);
  histo->Fit("fit","","S",2.7,3.4);
  fittu->ReleaseParameter(1);
  histo->Fit("fit","","S",2.7,3.4);
  histo->SetLineColor(1);
//  fittu->SetLineColor(2);
//  fittu->Draw();
  
//  histo->GetXaxis()->SetRangeUser(2,4.2);
  histo->SetTitle("Jpsi invariant mass");
  histo->SetName("JpsiMass");
  histo->GetXaxis()->SetTitle("GeV");
  canvas->Update();
  canvas->SaveAs("Mass.eps");


//-------------------------------------------------------
//endProgram

  canvas->Write();
  histo->Write();

  file->Close();
   outfile->Close();
   std::cout << "outfile closed" << std::endl;

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

