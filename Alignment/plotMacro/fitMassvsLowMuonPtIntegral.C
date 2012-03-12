#include <memory>

#include "Riostream.h" 
#include <math.h>
#include <sstream>
#include <vector>
#include <iomanip>
#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TPaveStats.h>
#include <TMath.h>

#define DEBUG 1
double fitf(Double_t *x, Double_t *par);
// vector < TProfile *>  hSlice, hSliceUpsilon; 
vector < TH1D *>  hSlice, hSliceUpsilon; 

using namespace std;
void fitMassvsLowMuonPtIntegral(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gROOT->ForceStyle();
  string file1;
//   file1 = "myUpsilonMuMu_noIso.root";
//   file1 = "myUpsilonMuMu_03.root";
//   file1 = "myUpsilonMuMu_05.root";
   file1 = "myUpsilonMuMu_all.root";

  cout << __LINE__ << " Opening file: " << file1 << endl;
  if (DEBUG) cout << __LINE__ << endl;
  TFile * f1 = new TFile();  
  f1 = TFile::Open(file1.c_str());
  if (!f1){
    cout << "FILE NOT FOUND" << endl;
    return;
  }

  if (DEBUG) cout << __LINE__ << endl;

//   stringstream outstring;
//   outstring.str("");
//   outstring << "plotMyMass.root";
//   if (DEBUG) cout << __LINE__ << endl;
  
//   TFile* outfile = new TFile(outstring.str().c_str(),"RECREATE");
//   outfile -> cd(); 

  TDirectory * dir1 = (TDirectory *) f1 ->Get("demo");
  if (DEBUG) cout << __LINE__ << endl;

  TH2D * h1 = (TH2D*)dir1->Get("hMassvsPtHigh");
  if (DEBUG) cout << __LINE__ << endl;
  TH2D * h2 = (TH2D*)dir1->Get("hMassvsPtLow");
  if (DEBUG) cout << __LINE__ << endl;
  TH2D * h3 = (TH2D*)dir1->Get("hMassvsPtUpsilon");
  if (DEBUG) cout << __LINE__ << endl;
  TH1D* hYield= new TH1D ("hYield","hYield",50,0,50);
  TH1D* hMean = new TH1D ("hMean","hMean",50,0,50);
  TH1D* hSigma= new TH1D ("hSigma","hSigma",50,0,50);
  TH1D* hBackground = new TH1D ("hBackground","hBackground",50,0,50);
  TH1D* hSnratio = new TH1D ("hSnratio","hSnratio",50,0,50);


  if (DEBUG) cout << __LINE__ << endl;

  stringstream hName, textstream; 
  int ihisto = 0 ; 
  int entot = 0;  
  TCanvas * c3 = new TCanvas ("c3", "c3", 1400, 800);
//   TCanvas * c3 = new TCanvas ("c3", "c3", 1200, 700);
  if (DEBUG) cout << __LINE__ << endl;
  c3->Draw();
  c3->Divide(5,4) ;
  //  TF1 *fittu = new TF1("fit",fitf,2.5,4,5);
  TF1 *fittu = new TF1("fit",fitf,8.9,9.9,4);
  fittu->SetParName(0,"Yieldgau");
  fittu->SetParName(1,"Meangau");
  fittu->SetParName(2,"Sigmagau");
  fittu->SetParName(3,"p0");
  //  fittu->SetParName(4,"p1");
  //  fittu->SetParName(5,"p2");

  for (int ibin =1 ; ibin <=20 ; ibin ++){
    fittu->SetParameter(0,10000.);
    fittu->SetParameter(1,9.4);
    fittu->SetParameter(2,0.02);
    fittu->SetParameter(3,100.);
    //    fittu->SetParameter(4,0.1);
    fittu->SetParLimits(2,0.002,0.2);
    fittu->FixParameter(1,9.4);
    fittu->SetLineColor(2);
    cout << "bin = " << ibin << endl;
    hName.str("");
    if (DEBUG) cout << __LINE__ << endl;
    hName << "Upsilon pt >" <<  ibin-1 ; 
    if (DEBUG) cout << __LINE__ << endl;
    hSliceUpsilon.push_back( h2->ProjectionY(hName.str().c_str(), ibin, h3->GetYaxis()->GetLast(), ""));
//     hSliceUpsilon.push_back( h1->ProfileY(hName.str().c_str(), ibin, h3->GetYaxis()->GetLast(), ""));
    if (DEBUG) cout << __LINE__ << endl;
    c3->cd(ihisto+1);
    if (DEBUG) cout << __LINE__ << endl;
    hSliceUpsilon[ihisto] -> SetTitle(hName.str().c_str());
    hSliceUpsilon[ihisto] -> Draw();
    hSliceUpsilon[ihisto] -> GetYaxis()->SetRangeUser(0,hSliceUpsilon[ihisto] -> GetBinContent(hSliceUpsilon[ihisto] ->GetMaximumBin())*1.1);
    if (DEBUG) cout << __LINE__ << endl;
    hSliceUpsilon[ihisto] ->Fit("fit","","",8.9,9.9);
    fittu->ReleaseParameter(1);
    hSliceUpsilon[ihisto] ->Fit("fit","","",8.9,9.9);
    hSliceUpsilon[ihisto] ->SetLineColor(1);
    hSliceUpsilon[ihisto] ->Draw("SAMES");
    
    

    double bkg = fittu->Integral(8.9,9.9)/ hSliceUpsilon[ihisto] ->GetBinWidth(1) - fittu->GetParameter(0);
    

    cout << "Integral  " << fittu->Integral(8.9,9.9)/ hSliceUpsilon[ihisto] ->GetBinWidth(1) <<", signal  "<<     fittu->GetParameter(0) << ", noise " << bkg << ", S/(S+N) " << fittu->GetParameter(0) / (fittu->GetParameter(0) + bkg) << endl; 
    textstream.str("");
    textstream << " S/S+N = " <<  std::setprecision(2)<< fittu->GetParameter(0) / (fittu->GetParameter(0) + bkg) ;
    
    TText * text = new TText (0.2,0.8, textstream.str().c_str());
    text -> SetTextSize( text ->GetTextSize()*1.5);
    text -> DrawText(8.93, hSliceUpsilon[ihisto]->GetBinContent( hSliceUpsilon[ihisto]-> GetMaximumBin())*0.9,textstream.str().c_str());
    
    hYield->Fill(ibin-0.5,fittu->GetParameter(0) );
    hMean->Fill(ibin-0.5,fittu->GetParameter(1) );
    hSigma->Fill(ibin-0.5,fittu->GetParameter(2) );
    hBackground -> Fill(ibin-0.5, bkg);
    hSnratio->Fill (ibin-0.5, fittu->GetParameter(0) / (fittu->GetParameter(0) + bkg));

    entot+=hSliceUpsilon[ihisto] ->Integral (); 
    cout << __LINE__ <<  " " << h3->Integral() <<  " " << hSliceUpsilon[ihisto] ->Integral () << " " << entot << endl;

    ihisto++;

  }
  c3->Update();
  c3->SaveAs("AllFit_PtLowMuIntegral.png");
  TCanvas* c2= new TCanvas("c2", "c2", 1400,800) ;
//   TCanvas* c2= new TCanvas("c2", "c2", 1200,700) ;
  c2->Draw(); 
  c2->Divide (3,2); 
  c2->cd(1); 
  hYield->Draw(); 
  hYield->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  c2->cd(2);
  hMean->Draw();
  hMean->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  hMean->GetYaxis()->SetRangeUser(9.3,9.6);
  c2->cd(3);
  hSigma->Draw(); 
  hSigma->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  c2->cd(4);
  hBackground ->Draw();
  hBackground ->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");
  c2->cd(5);
  hSnratio->Draw(); 
  hSnratio->GetXaxis()->SetTitle("p_{T}(#mu) (GeV)");

  c2->SaveAs("FitResults_PtLowMuIntegral.png");
  //  c3->SaveAs("higherMuonPt.png");
//   c->Write();
//   h1->Write();
//   h2->Write();
//   h3->Write();
//   h4->Write();
//   h5->Write();

//   outfile->Close();
//  f1->Close();

}
Double_t fitf(Double_t *x, Double_t *par)
   { 
//    Gaussian plus a polynomial
//     
//    Sqrt(2*pi)
      Double_t Rt2pi = 2.506628274922;
//            
      Double_t Binwid =hSliceUpsilon[0]->GetBinWidth(0);
      Double_t arg = 0;
      arg = -0.5*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]);
      if( arg < -70. ) arg = -70.;
      Double_t gaus = Binwid*par[0] / (par[2]*Rt2pi)*TMath::Exp(arg);
//       Double_t poli = par[3] + par[4]*(x[0]-par[1]) + par[5]*
//       (x[0]-par[1])*(x[0]-par[1]); 
//       poli += par[6]*(x[0]-par[1])*(x[0]-par[1])*(x[0]-par[1]) +
//        par[7]*(x[0]-par[1])*(x[0]-par[1])*(x[0]-par[1])*(x[0]-par[1]);
      Double_t poli = par[3] ;

      Double_t fitval = gaus + poli;
      return fitval;
   }

