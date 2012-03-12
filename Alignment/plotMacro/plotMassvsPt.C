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
vector < TProfile *>  hSlice, hSliceUpsilon; 

using namespace std;
void plotMassvsPt(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gROOT->ForceStyle();
  string file1;
  file1 = "myUpsilonMuMu_noIso.root";


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

  TCanvas * c = new TCanvas ("c", "c", 1200, 800);
  if (DEBUG) cout << __LINE__ << endl;
  c->Draw();
  c->Divide(5,4) ;
//   c->cd(1);
  TDirectory * dir1 = (TDirectory *) f1 ->Get("demo");
  if (DEBUG) cout << __LINE__ << endl;

  TH2D * h1 = (TH2D*)dir1->Get("hMassvsPtHigh");
  if (DEBUG) cout << __LINE__ << endl;
  TH2D * h2 = (TH2D*)dir1->Get("hMassvsPtLow");
  if (DEBUG) cout << __LINE__ << endl;
  TH2D * h3 = (TH2D*)dir1->Get("hMassvsPtUpsilon");
  if (DEBUG) cout << __LINE__ << endl;

  if (DEBUG) cout << __LINE__ << endl;

  stringstream hName, textstream; 
   int ihisto = 0 ; 
  for (int ibin = 4; ibin <=24 ; ibin ++){
    cout << "bin = " << ibin << endl;
    hName.str("");
    if (DEBUG) cout << __LINE__ << endl;
    hName <<  ibin -1 << " < pt < "<< ibin ; 
    if (DEBUG) cout << __LINE__ << endl;
//     hSlice.push_back( h1->ProjectionY(hName.str().c_str(), ibin, ibin, ""));
    hSlice.push_back( h1->ProfileY(hName.str().c_str(), ibin, ibin, ""));
    if (DEBUG) cout << __LINE__ << endl;
    c->cd(ihisto+1);
    if (DEBUG) cout << __LINE__ << endl;
    hSlice[ihisto] -> SetTitle(hName.str().c_str());
    hSlice[ihisto] -> Draw();
    if (DEBUG) cout << __LINE__ << endl;
    
    ihisto++;

  }
  ihisto = 0 ; 
  TCanvas * c3 = new TCanvas ("c3", "c3", 1400, 800);
  if (DEBUG) cout << __LINE__ << endl;
  c3->Draw();
  c3->Divide(7,5) ;
  //  TF1 *fittu = new TF1("fit",fitf,2.5,4,5);
  TF1 *fittu = new TF1("fit",fitf,8.9,9.9,4);
  fittu->SetParName(0,"Yieldgau");
  fittu->SetParName(1,"Meangau");
  fittu->SetParName(2,"Sigmagau");
  fittu->SetParName(3,"p0");
  //  fittu->SetParName(4,"p1");
  //  fittu->SetParName(5,"p2");

  for (int ibin = 1; ibin <=35 ; ibin ++){
    fittu->SetParameter(0,10000.);
    fittu->SetParameter(1,9.4);
    fittu->SetParameter(2,0.02);
    fittu->SetParameter(3,1.);
    //    fittu->SetParameter(4,0.1);
    fittu->SetParLimits(2,0.002,0.2);
    fittu->FixParameter(1,9.4);
    fittu->SetLineColor(2);
    cout << "bin = " << ibin << endl;
    hName.str("");
    if (DEBUG) cout << __LINE__ << endl;
    hName <<  ibin << " < Upsilon pt < "<< ibin+1 ; 
    if (DEBUG) cout << __LINE__ << endl;
//     hSliceUpsilon.push_back( h3->ProjectionY(hName.str().c_str(), ibin, ibin, ""));
    hSliceUpsilon.push_back( h3->ProfileY(hName.str().c_str(), ibin, ibin, ""));
    if (DEBUG) cout << __LINE__ << endl;
    c3->cd(ihisto+1);
    if (DEBUG) cout << __LINE__ << endl;
    hSliceUpsilon[ihisto] -> SetTitle(hName.str().c_str());
    hSliceUpsilon[ihisto] -> Draw();
    hSliceUpsilon[ihisto] -> GetYaxis()->SetRangeUser(0,hSliceUpsilon[ihisto] -> GetBinContent(hSliceUpsilon[ihisto] ->GetMaximumBin())*1.1);
    if (DEBUG) cout << __LINE__ << endl;
    hSliceUpsilon[ihisto] ->Fit("fit","","S",8.9,9.9);
    fittu->ReleaseParameter(1);
    hSliceUpsilon[ihisto] ->Fit("fit","","S",8.9,9.9);
    hSliceUpsilon[ihisto] ->SetLineColor(1);
    hSliceUpsilon[ihisto] ->Draw("SAMES");
    
    

    double bkg = fittu->Integral(8.9,9.9)/ hSliceUpsilon[ihisto] ->GetBinWidth(1) - fittu->GetParameter(0);
    

    cout << "Integral  " << fittu->Integral(8.9,9.9)/ hSliceUpsilon[ihisto] ->GetBinWidth(1) <<", signal  "<<     fittu->GetParameter(0) << ", noise " << bkg << ", S/(S+N) " << fittu->GetParameter(0) / (fittu->GetParameter(0) + bkg) << endl; 
    textstream.str("");
    textstream << " S/S+N = " <<  std::setprecision(2)<< fittu->GetParameter(0) / (fittu->GetParameter(0) + bkg) ;
    
    TText * text = new TText (0.2,0.8, textstream.str().c_str());
    text -> SetTextSize( text ->GetTextSize()*1.5);
    text -> DrawText(8.93, hSliceUpsilon[ihisto]->GetBinContent( hSliceUpsilon[ihisto]-> GetMaximumBin())*0.9,textstream.str().c_str());
    
    ihisto++;

  }
  c3->Update();
//   h1->SetLineWidth(2);
// //   h1->SetLineColor(1);
//   h1->Draw();
//   char bufferE[100];
//   char bufferM[100];
//   char bufferR[100];

//   TPaveStats *ptstatsDATA = new TPaveStats(0.8,0.9,.99,0.99,"brNDC");
//   ptstatsDATA->SetName("stats");
//   ptstatsDATA->SetBorderSize(2);
//   ptstatsDATA->SetFillColor(kWhite);
//   ptstatsDATA->SetTextAlign(12);
//   ptstatsDATA->SetTextColor(1);
//   sprintf(bufferE,"Entries = %1.0f",h1->GetEntries());
//   sprintf(bufferM,"Mean = %1.2f",h1->GetMean());
//   sprintf(bufferR,"RMS = %1.2f",h1->GetRMS());
//   ptstatsDATA->AddText(bufferE);
//   ptstatsDATA->AddText(bufferM);
//   ptstatsDATA->AddText(bufferR);
// //   ptstatsDATA->SetOptStat(10001);
// //   ptstatsDATA->SetOptFit(111);

//   ptstatsDATA->Draw();
//   if (DEBUG) cout << __LINE__ << endl;
//   h2->SetLineWidth(2);
//   h2->SetLineColor(2);
//   h2->Draw("SAMES");
//   TPaveStats * ps2 = new TPaveStats(0.8,0.8,.99,0.89,"brNDC");
//   ps2->SetName("stats");
//   ps2->SetBorderSize(2);
//   ps2->SetFillColor(kWhite);
//   ps2->SetTextAlign(12);
//   ps2->SetTextColor(2);
//   sprintf(bufferE,"Entries = %1.0f",h2->GetEntries());
//   sprintf(bufferM,"Mean = %1.2f",h2->GetMean());
//   sprintf(bufferR,"RMS = %1.2f",h2->GetRMS());
//   ps2->AddText(bufferE);
//   ps2->AddText(bufferM);
//   ps2->AddText(bufferR);
//   ps2->Draw();


//   TLegend * leg = new TLegend (0.15,0.7,0.35,0.85); 
//   leg->AddEntry(h1,"noIso");
//   leg->AddEntry(h2,"RelIso 0.3");
//   leg->AddEntry(h3,"RelIso (Em+Had) 0.3");
//   leg->AddEntry(h4,"RelIso (Em+Had) 0.3 + RelIso (Pt) 0.15");
//   leg->AddEntry(h5,"RelIso (Em+Had) 0.15 + RelIso (Pt) 0.15");
//   leg->AddEntry(h6,"RelIso (Em+Had) 0.15");
//   leg->AddEntry(h7,"RelIso (Pt) 0.15");
//   leg->AddEntry(h8,"RelIso (Had) 0.15 + RelIso (Em) 0.15 + RelIso (Pt) 0.15");
//   leg->AddEntry(h9,"4_4_X cfi");
//   leg -> SetFillColor(0); 
//   leg -> Draw(); 

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
      Double_t poli = par[3];

      Double_t fitval = gaus + poli;
      return fitval;
   }

