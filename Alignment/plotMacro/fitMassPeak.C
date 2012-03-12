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
// //vector < TProfile *>  hSlice, hSliceUpsilon; 
// vector < TH1D *>  hSlice, hSliceUpsilon; 
TH1D * h1;

using namespace std;
void fitMassPeak(){
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
//   file1 = "myUpsilonMuMu_all.root";
//   file1 = "myUpsilonMuMu_05_th3.root";
  file1 = "myUpsilonMuMu_03_pt3.root";

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

//   h1 = (TH1D*)dir1->Get("hInvMassUpsilon_thr2mu");
  h1 = (TH1D*)dir1->Get("hInvMassUpsilon");
  if (DEBUG) cout << __LINE__ << endl;
  h1->Draw();
  //  TF1 *fittu = new TF1("fit",fitf,2.5,4,5);
  TF1 *fittu = new TF1("fit",fitf,8.9,9.9,5);
  if (DEBUG) cout << __LINE__ << endl;
  fittu->SetParName(0,"Yieldgau");
  fittu->SetParName(1,"Meangau");
  fittu->SetParName(2,"Sigmagau");
  fittu->SetParName(3,"p0");
  if (DEBUG) cout << __LINE__ << endl;
  fittu->SetParName(4,"p1");
  //  fittu->SetParName(5,"p2");
  fittu->SetParameter(0,1000.);
  fittu->SetParameter(1,9.4);
  fittu->SetParameter(2,0.02);
  fittu->SetParameter(3,100.);
  if (DEBUG) cout << __LINE__ << endl;
  fittu->SetParameter(4,0.1);
  fittu->SetParLimits(2,0.002,0.2);
  if (DEBUG) cout << __LINE__ << endl;
  fittu->FixParameter(1,9.4);
  if (DEBUG) cout << __LINE__ << endl;
  fittu->SetLineColor(2);
  h1->Fit("fit","R");
  if (DEBUG) cout << __LINE__ << endl;
  fittu->ReleaseParameter(1);
  h1 ->Fit("fit","R","S",8.9,9.9);
  h1->SetLineColor(1);
  if (DEBUG) cout << __LINE__ << endl;
  h1->Draw("SAMES") ;
  
  double bkg = fittu->Integral(8.9,9.9)/ h1 ->GetBinWidth(1) - fittu->GetParameter(0);
    
  stringstream textstream;
  textstream.str("");
  textstream << " S/S+N = " <<  std::setprecision(2)<< fittu->GetParameter(0) / (fittu->GetParameter(0) + bkg) ;
  cout << "Yield "<<  fittu->GetParameter(0) << ", bkg "  << bkg << endl;
  
  TText * text = new TText (0.1,0.8, textstream.str().c_str());
  text -> SetTextSize( text ->GetTextSize()*1.5);
  text -> DrawText(8.1, h1->GetBinContent(h1-> GetMaximumBin())*0.9,textstream.str().c_str());
    


    

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
      Double_t Binwid =h1->GetBinWidth(0);
      Double_t arg = 0;
      arg = -0.5*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]);
      if( arg < -70. ) arg = -70.;
      Double_t gaus = Binwid*par[0] / (par[2]*Rt2pi)*TMath::Exp(arg);
//       Double_t poli = par[3] + par[4]*(x[0]-par[1]) + par[5]*
//       (x[0]-par[1])*(x[0]-par[1]); 
//       poli += par[6]*(x[0]-par[1])*(x[0]-par[1])*(x[0]-par[1]) +
//        par[7]*(x[0]-par[1])*(x[0]-par[1])*(x[0]-par[1])*(x[0]-par[1]);
      Double_t poli = par[3]+ par[4]*(x[0]-par[1]) ;

      Double_t fitval = gaus + poli;
      return fitval;
   }

