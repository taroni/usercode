#include "Riostream.h" 
#include <memory>
#include <math.h>
#include <sstream>

#include <TBranch.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TKey.h>
#include "THashList.h"

using namespace std;
void plotOutFile(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  vector <TH1F*> vHisto;

  TCanvas * c1= new TCanvas ("c1","c1", 1200, 800);
  c1->Draw();
  c1->SetGridy();
  c1->cd();
  TLegend * legend = new TLegend(0.6,0.55,0.85,0.8, "","brNDC");

  TFile * f = new TFile();
    f = TFile::Open("outfile.root");
  if (!f){
    cout << "FILE NOT FOUND" << endl;
    return;
  }
  TIter hIt(f->GetListOfKeys());
  TKey *keyHisto;
  vector <string> hNames;
  while ((keyHisto=(TKey*)hIt())) {
    hNames.push_back(keyHisto->GetName());
  }
  for(unsigned iHisto=0;iHisto<hNames.size();iHisto++){
    TH1F *h = (TH1F*) f->FindObjectAny(hNames[iHisto].c_str());
    vHisto.push_back(h);
    if (iHisto==0) {
      h->Draw();
      h->GetYaxis()->SetTitle("evt/LS");
      h->GetXaxis()->LabelsOption("v");
    }
    if (iHisto!=0) h->Draw("SAME");
    // h->SetMarkerStyle(20+iHisto);
    // h->SetMarkerSize(2);
    // h->SetMarkerColor(iHisto+1);
    // h->ShowMarker();
    h->SetLineColor(iHisto+1);
    h->SetLineWidth(2);
    legend ->AddEntry(h);

  }
  legend -> Draw();
  legend->SetFillColor(0);  

  c1->SaveAs("canv.root");
  
}
