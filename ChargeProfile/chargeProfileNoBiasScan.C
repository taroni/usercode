#include <memory>

#include "Riostream.h" 
#include <math.h>
#include <sstream>
#include <vector>
#include <string>

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TKey.h>

#define DEBUG 1

using namespace std;
void chargeProfileNoBiasScan(string filename="lorentzangleALCARECO.root",  bool saveFig= 0){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gROOT->ForceStyle();
  gStyle->SetPadGridY(1); 
  gStyle->SetPadGridX(1); 
  

  double pigreco = 3.141592;
  
  TFile * file = new TFile(filename.c_str(), "READ"); 
  stringstream outname;
  outname.str(""); 
  outname<< "chargeProfile_" << filename; 
  TFile *  outfile = new TFile(outname.str().c_str(), "RECREATE"); 
  
  vector <TH1F*> vHisto;
 
  TIter hIt(file->GetListOfKeys());
  TKey *keyHisto;
  vector <string> hNames;
  while ((keyHisto=(TKey*)hIt())) {
    hNames.push_back(keyHisto->GetName());
  }
  TCanvas * c= new TCanvas(); 
  c->Draw(); 
  for (int ilayer = 1; ilayer<4; ilayer ++){
    for (int imodule = 1; imodule<9; imodule++){
      file->cd(); 
      stringstream name; 
      name.str(""); 
      name << "h_drift_depth_adc_layer"<<ilayer<<"_module"<<imodule<< "_bias_150";
      TH2F *h = (TH2F*) file->Get(name.str().c_str());
      name.str(""); 
      name << "h_drift_depth_noadc_layer"<<ilayer<<"_module"<<imodule<< "_bias_150";
      TH2F *h2 = (TH2F*) file->Get(name.str().c_str());

      TH1D * hy = (TH1D*) h->ProjectionY(); 
      TH1D * h2y = (TH1D*) h2->ProjectionY(); 
      
      hy->Rebin(2);
      h2y->Rebin(2);
      
      hy->Sumw2(); 
      h2y->Sumw2(); 
      
      hy->Divide(h2y); 
      name.str("");
      name<< "chargeProfile_layer"<< ilayer<< "_module"<< imodule<< "_bias_150"; 
	
      hy->SetTitle(name.str().c_str());
      c->SetName(name.str().c_str()); 
      c->cd(); 
      hy->Draw("E"); 
      name << ".png";
      if (saveFig) c->SaveAs(name.str().c_str());
      outfile->cd(); 
      hy->Write();
      c->Write(); 
      
				  
//     vHisto.push_back(h);
  }//imodule
  }//ilayer
  outfile->Close(); 
  file ->Close();

}
