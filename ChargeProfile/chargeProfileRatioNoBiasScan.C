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
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TKey.h>

#define DEBUG 1

using namespace std;
void chargeProfileRatioNoBiasScan(string filename="lorentzangleALCARECO.root",  bool saveFig= 0){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  gStyle->SetPadGridY(1); 
  gStyle->SetPadGridX(1); 
  
  map <int, pair<double, double> > myRatioMap;

  double pigreco = 3.141592;
  int bias = 150 ;

  TFile * file = new TFile(filename.c_str(), "READ"); 
  stringstream outname;
  outname.str(""); 
  outname<< "chargeRatio_" << filename; 
  TFile *  outfile = new TFile(outname.str().c_str(), "RECREATE"); 
  outname.str(outname.str().substr(0,outname.str().find(".root")));
  outname <<".txt"; 
  FILE * pFile;
  pFile = fopen (outname.str().c_str(),"w");
  fprintf (pFile,"layer\t module\t bias\t ratio\t errratio\n");
  
 
  TIter hIt(file->GetListOfKeys());
  TKey *keyHisto;
  vector <string> hNames;
  while ((keyHisto=(TKey*)hIt())) {
    hNames.push_back(keyHisto->GetName());
  }
  TCanvas * c= new TCanvas(); 
  c->Draw(); 
  stringstream name; 

  for (int ilayer = 1; ilayer<4; ilayer ++){
    vector <TH1D*> vHisto;

    name.str(""); 
    name << "h_ratiovsmodule_layer"<< ilayer<<"_bias" << bias ; 
    TH1F * hratio = new TH1F (name.str().c_str(), name.str().c_str(), 400, 0.5, 4.5); 
    for (int imodule = 1; imodule<9; imodule++){
      file->cd(); 
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
      
      vHisto.push_back(hy); 

      
				  
//     vHisto.push_back(h);
    }//imodule
    if (DEBUG) cout << __LINE__ << " " << vHisto.size() << endl;
    if (DEBUG) cout <<__LINE__ << endl;
    for (int imodule=1; imodule < 5; imodule++) {
      if (DEBUG) cout <<__LINE__ << endl;
      name.str("");
      name << "hratio_module"<< 9-imodule << "_module" << imodule << "_bias"<< bias; 
      if (DEBUG) cout <<__LINE__ << endl;

      TH1D * hh = vHisto[8-imodule];
      if (DEBUG) cout <<__LINE__ << endl;
      
      hh->SetName(name.str().c_str()); 
      hh->SetTitle(name.str().c_str()); 
      stringstream name2 ; 
      name2.str() ;
      name2 << "c_"<< name.str() ; 
      c->SetName(name2.str().c_str()); 
      if (DEBUG) cout <<__LINE__ << endl;

      hh->Divide((TH1F*)vHisto[imodule-1]->Clone());
      if (DEBUG) cout <<__LINE__ << endl;
    
      TF1 * f = new TF1 ("f", "pol0", -500, 500); 
      hh->Fit("f", "RS", "S", 50, 250); 
      int histoID = 10000* ilayer + 1000*imodule+ bias; 
      if (DEBUG) cout <<__LINE__ << endl;
      pair <double, double> myPair (f->GetParameter(0), f->GetParError(0)); 
      myRatioMap[histoID] = myPair; 
      fprintf (pFile,"%i\t %i\t %i\t %f\t %f\n", ilayer, imodule, bias, f->GetParameter(0), f->GetParError(0)); 
      hratio->Fill(imodule,  f->GetParameter(0)); 
      hratio->SetBinError(hratio->GetXaxis()->FindBin(imodule), f->GetParError(0)); 
      hratio->SetMarkerStyle(20); 
      hratio->GetYaxis()->SetRangeUser(0.5,1.5); 
      name << ".png";
      if (DEBUG) cout <<__LINE__ << endl;
      if (saveFig) c->SaveAs(name.str().c_str());
      outfile->cd(); 
      hh->Write();
      c->Write(); 
    }
    outfile->cd(); 
    if (DEBUG) cout <<__LINE__ << endl;
    c->cd(); 
    hratio->Draw("E"); 
    TPaveText *pt = new TPaveText(.66,.53,1.6,.8);
    pt->SetFillColor(0); 
    pt->AddText("1=ring8/ring1");
    pt->AddText("2=ring7/ring2");
    pt->AddText("3=ring6/ring3");
    pt->AddText("4=ring5/ring4");
      
    pt->Draw("SAME");
    name.str(""); 
    name << "c_" << hratio->GetName();
    c->SetName(name.str().c_str());
    name<< ".png"; 
    if (saveFig)c->SaveAs(name.str().c_str()); 
    c->Write();
    hratio->Write(); 

  }//ilayer
  fclose (pFile);
  outfile->Close(); 
  file ->Close();

}
