#include <memory>

#include "Riostream.h" 
#include <math.h>
#include <sstream>
#include <vector>
#include <string>

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TKey.h>

#define DEBUG 1

using namespace std;
void LorentzFitNoBiasScan(string filename="lorentzangleALCARECO.root",  bool saveFig= 0){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gROOT->ForceStyle();
  gStyle->SetPadGridY(1); 
  gStyle->SetPadGridX(1); 
  gStyle->SetFuncColor(2);
  gStyle->SetFuncWidth(2);
  gStyle->SetHistLineWidth(2);

  //  int bias = 150; 
  int mybias[]  = {300, 100, 70};
  //int mybias[]  = {150};
  std::vector<int> vbias (mybias, (mybias + sizeof(mybias)/sizeof(int)));
  double pigreco = 3.141592;
  
  TFile * file = new TFile(filename.c_str(), "READ"); 
  stringstream outname;
  outname.str(""); 
  outname<< "LorentzFit_" << filename; 
  TFile *  outfile = new TFile(outname.str().c_str(), "RECREATE"); 
  
  vector <TH1F*> vHisto;
 
  TIter hIt(file->GetListOfKeys());
  TKey *keyHisto;
  vector <string> hNames;
  while ((keyHisto=(TKey*)hIt())) {
    hNames.push_back(keyHisto->GetName());
  }
  TCanvas * c= new TCanvas("c","c"); 
  c->Draw(); 
  stringstream name; 
  TLegend * leg = new TLegend (0.4, 0.12, 0.6, 0.22); 
  leg->SetFillColor(0); 
  for (unsigned int ibias = 0; ibias < vbias.size(); ibias++){
    for (int ilayer = 1; ilayer<4; ilayer ++){
      name.str(""); 
      name << "tanLAvsModule_layer" << ilayer <<"_bias"<< mybias[ibias] ; 
      TH1F * hAngle = new TH1F(name.str().c_str(), name.str().c_str(), 800, 0.5, 8.5);
      for (int imodule = 1; imodule<9; imodule++){
	file->cd(); 
	name.str(""); 
	name << "h_drift_depth_noadc_layer"<<ilayer<<"_module"<<imodule<< "_bias_"<< mybias[ibias];
	TH2F *h2 = (TH2F*) file->Get(name.str().c_str());
	h2->Rebin2D(2,2); 
	TProfile * h2y = (TProfile*) h2->ProfileY(); 
	h2y->Sumw2(); 
	name.str("");
	name<< "lorentzFit_layer"<< ilayer<< "_module"<< imodule<< "_bias_"<<  mybias[ibias]; 
	h2y->SetTitle(name.str().c_str());
	h2y->SetName(name.str().c_str());
	TF1 * f = new TF1 ("f", "pol1",-500, 500);
	h2y->Fit("f", "R", "", 100, 200);
	h2y->GetYaxis()->SetRangeUser(-50, 150); 
	hAngle->Fill (imodule, f->GetParameter(1)); 
	hAngle->SetBinError (hAngle->GetXaxis()->FindBin(imodule), f->GetParError(1)); 
	outfile->cd();
	h2y->Write();
	name.str(""); 
	name<< "canvas/lorentzFit_layer"<< ilayer<< "_module"<< imodule<< "_bias_"<<  mybias[ibias]<<".png"; 

	c->SaveAs(name.str().c_str());

      }//imodule
      outfile->cd(); 
      hAngle->Write(); 
      vHisto.push_back(hAngle); 
    }//ilayer
  }//ibias
  for (unsigned int ibias = 0; ibias < vbias.size(); ibias++){
    for (int ilayer = 1; ilayer < 4; ilayer ++){
      
      if (ilayer==1) {
	c->Clear(); 
	vHisto[ibias*3+ilayer-1]->Draw("E"); 
	vHisto[ibias*3+ilayer-1]->GetYaxis()->SetRangeUser(0.2, 0.6); 
	vHisto[ibias*3+ilayer-1]->SetMarkerStyle(20);
	vHisto[ibias*3+ilayer-1]->SetTitle("tanLAvsRing"); 
	
      } else {
	vHisto[ibias*3+ilayer-1]->Draw("ESAME");
	vHisto[ibias*3+ilayer-1]->SetMarkerStyle(20); 
	vHisto[ibias*3+ilayer-1]->SetMarkerColor(ilayer);	  
      }
      if (ibias == 0 ){
	name.str("");
	name << "layer "<< ilayer; 
	leg->AddEntry(vHisto[ibias*3+ilayer-1], name.str().c_str(), "lpf");
      }

      c->Update();
      
      if (ilayer ==3)   leg->Draw();
      
    }//ilayer
    name.str(""); 
    name << "canvas/LorentzAnglevsRing_bias"<< vbias[ibias]<< ".png";
    c->SaveAs(name.str().c_str());
    c->Write(); 
  }//ibias
//   outfile->Close(); 
//   file ->Close();

}
