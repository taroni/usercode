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
void plotProfileBiasScan(string file1="driftdepthBiasScan.root"){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  gStyle->SetPadGridY(1); 
  gStyle->SetPadGridX(1); 

  stringstream name, label; 


  if (DEBUG) cout << __LINE__ << endl;
  int mybias[]  = {300,100,70};
  std::vector<int> vbias (mybias, (mybias + sizeof(mybias)/sizeof(int)));
  
  TFile *_file0 = TFile::Open(file1.c_str());

  TCanvas * c = new TCanvas(); 

  c->Draw(); 
  TLegend *leg = new TLegend(0.6, 0.12, 0.85, 0.27) ; 
  leg->SetFillColor(0);

  if (DEBUG) cout << __LINE__ << endl;
  for (int ilayer=1; ilayer < 4; ilayer++){
    for (int imodule =1; imodule<9; imodule++){
      for (int ibias =0; ibias < vbias.size(); ibias++){
	if (DEBUG) cout << __LINE__ << endl;
	name.str("");
	name << "h_drift_depth_adc_layer"<<ilayer<< "_module"<< imodule << "_bias_"<< vbias[ibias]; 
	cout << __LINE__ << " " << name.str().c_str()<< endl;
	TH2F * h = (TH2F *) _file0->Get(name.str().c_str());
	name.str("");
	name << "h_drift_depth_noadc_layer"<<ilayer<< "_module"<< imodule << "_bias_"<< vbias[ibias]; 
	cout << __LINE__ << " " << name.str().c_str()<< endl;
	TH2F * h2 = (TH2F *) _file0->Get(name.str().c_str());
	if (DEBUG) cout << __LINE__ << endl;
	TH1F * hy = (TH1F*) h->ProjectionY(); 
	if (DEBUG) cout << __LINE__ << endl;
	TH1F * h2y = (TH1F*) h2->ProjectionY(); 
	if (DEBUG) cout << __LINE__ << endl;
	hy->Sumw2();
	h2y->Sumw2(); 

	hy->Divide(h2y) ;
	
	if (DEBUG) cout << __LINE__ << endl;
	c->cd(); 
	hy->Sumw2();
	//hy->Scale(1./hy->Integral());
	//	h->SetMarkerStyle(20);
	hy->SetLineWidth(2);
	hy->SetLineColor(ibias+1);
	if (ilayer ==1 && imodule ==1){
	  label.str("");
	  label << "V = " << vbias[ibias]; 
	  leg->AddEntry(hy, label.str().c_str(), "lpf"); 
	}
	if (ibias == 0){
	  hy->Draw("E");
	}else{
	  hy->Draw("ESAME");
	}
	
	hy->GetYaxis()->SetRangeUser(0,15000);
	if (ibias == vbias.size()-1){
	  leg->Draw(); 
	  name.str("");
	  name<< "canvas/profile_layer"<<ilayer<< "_module"<< imodule << ".png";
	  c->SaveAs(name.str().c_str());
	}
	  

      }//bias

    }//imodule
  }//ilayer

}
