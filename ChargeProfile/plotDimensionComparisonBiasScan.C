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
void plotDimensionComparisonBiasScan(string file1){
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
  TCanvas * c2 = new TCanvas(); 
  TCanvas * c3 = new TCanvas(); 
  c->Draw(); 
  c->SetLogy(1);
  c2->Draw(); 
  c2->SetLogy(1);
  c3->Draw(); 
  c3->SetLogy(1);
  TLegend *leg = new TLegend(0.6, 0.12, 0.85, 0.27) ; 
  leg->SetFillColor(0);

  if (DEBUG) cout << __LINE__ << endl;
  for (int ilayer=1; ilayer < 4; ilayer++){
    for (int imodule =1; imodule<9; imodule++){
      for (int ibias =0; ibias < vbias.size(); ibias++){
	if (DEBUG) cout << __LINE__ << endl;
	name.str("");
	name << "clusterSize_layer" << ilayer<<"_module"<< imodule << "_bias"<< vbias[ibias]; 
	cout << __LINE__ << " " << name.str().c_str()<< endl;
	TH1F * h = (TH1F *) _file0->Get(name.str().c_str());
	if (DEBUG) cout << __LINE__ << endl;
	c->cd(); 
	h->Sumw2();
	h->Scale(1./h->Integral());
	h->SetMarkerStyle(20);
	h->SetMarkerColor(ibias+1);
	if (ilayer ==1 && imodule ==1){
	  label.str("");
	  label << "V = " << vbias[ibias]; 
	  leg->AddEntry(h, label.str().c_str(), "lpf"); 
	}
	if (ibias == 0){
	  h->Draw("EP");
	}else{
	  h->Draw("EPSAME");
	}
	
	h->GetXaxis()->SetRangeUser(0,50);
	if (ibias == vbias.size()-1){
	  leg->Draw(); 
	  name.str("");
	  name<< "clusterSize_layer" << ilayer<<"_module"<< imodule  << "_comparison.png";
	  c->SaveAs(name.str().c_str());
	}
	  
	name.str("");
	name << "clusterSizeX_layer" << ilayer<<"_module"<< imodule << "_bias"<< vbias[ibias]; 
	cout << __LINE__ << " " << name.str().c_str()<< endl;
	TH1F * hx = (TH1F *) _file0->Get(name.str().c_str());
	if (DEBUG) cout << __LINE__ << endl;
	c2->cd(); 
	hx->Sumw2();
	hx->Scale(1./hx->Integral());
	hx->GetXaxis()->SetRangeUser(0,50);
	hx->SetMarkerStyle(20);
	hx->SetMarkerColor(ibias+1);
	if (ibias == 0){
	  hx->Draw("EP");
	}else{
	  hx->Draw("EPSAME");
	}

 
	if (ibias == vbias.size()-1){
	  leg->Draw(); 
	  name.str("");
	  name<< "clusterSizeX_layer" << ilayer<<"_module"<< imodule  << "_comparison.png";
	  c->SaveAs(name.str().c_str());
	}

	name.str("");
	name << "clusterSizeY_layer" << ilayer<<"_module"<< imodule << "_bias"<< vbias[ibias]; 
	cout << __LINE__ << " " << name.str().c_str()<< endl;
	TH1F * hy = (TH1F *) _file0->Get(name.str().c_str());
	if (DEBUG) cout << __LINE__ << endl;
	c3->cd(); 
	hy->Sumw2();
	hy->Scale(1./hy->Integral());
	hy->GetXaxis()->SetRangeUser(0,50);
	hy->SetMarkerStyle(20);
	hy->SetMarkerColor(ibias+1);
	if (ibias == 0){
	  hy->Draw("EP");
	}else{
	  hy->Draw("EPSAME");
	}

 
	if (ibias == vbias.size()-1){
	  leg->Draw(); 
	  name.str("");
	  name<< "clusterSizeY_layer" << ilayer<<"_module"<< imodule  << "_comparison.png";
	  c->SaveAs(name.str().c_str());
	}

      }//bias

    }//imodule
  }//ilayer

}
