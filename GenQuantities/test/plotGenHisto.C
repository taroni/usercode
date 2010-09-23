#include <TROOT.h>
#include <TStyle.h>
#include <TVector.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TKey.h>
#include <TList.h>
#include <TLegend.h>
#include <iostream>
#include <sstream>
using namespace std;
void plotGenHisto () {
  //gROOT->Reset();
  gROOT->SetStyle("Plain");
  //gStyle->SetOptFit(1);
  gStyle->SetHistLineWidth(2);
  
  
  TFile * f1 = TFile::Open("H140toWWto2Mu2NuHistos-100000.root");
  
  TH1F * genHiggsPt  = (TH1F *) f1->Get("genHiggsPt");
  TH1F * genHiggsP   = (TH1F *) f1->Get("genHiggsP");
  TH1F * genHiggsE   = (TH1F *) f1->Get("genHiggsE");
  TH1F * genHiggsEta = (TH1F *) f1->Get("genHiggsEta");
  TH1F * genHiggsPhi = (TH1F *) f1->Get("genHiggsPhi");

  TH1F * genWPt  = (TH1F *) f1->Get("genWPt");
  TH1F * genWP   = (TH1F *) f1->Get("genWP");
  TH1F * genWE   = (TH1F *) f1->Get("genWE");
  TH1F * genWEta = (TH1F *) f1->Get("genWEta");
  TH1F * genWPhi = (TH1F *) f1->Get("genWPhi");

  TH1F * genMuonPt  = (TH1F *) f1->Get("genMuonPt");
  TH1F * genMuonP   = (TH1F *) f1->Get("genMuonP");
  TH1F * genMuonE   = (TH1F *) f1->Get("genMuonE");
  TH1F * genMuonEta = (TH1F *) f1->Get("genMuonEta");
  TH1F * genMuonPhi = (TH1F *) f1->Get("genMuonPhi");

  TCanvas * c1 = new TCanvas ("c1","c1"); 
  c1->Draw();
  c1->cd();
  
  c1->SetLogy();
  c1->SetLogx(0);
  genHiggsPt->Draw();
  genHiggsPt->GetXaxis()->SetTitle("GeV/c");
  c1->SaveAs("genHiggsPt.eps");
  
  genHiggsP->Draw();
  genHiggsP->GetXaxis()->SetTitle("GeV/c");
  c1->Update();
  c1->SaveAs("genHiggsP.eps");
 
  genHiggsE->Draw();
  genHiggsE->GetXaxis()->SetTitle("GeV");
  c1->Update();
  c1->SaveAs("genHiggsE.eps");

  genHiggsEta->Draw();
  genHiggsEta->GetXaxis()->SetTitle("#eta");
  c1->Update();
  c1->SaveAs("genHiggsEta.eps");

  genHiggsPhi->Draw();
  genHiggsPhi->GetXaxis()->SetTitle("#phi");
  c1->Update();
  c1->SaveAs("genHiggsPhi.eps");

  genWPt->Draw();
  genWPt->GetXaxis()->SetTitle("GeV/c");
  c1->SaveAs("genWPt.eps");
  
  genWP->Draw();
  genWP->GetXaxis()->SetTitle("GeV/c");
  c1->Update();
  c1->SaveAs("genWP.eps");
 
  genWE->Draw();
  genWE->GetXaxis()->SetTitle("GeV");
  c1->Update();
  c1->SaveAs("genWE.eps");

  genWEta->Draw();
  genWEta->GetXaxis()->SetTitle("#eta");
  c1->Update();
  c1->SaveAs("genWEta.eps");

  c1->SetLogy(0);
  genWPhi->Draw();
  genWPhi->GetXaxis()->SetTitle("#phi");
  c1->Update();
  c1->SaveAs("genWPhi.eps");


  c1->SetLogy(1);

  genMuonPt->Draw();
  genMuonPt->GetXaxis()->SetTitle("GeV/c");
  genMuonPt->GetXaxis()->SetRangeUser(0,500);
  c1->Update();
  c1->SaveAs("genMuonPt.eps");
  
  genMuonP->Draw();
  genMuonP->GetXaxis()->SetTitle("GeV/c");
  c1->Update();
  c1->SaveAs("genMuonP.eps");
 
  genMuonE->Draw();
  genMuonE->GetXaxis()->SetTitle("GeV");
  c1->Update();
  c1->SaveAs("genMuonE.eps");

  genMuonEta->Draw();
  genMuonEta->GetXaxis()->SetTitle("#eta");
  c1->Update();
  c1->SaveAs("genMuonEta.eps");

  c1->SetLogy(0);
  genMuonPhi->Draw();
  genMuonPhi->GetXaxis()->SetTitle("#phi");
  c1->Update();
  c1->SaveAs("genMuonPhi.eps");

  f1->Close();

}
