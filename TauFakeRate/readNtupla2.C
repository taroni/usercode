//complining in root: root [0] .L readNtupla.C+
//executing in root : root [1] readNtupla("/data/taroni/outfile_fakerRate_ntuple_jetID_QCD50.root")
#include <TROOT.h>
#include <TStyle.h>
#include <TVector.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TH2.h"
#include "TProfile.h"
#include <TTree.h>
#include <TKey.h>
#include <TList.h>
#include <TLegend.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <memory>

using namespace std;
void readNtupla2 ( std::string & strInFile) {

  //gROOT->Reset();
  gROOT->SetStyle("Plain");
  //gStyle->SetOptFit(1);
  gStyle->SetHistLineWidth(2);

  std::string argStr;
  std::ifstream ifs(strInFile.c_str());
  cout << "readFile " << strInFile << endl;
  std::vector<std::string> fileList;

  while (ifs.good()){
    std::getline(ifs,argStr);
    string tmp = argStr;
    fileList.push_back(tmp);
    cout << "reading file "  << tmp << endl;
    if(fileList.size()>0){
      string tmp2 (".root");
      size_t foundFile;
      foundFile=fileList[fileList.size()-1].find(tmp2);
      if (foundFile==string::npos) fileList.pop_back();
  }
  }
  if (ifs.is_open()) ifs.close();
  cout << "I'm going to read " << fileList.size() << " files" << endl;
  std::string str(".list");
  size_t found;
  found = strInFile.find(str);
  cout << "found " << found << endl;
  stringstream  outfilename;
  outfilename <<"jetID_"<<  strInFile.substr(0,int(found))<< ".root";
  
  cout << " outfile name " << outfilename.str().c_str() << endl;
  TFile* outfile = new TFile(outfilename.str().c_str(),"recreate");
  outfile->cd();
  TH2D * h_TauLoose_1 = new TH2D ("h_TauLoose_1","Pt vs Eta, Loose1", 3, 0, 2.4, 20, 0, 200);
  TH2D * h_TauTight_1 = new TH2D ("h_TauTight_1","Pt vs Eta, Tight1", 3, 0, 2.4, 20, 0, 200);
  TH2D * h_TauTight_2 = new TH2D ("h_TauTight_2","Pt vs Eta, Tight2", 3, 0, 2.4, 20, 0, 200);

  for (int nFile=0; nFile < fileList.size();nFile++){
    string filename= fileList[nFile];
    TFile * f1 = TFile::Open(filename.c_str());
    if (!f1) {
      cout << "file not found"<< endl;
      return;
    }
    f1->cd();
    TTree *t1 = (TTree*)f1->FindObjectAny("ntuple");
    cout << "Ntuple read" << endl;
    if (!t1) {
      cout << "Tree not found" << endl;
      return;
    }
    
    cout << " reading file: " << filename << endl;

    std::vector<double>  *tauPt=0, *tauEta=0, *bTag=0, *mcTruth=0, *hT=0;
    std::vector<double>  *looseIso=0, *mediumIso=0, *tightIso=0, *againstEle=0, *againstMu=0, *byDecay=0;
    
    TBranch * tauPtB      = (TBranch *) t1-> GetBranch("tauPt");
    TBranch * tauEtaB     = (TBranch *) t1-> GetBranch("tauEta");
    TBranch * bTagB       = (TBranch *) t1-> GetBranch("bTag");
    TBranch * mcTruthB    = (TBranch *) t1-> GetBranch("mcTruth");
    TBranch * hTB         = (TBranch *) t1-> GetBranch("hT");
    TBranch * looseIsoB   = (TBranch *) t1-> GetBranch("looseIso");
    TBranch * mediumIsoB  = (TBranch *) t1-> GetBranch("mediumIso");
    TBranch * tightIsoB   = (TBranch *) t1-> GetBranch("tightIso");
    TBranch * againstEleB = (TBranch *) t1-> GetBranch("againstEle");
    TBranch * againstMuB  = (TBranch *) t1-> GetBranch("againstMu");
    TBranch * byDecayB    = (TBranch *) t1-> GetBranch("byDecay");


    tauPtB      ->SetAddress(&(tauPt));
    tauEtaB     ->SetAddress(&(tauEta));
    bTagB       ->SetAddress(&(bTag));
    mcTruthB    ->SetAddress(&(mcTruth));
    hTB         ->SetAddress(&(hT));
    looseIsoB   ->SetAddress(&(looseIso));
    mediumIsoB  ->SetAddress(&(mediumIso));
    tightIsoB   ->SetAddress(&(tightIso));
    againstEleB ->SetAddress(&(againstEle));
    againstMuB  ->SetAddress(&(againstMu));
    byDecayB    ->SetAddress(&(byDecay));

    int n = (int) t1->GetEntries();
    cout << " Tree Entries " << n << endl; 
    
    for (int evt = 0; evt < n; ++evt) { 
    t1->GetEntry(evt);	
    if(hT->size()<=0 ) continue;
    if (hT->at(0) < 300) continue ;
    
    for (unsigned int itau =0 ; itau <  tauPt->size(); itau++){
      if (againstEle->at(itau)>0.5 && againstMu->at(itau)>0.5 && byDecay->at(itau)>0.5){
	h_TauLoose_1->Fill(fabs(tauEta->at(itau)),tauPt->at(itau));
	if (looseIso->at(itau)>0.5){
	  h_TauTight_1->Fill (fabs(tauEta->at(itau)),tauPt->at(itau));
	  if (mediumIso->at(itau)>0.5){
	    h_TauTight_2->Fill (fabs(tauEta->at(itau)),tauPt->at(itau));
	  }//tight2
	}//tight1
      }//loose1
    }
    
    }//evt
    
    tauPtB      ->~TBranch();
    tauEtaB     ->~TBranch();
    bTagB       ->~TBranch();
    mcTruthB    ->~TBranch();
    hTB         ->~TBranch();
    looseIsoB   ->~TBranch();
    mediumIsoB  ->~TBranch();
    tightIsoB   ->~TBranch();
    againstEleB ->~TBranch();
    againstMuB  ->~TBranch();
    byDecayB    ->~TBranch();
    t1->~TTree();
    f1->Close();
    f1->~TFile();

  }//nFile

      cout << __LINE__<<" readFile " << strInFile << endl;

  TH2D * eps_1 = (TH2D*) h_TauTight_1->Clone();
  TH2D * eps_2 = (TH2D*) h_TauTight_2->Clone();
  eps_1->SetName("h_tightToLoose1");
  eps_2->SetName("h_tightToLoose2");
  eps_1->SetTitle("h_tightToLoose1");
  eps_2->SetTitle("h_tightToLoose2");

  eps_1->Sumw2();
  eps_2->Sumw2();
  eps_1->Divide(h_TauLoose_1);
  eps_2->Divide(h_TauTight_1);

  cout<<__LINE__<< endl;
  //binomial errors
  cout << "Pt bins: " << h_TauTight_1->GetYaxis()->GetNbins() << "; Eta bins " << h_TauTight_1->GetXaxis()->GetNbins() << endl;
  for (unsigned int ibinPt = 1; ibinPt!=h_TauTight_1->GetYaxis()->GetNbins()+1; ibinPt++){
  cout<<__LINE__<< endl;
    for (unsigned int ibinEta = 1; ibinEta!=h_TauTight_1->GetXaxis()->GetNbins()+1; ibinEta++){
      if (h_TauLoose_1->GetBinContent(ibinEta, ibinPt) !=0){
	  double error_1 = 0;
	if (h_TauTight_1->GetBinContent(ibinEta, ibinPt) != 0 ){
	  if (h_TauTight_1->GetBinContent(ibinEta, ibinPt) ==1 && h_TauLoose_1->GetBinContent(ibinEta, ibinPt)==1 ){
	    error_1=1;
	  }else {
	    error_1 = (1./h_TauLoose_1->GetBinContent(ibinEta, ibinPt))*sqrt(h_TauTight_1->GetBinContent(ibinEta, ibinPt) * (1- h_TauTight_1->GetBinContent(ibinEta, ibinPt) /h_TauLoose_1->GetBinContent(ibinEta, ibinPt)));
	  }
	  eps_1->SetBinError(ibinEta, ibinPt, error_1);
	}
      }//ifbincontent1
      cout<<__LINE__<< endl;
      if (h_TauTight_1->GetBinContent(ibinEta, ibinPt) !=0){
	double error_2 = 0; 
	if (h_TauTight_2->GetBinContent(ibinEta, ibinPt)!=0){
	  if (h_TauTight_2->GetBinContent(ibinEta, ibinPt) ==1 && h_TauTight_1->GetBinContent(ibinEta, ibinPt)==1 ){
	    error_2=1;
	  }else {
	    error_2 = (1./h_TauTight_1->GetBinContent(ibinEta, ibinPt))*sqrt(h_TauTight_2->GetBinContent(ibinEta, ibinPt) * (1- h_TauTight_2->GetBinContent(ibinEta, ibinPt) /h_TauTight_1->GetBinContent(ibinEta, ibinPt)));
	  }
	eps_2->SetBinError(ibinEta, ibinPt, error_2);
	}
      }//ifbincontent2
    }//ibinEta
    cout << "pt bin " << ibinPt << endl;
  }//ibinPt
  cout<<__LINE__<< endl;

  outfile->cd();
  h_TauLoose_1->Write();
  h_TauTight_1->Write();
  h_TauTight_2->Write();
  eps_1->Write();
  eps_2->Write();
//   profy_eta5->Write();
  TH1D *  proj_tight1_eta1 = eps_1->ProjectionY("proj_tight1_eta1",1,1);
  TH1D *  proj_tight1_eta2 = eps_1->ProjectionY("proj_tight1_eta2",2,2);
  TH1D *  proj_tight1_eta3 = eps_1->ProjectionY("proj_tight1_eta3",3,3);
  TH1D *  proj_tight2_eta1 = eps_2->ProjectionY("proj_tight2_eta1",1,1);
  TH1D *  proj_tight2_eta2 = eps_2->ProjectionY("proj_tight2_eta2",2,2);
  TH1D *  proj_tight2_eta3 = eps_2->ProjectionY("proj_tight2_eta3",3,3);

  stringstream histoName;
  histoName << outfilename.str().substr(0, outfilename.str().length()-5) << "proj_tight1_eta08";
  cout << histoName.str()<< endl;

  proj_tight1_eta1->SetName (histoName.str().c_str());
  proj_tight1_eta1->SetTitle(histoName.str().c_str());
  proj_tight1_eta1->Draw();

  histoName.str("");
  histoName << outfilename.str().substr(0, outfilename.str().length()-5) << "proj_tight1_eta16";
  proj_tight1_eta2->SetName(histoName.str().c_str());
  proj_tight1_eta2->SetTitle(histoName.str().c_str());
  proj_tight1_eta2->Draw();

  histoName.str("");
  histoName << outfilename.str().substr(0, outfilename.str().length()-5) << "proj_tight1_eta24";
  proj_tight1_eta3->SetName(histoName.str().c_str());
  proj_tight1_eta3->SetTitle(histoName.str().c_str());
  proj_tight1_eta3->Draw();

  histoName.str("");
  histoName << outfilename.str().substr(0, outfilename.str().length()-5) << "proj_tight2_eta08";
  proj_tight2_eta1->SetName(histoName.str().c_str());
  proj_tight2_eta1->SetTitle(histoName.str().c_str());
  proj_tight2_eta1->Draw();


  histoName.str("");
  histoName << outfilename.str().substr(0, outfilename.str().length()-5) << "proj_tight2_eta16";
  proj_tight2_eta2->SetName(histoName.str().c_str());
  proj_tight2_eta2->SetTitle(histoName.str().c_str());
  proj_tight2_eta2->Draw();


  histoName.str("");
  histoName << outfilename.str().substr(0, outfilename.str().length()-5) << "proj_tight2_eta24";
  proj_tight2_eta3->SetName(histoName.str().c_str());
  proj_tight2_eta3->SetTitle(histoName.str().c_str()); 
  proj_tight2_eta3->Draw();

  TCanvas *c = new TCanvas();
  c->Draw();
  c->cd();
  proj_tight1_eta1->Draw();
  proj_tight1_eta2->Draw("SAME");
  proj_tight1_eta2->SetLineColor(2);
  proj_tight1_eta3->Draw("SAME");
  proj_tight1_eta3->SetLineColor(3);
  proj_tight2_eta1->Draw("SAME");
  proj_tight2_eta1->SetLineColor(4);
  proj_tight2_eta2->Draw("SAME");
  proj_tight2_eta2->SetLineColor(5);
  proj_tight2_eta3->Draw("SAME");
  proj_tight2_eta3->SetLineColor(6);
  TLegend * leg = new TLegend (0.2,0.2,0.6,0.6,"", "brNDC");
  leg->AddEntry(proj_tight1_eta1, "|#eta| < 0.8, tight1/loose" , "lpf" );
  leg->AddEntry(proj_tight2_eta1, "|#eta| < 0.8, tight2/tight1", "lpf" );
  leg->AddEntry(proj_tight1_eta2, "0.8 < |#eta| < 1.6, tight1/loose" , "lpf" );
  leg->AddEntry(proj_tight2_eta2, "0.8 < |#eta| < 1.6, tight2/tight1", "lpf" );
  leg->AddEntry(proj_tight1_eta3, "1.6 < |#eta| < 2.4, tight1/loose" , "lpf" );
  leg->AddEntry(proj_tight2_eta3, "1.6 < |#eta| < 2.4, tight2/tight1", "lpf" );
  leg->Draw();
  c->Update();
  stringstream canvasName;
  canvasName.str("");
  canvasName <<"canvas" <<outfilename.str().substr(0, outfilename.str().length()-5) << ".root";
  c->SaveAs(canvasName.str().c_str());

  outfile->cd();
  c->Write();
  proj_tight1_eta1->Write();
  proj_tight1_eta2->Write();
  proj_tight1_eta3->Write();
  proj_tight2_eta1->Write();
  proj_tight2_eta2->Write();
  proj_tight2_eta3->Write();


  outfile->Close();
}
