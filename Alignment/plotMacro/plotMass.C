#include <memory>

#include "Riostream.h" 
#include <math.h>
#include <sstream>
#include <vector>

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

using namespace std;
void plotMass(string decay){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gROOT->ForceStyle();
  stringstream file,file2;
  file.str("");
  file << "outfileMyMass"<<decay<<".root";


  cout << __LINE__ << " Opening file: " << file.str().c_str() << endl;
  
  TFile * f = new TFile();  
  f = TFile::Open(file.str().c_str());
  if (!f){
    cout << "FILE NOT FOUND" << endl;
    return;
  }
  cout << __LINE__ << " file " << file.str().c_str() <<  " opened" << endl;
  stringstream outstring;
  outstring.str("");
  outstring << "plotMyMass" << decay << ".root";
  
  TFile* outfile = new TFile(outstring.str().c_str(),"RECREATE");
  
  TCanvas * c = new TCanvas ("c", "c", 1400, 800);
  c->Divide(3,2) ;
  c->cd(1);
  TH1D * h1 = (TH1D*)f->Get("hMassOneCand");
  TH1D * h2 = (TH1D*)f->Get("hMassTwoCand");
  TH1D * h3 = (TH1D*)f->Get("hMassThreeCand");
  TH1D * h4 = (TH1D*)f->Get("hMassFourCand");
  TH1D * h5 = (TH1D*)f->Get("hMassFiveCand");

  TH1D * hh1 = (TH1D*)f->Get("hMassFirstCand");
  TH1D * hh2 = (TH1D*)f->Get("hMassSecondCand");
  TH1D * hh3 = (TH1D*)f->Get("hMassThirdCand");
  TH1D * hh4 = (TH1D*)f->Get("hMassForthCand");
  TH1D * hh5 = (TH1D*)f->Get("hMassFifthCand");


  h1->SetLineWidth(2);
  h1->Draw();
  c->cd(2);
  h2->SetLineWidth(2);
  h2->Draw();
  c->cd(3);
  h3->SetLineWidth(2);
  h3->Draw();
  c->cd(4);
  h4->SetLineWidth(2);
  h4->Draw();
  c->cd(5);
  h5->SetLineWidth(2);
  h5->Draw();


  TCanvas * c1 = new TCanvas ("c1", "c1", 1400, 800);
  c1->Divide(3,2) ;
  c1->cd(1);
  h1->Draw();

  c1->cd(2);
  h2->Draw();

  c1->cd(3);
  h3->Draw();

  c1->cd(4);
  h4->Draw();

  c1->cd(5);
  h5->Draw();


  
  TCanvas * c2 = new TCanvas ("c2", "c2", 1400, 800);
  c2->Divide(3,2) ;
  c2->cd(1);
  hh1->Draw();

  c2->cd(2);
  hh2->Draw();

  c2->cd(3);
  hh3->Draw();

  c2->cd(4);
  hh4->Draw();

  c2->cd(5);
  hh5->Draw();

  outfile->cd (); 

  c->Write();
  c1->Write();
  c2->Write();
//   h1->Write();
//   h2->Write();
//   h3->Write();
//   h4->Write();
//   h5->Write();

  outfile->Close();
  f->Close();

}
