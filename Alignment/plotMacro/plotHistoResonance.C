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
void plotHistoResonance(string resonance){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptStat(1);
  gROOT->ForceStyle();
  stringstream file,file2;
  file.str("");
  file << "my"<<resonance<<"MuMu.root";
  cout << __LINE__ << " Opening file: " << file.str().c_str() << endl;
  TFile * f = new TFile();  
  f = TFile::Open(file.str().c_str());
  if (!f){
    cout << "FILE NOT FOUND" << endl;
    return;
  }
  cout << __LINE__ << " file " << file.str().c_str() <<  " opened" << endl;
//   stringstream outstring;
//   outstring.str("");
//   outstring << "plotMyMass" << decay << ".root";
//   TFile* outfile = new TFile(outstring.str().c_str(),"RECREATE");

  TDirectory * dir = (TDirectory *) f->Get("myanalysis"); 
 
  TCanvas * c = new TCanvas ("c", "c");

  stringstream masshisto; 
  masshisto.str(""); 
  masshisto << "hInvMass" << resonance; 
  TH1D * h1 = (TH1D*) dir->Get(masshisto.str().c_str());
  TH1D * h2 = (TH1D*) dir->Get("hchi2");
  TH1D * h3 = (TH1D*) dir->Get("hNtrk");
  TH1D * h4 = (TH1D*) dir->Get("hP");   
  TH1D * h5 = (TH1D*) dir->Get("hPt");  
  TH1D * h6 = (TH1D*) dir->Get("hHit");
  TH1D * h7 = (TH1D*) dir->Get("hEta");
  TH1D * h8 = (TH1D*) dir->Get("hPhi");
  TH1D * h9 = (TH1D*) dir->Get("hvx"); 
  TH1D * h10 = (TH1D*) dir->Get("hvy"); 
  TH1D * h11 = (TH1D*) dir->Get("hvz"); 
  TH1D * h12 = (TH1D*) dir->Get("hd0"); 
  TH1D * h13 = (TH1D*) dir->Get("nhpxb");
  TH1D * h14 = (TH1D*) dir->Get("nhpxe");
  TH1D * h15 = (TH1D*) dir->Get("nhTIB");
  TH1D * h16 = (TH1D*) dir->Get("nhTID");
  TH1D * h17 = (TH1D*) dir->Get("nhTOB");
  TH1D * h18 = (TH1D*) dir->Get("nhTEC");
  TH1D * h19 = (TH1D*) dir->Get("hPx"); 
  TH1D * h20 = (TH1D*) dir->Get("hPy"); 
  TH1D * h21 = (TH1D*) dir->Get("hPz"); 
  TH1D * h22 = (TH1D*) dir->Get("hmuQ");
  TH1D * h23 = (TH1D*) dir->Get("hdxy");
  TH1D * h24 = (TH1D*) dir->Get("hdz"); 
  TH1D * h25 = (TH1D*) dir->Get("hEtaMother"); 
  TH1D * h26 = (TH1D*) dir->Get("hPtMother");  
  TH1D * h27 = (TH1D*) dir->Get("hPhiMother"); 
  TH1D * h28 = (TH1D*) dir->Get("hDeltaPhi");  
  TH1D * h29 = (TH1D*) dir->Get("hDeltaEta");  
  TH1D * h30 = (TH1D*) dir->Get("R"); 

  vector <TH1D *> vHisto; 
  vHisto.push_back(h1); 
  vHisto.push_back(h2); 
  vHisto.push_back(h3); 
  vHisto.push_back(h4); 
  vHisto.push_back(h5); 
  vHisto.push_back(h6); 
  vHisto.push_back(h7); 
  vHisto.push_back(h8); 
  vHisto.push_back(h9); 
  vHisto.push_back(h10); 
  vHisto.push_back(h11); 
  vHisto.push_back(h12); 
  vHisto.push_back(h13); 
  vHisto.push_back(h14); 
  vHisto.push_back(h15); 
  vHisto.push_back(h16); 
  vHisto.push_back(h17); 
  vHisto.push_back(h18); 
  vHisto.push_back(h19); 
  vHisto.push_back(h20); 
  vHisto.push_back(h21); 
  vHisto.push_back(h22); 
  vHisto.push_back(h23); 
  vHisto.push_back(h24); 
  vHisto.push_back(h25); 
  vHisto.push_back(h26); 
  vHisto.push_back(h27); 
  vHisto.push_back(h28); 
  vHisto.push_back(h29); 
  vHisto.push_back(h30); 
  

  stringstream filename; 
  for (int i=0 ; i <  vHisto.size() ; i++ ) {
    vHisto[i] ->Draw();
    filename.str("");
    filename << vHisto[i]->GetName() << ".png"; 
    c->SaveAs(filename.str().c_str()); 
  }

//   c->~TCanvas();
  
//   for (int i=1 ; i < 31 ; i++ ) {
//     vHisto[i] ->~TH1D();
//   }
  f->Close(); 

}
