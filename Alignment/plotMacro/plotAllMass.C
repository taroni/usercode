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
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#define DEBUG 1

using namespace std;
void plotAllMass(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gROOT->ForceStyle();
  string file1,file2, file3,file4, file5,file6,file7, file8,file9;
  
  vector < string > filename; 

  filename.push_back("myUpsilonMuMu_noIso.root");
  filename.push_back("myUpsilonMuMu_03.root");
  filename.push_back("myUpsilonMuMu_03NoTrk.root");
  filename.push_back("myUpsilonMuMu_03NoTrk_Tk015.root");
  filename.push_back("myUpsilonMuMu_015NoTrk_Tk015.root");
  filename.push_back("myUpsilonMuMu_015NoTrk_Tk005.root");
  filename.push_back("myUpsilonMuMu_015NoTrk_.root");
  filename.push_back("myUpsilonMuMu_015OnlyTrk_.root");
  filename.push_back("myUpsilonMuMu_separated.root");
  filename.push_back("myUpsilonMuMu_all.root");
  filename.push_back("myUpsilonMuMu_005OnlyTrk_.root");
  filename.push_back("myUpsilonMuMu_005OnlyEm.root");
  filename.push_back("myUpsilonMuMu_005OnlyHad.root");
  filename.push_back("myUpsilonMuMu_separated005.root");
  filename.push_back("myUpsilonMuMu_005NoTrk.root");
  filename.push_back("myUpsilonMuMu_03NoTrk_Tk005.root");
  filename.push_back("myUpsilonMuMu_05NoTrk_Tk003.root");

  stringstream outstring;
  outstring.str("");
  outstring << "plotAllMyMass.root";
  if (DEBUG) cout << __LINE__ << endl;
  
  TFile* outfile = new TFile(outstring.str().c_str(),"RECREATE");
  outfile -> cd(); 
//   TCanvas * c = new TCanvas ("c", "c", 1400, 800);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gROOT->ForceStyle();
  TCanvas * c = new TCanvas ("c", "c", 1200, 800);
  if (DEBUG) cout << __LINE__ << endl;
  c->Draw();
  //  c->Divide(4,4);
  vector < TFile * > vFile(filename.size(), new TFile());
  int thePad = 1 ;
  for (int ifile = 0 ; ifile<filename.size(); ifile++){
    cout << __LINE__ << " Opening file: " << filename[ifile] << endl;
    if (DEBUG) cout << __LINE__ << endl;

   vFile [ifile] = TFile::Open(filename[ifile].c_str());
    if (DEBUG) cout << __LINE__ << endl;

//   c->Divide(3,2) ;
//   c->cd(1);
    stringstream cName; 
    cName.str("");
    cName << filename[ifile].substr(0,filename[ifile].length()-5) << ".png";
    cout << cName.str().c_str() << endl;
//     c->SetName(cName.str().c_str());
//     c->SetTitle(cName.str().c_str());
    TDirectory * dir1 = (TDirectory *) vFile[ifile] ->Get("demo");
    if (DEBUG) cout << __LINE__ << endl;
    TH1D * h1 = (TH1D*)dir1->Get("hInvMassUpsilon");
    if (DEBUG) cout << __LINE__ << endl;
    outfile -> cd(); 
    c->cd(thePad);
    h1->SetLineWidth(2);
    h1->SetLineColor(1);
    h1->Draw();
//     thePad++;
    c->Update(); 
    c-> Write();
    c->SaveAs(cName.str().c_str());

  }
 

}
