//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug  9 10:50:31 2013 by ROOT version 5.32/00
// from TTree SiPixelLorentzAngleTree_/SiPixel LorentzAngle tree
// found on file: lorentzangleTree_all.root
//////////////////////////////////////////////////////////

#ifndef plotAngle_h
#define plotAngle_h
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "Riostream.h" 
#include <math.h>
#include <sstream>
#include <vector>
#include <string>
#include <map> 

#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxgamma = 1;


class plotChargevsBeta {
public :
  //   TFile          *f; 
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           orbit;
   Int_t           event;
   Int_t           module;
   Int_t           ladder;
   Int_t           layer;
   Int_t           isflipped;
   Float_t         pt;
   Float_t         eta;
   Float_t         PV_vx;
   Float_t         PV_vy;
   Float_t         PV_vz;
   Float_t         phi;
   Double_t        chi2;
   Double_t        ndof;
   Float_t         trackhit_x;
   Float_t         trackhit_y;
   Double_t        trackhit_alpha;
   Double_t        trackhit_beta;
   Double_t        trackhit_gamma_;
   Float_t         simhit_x;
   Float_t         simhit_y;
   Double_t        simhit_alpha;
   Double_t        simhit_beta;
   Double_t        simhit_gamma_;
   Int_t           npix;
   Float_t         rowpix[208];   //[npix]
   Float_t         colpix[208];   //[npix]
   Float_t         adc[208];   //[npix]
   Float_t         norm_charge[208];   //[npix]
   Float_t         xpix[208];   //[npix]
   Float_t         ypix[208];   //[npix]
   Float_t         clust_x;
   Float_t         clust_y;
   Float_t         clust_charge;
   Int_t           clust_size_x;
   Int_t           clust_size_y;
   Int_t           clust_maxPixelCol;
   Int_t           clust_maxPixelRow;
   Int_t           clust_minPixelCol;
   Int_t           clust_minPixelRow;
   Float_t         rechit_x;
   Float_t         rechit_y;
   Int_t           bias;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_event;   //!
   TBranch        *b_module;   //!
   TBranch        *b_ladder;   //!
   TBranch        *b_layer;   //!
   TBranch        *b_isflipped;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_PV;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_ndof;   //!
   TBranch        *b_trackhit;   //!
   TBranch        *b_simhit;   //!
   TBranch        *b_npix;   //!
   TBranch        *b_rowpix;   //!
   TBranch        *b_colpix;   //!
   TBranch        *b_adc;   //!
   TBranch        *b_norm_charge;   //!
   TBranch        *b_xpix;   //!
   TBranch        *b_ypix;   //!
   TBranch        *b_clust;   //!
   TBranch        *b_rechit;   //!
   TBranch        *b_bias;   //!

   plotChargevsBeta(TTree *tree=0);
   virtual ~plotChargevsBeta();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef plotChargevsBeta_cxx
plotChargevsBeta::plotChargevsBeta(TTree *tree) : fChain(0) 
{
  string filename = "lorentzangleMC.root";
  
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile * f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.c_str());
      if (!f || !f->IsOpen()) {
	f = new TFile(filename.c_str());
      }
      f->GetObject("SiPixelLorentzAngleTree_",tree);

   }
   Init(tree);
}

plotChargevsBeta::~plotChargevsBeta()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t plotChargevsBeta::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t plotChargevsBeta::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void plotChargevsBeta::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("module", &module, &b_module);
   fChain->SetBranchAddress("ladder", &ladder, &b_ladder);
   fChain->SetBranchAddress("layer", &layer, &b_layer);
   fChain->SetBranchAddress("isflipped", &isflipped, &b_isflipped);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("PV", &PV_vx, &b_PV);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("ndof", &ndof, &b_ndof);
   fChain->SetBranchAddress("trackhit", &trackhit_x, &b_trackhit);
   fChain->SetBranchAddress("simhit", &simhit_x, &b_simhit);
   fChain->SetBranchAddress("npix", &npix, &b_npix);
   fChain->SetBranchAddress("rowpix", rowpix, &b_rowpix);
   fChain->SetBranchAddress("colpix", colpix, &b_colpix);
   fChain->SetBranchAddress("adc", adc, &b_adc);
   fChain->SetBranchAddress("norm_charge", norm_charge, &b_norm_charge);
   fChain->SetBranchAddress("xpix", xpix, &b_xpix);
   fChain->SetBranchAddress("ypix", ypix, &b_ypix);
   fChain->SetBranchAddress("clust", &clust_x, &b_clust);
   fChain->SetBranchAddress("rechit", &rechit_x, &b_rechit);
   fChain->SetBranchAddress("bias", &bias, &b_bias);
   Notify();
}

Bool_t plotChargevsBeta::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void plotChargevsBeta::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t plotChargevsBeta::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  if (entry) {}
   return 1;
}
#endif // #ifdef readTree_cxx
