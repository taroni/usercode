//Plot the validation histograms
// last update 9th Sep 2010
//----------------------------------------------
// to run  in ROOT 
// root [1] .L compareZ.C+
// root [2] compareZ(saveFig)
// where decay is ZMuMu, JpsiMuMu, UpsilonMuMU
// and saveFig = 1 or 0, if 1 the canvases will 
// be saved in png and eps format
//-----------------------------------------------
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
void compareZ( bool saveFig= 0){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gROOT->ForceStyle();

  double pigreco = 3.141592;
  
  TFile * f = new TFile("");
  TFile * f1 = new TFile("");
  
  f = TFile::Open("outfileZMuMu.root");
  f1 = TFile::Open("outfileZMuMu_newHLT.root");
  if (!f){
    cout << "FILE NOT FOUND" << endl;
    return;
  }
  if (!f1){
    cout << "FILE NOT FOUND" << endl;
    return;
  }
  TFile* outfile = new TFile("compareZ.root","RECREATE");
  
  f->cd();
 
  TH1D *hchi2 = (TH1D *) f->Get("hchi2");
  TH1D *hNtrk = (TH1D *) f->Get("hNtrk");
  TH1D *hMass = (TH1D *) f->Get("hMass");
  TH1D *hP    = (TH1D *) f->Get("hP");
  TH1D *hPPlus= (TH1D *) f->Get("hPPlus");
  TH1D *hPMinus=(TH1D *) f->Get("hPMinus");
  TH1D *hPt   = (TH1D *) f->Get("hPt");
  TH1D *hMinPt= (TH1D *) f->Get("hMinPt");
  TH1D *hPtPlus = (TH1D *) f->Get("hPtPPlus");
  TH1D *hPtMinus= (TH1D *) f->Get("hPtMinus");
  TH1D *hHit  = (TH1D *) f->Get("hHit");
  TH1I *hHitPlus= (TH1I *) f->Get("hHitPlus");
  TH1I *hHitMinus=(TH1I *) f->Get("hHitMinus");


  TH1D *hEta     = (TH1D *) f->Get("hEta");
  TH1D *hEtaPlus = (TH1D *) f->Get("hEtaPlus");
  TH1D *hEtaMinus= (TH1D *) f->Get("hEtaMinus");
  TH1D *hPhi     = (TH1D *) f->Get("hPhi");
  TH1D *hPhiPlus = (TH1D *) f->Get("hPhiPlus");
  TH1D *hPhiMinus= (TH1D *) f->Get("hPhiMinus");

  TH1D *hDeltaPhi= (TH1D *) f->Get("hDeltaPhi");
  TH1D *hDeltaEta= (TH1D *) f->Get("hDeltaEta");
  TH1D *hDeltaR  = (TH1D *) f->Get("R");

//   TH2D *hd0DeltaPhivsPtMu1  = new TH2D("hd0DeltaPhivsPtMu1"  ,"hd0DeltaPhivsPtMu1"  ,100,0.,100.,100,-10.,10.);
//   TH2D *hd0DeltaEtavsPtMu1  = new TH2D("hd0DeltaEtavsPtMu1"  ,"hd0DeltaEtavsPtMu1"  ,100,0.,100.,100,-5.,5.);
//   TH2D *hd0DeltaPhivsPtMu2  = new TH2D("hd0DeltaPhivsPtMu2"  ,"hd0DeltaPhivsPtMu2"  ,100,0.,100.,100,-10.,10.);
//   TH2D *hd0DeltaEtavsPtMu2  = new TH2D("hd0DeltaEtavsPtMu2"  ,"hd0DeltaEtavsPtMu2"  ,100,0.,100.,100,-5.,5.);


//   TH1D *hvx  = new TH1D("hvx"  ,"hvx"  ,100,0.,.5);
//   TH1D *hvy  = new TH1D("hvy"  ,"hvy"  ,100,0.,.5);
//   TH1D *hvz  = new TH1D("hvz"  ,"hvz"  ,100,-10.,10.);
  TH1D *hd0  = (TH1D *) f->Get("hd0");

  TH2D *hd0vsphi  = (TH2D *) f->Get("hd0vsphi");
  TH2D *hd0vseta  = (TH2D *) f->Get("hd0vseta");
  TH2D *hd0vspt   = (TH2D *) f->Get("hd0vspt" );

  TH1D *hnhpxb  = (TH1D *) f->Get("nhpxb" );
  TH1D *hnhpxe  = (TH1D *) f->Get("nhpxe" );
  TH1D *hnhTIB  = (TH1D *) f->Get("nhTIB" );
  TH1D *hnhTID  = (TH1D *) f->Get("nhTID" );
  TH1D *hnhTOB  = (TH1D *) f->Get("nhTOB" );
  TH1D *hnhTEC  = (TH1D *) f->Get("nhTEC" );

  TH1D *hnhpxbe  = (TH1D *) f->Get("nhpxbe" );
  TH1D *hnhTIBTOB  = (TH1D *) f->Get("nhTIBTOB" );
  TH1D *hnhTIBTOBTECTID= (TH1D *) f->Get("hnhTIBTOBTECTID");

//   TH1D *hnhpxbMuSameCharge  = new TH1D("nhpxbMuSameCharge"  ,"nhpxbMuSameCharge"  ,10,0.,10.);
//   TH1D *hnhpxeMuSameCharge  = new TH1D("nhpxeMuSameCharge"  ,"nhpxeMuSameCharge"  ,10,0.,10.);
//   TH1D *hnhTIBMuSameCharge  = new TH1D("nhTIBMuSameCharge"  ,"nhTIBMuSameCharge"  ,20,0.,20.);
//   TH1D *hnhTIDMuSameCharge  = new TH1D("nhTIDMuSameCharge"  ,"nhTIDMuSameCharge"  ,20,0.,20.);
//   TH1D *hnhTOBMuSameCharge  = new TH1D("nhTOBMuSameCharge"  ,"nhTOBMuSameCharge"  ,20,0.,20.);
//   TH1D *hnhTECMuSameCharge  = new TH1D("nhTECMuSameCharge"  ,"nhTECMuSameCharge"  ,20,0.,20.);
//   TH1D *hEtaMuSameCharge = new TH1D("hEtaMuSameCharge","EtaMuSameCharge",100,-4.,4.);
//   TH1D *hPhiMuSameCharge = new TH1D("hPhiMuSameCharge","PhiMuSameCharge",100,-4.,4.);
//   TH1D *hPtMuSameCharge = new TH1D("hPtMuSameCharge","PtMuSameCharge",100,0.,100.);
//   TH1D *hEtaBothMuP = new TH1D("hEtaBothMuP","EtaBothMuP",100,-4.,4.);
//   TH1D *hPhiBothMuP = new TH1D("hPhiBothMuP","PhiBothMuP",100,-4.,4.);
//   TH1D *hEtaBothMuM = new TH1D("hEtaBothMuM","EtaBothMuM",100,-4.,4.);
//   TH1D *hPhiBothMuM = new TH1D("hPhiBothMuM","PhiBothMuM",100,-4.,4.);

//   TH1I *hMultCand = new TH1I("hMultCand","MultipleCandidate",50,-.5,49.5);
  
  f1->cd();
  
  TH1D *h2chi2 = (TH1D *) f1->Get("hchi2");
  TH1D *h2Ntrk = (TH1D *) f1->Get("hNtrk");
  TH1D *h2Mass = (TH1D *) f1->Get("hMass");
  TH1D *h2P    = (TH1D *) f1->Get("hP");
  TH1D *h2PPlus= (TH1D *) f1->Get("hPPlus");
  TH1D *h2PMinus=(TH1D *) f1->Get("hPMinus");
  TH1D *h2Pt   = (TH1D *) f1->Get("hPt");
  TH1D *h2MinPt= (TH1D *) f1->Get("hMinPt");
  TH1D *h2PtPlus = (TH1D *) f1->Get("hPtPPlus");
  TH1D *h2PtMinus= (TH1D *) f1->Get("hPtMinus");
  TH1D *h2Hit  = (TH1D *) f1->Get("hHit");
  TH1I *h2HitPlus= (TH1I *) f1->Get("hHitPlus");
  TH1I *h2HitMinus=(TH1I *) f1->Get("hHitMinus");
  
  
  TH1D *h2Eta     = (TH1D *) f1->Get("hEta");
  TH1D *h2EtaPlus = (TH1D *) f1->Get("hEtaPlus");
  TH1D *h2EtaMinus= (TH1D *) f1->Get("hEtaMinus");
  TH1D *h2Phi     = (TH1D *) f1->Get("hPhi");
  TH1D *h2PhiPlus = (TH1D *) f1->Get("hPhiPlus");
  TH1D *h2PhiMinus= (TH1D *) f1->Get("hPhiMinus");

  TH1D *h2DeltaPhi= (TH1D *) f1->Get("hDeltaPhi");
  TH1D *h2DeltaEta= (TH1D *) f1->Get("hDeltaEta");
  TH1D *h2DeltaR  = (TH1D *) f1->Get("R");

  TH1D *h2d0  = (TH1D *) f1->Get("hd0");

  TH2D *h2d0vsphi  = (TH2D *) f1->Get("hd0vsphi");
  TH2D *h2d0vseta  = (TH2D *) f1->Get("hd0vseta");
  TH2D *h2d0vspt   = (TH2D *) f1->Get("hd0vspt" );

  TH1D *h2nhpxb  = (TH1D *) f1->Get("nhpxb" );
  TH1D *h2nhpxe  = (TH1D *) f1->Get("nhpxe" );
  TH1D *h2nhTIB  = (TH1D *) f1->Get("nhTIB" );
  TH1D *h2nhTID  = (TH1D *) f1->Get("nhTID" );
  TH1D *h2nhTOB  = (TH1D *) f1->Get("nhTOB" );
  TH1D *h2nhTEC  = (TH1D *) f1->Get("nhTEC" );

  TH1D *h2nhpxbe  = (TH1D *) f1->Get("nhpxbe" );
  TH1D *h2nhTIBTOB  = (TH1D *) f1->Get("nhTIBTOB" );
  TH1D *h2nhTIBTOBTECTID= (TH1D *) f1->Get("hnhTIBTOBTECTID");

  hchi2 ->Scale(1/hchi2 ->Integral());
  h2chi2->Scale(1/h2chi2->Integral());
  hMass ->Scale(1/hMass ->Integral());
  h2Mass->Scale(1/h2Mass->Integral());
  hP ->Scale(1/hP ->Integral());
  h2P->Scale(1/h2P->Integral());
  hPt ->Scale(1/hPt ->Integral());
  h2Pt->Scale(1/h2Pt->Integral());
  hHit ->Scale(1/hHit ->Integral());
  h2Hit->Scale(1/h2Hit->Integral());
  hEta ->Scale(1/hEta ->Integral());
  h2Eta->Scale(1/h2Eta->Integral());
  hPhi ->Scale(1/hPhi ->Integral());
  h2Phi->Scale(1/h2Phi->Integral());
  hd0 ->Scale(1/hd0 ->Integral());
  h2d0->Scale(1/h2d0->Integral());
  hDeltaPhi ->Scale(1/hDeltaPhi ->Integral());
  h2DeltaPhi->Scale(1/h2DeltaPhi->Integral());
  hDeltaEta ->Scale(1/hDeltaEta ->Integral());
  h2DeltaEta->Scale(1/h2DeltaEta->Integral());
  hDeltaR ->Scale(1/hDeltaR ->Integral());
  h2DeltaR->Scale(1/h2DeltaR->Integral());

  outfile->cd();

  TCanvas *c1= new TCanvas("c1","c1", 1000,700);
  c1->Draw();
  c1->Divide(2,2);
  c1->cd(1);
  hchi2->Draw();
  hchi2->SetFillColor(0);
  hchi2->SetLineColor(1);
  hchi2->SetLineWidth(2);
  h2chi2->Draw("SAMES");
  h2chi2->SetFillColor(0);
  h2chi2->SetLineColor(2);
  h2chi2->SetLineWidth(2);
  c1->cd(2);
  hMass->Draw();
  hMass->SetFillColor(0);
  hMass->SetLineColor(1);
  hMass->SetLineWidth(2);
  h2Mass->Draw("SAMES");
  h2Mass->SetFillColor(0);
  h2Mass->SetLineColor(2);
  h2Mass->SetLineWidth(2);
  c1->cd(3);
  hP->Draw();
  hP->SetFillColor(0);
  hP->SetLineColor(1);
  hP->SetLineWidth(2);
  h2P->Draw("SAMES");
  h2P->SetFillColor(0);
  h2P->SetLineColor(2);
  h2P->SetLineWidth(2);
  c1->cd(4);
  hPt->Draw();
  hPt->SetFillColor(0);
  hPt->SetLineColor(1);
  hPt->SetLineWidth(2);
  h2Pt->Draw("SAMES");
  h2Pt->SetFillColor(0);
  h2Pt->SetLineColor(2);
  h2Pt->SetLineWidth(2);

  TCanvas *c2= new TCanvas("c2","c2", 1000,700);
  c2->Draw();
  c2->Divide(2,2);
  c2->cd(1);
  hHit->Draw();
  hHit->SetFillColor(0);
  hHit->SetLineColor(1);
  hHit->SetLineWidth(2);
  h2Hit->Draw("SAMES");
  h2Hit->SetFillColor(0);
  h2Hit->SetLineColor(2);
  h2Hit->SetLineWidth(2);
  c2->cd(2);
  hEta->Draw();
  hEta->SetFillColor(0);
  hEta->SetLineColor(1);
  hEta->SetLineWidth(2);
  h2Eta->Draw("SAMES");
  h2Eta->SetFillColor(0);
  h2Eta->SetLineColor(2);
  h2Eta->SetLineWidth(2);
  c2->cd(3);
  hPhi->Draw();
  hPhi->SetFillColor(0);
  hPhi->SetLineColor(1);
  hPhi->SetLineWidth(2);
  h2Phi->Draw("SAMES");
  h2Phi->SetFillColor(0);
  h2Phi->SetLineColor(2);
  h2Phi->SetLineWidth(2);
  c2->cd(4);
  hd0->Draw();
  hd0->GetXaxis()->SetTitle("cm");
  hd0->SetFillColor(0);
  hd0->SetLineColor(1);
  hd0->SetLineWidth(2);
  h2d0->Draw("SAMES");
  h2d0->SetFillColor(0);
  h2d0->SetLineColor(2);
  h2d0->SetLineWidth(2);
  TCanvas *c3= new TCanvas("c3","c3", 1000,700);
  c3->Draw();
  c3->Divide(2,2);
  c3->cd(1);
  h2d0vsphi->Draw();
  h2d0vsphi->SetFillColor(0);
  h2d0vsphi->SetMarkerColor(2);
  h2d0vsphi->SetLineColor(1);
  h2d0vsphi->SetLineWidth(2);
  hd0vsphi->Draw("SAMES");
  hd0vsphi->SetFillColor(0);
  hd0vsphi->SetLineColor(2);
  hd0vsphi->SetLineWidth(2);
  c3->cd(2);
  h2d0vseta->Draw();
  h2d0vseta->SetFillColor(0);
  h2d0vseta->SetMarkerColor(2);
  h2d0vseta->SetLineColor(1);
  h2d0vseta->SetLineWidth(2);
  hd0vseta->Draw("SAMES");
  hd0vseta->SetFillColor(0);
  hd0vseta->SetLineColor(2);
  hd0vseta->SetLineWidth(2);
  c3->cd(3);
  h2d0vspt->Draw();
  h2d0vspt->SetFillColor(0);
  h2d0vspt->SetMarkerColor(2);
  h2d0vspt->SetLineColor(1);
  h2d0vspt->SetLineWidth(2);
  hd0vspt->Draw("SAMES");
  hd0vspt->SetFillColor(0);
  hd0vspt->SetLineColor(2);
  hd0vspt->SetLineWidth(2);

  TCanvas *c4= new TCanvas("c4","c4", 1000,700);
  c4->Draw();
  c4->Divide(2,2);
  c4->cd(1);
  hDeltaEta->Draw();
  hDeltaEta->SetFillColor(0);
  hDeltaEta->SetLineColor(1);
  hDeltaEta->SetLineWidth(2);
  h2DeltaEta->Draw("SAMES");
  h2DeltaEta->SetFillColor(0);
  h2DeltaEta->SetLineColor(2);
  h2DeltaEta->SetLineWidth(2);
  c4->cd(2);
  hDeltaPhi->Draw();
  hDeltaPhi->SetFillColor(0);
  hDeltaPhi->SetLineColor(1);
  hDeltaPhi->SetLineWidth(2);
  h2DeltaPhi->Draw("SAMES");
  h2DeltaPhi->SetFillColor(0);
  h2DeltaPhi->SetLineColor(2);
  h2DeltaPhi->SetLineWidth(2);
  c4->cd(3);
  hDeltaR->Draw();
  hDeltaR->SetFillColor(0);
  hDeltaR->SetLineColor(1);
  hDeltaR->SetLineWidth(2);
  h2DeltaR->Draw("SAMES");
  h2DeltaR->SetFillColor(0);
  h2DeltaR->SetLineColor(2);
  h2DeltaR->SetLineWidth(2);

//    f->Close();
//    f2->Close();

//    outfile->Close();



 /////////////////////////////////////////////////////////////////////////////


}
