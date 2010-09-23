#include <TROOT.h>
#include <TStyle.h>
#include <TVector.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TKey.h>
#include <TList.h>
#include <TLegend.h>
#include <iostream>
using namespace std;
void readIso () {
  //gROOT->Reset();
  gROOT->SetStyle("Plain");
  //gStyle->SetOptFit(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptLogy();
  
  TFile * f1 = TFile::Open("rfio:/castor/cern.ch/user/t/taroni/HtoWWOutput/HtoWWGenTree.root");
  if (!f1){
    cout << "FILE NOT FOUND" << endl;
    return;
  }

  TTree *t1 = (TTree*)f1->FindObjectAny("ntuple");
  if (!t1){
    cout << "Tree NOT FOUND" << endl;
    return;
  }
  
  std::vector<float> *muPt=0, *muP=0, *muEta=0, *muPhi=0, *muPx=0, *muPy=0, *muPz=0, *muE=0, 
    *muIsoTRK03=0, *muIsoECAL03=0, *muIsoHCAL03=0, *muIsoHO03=0, *muIsoTRK05=0, *muIsoECAL05=0, *muIsoHCAL05=0, *muIsoHO05=0;
  std::vector<int> *muQ=0;
  std::vector<double> *muConeChPt=0, *muConeChP=0, *muConeChEt=0, *muConeChE=0, *muConeChNPart=0, *muConeNePt=0, *muConeNeP=0, *muConeNeEt=0, *muConeNeE=0, *muConeNeNPart=0;

  TBranch * muPtB    = (TBranch*) t1-> GetBranch("muPt");
  TBranch * muPB     = (TBranch*) t1-> GetBranch("muP");
  TBranch * muPxB    = (TBranch*) t1-> GetBranch("muPx");
  TBranch * muPyB    = (TBranch*) t1-> GetBranch("muPy");
  TBranch * muPzB    = (TBranch*) t1-> GetBranch("muPz");
  TBranch * muEB     = (TBranch*) t1-> GetBranch("muE");
  TBranch * muEtaB   = (TBranch*) t1-> GetBranch("muEta");
  TBranch * muPhiB   = (TBranch*) t1-> GetBranch("muPhi");

  TBranch * muQB   = (TBranch*) t1-> GetBranch("muQ");

  muPtB ->SetAddress(&muPt);
  muPB  ->SetAddress(&muP);
  muPxB ->SetAddress(&muPx);
  muPyB ->SetAddress(&muPy);
  muPzB ->SetAddress(&muPz);
  muEB  ->SetAddress(&muE);
  muEtaB->SetAddress(&muEta);
  muPhiB->SetAddress(&muPhi);

  muQB->SetAddress(&muQ);

  TBranch * muConeChPtB    = (TBranch*) t1-> GetBranch("muConeChPt");
  TBranch * muConeChPB     = (TBranch*) t1-> GetBranch("muConeChP");
  TBranch * muConeChEtB    = (TBranch*) t1-> GetBranch("muConeChEt");
  TBranch * muConeChEB     = (TBranch*) t1-> GetBranch("muConeChE");
  TBranch * muConeChNPartB = (TBranch*) t1-> GetBranch("muConeChNPart");
  TBranch * muConeNePtB    = (TBranch*) t1-> GetBranch("muConeNePt");
  TBranch * muConeNePB     = (TBranch*) t1-> GetBranch("muConeNeP");
  TBranch * muConeNeEtB    = (TBranch*) t1-> GetBranch("muConeNeEt");
  TBranch * muConeNeEB     = (TBranch*) t1-> GetBranch("muConeNeE");
  TBranch * muConeNeNPartB = (TBranch*) t1-> GetBranch("muConeNeNPart");

  muConeChPtB   -> SetAddress(&(muConeChPt));
  muConeChPB    -> SetAddress(&(muConeChP));
  muConeChEtB   -> SetAddress(&(muConeChEt));
  muConeChEB    -> SetAddress(&(muConeChE));
  muConeChNPartB-> SetAddress(&(muConeChNPart));
  muConeNePtB   -> SetAddress(&(muConeNePt));
  muConeNePB    -> SetAddress(&(muConeNeP));
  muConeNeEtB   -> SetAddress(&(muConeNeEt));
  muConeNeEB    -> SetAddress(&(muConeNeE));
  muConeNeNPartB-> SetAddress(&(muConeNeNPart));


  int n = (int) t1->GetEntries();
  cout << " Tree Entries " << n << endl; 

  vector<double> doubleMuPt (1000,0);
  vector<double> singleMuPt (1000,0);
//   for (int i=0; i<1000; i++){
//     doubleMuPt[i].assign(0);
//     singleMuPt[i].assign(0);
//   }

  TH1F * singlemuRChPt   =new TH1F ("singlemuRChPt","singlemuRChPt",500,0,500);
  TH1F * singlemuRChP    =new TH1F ("singlemuRChP" ,"singlemuRChP" ,500,0,500);
  TH1F * singlemuRChE    =new TH1F ("singlemuRChE" ,"singlemuRChE" ,500,0,500);
  TH1F * singlemuRChEt   =new TH1F ("singlemuRChEt","singlemuRChEt",500,0,500);
  TH1F * singlemuRChNPart=new TH1F ("singlemuRChNPart","singlemuRChNPart",500,0,500);

  TH1F * singlemuRNePt   =new TH1F ("singlemuRNePt","singlemuRNePt",500,0,500);
  TH1F * singlemuRNeP    =new TH1F ("singlemuRNeP" ,"singlemuRNeP" ,500,0,500);
  TH1F * singlemuRNeE    =new TH1F ("singlemuRNeE" ,"singlemuRNeE" ,500,0,500);
  TH1F * singlemuRNeEt   =new TH1F ("singlemuRNeEt","singlemuRNeEt",500,0,500);
  TH1F * singlemuRNeNPart=new TH1F ("singlemuRNeNPart","singlemuRNeNPart",500,0,500);

  TH1F * doublemuRChPt=new TH1F ("doublemuRChPt","doublemuRChPt",500,0,500);
  TH1F * doublemuRChP =new TH1F ("doublemuRChP" ,"doublemuRChP" ,500,0,500);
  TH1F * doublemuRChE =new TH1F ("doublemuRChE" ,"doublemuRChE" ,500,0,500);
  TH1F * doublemuRChEt=new TH1F ("doublemuRChEt","doublemuRChEt",500,0,500);
  TH1F * doublemuRChNPart=new TH1F ("doublemuRChNPart","doublemuRChNPart",500,0,500);

  TH1F * doublemuRNePt=new TH1F ("doublemuRNePt","doublemuRNePt",500,0,500);
  TH1F * doublemuRNeP =new TH1F ("doublemuRNeP" ,"doublemuRNeP" ,500,0,500);
  TH1F * doublemuRNeE =new TH1F ("doublemuRNeE" ,"doublemuRNeE" ,500,0,500);
  TH1F * doublemuRNeEt=new TH1F ("doublemuRNeEt","doublemuRNeEt",500,0,500);
  TH1F * doublemuRNeNPart=new TH1F ("doublemuRNeNPart","doublemuRNeNPart",500,0,500);

   for (int evt = 0; evt < n; ++evt) { //event loop
 
     t1->GetEntry(evt);
//      muConeChPtB   ->GetEntry(evt);
//      muConeChPB    ->GetEntry(evt);
//      muConeChEtB   ->GetEntry(evt);
//      muConeChEB    ->GetEntry(evt);
//      muConeChNPartB   ->GetEntry(evt);
//      muConeNePtB ->GetEntry(evt);
//      muConeNePB  ->GetEntry(evt);
//      muConeNeEtB ->GetEntry(evt);
//      muConeNeEB  ->GetEntry(evt);
//      muConeNeNPartB  ->GetEntry(evt);
          
     
     for (int imu = 0; imu <2; imu++){
       doublemuRChPt->Fill(muConeChPt->at(imu));
       doublemuRChP ->Fill(muConeChP ->at(imu));
       doublemuRChE ->Fill(muConeChE ->at(imu));
       doublemuRChEt->Fill(muConeChEt->at(imu));
       doublemuRChNPart->Fill(muConeChNPart->at(imu));
       doublemuRNePt->Fill(muConeNePt->at(imu));
       doublemuRNeP ->Fill(muConeNeP ->at(imu));
       doublemuRNeE ->Fill(muConeNeE ->at(imu));
       doublemuRNeEt->Fill(muConeNeEt->at(imu));
       doublemuRNeNPart->Fill(muConeNeNPart->at(imu));
       
       
     } 
      

     

   }//evt
   TFile* outfile = new TFile("ISO.root","recreate");
   outfile -> cd();
   
   doublemuRChPt->Write();
   doublemuRChP ->Write();
   doublemuRChE ->Write();
   doublemuRChEt->Write();
   doublemuRChNPart->Write();
   doublemuRNePt->Write();
   doublemuRNeP ->Write();
   doublemuRNeE ->Write();
   doublemuRNeEt->Write();
   doublemuRNeNPart->Write();
   

   outfile->Close();
 //   c2->SaveAs("effSingle.eps");
   TCanvas *c2 = new TCanvas();
   c2 -> Divide(2,2);
   c2->cd(1);
   doublemuRChPt->Draw();
   doublemuRChPt->GetXaxis()->SetRangeUser(0,50);
   c2->cd(2);
   doublemuRChEt->Draw();
   doublemuRChEt->GetXaxis()->SetRangeUser(0,50);
   c2->cd(3);
   doublemuRNePt->Draw();
   doublemuRNePt->GetXaxis()->SetRangeUser(0,50);
   c2->cd(4);
   doublemuRNeEt->Draw();
   doublemuRNeEt->GetXaxis()->SetRangeUser(0,50);
   
   c2->Update();
   c2->SaveAs("MuPtEt.eps");

   c2->cd(1);
   doublemuRChP->Draw();
   doublemuRChP->GetXaxis()->SetRangeUser(0,50);
   c2->cd(2);
   doublemuRChE->Draw();
   doublemuRChE->GetXaxis()->SetRangeUser(0,50);
   c2->cd(3);
   doublemuRNeP->Draw();
   doublemuRNeP->GetXaxis()->SetRangeUser(0,50);
   c2->cd(4);
   doublemuRNeE->Draw();
   doublemuRNeE->GetXaxis()->SetRangeUser(0,50);
   c2->Update();
   c2->SaveAs("MuPE.eps");
   
   c2->Clear();
   c2->Divide(2,2);
   c2->cd(1);
   doublemuRChNPart->Draw();
   doublemuRChNPart->GetXaxis()->SetRangeUser(0,50);
   c2->cd(2);
   doublemuRNeNPart->Draw();
   doublemuRNeNPart->GetXaxis()->SetRangeUser(0,50);
   c2->Update();
   c2->SaveAs("MuNPart.eps");



   f1->Close();


}
