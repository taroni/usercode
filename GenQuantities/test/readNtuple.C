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
void readNtuple () {
  //gROOT->Reset();
  gROOT->SetStyle("Plain");
  //gStyle->SetOptFit(1);
  gStyle->SetHistLineWidth(2);
  
  
//   TFile * f1 = TFile::Open("/data06/users/taroni/trigger/HtoWW/MuonsTree.root");
  TFile * f1 = TFile::Open("HtoWWGenTree.root");

  TTree *t1 = (TTree*)f1->FindObjectAny("ntuple");
  cout << "Ntuple read" << endl;

//   int run, event;
//   TBranch * runB   = (TBranch*) t1-> GetBranch("run");
//   TBranch * eventB = (TBranch*) t1-> GetBranch("event");
//   t1->SetBranchAddress("run",&run);
//   t1->SetBranchAddress("event",&event);

  
  std::vector<float> *muPt=0, *muP=0, *muEta=0, *muPhi=0, *muPx=0, *muPy=0, *muPz=0, *muE=0, 
    *muIsoTRK03=0, *muIsoECAL03=0, *muIsoHCAL03=0, *muIsoHO03=0, *muIsoTRK05=0, *muIsoECAL05=0, *muIsoHCAL05=0, *muIsoHO05=0;
  std::vector<int> *muQ=0;
  std::vector<double> *muConeChPt=0, *muConeChP=0, *muConeChEt=0, *muConeChE=0, *muConeChNPart=0, *muConeNePt=0, *muConeNeP=0, *muConeNeEt=0, *muConeNeE=0, *muConeNeNPart=0;

  TBranch * muPtB   = (TBranch*) t1-> GetBranch("muPt");
  TBranch * muPB    = (TBranch*) t1-> GetBranch("muP");
  TBranch * muPxB   = (TBranch*) t1-> GetBranch("muPx");
  TBranch * muPyB   = (TBranch*) t1-> GetBranch("muPy");
  TBranch * muPzB   = (TBranch*) t1-> GetBranch("muPz");
  TBranch * muEB    = (TBranch*) t1-> GetBranch("muE");
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

  muQB  ->SetAddress(&muQ);

  TBranch * muConeChPtB   = (TBranch*) t1-> GetBranch("muConeChPt");
  TBranch * muConeChPB    = (TBranch*) t1-> GetBranch("muConeChP");
  TBranch * muConeChEtB   = (TBranch*) t1-> GetBranch("muConeChEt");
  TBranch * muConeChEB    = (TBranch*) t1-> GetBranch("muConeChE");
  TBranch * muConeChNPartB   = (TBranch*) t1-> GetBranch("muConeChNPart");
  TBranch * muConeNePtB   = (TBranch*) t1-> GetBranch("muConeNePt");
  TBranch * muConeNePB    = (TBranch*) t1-> GetBranch("muConeNeP");
  TBranch * muConeNeEtB   = (TBranch*) t1-> GetBranch("muConeNeEt");
  TBranch * muConeNeEB    = (TBranch*) t1-> GetBranch("muConeNeE");
  TBranch * muConeNeNPartB   = (TBranch*) t1-> GetBranch("muConeNeNPart");

   muConeChPtB   -> SetAddress(&(muConeChPt));
   muConeChPB	 -> SetAddress(&(muConeChP));
   muConeChEtB   -> SetAddress(&(muConeChEt));
   muConeChEB	 -> SetAddress(&(muConeChE));
   muConeChNPartB-> SetAddress(&(muConeChNPart));
   muConeNePtB   -> SetAddress(&(muConeNePt));
   muConeNePB	 -> SetAddress(&(muConeNeP));
   muConeNeEtB   -> SetAddress(&(muConeNeEt));
   muConeNeEB	 -> SetAddress(&(muConeNeE));
   muConeNeNPartB-> SetAddress(&(muConeNeNPart));

  int n = (int) t1->GetEntries();
  cout << " Tree Entries " << n << endl; 

  vector<double> doubleMuPt (1000,0);
  vector<double> singleMuPt (1000,0);
//   for (int i=0; i<1000; i++){
//     doubleMuPt[i].assign(0);
//     singleMuPt[i].assign(0);
//   }

  TH2F * mu2Dpt = new TH2F ("mu2Dpt","mu2Dpt",500,0, 500, 500, 0, 500);

  TH2F * singlemuRChPt=new TH2F ("singlemuRChPt","singlemuRChPt",1000,0,1000,500,0,500);
  TH2F * singlemuRChP =new TH2F ("singlemuRChP" ,"singlemuRChP" ,1000,0,1000,500,0,500);
  TH2F * singlemuRChE =new TH2F ("singlemuRChE" ,"singlemuRChE" ,1000,0,1000,500,0,500);
  TH2F * singlemuRChEt=new TH2F ("singlemuRChEt","singlemuRChEt",1000,0,1000,500,0,500);
  TH2F * singlemuRChNPart=new TH2F ("singlemuRChNPart","singlemuRChNPart",1000,0,1000,500,0,500);

  TH2F * singlemuRNePt=new TH2F ("singlemuRNePt","singlemuRNePt",1000,0,1000,500,0,500);
  TH2F * singlemuRNeP =new TH2F ("singlemuRNeP" ,"singlemuRNeP" ,1000,0,1000,500,0,500);
  TH2F * singlemuRNeE =new TH2F ("singlemuRNeE" ,"singlemuRNeE" ,1000,0,1000,500,0,500);
  TH2F * singlemuRNeEt=new TH2F ("singlemuRNeEt","singlemuRNeEt",1000,0,1000,500,0,500);
  TH2F * singlemuRNeNPart=new TH2F ("singlemuRNeNPart","singlemuRNeNPart",1000,0,1000,500,0,500);

  TH2F * doublemuRChPt=new TH2F ("doublemuRChPt","doublemuRChPt",1000,0,1000,500,0,500);
  TH2F * doublemuRChP =new TH2F ("doublemuRChP" ,"doublemuRChP" ,1000,0,1000,500,0,500);
  TH2F * doublemuRChE =new TH2F ("doublemuRChE" ,"doublemuRChE" ,1000,0,1000,500,0,500);
  TH2F * doublemuRChEt=new TH2F ("doublemuRChEt","doublemuRChEt",1000,0,1000,500,0,500);
  TH2F * doublemuRChNPart=new TH2F ("doublemuRChNPart","doublemuRChNPart",1000,0,1000,500,0,500);

  TH2F * doublemuRNePt=new TH2F ("doublemuRNePt","doublemuRNePt",1000,0,1000,500,0,500);
  TH2F * doublemuRNeP =new TH2F ("doublemuRNeP" ,"doublemuRNeP" ,1000,0,1000,500,0,500);
  TH2F * doublemuRNeE =new TH2F ("doublemuRNeE" ,"doublemuRNeE" ,1000,0,1000,500,0,500);
  TH2F * doublemuRNeEt=new TH2F ("doublemuRNeEt","doublemuRNeEt",1000,0,1000,500,0,500);
  TH2F * doublemuRNeNPart=new TH2F ("doublemuRNeNPart","doublemuRNeNPart",1000,0,1000,500,0,500);

  TH1F * ptRel = new TH1F ("ptChRel","ptChRel",100, 0,50);
  TH1F * EtChRel = new TH1F ("EtChRel","EtChRel",100, 0,50);
  TH1F * EtNeRel = new TH1F ("EtNeRel","EtNeRel",100, 0,50);
  TH1F * EtRel = new TH1F ("EtTotRel","EtTotRel",100, 0,50);

  TH2F * ptRelVsPt = new TH2F ("ptChRelvsPt","ptChRelvsPt",100, 0,50,100,0,500);
  TH2F * EtChRelVsPt = new TH2F ("EtChRelvsPt","EtChRelvsPt",100, 0,50,100,0,500);
  TH2F * EtNeRelVsPt = new TH2F ("EtNeRelvsPt","EtNeRelvsPt",100, 0,50,100,0,500);
  TH2F * EtRelVsPt = new TH2F ("EtTotRelvsPt","EtTotRelvsPt",100, 0,50,100,0,500);
  

   for (int evt = 0; evt < n; ++evt) { //event loop
//    runB->GetEntry(evt); //fill the variables with the i-th event
//    eventB->GetEntry(evt); 
//      cout << "run=" << run 
// 	  << " evt=" << event 
// 	  << endl;
     t1->GetEntry(evt);	

//   muPtB->GetEntry(evt);
//   muConeChPtB   ->GetEntry(evt);
//   muConeChPB    ->GetEntry(evt);
//   muConeChEtB   ->GetEntry(evt);
//   muConeChEB    ->GetEntry(evt);
//   muConeChNPartB   ->GetEntry(evt);
//   muConeNePtB ->GetEntry(evt);
//   muConeNePB  ->GetEntry(evt);
//   muConeNeEtB ->GetEntry(evt);
//   muConeNeEB  ->GetEntry(evt);
//   muConeNeNPartB  ->GetEntry(evt);
     
//      cout << __LINE__ << " " << muConeChPt->at(0)<< " "<<muConeChP->at(0)
// // 	  << " "<<muConeChE->at(0)<< " "<<muConeChEt->at(0) << " " << muConeChNPart->at(0)
// 	  << endl;
     
     
     for (int i=0; i<1000; i++){
       if (muPt->at(0) > i && muPt->at(1) > i) {
	 doubleMuPt[i]++;
	 for (int imu = 0; imu <2; imu++){
	   doublemuRChPt->Fill(i,muConeChPt->at(imu));
	   doublemuRChP ->Fill(i,muConeChP ->at(imu));
	   doublemuRChE ->Fill(i,muConeChE ->at(imu));
	   doublemuRChEt->Fill(i,muConeChEt->at(imu));
	   doublemuRChNPart->Fill(i,muConeChNPart->at(imu));
	   doublemuRNePt->Fill(i,muConeNePt->at(imu));
	   doublemuRNeP ->Fill(i,muConeNeP ->at(imu));
	   doublemuRNeE ->Fill(i,muConeNeE ->at(imu));
	   doublemuRNeEt->Fill(i,muConeNeEt->at(imu));
	   doublemuRNeNPart->Fill(i,muConeNeNPart->at(imu));
 	 }
	 //2Dpt
	 mu2Dpt -> Fill(muPt->at(0),muPt->at(1));
       } 
       if  (muPt->at(0) > i || muPt->at(1) > i) {
	 singleMuPt[i]++;
	 for (int imu = 0; imu <2; imu++){
	   singlemuRChPt->Fill(i,muConeChPt->at(imu));
	   singlemuRChP ->Fill(i,muConeChP ->at(imu));
	   singlemuRChE ->Fill(i,muConeChE ->at(imu));
	   singlemuRChEt->Fill(i,muConeChEt->at(imu));
	   singlemuRChNPart->Fill(i,muConeChNPart->at(imu));
	   singlemuRNePt->Fill(i,muConeNePt->at(imu));
	   singlemuRNeP ->Fill(i,muConeNeP ->at(imu));
	   singlemuRNeE ->Fill(i,muConeNeE ->at(imu));
	   singlemuRNeEt->Fill(i,muConeNeEt->at(imu));
	   singlemuRNeNPart->Fill(i,muConeNeNPart->at(imu));
 	 }

       }
     }


     for (int imu=0; imu<2; imu++){
//        cout << evt << " "<<  muConeChPt->at(imu) << " " << muPt->at(imu) << " " << muConeChPt->at(imu)/muPt->at(imu) << endl;
       ptRel ->Fill(muConeChPt->at(imu)/muPt->at(imu));

       EtRel ->Fill((muConeChEt->at(imu)+muConeNeEt->at(imu))/muPt->at(imu));
       EtChRel ->Fill(muConeChEt->at(imu)/muPt->at(imu));
       EtNeRel ->Fill(muConeNeEt->at(imu)/muPt->at(imu));

       ptRelVsPt -> Fill (muConeChPt->at(imu)/muPt->at(imu),muPt->at(imu));
       EtRelVsPt ->Fill((muConeChEt->at(imu)+muConeNeEt->at(imu))/muPt->at(imu),muPt->at(imu));
       EtChRelVsPt ->Fill(muConeChEt->at(imu)/muPt->at(imu),muPt->at(imu));
       EtNeRelVsPt ->Fill(muConeNeEt->at(imu)/muPt->at(imu),muPt->at(imu));
 

     }



     
   }//evt
//    for (int i=0; i<10; i++){
//      cout <<  doubleMuPt[i] << " " << singleMuPt[i] << endl;
//    }
   TH1F * effDoubleHisto = new TH1F("effDouble", "effDouble", 1000,0,1000);
   TH1F * effSingleHisto = new TH1F("effSingle", "effSingle", 1000,0,1000);

   for (int i = 0; i < 1000; i++){
     double effDouble =doubleMuPt[i]/((double) n);
     effDoubleHisto -> Fill(i,effDouble);
     double effSingle =singleMuPt[i]/((double) n);
     effSingleHisto -> Fill(i,effSingle);
   }
   TCanvas * c1 = new TCanvas ("c1","c1");
   TCanvas * c2 = new TCanvas ("c2","c2");
   TCanvas * c3 = new TCanvas ("c3","c3");
   TCanvas * c4 = new TCanvas ("c4","c4");
   TCanvas * c5 = new TCanvas ("c5","c5");

   c1->Draw();
   c1->cd();
   c1->SetLogy();
   effDoubleHisto->Draw();
   effDoubleHisto->GetXaxis()->SetTitle("GeV");
   effDoubleHisto->GetXaxis()->SetRangeUser(0,100);
   effSingleHisto->Draw("SAMES");
   effSingleHisto->SetLineColor(2);
   TLegend* legend = new TLegend(0.6,0.6,0.85,0.8, "","brNDC");
   legend ->AddEntry(effDoubleHisto, "Double Mu");
   legend ->AddEntry(effSingleHisto, "Single Mu");
   legend -> Draw();
   legend->SetFillColor(0);

   c1->SaveAs("plots/efficiency.eps");
   

   c2->Draw();
   c2->Divide(1,2);
   c2->cd(1);
   c2->GetPad(1)->SetLogy();
   ptRel->Draw();
//    ptRel->GetXaxis()->SetRangeUser(0,10);
   c2->cd(2);
   c2->GetPad(2)->SetLogy();
   EtRel->Draw();
//    EtRel->GetXaxis()->SetRangeUser(0,10);
   c2->Update();
   c2->SaveAs("plots/PtEtRel.eps");

   c3->Draw();
   c3->Divide(1,2);
   c3->cd(1);
   c3->GetPad(1)->SetLogy();
   EtChRel->Draw();
//    EtChRel->GetXaxis()->SetRangeUser(0,10);
   c3->cd(2);
   c3->GetPad(2)->SetLogy();
   EtNeRel->Draw();
//    EtNeRel->GetXaxis()->SetRangeUser(0,10);
   c3->Update();
   c3->SaveAs("plots/EtRel.eps");

   c4->Draw();
   c4->Divide(1,2);
   c4->cd(1);
//    c4->GetPad(1)->SetLogy();
   ptRelVsPt->Draw();
   c4->cd(2);
//    c4->GetPad(2)->SetLogy();
   EtRelVsPt->Draw();
   c4->SaveAs("plots/PtEtRelvsPt.eps");

   c5->Draw();
   c5->Divide(1,2);
   c5->cd(1);
//    c5->GetPad(1)->SetLogy();
   EtChRelVsPt->Draw();
   c5->cd(2);
//    c5->GetPad(2)->SetLogy();
   EtNeRelVsPt->Draw();
   c5->SaveAs("plots/EtRelvsPt.eps");

   
   
//   c2->Draw();
//   c2->cd();
   TFile* outfile = new TFile("ISOvsTrigger1M.root","recreate");
   outfile -> cd();

   mu2Dpt->Write();
   
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
   
   singlemuRChPt->Write();
   singlemuRChP ->Write();
   singlemuRChE ->Write();
   singlemuRChEt->Write();
   singlemuRChNPart->Write();
   singlemuRNePt->Write();
   singlemuRNeP ->Write();
   singlemuRNeE ->Write();
   singlemuRNeEt->Write();
   singlemuRNeNPart->Write();


   EtRel   ->Write();
   EtChRel ->Write();
   EtNeRel ->Write();
   
   ptRelVsPt   ->Write();
   EtRelVsPt   ->Write();
   EtChRelVsPt ->Write();
   EtNeRelVsPt ->Write();
   


   outfile->Close();
 //   c2->SaveAs("effSingle.eps");
   f1->Close();


}
