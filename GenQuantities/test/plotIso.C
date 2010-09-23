#include <TROOT.h>
#include <TStyle.h>
#include <TVector.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TTree.h>
#include <TKey.h>
#include <TList.h>
#include <TLegend.h>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;
void plotIso () {
  //gROOT->Reset();
  gROOT->SetStyle("Plain");
  //gStyle->SetOptFit(1);
  gStyle->SetHistLineWidth(2);
  
  TFile * f1 = TFile::Open("ISOvsTrigger1M.root");

  TH2F * doublemuRChPt  = (TH2F *) f1->Get("doublemuRChPt");
  TH2F * doublemuRChP   = (TH2F *) f1->Get("doublemuRChP");
  TH2F * doublemuRChE = (TH2F *)f1->Get("doublemuRChE");
  TH2F * doublemuRChEt = (TH2F *)f1->Get("doublemuRChEt");
  TH2F * doublemuRChNPart = (TH2F *)f1->Get("doublemuRNeNPart");
  TH2F * doublemuRNeNPart = (TH2F *)f1->Get("doublemuRNeNPart");
  TH2F * doublemuRNePt = (TH2F *)f1->Get("doublemuRNePt");
  TH2F * doublemuRNeP = (TH2F *)f1->Get("doublemuRNeP");
  TH2F * doublemuRNeE = (TH2F *)f1->Get("doublemuRNeE");
  TH2F * doublemuRNeEt = (TH2F *)f1->Get("doublemuRNeEt");
  TH2F * singlemuRChPt = (TH2F *)f1->Get("singlemuRChPt");
  TH2F * singlemuRChP = (TH2F *)f1->Get("singlemuRChP");
  TH2F * singlemuRChE = (TH2F *)f1->Get("singlemuRChE");
  TH2F * singlemuRChEt = (TH2F *)f1->Get("singlemuRChEt");
  TH2F * singlemuRChNPart= (TH2F *)f1->Get("singlemuRNeNPart");
  TH2F * singlemuRNeNPart= (TH2F *)f1->Get("singlemuRNeNPart");
  TH2F * singlemuRNePt = (TH2F *)f1->Get("singlemuRNePt");
  TH2F * singlemuRNeP = (TH2F *)f1->Get("singlemuRNeP");
  TH2F * singlemuRNeE = (TH2F *)f1->Get("singlemuRNeE");
  TH2F * singlemuRNeEt = (TH2F *)f1->Get("singlemuRNeEt");
  
  
//  vector <TH1D *> singleH, doubleH;
  TCanvas *c1 = new TCanvas("c1","ChargedPtEt",1200,800);
  TCanvas *c2 = new TCanvas("c2","NeutralPtEt",1200,800);
  TCanvas *c3 = new TCanvas("c3","c3",1200,800);
  TCanvas *c4 = new TCanvas("c4","c4",1200,800);
  TCanvas *c5 = new TCanvas("c5","c5",1200,800);

  
  c1->Draw();
  c1->Divide(2,2);
  c2->Draw();
  c2->Divide(2,2);
  c3->Draw();
  c3->Divide(2,2);
  c4->Draw();
  c4->Divide(2,2);
  c5->Draw();
  c5->Divide(2,2);

  stringstream plotName;
  stringstream binName;
  

  
  for (int ibin = 0; ibin<6; ibin++){
     int newbin= 1+ibin*5;
     cout << newbin << " " << ibin << endl;
     binName.str("");
     binName<< "doubleMuChPt" << newbin;
     TH1D * dh1 = doublemuRChPt->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "doubleMuChP" << newbin;     
     TH1D * dh2 = doublemuRChP ->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "doubleMuChEt" << newbin;          
     TH1D * dh3 = doublemuRChEt->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "doubleMuChE" << newbin;          
     TH1D * dh4 = doublemuRChE ->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "doubleMuChNPart" << newbin;          
     TH1D * dh5 = doublemuRChNPart->ProjectionY(binName.str().c_str(),newbin,newbin);
     // dh1 ->SetName(doublemuRChPt->GetName());
     dh1 ->SetTitle(doublemuRChPt->GetName());
     // dh2 ->SetName(doublemuRChP->GetName());
     dh2 ->SetTitle(doublemuRChP->GetName());
     // dh3 ->SetName(doublemuRChEt->GetName());
     dh3 ->SetTitle(doublemuRChEt->GetName());
     // dh4 ->SetName(doublemuRChE->GetName());
     dh4 ->SetTitle(doublemuRChE->GetName());
     // dh5 ->SetName(doublemuRChNPart->GetName());
     dh5 ->SetTitle(doublemuRChNPart->GetName());

     binName.str("");
     binName<< "singleMuChPt" << newbin;
     TH1D * sh1 = singlemuRChPt->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "singleMuChP" << newbin;
     TH1D * sh2 = singlemuRChP ->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "singleMuChEt" << newbin;
     TH1D * sh3 = singlemuRChEt->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "singleMuChE" << newbin;
     TH1D * sh4 = singlemuRChE ->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "singleChNPart" << newbin;
     TH1D * sh5 = singlemuRChNPart->ProjectionY(binName.str().c_str(),newbin,newbin);
     // sh1 ->SetName(singlemuRChPt->GetName());
     sh1 ->SetTitle(singlemuRChPt->GetName());
     // sh2 ->SetName(singlemuRChP->GetName());
     sh2 ->SetTitle(singlemuRChP->GetName());
     // sh3 ->SetName(singlemuRChEt->GetName());
     sh3 ->SetTitle(singlemuRChEt->GetName());
     // sh4 ->SetName(singlemuRChE->GetName());
     sh4 ->SetTitle(singlemuRChE->GetName());
     // sh5 ->SetName(singlemuRChNPart->GetName());
     sh5 ->SetTitle(singlemuRChNPart->GetName());

     binName.str("");
     binName<< "doubleMuNePt" << newbin;
     TH1D * dn1 = doublemuRNePt->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "doubleMuNeP" << newbin;
     TH1D * dn2 = doublemuRNeP ->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "doubleMuNeEt" << newbin;
     TH1D * dn3 = doublemuRNeEt->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "doubleMuNeE" << newbin;
     TH1D * dn4 = doublemuRNeE ->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "doubleMuNeNPart" << newbin;
     TH1D * dn5 = doublemuRNeNPart->ProjectionY(binName.str().c_str(),newbin,newbin);
     // dn1 ->SetName(doublemuRNePt->GetName());
     dn1 ->SetTitle(doublemuRNePt->GetName());
     // dn2 ->SetName(doublemuRNeP->GetName());
     dn2 ->SetTitle(doublemuRNeP->GetName());
     // dn3 ->SetName(doublemuRNeEt->GetName());
     dn3 ->SetTitle(doublemuRNeEt->GetName());
     // dn4 ->SetName(doublemuRNeE->GetName());
     dn4 ->SetTitle(doublemuRNeE->GetName());
     // dn5 ->SetName(doublemuRNeNPart->GetName());
     dn5 ->SetTitle(doublemuRNeNPart->GetName());

     binName.str("");
     binName<< "singleMuNePt" << newbin;
     TH1D * sn1 = singlemuRNePt->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "singleMuNeP" << newbin;
     TH1D * sn2 = singlemuRNeP ->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "singleMuNeEt" << newbin;
     TH1D * sn3 = singlemuRNeEt->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "singleMuNeE" << newbin;
     TH1D * sn4 = singlemuRNeE ->ProjectionY(binName.str().c_str(),newbin,newbin);
     binName.str("");
     binName<< "singleMuNeNPart" << newbin;
     TH1D * sn5 = singlemuRNeNPart->ProjectionY(binName.str().c_str(),newbin,newbin);
     // sn1 ->SetName(singlemuRNePt->GetName());
     sn1 ->SetTitle(singlemuRNePt->GetName());
     // sn2 ->SetName(singlemuRNeP->GetName());
     sn2 ->SetTitle(singlemuRNeP->GetName());
     // sn3 ->SetName(singlemuRNeEt->GetName());
     sn3 ->SetTitle(singlemuRNeEt->GetName());
     // sn4 ->SetName(singlemuRNeE->GetName());
     sn4 ->SetTitle(singlemuRNeE->GetName());
     // sn5 ->SetName(singlemuRNeNPart->GetName());
     sn5 ->SetTitle(singlemuRNeNPart->GetName());
   
     c1->cd(1);
     // if (ibin == 0) sh1->Draw();
     // if (ibin!= 0) sh1->Draw("SAMES");
     if (ibin == 0){
       sh1->Draw();
       // sh1->GetXaxis()->SetRangeUser(0,100);
     }
     if (ibin!= 0) sh1->Draw("SAMES");
     sh1->SetLineColor(ibin+1);
     c1->GetPad(1)->SetLogy();
     // c1->Update();

     c1->cd(2);
     if (ibin == 0) dh1->Draw();
     if (ibin!= 0) dh1->Draw("SAMES");
     dh1->SetLineColor(ibin+1);
     c1->GetPad(2)->SetLogy();
     // c1->Update();

     c1->cd(3);
     if (ibin == 0) sh3->Draw();
     if (ibin!= 0) sh3->Draw("SAMES");
     sh3->SetLineColor(ibin+1);
     c1->GetPad(3)->SetLogy();
     // c1->Update();

     c1->cd(4);
     if (ibin == 0) dh3->Draw();
     if (ibin!= 0) dh3->Draw("SAMES");
     dh3->SetLineColor(ibin+1);
     c1->GetPad(4)->SetLogy();
     c1->Update();

     c2->cd(1);
     if (ibin == 0) sn1->Draw();
     if (ibin!= 0) sn1->Draw("SAMES");
     sn1->SetLineColor(ibin+1);
     c2->GetPad(1)->SetLogy();
     c2->Update();

     c2->cd(2);
     if (ibin == 0) dn1->Draw();
     if (ibin!= 0) dn1->Draw("SAMES");
     dn1->SetLineColor(ibin+1);
     c2->GetPad(2)->SetLogy();
     c2->Update();

     c2->cd(3);
     if (ibin == 0) sn3->Draw();
     if (ibin!= 0) sn3->Draw("SAMES");
     sn3->SetLineColor(ibin+1);
     c2->GetPad(3)->SetLogy();
     c2->Update();

     c2->cd(4);
     if (ibin == 0) dn3->Draw();
     if (ibin!= 0) dn3->Draw("SAMES");
     dn3->SetLineColor(ibin+1);
     c2->GetPad(4)->SetLogy();
     c2->Update();
     //P e E
     c3->cd(1);
     if (ibin == 0) sh2->Draw();
     if (ibin!= 0) sh2->Draw("SAMES");
     sh2->SetLineColor(ibin+1);
     c3->GetPad(1)->SetLogy();
     c3->Update();

     c3->cd(2);
     if (ibin == 0) dh2->Draw();
     if (ibin!= 0) dh2->Draw("SAMES");
     dh2->SetLineColor(ibin+1);
     c3->GetPad(2)->SetLogy();
     c3->Update();

     c3->cd(3);
     if (ibin == 0) sh4->Draw();
     if (ibin!= 0) sh4->Draw("SAMES");
     sh4->SetLineColor(ibin+1);
     c3->GetPad(3)->SetLogy();
     c3->Update();

     c3->cd(4);
     if (ibin == 0) dh4->Draw();
     if (ibin!= 0) dh4->Draw("SAMES");
     dh4->SetLineColor(ibin+1);
     c3->GetPad(4)->SetLogy();
     c3->Update();

     
     c4->cd(1);
     if (ibin == 0) sn2->Draw();
     if (ibin!= 0) sn2->Draw("SAMES");
     sn2->SetLineColor(ibin+1);
     c4->GetPad(1)->SetLogy();
     c4->Update();

     c4->cd(2);
     if (ibin == 0) dn2->Draw();
     if (ibin!= 0) dn2->Draw("SAMES");
     dn2->SetLineColor(ibin+1);
     c4->GetPad(2)->SetLogy();
     c4->Update();

     c4->cd(3);
     if (ibin == 0) sn4->Draw();
     if (ibin!= 0) sn4->Draw("SAMES");
     sn4->SetLineColor(ibin+1);
     c4->GetPad(3)->SetLogy();
     c4->Update();

     c4->cd(4);
     if (ibin == 0) dn4->Draw();
     if (ibin!= 0) dn4->Draw("SAMES");
     dn4->SetLineColor(ibin+1);
     c4->GetPad(4)->SetLogy();
     c4->Update();



     //NPart
     c5->cd(1);
     if (ibin == 0) sh5->Draw();
     if (ibin!= 0) sh5->Draw("SAMES");
     sh5->SetLineColor(ibin+1);
     sh5->GetXaxis()->SetRangeUser(0,50);
     c5->GetPad(1)->SetLogy();
     c5->Update();

     c5->cd(2);
     if (ibin == 0) dh5->Draw();
     if (ibin!= 0) dh5->Draw("SAMES");
     dh5->SetLineColor(ibin+1);
     dh5->GetXaxis()->SetRangeUser(0,50);
     c5->GetPad(2)->SetLogy();
     c5->Update();

     c5->cd(3);
     if (ibin == 0) sn5->Draw();
     if (ibin!= 0) sn5->Draw("SAMES");
     sn5->SetLineColor(ibin+1);
     sn5->GetXaxis()->SetRangeUser(0,50);
     c5->GetPad(3)->SetLogy();
     c5->Update();

     c5->cd(4);
     if (ibin == 0){
       dn5->Draw();
       dn5->GetXaxis()->SetRangeUser(0,50);
     }
     if (ibin!= 0) dn5->Draw("SAMES");
     dn5->SetLineColor(ibin+1);
     c5->GetPad(4)->SetLogy();
     c5->Update();
     

     // plotName.str(" ");
     // plotName<<doublemuRChPt->GetName()<<newbin<<".eps";
     // c1->SaveAs(plotName.str().c_str());


     }
  c1->SaveAs("muIsoPtEt.eps");
  // c1->SaveAs("muIsoPtEt.pdf");
  c2->SaveAs("muIsoPE.eps");
  // c2->SaveAs("muIsoPE.pdf");
  c3->SaveAs("muIsoNePtEt.eps");
  // c3->SaveAs("muIsoNePtEt.pdf");
  c4->SaveAs("muIsoNePE.eps");
  // c4->SaveAs("muIsoNePE.pdf");
  c5->SaveAs("muIsoNPart.eps");
  // c5->SaveAs("muIsoNPart.pdf");
  	
 
}
