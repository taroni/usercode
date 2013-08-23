#include <memory>

#include "Riostream.h" 
#include <math.h>
#include <sstream>
#include <vector>
#include <string>

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TProfile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TMath.h>
#define DEBUG 1

using namespace std;
void LAhalfLayers(string filename="lorentzangleALCARECO.root",  bool saveFig= 0){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gROOT->ForceStyle();
  gStyle->SetPadGridY(1); 
  gStyle->SetPadGridX(1); 
  gStyle->SetHistLineWidth(2); 

  double pigreco = 3.141592;
  
  TFile * file = new TFile(filename.c_str(), "READ"); 
  stringstream outname, name;
  outname.str(""); 
  outname<< "compareForwBack_" << filename; 
  TFile *  outfile = new TFile(outname.str().c_str(), "RECREATE"); 
  
  int mybias[]  = {300, 100, 70};
  std::vector<int> vbias (mybias, (mybias + sizeof(mybias)/sizeof(int)));
  vector <TH1D*> vHistoForw, vHistoBack;
 
  TIter hIt(file->GetListOfKeys());
  TKey *keyHisto;
  vector <string> hNames;
  while ((keyHisto=(TKey*)hIt())) {
    hNames.push_back(keyHisto->GetName());
  }
  TCanvas * c= new TCanvas("c", "c"); 
  c->Draw(); 
  TCanvas * c1= new TCanvas("c1", "c1"); 
  c1->Draw(); 
  TLegend *leg = new TLegend(0.4,0.1,0.6,0.2); 
  leg->SetFillColor(0); 
  TH2F * hForw;
  TH2F * hBack;
  TH2F * hForw_noAdc;
  TH2F * hBack_noAdc;
  vector<TGraphErrors*> myGraph; 
  for (int ilayer = 1; ilayer<4; ilayer ++){
    vHistoForw.clear();
    vHistoBack.clear();
    myGraph.clear();
    for (int ibias = 0; ibias< vbias.size(); ibias++){
      for (int imodule = 1; imodule<8; imodule++){
	file->cd(); 
	name.str(""); 
	name << "h_drift_depth_adc_layer"<<ilayer<<"_module"<<imodule<< "_bias_"<<vbias[ibias];
	TH2F *h = (TH2F*) file->Get(name.str().c_str());
	name.str(""); 
	name << "h_drift_depth_noadc_layer"<<ilayer<<"_module"<<imodule<< "_bias_"<<vbias[ibias];
	TH2F *h2 = (TH2F*) file->Get(name.str().c_str());
	
//       TH1D * hy = (TH1D*) h->ProjectionY(); 
//       TH1D * h2y = (TH1D*) h2->ProjectionY(); 
	
	if (imodule < 5){
	  if (imodule ==1 ) {
	    hBack = (TH2F *)h->Clone();
	    hBack_noAdc = (TH2F *)h2->Clone();
	  }else {
	    hBack -> Add(h) ; 
	    hBack_noAdc -> Add(h2) ; 
	  } 
	} else {
	  if (imodule == 5){
	    hForw = (TH2F*) h->Clone(); 
	    hForw_noAdc = (TH2F*) h2->Clone(); 
	  }else{
	    hForw ->Add(h) ; 
	    hForw_noAdc ->Add(h2) ; 
	  }
	  
	}
	
	//    
      }//imodule
      
       hForw->Sumw2(); 
       hBack->Sumw2(); 
//       hForw->Divide(hForw_noAdc);
//       hBack->Divide(hBack_noAdc); 
      

       TCanvas * c2 = new TCanvas("c2", "c2");
      c2->Draw(); 
    
      TH1D * hyForw =hForw->ProjectionY();
      hyForw->Sumw2(); 
      hyForw->Rebin(2); 
      hyForw->Divide(hForw_noAdc->ProjectionY()->Rebin(2)); 
      hyForw->Draw("E"); 
      TH1D * hyBack =hBack->ProjectionY();
      hyBack->Sumw2(); 
      hyBack->Rebin(2); 
      hyBack->Divide(hBack_noAdc->ProjectionY()->Rebin(2)); 

//       hForw_noAdc->Rebin2D(2,2); 
//       hBack_noAdc->Rebin2D(2,2); 

      TProfile  * hpfForw =hForw_noAdc->ProfileY();
      TProfile  * hpfBack =hBack_noAdc->ProfileY();
     

//       TF1 * f = new TF1 ("f", "[3]+[2]*atan((x-[0])/[1])", -1000, 1000);
      TF1 * f = new TF1 ("f", "pol1", -1000, 1000); 
//       f->SetParameter(1, 100); 
//       f->SetParameter(0, 150) ; 
//       f->SetParameter(3, 50) ; 
      
      TGraphErrors * gr = new TGraphErrors (2); 
//       hpfBack->Fit ("f", "R", "", 0, 300); 
//       cout << "Backward: " << f->GetParameter(2)/f->GetParameter(1)  << "+/-" << sqrt(pow(f->GetParError(2)/f->GetParameter(1),2) +  pow(f->GetParameter(2)*f->GetParError(1)/pow(f->GetParameter(1),2),2))<< endl; 

//       gr->SetPoint(0, -1,  f->GetParameter(2)/f->GetParameter(1) );
//       gr->SetPointError(0, 1,sqrt(pow(f->GetParError(2)/f->GetParameter(1),2) +  pow(f->GetParameter(2)*f->GetParError(1)/pow(f->GetParameter(1),2),2)));
//       hpfBack->Fit ("f", "R", "", 120, 180);
      hpfBack->Fit ("f", "R", "", 100, 200);
      gr->SetPoint(0, -1, f->GetParameter(1)); 
      gr->SetPointError(0, 1, f->GetParError(1));

//      hpfForw->Fit ("f", "R", "", -100, 400); 
//       cout << "Forward: " << f->GetParameter(2)/f->GetParameter(1) << "+/-" << sqrt(pow(f->GetParError(2)/f->GetParameter(1),2) +  pow(f->GetParameter(2)*f->GetParError(1)/pow(f->GetParameter(1),2),2)) << endl;
//       gr->SetPoint(1, 1,  f->GetParameter(2)/f->GetParameter(1) );
//       gr->SetPointError(1, 1,sqrt(pow(f->GetParError(2)/f->GetParameter(1),2) +  pow(f->GetParameter(2)*f->GetParError(1)/pow(f->GetParameter(1),2),2)));
			
//       hpfForw->Fit ("f", "R", "", 120, 180); 
      hpfForw->Fit ("f", "R", "", 100, 200); 
      gr->SetPoint(1, 1, f->GetParameter(1)); 
      gr->SetPointError(1, 1, f->GetParError(1));
      
      
      myGraph.push_back(gr); 

      

      outfile ->cd();
      hpfForw ->Write();
      hpfBack ->Write(); 
     


      hyBack->Draw("ESAME"); 
      hyBack->SetLineColor(2); 
      
      name.str("");
      name<< "canvas/rawProfile_layer"<< ilayer<<"_Forw-Back_bias"<< vbias[ibias]<<".png";
      //c2->SaveAs(name.str().c_str());
      vHistoForw.push_back(hyForw);
      vHistoBack.push_back(hyBack);
    }//ibias
    c1->Clear(); 
    c1->Divide(1,vbias.size()); 
    for(int ibias=0; ibias<vbias.size(); ibias++){
      c1->cd(ibias+1) ; 
      //c1->cd(); 
      TH1D* h1= (TH1D*)vHistoForw[ibias]->Clone(); 
      TH1D* h2= (TH1D*)vHistoBack[ibias]->Clone(); 
      h1->Divide(h2); 
      h1->Draw();

      h1->GetYaxis()->SetRangeUser(0.9, 1.1); 
      TLine * line = new TLine(-100, 1, 400, 1); 
      line->SetLineStyle(2);
      line->SetLineColor(4); 
      line->SetLineWidth(2); 
      line->Draw(); 
      name.str(""); 
      name << "chargeRatioForwOverBackvsdepth_layer"<< ilayer<<"_bias"<< vbias[ibias];
      h1->SetTitle(name.str().c_str()); 
      name.str(""); 
      name << "canvas/chargeRatiovsdepth_layer"<< ilayer<<".png";
      if (ibias == vbias.size()-1){
      //	c1->SaveAs(name.str().c_str()); 
	outfile->cd(); 
	c1->Write(); 
      }
      c->cd();
      if (ibias==0) {
	vHistoForw[ibias]->Draw("");
      }else{
	vHistoForw[ibias]->Draw("SAME");
      }	
      vHistoBack[ibias]->Draw("SAME"); 
      vHistoForw[ibias]->SetLineColor(ibias+1);
      vHistoBack[ibias]->SetLineColor(ibias+1); 
      vHistoBack[ibias]->SetLineStyle(2); 
      if (ilayer ==1) {
	name.str(""); 
	name << "V_{b} = "<< vbias[ibias]<< ", Forward half Layer "<<ilayer; 
	leg->AddEntry(vHistoForw[ibias], name.str().c_str(), "lpf"); 
	name.str(""); 
	name << "V_{b} = "<< vbias[ibias]<< ", Backward half Layer "<< ilayer; 
	leg->AddEntry(vHistoBack[ibias], name.str().c_str(), "lpf"); 
      }
      

    }//ibias
    name.str("");
    name<< "canvas/rawProfile_layer"<< ilayer<<"_Forw-Back.png";
    leg->Draw(); 
    c->SaveAs(name.str().c_str()); 
    outfile->cd(); 
    c->Write();
    c->Clear(); 
    TMultiGraph * mg = new TMultiGraph() ; 
    leg->Clear();
    for (int i = 0; i < myGraph.size(); i++){
      //     myGraph[i]->GetXaxis()->SetRangeUser(-2,2);
      myGraph[i]->SetMarkerStyle(20);
      myGraph[i]->SetMarkerColor(i+1); 
      myGraph[i]->SetLineColor(i+1); 
      if (i ==2){
	myGraph[i]->SetMarkerColor(i+2); 
	myGraph[i]->SetLineColor(i+2); 
      }
     
      name.str("");
      name<< "V_{b} = " << vbias[i];
      leg ->AddEntry (myGraph[i], name.str().c_str(), "lp");
      mg->Add(myGraph[i]);
    }
    mg->Draw("AP");
    leg->Draw(); 
    mg->GetXaxis()->SetRangeUser(-2,2);
    c->Update(); 
    name.str("");
    name << "multigraph_layer" << ilayer; 
    c->SetName(name.str().c_str());
    c->Write();
    name << ".png"; 
    c->SaveAs(name.str().c_str()); 
  }//ilayer
  outfile->Close(); 
  file ->Close();

}
