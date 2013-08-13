//      Root > .L readTree.C
//      Root > readTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
#define plotChargevsBeta_cxx
#include "plotChargevsBeta.h"



void plotChargevsBeta::Loop(){
//   gROOT->Reset();
//   gROOT->SetStyle("Plain");
//   gStyle->SetOptStat(1);
//   gROOT->ForceStyle();
//   gStyle->SetPadGridY(1); 
//   gStyle->SetPadGridX(1); 

//  cout << __LINE__ << endl; 
  
  map < int , TH2F * > mapDepthBeta;
  map < int , TH2F * > mapDepthAbsBeta;
  map < int , TH2F * > mapDepthBetaNoAdc;
  map < int , TH2F * > mapDepthAbsBetaNoAdc;

  double pigreco = 3.141592;
  int histoID = 0; 

  float width_ = 0.0285;

  if (fChain == 0) return;
  stringstream name; 
  //  cout << __LINE__ << endl; 
  Long64_t nentries = fChain->GetEntriesFast();

  for (int ilayer =1; ilayer < 4; ilayer++) {
    for (int imodule =1; imodule < 9; imodule++){
      //  cout << __LINE__ << endl; 

      histoID = 10*ilayer + imodule;
      name.str(""); 
      name << "depthbeta_adc_layer" << ilayer <<"_module"<< imodule;  
      TH2F * hbeta = new TH2F (name.str().c_str(), name.str().c_str(), 320, -3.2, 3.2, 50, -100, 400);
      name.str(""); 
      name << "depthbeta_noadc_layer" << ilayer <<"_module"<< imodule;  
      TH2F * hbetanoadc = new TH2F (name.str().c_str(), name.str().c_str(), 320, -3.2, 3.2, 50, -100, 400);
      mapDepthBeta[histoID] = hbeta; 
      mapDepthBetaNoAdc[histoID] = hbetanoadc; 

      name.str(""); 
      name << "depthabsbeta_adc_layer" << ilayer <<"_module"<< imodule;  
      TH2F * habsbeta = new TH2F (name.str().c_str(), name.str().c_str(), 160, 0, 3.2, 50, -100, 400);
      name.str(""); 
      name << "depthabsbeta_noadc_layer" << ilayer <<"_module"<< imodule;  
      TH2F * habsbetanoadc = new TH2F (name.str().c_str(), name.str().c_str(), 160, 0, 3.2, 50, -100, 400);

      mapDepthAbsBeta[histoID] = habsbeta; 
      mapDepthAbsBetaNoAdc[histoID] = habsbetanoadc;


      
      //  cout << __LINE__ << endl; 

    }//module
  }//layer

  Long64_t nbytes = 0, nb = 0;
  //  cout << __LINE__ << endl; 

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (chi2/ndof > 2) continue;
    if (clust_size_y < 4) continue;
    if (clust_charge> 120) continue; 
    double residual = TMath::Sqrt( (trackhit_x - rechit_x) * (trackhit_x - rechit_x) + (trackhit_y - rechit_y) * (trackhit_y - rechit_y) );
    if (residual > 0.005) continue; 
    bool large_pix = false;
    for (int j = 0; j <  npix; j++){
      int colpos = static_cast<int>(colpix[j]);
      if (rowpix[j] == 0 || rowpix[j] == 79 || rowpix[j] == 80 || rowpix[j] == 159 || colpos % 52 == 0 || colpos % 52 == 51 ){
	      large_pix = true;	
      }
    }//end of pix loop
    if (large_pix==true) continue; 
      // if (Cut(ientry) < 0) continue;
    //  cout << nb << endl; 
    if ( (int) jentry/10000. == jentry/10000.) cout <<"Processing "<< jentry+1 <<"th entry" <<endl;
    histoID = 10*layer + module;
    for (int i = 0 ; i < npix ; i++){

      float dy = (xpix[i]  - (trackhit_y - width_/2. / TMath::Tan(trackhit_beta))) * 10000.;
      float depth = dy * tan(trackhit_beta);
      
      mapDepthBeta[histoID] ->Fill(trackhit_beta, depth, norm_charge[i]);
      mapDepthBetaNoAdc[histoID] ->Fill(trackhit_beta, depth);

      if (fabs(trackhit_beta)<pigreco/2.){
	mapDepthAbsBeta[histoID] ->Fill(fabs(trackhit_beta), depth, norm_charge[i]);
	mapDepthAbsBetaNoAdc[histoID] ->Fill(fabs(trackhit_beta), depth);
      }
      if (fabs(trackhit_beta)>pigreco/2.) {
	mapDepthAbsBeta[histoID] ->Fill(pigreco - fabs(trackhit_beta), depth, norm_charge[i]);
	mapDepthAbsBetaNoAdc[histoID] ->Fill(pigreco - fabs(trackhit_beta), depth);
      }
    }

       
    
  }//entry

  name.str("");
  //  name << "angles_" << f->GetName(); 
  name << "ratiovsbeta.root"; 
  TFile * outfile = new TFile(name.str().c_str(), "RECREATE"); 
  outfile ->cd();
  TCanvas * c = new TCanvas ();
  c->Draw(); 
  gStyle ->SetOptStat(0); 
  gStyle ->SetHistLineWidth(2);
  TLegend * leg = new TLegend(0.1 ,0.8,0.3,0.9);
  leg->SetFillColor(0); 
  for (int ilayer =1; ilayer < 4; ilayer++) {
    for (int imodule =1; imodule < 9; imodule++){
      histoID = 10*ilayer + imodule;
      mapDepthBeta[histoID]->Divide(mapDepthBetaNoAdc[histoID]);
      mapDepthAbsBeta[histoID]->Divide(mapDepthAbsBetaNoAdc[histoID]);
      mapDepthBeta[histoID]->Write();
      mapDepthAbsBeta[histoID]->Write();

      mapDepthBeta[histoID]->Draw("COLZ") ;
      name.str("");
      name << "depthvsbeta_layer"<< ilayer << "_imodule" << imodule <<".png";
      c->SaveAs(name.str().c_str());
      mapDepthAbsBeta[histoID]->Draw("COLZ") ;
      name.str("");
      name << "depthvsabsbeta_layer"<< ilayer << "_imodule" << imodule <<".png";
      c->SaveAs(name.str().c_str());

   }
  }
  outfile->Close(); 
  
}
