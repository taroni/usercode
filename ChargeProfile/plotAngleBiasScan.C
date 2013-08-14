//      Root > .L readTree.C
//      Root > readTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
#define plotAngleBiasScan_cxx
#include "plotAngleBiasScan.h"




void plotAngleBiasScan::Loop(int ievent){
//   gROOT->Reset();
//   gROOT->SetStyle("Plain");
//   gStyle->SetOptStat(1);
//   gROOT->ForceStyle();
//   gStyle->SetPadGridY(1); 
//   gStyle->SetPadGridX(1); 

  if (DEBUG) cout << __LINE__ << endl; 
  
  map < int , TH1F * > mapAlpha;
  map < int , TH1F * > mapBeta;
  map < int , TH1F * > mapGamma;
  map < int , TH1F * > mapAbsAlpha;
  map < int , TH1F * > mapAbsBeta;
  map < int , TH1F * > mapAbsGamma;
  double pigreco = 3.141592;
  int histoID = 0; 

  if (fChain == 0) return;
  stringstream name; 
  if (DEBUG)  cout << __LINE__ << endl; 
  
  Long64_t nentries = 0; 
  if (ievent ==-1) {
    nentries = fChain->GetEntriesFast();
  }else{
    nentries = ievent;
  }
  
  int mybias[]  = {300, 100, 70};
  std::vector<int> vbias (mybias, (mybias + sizeof(mybias)/sizeof(int)));

  for (int ilayer =1; ilayer < 4; ilayer++) {
    for (int imodule =1; imodule < 9; imodule++){
      for (unsigned int ibias = 0; ibias < vbias.size(); ibias++){
       if (DEBUG)  cout << __LINE__ << endl; 

	histoID = 10000*ilayer + 1000*imodule + vbias[ibias];
	name.str(""); 
	name << "tkalpha_layer" << ilayer <<"_module"<< imodule << "_bias" << vbias[ibias];  
	TH1F * halpha = new TH1F (name.str().c_str(), name.str().c_str(), 320, -3.2, 3.2);
	//  cout << __LINE__ << endl; 
	name.str(""); 
	name << "tkbeta_layer" << ilayer <<"_module"<< imodule << "_bias" << vbias[ibias];  
	TH1F * hbeta = new TH1F (name.str().c_str(), name.str().c_str(), 320, -3.2, 3.2);
	name.str(""); 
	name << "tkgamma_layer" << ilayer <<"_module"<< imodule << "_bias" << vbias[ibias];  
	TH1F * hgamma = new TH1F (name.str().c_str(), name.str().c_str(), 320, -3.2, 3.2);
	mapAlpha[histoID]= halpha; 
	mapBeta[histoID] = hbeta; 
	mapGamma[histoID]= hgamma;


	name.str(""); 
	name << "abstkalpha_layer" << ilayer <<"_module"<< imodule << "_bias" << vbias[ibias];  
	TH1F * halphaAbs = new TH1F (name.str().c_str(), name.str().c_str(), 320, 0, 3.2);
	//  cout << __LINE__ << endl; 
	name.str(""); 
	name << "abstkbeta_layer" << ilayer <<"_module"<< imodule << "_bias" << vbias[ibias];  
	TH1F * hbetaAbs = new TH1F (name.str().c_str(), name.str().c_str(), 320, 0, 3.2);
	name.str(""); 
	name << "abstkgamma_layer" << ilayer <<"_module"<< imodule << "_bias" << vbias[ibias];  
	TH1F * hgammaAbs = new TH1F (name.str().c_str(), name.str().c_str(), 320, 0, 3.2);
	mapAbsAlpha[histoID]= halphaAbs; 
	mapAbsBeta[histoID] = hbetaAbs; 
	mapAbsGamma[histoID]= hgammaAbs;

      
	if (DEBUG)cout << __LINE__ << endl; 
      }//ibias
    }//module
  }//layer

  Long64_t nbytes = 0, nb = 0;
  if (DEBUG)  cout << __LINE__ << endl; 
  int oldbias = -999; 
  bool biasbool=false;
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
    if ( floor(jentry/10000.) == jentry/10000.) cout <<"Processing "<< jentry+1 <<"th entry" <<endl;
    if (bias != oldbias) {
      biasbool = false;
      for (unsigned int ibias=0; ibias < vbias.size() ; ibias++){
	if (bias == vbias[ibias]) biasbool = true;
	oldbias = bias; 
      }
    }
    if (biasbool==false) continue;
    histoID = 10000*layer + 1000*module + bias;
    if (DEBUG) cout << __LINE__<< endl;
     mapAlpha[histoID]->Fill(trackhit_alpha);
    if (DEBUG) cout << __LINE__<< endl;
     mapBeta[histoID] ->Fill(trackhit_beta);
     if (DEBUG) cout << __LINE__<< endl;
    mapGamma[histoID] ->Fill(trackhit_gamma_);
     if (DEBUG) cout << __LINE__<< endl;

     mapAbsAlpha[histoID]->Fill(fabs(trackhit_alpha));
     if (fabs(trackhit_beta)<pigreco/2.) mapAbsBeta[histoID] ->Fill(fabs(trackhit_beta));
     if (fabs(trackhit_beta)>pigreco/2.) mapAbsBeta[histoID] ->Fill(pigreco - fabs(trackhit_beta));
     if (fabs(trackhit_gamma_) > pigreco/2.) mapAbsGamma[histoID] ->Fill(fabs(trackhit_gamma_));
     if (fabs(trackhit_gamma_) < pigreco/2.) mapAbsGamma[histoID] ->Fill(pigreco - fabs(trackhit_gamma_));
       
     if (DEBUG) cout << __LINE__<< endl;
  }//entry
  if (DEBUG) cout << __LINE__<< endl;

  name.str("");
  //  name << "angles_" << f->GetName(); 
  name << "angles.root"; 
  TFile * outfile = new TFile(name.str().c_str(), "RECREATE"); 
  outfile ->cd();
  TCanvas * c = new TCanvas ();
  c->Draw(); 
  gStyle ->SetOptStat(0); 
  gStyle ->SetHistLineWidth(2);
  TLegend * leg = new TLegend(0.1 ,0.8,0.3,0.9);
  if (DEBUG)  cout << __LINE__ << endl; 
  leg->SetFillColor(0); 
  for (int ilayer =1; ilayer < 4; ilayer++) {
    for (int imodule =1; imodule < 9; imodule++){
      for (unsigned int ibias = 0; ibias < vbias.size(); ibias++){
      //  cout << __LINE__ << endl; 
	  if (DEBUG)  cout << __LINE__ << endl; 
	histoID = 10000*ilayer + 1000*imodule + vbias[ibias];
	if (DEBUG)  cout << __LINE__ << endl; 	
	mapAlpha[histoID]->Write(); 
	mapBeta[histoID]->Write(); 
	if (DEBUG)  cout << __LINE__ << endl; 
	mapGamma[histoID]->Write();  
	mapAbsAlpha[histoID]->Write(); 
	mapAbsBeta[histoID]->Write(); 
	mapAbsGamma[histoID]->Write();  
      }
      if (DEBUG)  cout << __LINE__ << endl; 
    }
    for (int imodule = 1; imodule <5; imodule ++){
      for (unsigned int ibias = 0; ibias < vbias.size(); ibias++){
	if (DEBUG)  cout << __LINE__ << endl; 	c->cd() ; 
	c->SetLogy(0); 
	histoID = 10000*ilayer + 1000*imodule + vbias[ibias];

	int histoID2= 10000*ilayer + 1000*(9 - imodule) + vbias[ibias];
	mapAbsAlpha[histoID]->Scale(1./(mapAbsAlpha[histoID]->Integral())); 
	if (DEBUG)  cout << __LINE__ << endl; 	
	mapAbsAlpha[histoID2]->Scale(1./(mapAbsAlpha[histoID2]->Integral())); 
	if (DEBUG)  cout << __LINE__ << endl; 
	mapAbsAlpha[histoID]->Draw(); 
	mapAbsAlpha[histoID2]->Draw("SAME") ; 
	name.str(""); 
	name << "module "<< imodule; 
	leg->Clear(); 
	leg->AddEntry(mapAbsAlpha[histoID],name.str().c_str(),"lpf");
	name.str(""); 
	name << "module "<< 9-imodule; 
	leg->AddEntry(mapAbsAlpha[histoID2],name.str().c_str(),"lpf");
	leg->Draw(); 
	mapAbsAlpha[histoID2]->SetLineColor(2); 
	mapAbsAlpha[histoID]->SetLineWidth(2); 
	mapAbsAlpha[histoID2]->SetLineWidth(2); 
	if (mapAbsAlpha[histoID2]->GetBinContent(mapAbsAlpha[histoID2]->GetMaximumBin())> mapAbsAlpha[histoID]->GetBinContent(mapAbsAlpha[histoID]->GetMaximumBin()))
	  mapAbsAlpha[histoID]->GetYaxis()->SetRangeUser(0, 1.2*mapAbsAlpha[histoID2]->GetBinContent(mapAbsAlpha[histoID2]->GetMaximumBin()));
	if (DEBUG)  cout << __LINE__ << endl; 
	name.str("");
	name << "abstkalpha_layer" << ilayer << "_module" << imodule<< "_module" << 9-imodule<< "_bias" << vbias[ibias];
	mapAbsAlpha[histoID] ->SetTitle(name.str().c_str()); 
	name.str("");
	name << "canvas/abstkalpha_layer" << ilayer << "_module" << imodule<< "_module" << 9-imodule<< "_bias" << vbias[ibias]<< ".png"; 
	c->SaveAs(name.str().c_str());
	//      c->Write();
	double scalefactor= 1./ mapAbsBeta[histoID]->Integral(); 
	double scalefactor2= 1./ mapAbsBeta[histoID2]->Integral(); 
	mapAbsBeta[histoID] ->Scale(1./mapAbsBeta[histoID]->GetEntries()); 
	mapAbsBeta[histoID2] ->Scale(1./mapAbsBeta[histoID2]->GetEntries()); 
	//       cout << __LINE__ << " " << mapAbsBeta[histoID]->GetEntries() <<" " << mapAbsBeta[histoID]->Integral() << endl; 
	//       cout << __LINE__ << " " << mapAbsBeta[histoID2]->GetEntries() <<" " << mapAbsBeta[histoID2]->Integral() << endl; 
	mapAbsBeta[histoID] ->Draw() ; 
	mapAbsBeta[histoID2]->Draw("SAME") ; 
	mapAbsBeta[histoID2]->SetLineColor(2); 
	mapAbsBeta[histoID]->SetLineWidth(2); 
	mapAbsBeta[histoID2]->SetLineWidth(2); 
	if (DEBUG)  cout << __LINE__ << endl; 
	if (mapAbsBeta[histoID2]->GetBinContent(mapAbsBeta[histoID2]->GetMaximumBin())> mapAbsBeta[histoID]->GetBinContent(mapAbsBeta[histoID]->GetMaximumBin()))
	mapAbsBeta[histoID]->GetYaxis()->SetRangeUser(0, 1.2*mapAbsBeta[histoID2]->GetBinContent(mapAbsBeta[histoID2]->GetMaximumBin()));
	name.str("");
	name << "abstkbeta_layer" << ilayer << "_module" << imodule<< "_module" << 9-imodule<< "_bias" << vbias[ibias];
	mapAbsBeta[histoID] ->SetTitle(name.str().c_str()); 
	name.str("");
	name  << "canvas/abstkbeta_layer" << ilayer << "_module" << imodule<< "_module" << 9-imodule<< "_bias" << vbias[ibias]<< ".png"; 
	leg->Draw(); 
	c->SaveAs(name.str().c_str()); 
	//      c->Write();
	//c->SetLogy(1); 
	mapAbsGamma[histoID] ->Scale(1./mapAbsGamma[histoID] ->Integral()); 
	mapAbsGamma[histoID2] ->Scale(1./mapAbsGamma[histoID2] ->Integral()); 
	mapAbsGamma[histoID] ->Draw() ; 
	mapAbsGamma[histoID2]->Draw("SAME") ; 
	mapAbsGamma[histoID2]->SetLineColor(2); 
	mapAbsGamma[histoID]->SetLineWidth(2); 
	mapAbsGamma[histoID2]->SetLineWidth(2); 
	  if (DEBUG)  cout << __LINE__ << endl; 
	if (mapAbsGamma[histoID2]->GetBinContent(mapAbsGamma[histoID2]->GetMaximumBin())> mapAbsGamma[histoID]->GetBinContent(mapAbsGamma[histoID]->GetMaximumBin()))
	  mapAbsGamma[histoID]->GetYaxis()->SetRangeUser(0, 1.2*mapAbsGamma[histoID2]->GetBinContent(mapAbsGamma[histoID2]->GetMaximumBin()));
	if (DEBUG)  cout << __LINE__ << endl; 
	name.str("");
	name << "abstkgamma_layer" << ilayer << "_module" << imodule<< "_module" << 9-imodule<< "_bias" << vbias[ibias];
	mapAbsGamma[histoID]->SetTitle(name.str().c_str()); 
	name.str("");
	name<< "canvas/abstkgamma_layer" << ilayer << "_module" << imodule<< "_module" << 9-imodule << "_bias" << vbias[ibias]<< ".png"; 
	leg->Draw(); 
	c->SaveAs(name.str().c_str()); 
	if (DEBUG)  cout << __LINE__ << endl; 
      //      c->Write();
      }//ibias
    }//imodule      
  }//ilayer
  outfile->Close(); 
  
}
