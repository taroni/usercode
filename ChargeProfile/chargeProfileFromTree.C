//      Root > .L readTree.C
//      Root > readTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
#define chargeProfileFromTree_cxx
#include "chargeProfileFromTree.h"



void chargeProfileFromTree::Loop(){
//   gROOT->Reset();
//   gROOT->SetStyle("Plain");
//   gStyle->SetOptStat(1);
//   gROOT->ForceStyle();
//   gStyle->SetPadGridY(1); 
//   gStyle->SetPadGridX(1); 

//  cout << __LINE__ << endl; 
  
  map < int , TH2F * > mapDepthBeta;
  //map < int , TH2F * > mapDepthAbsBeta;
  map < int , TH2F * > mapDepthBetaNoAdc;
//   map < int , TH2F * > mapDepthAbsBetaNoAdc;

  double pigreco = 3.141592;
  int histoID = 0; 

  float width_ = 0.0285;
  // int mybias[]  = {150};
  int mybias[]  = {300, 100, 70};
  std::vector<int> vbias (mybias, (mybias + sizeof(mybias)/sizeof(int)));
  int oldbias = -99;
  stringstream name ; 
  Long64_t nentries = fChain->GetEntriesFast();
  bool biasbool = false; 
  if (fChain == 0) return;

  for (int ilayer =1; ilayer < 4; ilayer++) {
    for (int imodule =1; imodule < 9; imodule++){
      for (unsigned int ibias = 0; ibias < vbias.size(); ibias++){
	if (DEBUG)  cout << __LINE__ << endl; 

	histoID = 10000*ilayer + 1000*imodule + vbias[ibias];
	
      //  cout << __LINE__ << endl; 

	name.str(""); 
	name << "h_drift_depth_adc_layer" << ilayer <<"_module"<< imodule<<"_bias_"<< vbias[ibias];  
	TH2F * hbeta = new TH2F (name.str().c_str(), name.str().c_str(), 200, -1000, 1000, 50, -100, 400);
	name.str(""); 
	name << "h_drift_depth_noadc_layer"  << ilayer <<"_module"<< imodule<<"_bias_"<< vbias[ibias];  
	TH2F * hbetanoadc = new TH2F (name.str().c_str(), name.str().c_str(), 200, -1000, 1000, 50, -100, 400);
	mapDepthBeta[histoID] = hbeta; 
	mapDepthBetaNoAdc[histoID] = hbetanoadc; 

//       name.str(""); 
//       name << "h_drift_depth_adc_layer" << ilayer <<"_module"<< imodule;  
//       TH2F * habsbeta = new TH2F (name.str().c_str(), name.str().c_str(), 160, 0, 3.2, 50, -100, 400);
//       name.str(""); 
//       name << "h_drift_depth_adc_layer"  << ilayer <<"_module"<< imodule;  
//       TH2F * habsbetanoadc = new TH2F (name.str().c_str(), name.str().c_str(), 160, 0, 3.2, 50, -100, 400);

//       mapDepthAbsBeta[histoID] = habsbeta; 
//       mapDepthAbsBetaNoAdc[histoID] = habsbetanoadc;


      
      //  cout << __LINE__ << endl; 
      }//vbias
    }//module
  }//layer

  Long64_t nbytes = 0, nb = 0;
  //  cout << __LINE__ << endl; 

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (DEBUG) cout << __LINE__ << endl; 
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( (int)(jentry/10000.) == jentry/10000.) cout <<"Processing "<< jentry+1 <<"th entry" <<endl;
    
    if (DEBUG)cout << __LINE__ << " chi2/ndof " << chi2 << " " << ndof << " " << chi2/ndof<< endl;     
    if (pt<1.) continue;
    if (chi2/ndof > 2) continue;
    if (isHighPurity!=1) continue;
    if (nhits<7) continue;
    //  if (npixhits < 3 ) continue;
    if (clust_size_y < 4) continue;
    if (DEBUG) cout << __LINE__ << " clust_size_y "<< clust_size_y<< endl; 
    if (clust_charge > 120) continue;
    if (DEBUG)cout << __LINE__ << endl;  
    double residual = TMath::Sqrt( (trackhit_x - rechit_x) * (trackhit_x - rechit_x) + (trackhit_y - rechit_y) * (trackhit_y - rechit_y) );
    if (residual > 0.005) continue; 
   if (DEBUG) cout << __LINE__ << endl; 
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
    bool layerbool[3] = {false, false, false}; 
    if (vbias.size()>1) {
      if (run == 208392){
	layerbool[0] = true; 
	layerbool[1] = false; 
	layerbool[2] = false; 
      }
      if (DEBUG) cout << __LINE__ << " " << __PRETTY_FUNCTION__ << endl;
    
      if (run == 208393){
	if ((orbit > 1511670 && orbit < 4102155) || (orbit > 4102155 && orbit < 6119516) || (orbit> 6119516 && orbit < 8504600) || 
	  (orbit > 8504600 && orbit < 11338352) ||(orbit > 11338352 && orbit < 14205283)||  (orbit > 14205283 && orbit < 16640930)||(orbit > 16640930 && orbit < 18756493)||(orbit > 18756493 && orbit < 20815870)|| orbit < 20815870) {
          layerbool[0] = true; 
	  layerbool[1] = false; 
	  layerbool[2] = false;
        }else if (orbit >21140745&& orbit <50577268) {
	  layerbool[0] = false; 
	  layerbool[1] = true; 
	  layerbool[2] = false;
	} else if (orbit > 50825857) {
	layerbool[0] = false; 
	layerbool[1] = false; 
	layerbool[2] = true;
	}
      }
    
      if (run == 208394 ||( run == 208395 && orbit < 19753949 )) {
	layerbool[0] = false; 
	layerbool[1] = false; 
	layerbool[2] = true;
      }
    }else{
      layerbool[0] = true; 
      layerbool[1] = true; 
      layerbool[2] = true;
    }
    
  //  cout << __LINE__ << endl; 

    if (bias != oldbias) {
      biasbool = false;
      for (unsigned int ibias=0; ibias < vbias.size() ; ibias++){
	if (bias == vbias[ibias]) biasbool = true;
	oldbias = bias; 
      }
    }
    if (biasbool==false) continue;
    histoID = 10000*layer + 1000*module + bias;
    for (int i = 0 ; i < npix ; i++){
      float dx = (xpix[i]  - (trackhit_x - width_/2. / TMath::Tan(trackhit_alpha))) * 10000.;
      float dy = (ypix[i]  - (trackhit_y - width_/2. / TMath::Tan(trackhit_beta))) * 10000.;
 //      cout <<__LINE__<< " "<<  dx << " " <<xpix[i]<< " "  << trackhit_x  << " " << width_/2 << " " <<  (trackhit_x - width_/2. / TMath::Tan(trackhit_alpha))<< endl; 
//       cout <<__LINE__<< " "<< dy << " " << ypix[i]<< " " << trackhit_y  << " " << width_/2 <<" " <<  (trackhit_y - width_/2. / TMath::Tan(trackhit_beta)) << endl; 
      float depth = dy * tan(trackhit_beta);
      float drift = dx - dy * tan(trackhit_gamma_);
      
      
      //     mapDepthBeta[histoID] ->Fill(drift, depth, adc[i]);
      if (layerbool[layer -1]==true ){
	double mynorm_charge =  adc[i] *sqrt(1.0/(1.0/pow(tan(trackhit_alpha),2)+1.0/pow(tan(trackhit_beta),2)+1.0));  // pixeltree corrections
	//	mapDepthBeta[histoID] ->Fill(drift, depth,mynorm_charge);
	mapDepthBeta[histoID] ->Fill(drift, depth,adc[i]);
	mapDepthBetaNoAdc[histoID] ->Fill(drift, depth);
	//	cout << bias << " "<<layerbool[0] << " " << layerbool[1]<< " " << layerbool[2]<< " "  << mapDepthBetaNoAdc[histoID] ->GetName() << endl; 
    
      }
    }

       
    
  }//entry

  name.str("");
  //  name << "angles_" << f->GetName(); 
  name << "driftdepth.root"; 
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
      for (unsigned int ibias = 0; ibias < vbias.size(); ibias++){

	histoID = 10000*ilayer + 1000*imodule + vbias[ibias];
//       mapDepthBeta[histoID]->Divide(mapDepthBetaNoAdc[histoID]);
//       mapDepthAbsBeta[histoID]->Divide(mapDepthAbsBetaNoAdc[histoID]);
      mapDepthBeta[histoID]->Write();
      mapDepthBetaNoAdc[histoID]->Write();
//       mapDepthAbsBeta[histoID]->Write();
      
//       mapDepthBeta[histoID]->Draw("COLZ") ;
//       name.str("");
//       name << "depthvsbeta_layer"<< ilayer << "_imodule" << imodule <<".png";
//       c->SaveAs(name.str().c_str());
//       mapDepthAbsBeta[histoID]->Draw("COLZ") ;
//       name.str("");
//       name << "depthvsabsbeta_layer"<< ilayer << "_imodule" << imodule <<".png";
//       c->SaveAs(name.str().c_str());
      }//bias
    }//module
  }//layer
  outfile->Close(); 
  
}
