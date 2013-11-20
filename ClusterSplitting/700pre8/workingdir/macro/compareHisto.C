//g++ -o compareHistoVal compareHistoVal.C `root-config --cflags --libs`
#include <memory>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "Riostream.h" 
#include <math.h>
#include <map>
#include <vector>
#include <TApplication.h>
#include <TBranch.h>
#include <TFile.h>
#include <TLegend.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TPaveStats.h>
#include <TMath.h>
#include <iostream> 
using namespace std;


#define DEBUG 1



void validation(string type="ByChi2"){

//   gROOT->Reset();
//   gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  
  gROOT->ForceStyle();
  
  
  stringstream outfilename ;
  outfilename.str("");
  outfilename << "validation" << type << ".root"; 
  TFile * outfile = new TFile (outfilename.str().c_str(), "RECREATE");

  vector<string> filename;
  filename.clear(); 
  if (type=="ByChi2"){
    filename.push_back("../valid_NS_MC_chi2.root");
    filename.push_back("../valid_TS_MC_chi2.root");
//     filename.push_back("../valid_SS_MC_chi2.root");
  }else if (type=="ByHits"){
    filename.push_back("../valid_NS_MC_hits.root");
    filename.push_back("../valid_TS_MC_hits.root");
//     filename.push_back("../valid_SS_MC_hits.root");
   } else {
    cout << "ERROR type not known" << endl;
  }
  

  map < unsigned int, vector<TH1F* > > histoMap;
  vector <string> vHistoName;
  if (DEBUG) cout << __LINE__<< endl;

  vHistoName.push_back("effic"                );
  vHistoName.push_back("efficPt"              );
//   vHistoName.push_back("effic_PU_eta"         );
//   vHistoName.push_back("effic_PU_vertcount"   );
//   vHistoName.push_back("effic_PU_zpos"	      );
  vHistoName.push_back("effic_vertz_barrel"   );
  vHistoName.push_back("effic_vertz_entire"   );
  vHistoName.push_back("effic_vertz_fwdneg"   );
  vHistoName.push_back("effic_vertz_fwdpos"   );
  vHistoName.push_back("effic_vs_dxy"	      );
  vHistoName.push_back("effic_vs_dz"	      );
  vHistoName.push_back("effic_vs_hit"	      );
  vHistoName.push_back("effic_vs_phi"	      );
  vHistoName.push_back("effic_vs_vertpos"     );
  vHistoName.push_back("effic_vs_zpos"	      );
  vHistoName.push_back("efficdR"       	      );
//   vHistoName.push_back("effic_vertcount_barrel");
//   vHistoName.push_back("effic_vertcount_entire");
//   vHistoName.push_back("effic_vertcount_fwdneg");
//   vHistoName.push_back("effic_vertcount_fwdpos");

  vHistoName.push_back("fakerate"                 );
  vHistoName.push_back("fakeratePt"               );
//   vHistoName.push_back("fakerate_ootpu_barrel"    );
//   vHistoName.push_back("fakerate_ootpu_eta_entire");
//   vHistoName.push_back("fakerate_ootpu_fwdneg"    );
//   vHistoName.push_back("fakerate_ootpu_fwdpos"    );
//   vHistoName.push_back("fakerate_vertcount_barrel");
//   vHistoName.push_back("fakerate_vertcount_entire");
//   vHistoName.push_back("fakerate_vertcount_fwdneg");
//   vHistoName.push_back("fakerate_vertcount_fwdpos");
  vHistoName.push_back("fakerate_vs_dxy"          );
  vHistoName.push_back("fakerate_vs_dz"           );
  vHistoName.push_back("fakerate_vs_hit"          );
  vHistoName.push_back("fakerate_vs_phi"	  );
  vHistoName.push_back("fakeratedR"		  );
  
  if (DEBUG ) cout << __LINE__ << " " << filename.size() << endl;
 
    
  for (unsigned int ifile = 0 ; ifile < filename.size(); ifile++){
    TFile * file; 
    file = TFile::Open(filename[ifile].c_str());
    vector <TH1F*> vHisto;
    if (!file) {
       if (DEBUG ) cout << "FILE "<< filename[ifile] << "NOT FOUND" << endl;
      return;
    } 
    if (DEBUG) cout << __LINE__<< " " << filename[ifile] << endl;
    
    
    for (unsigned int iname=0; iname<vHistoName.size(); iname++){
      TDirectory * dir;
      if (type == "ByChi2" ){
	dir= (TDirectory * ) file->Get("DQMData/Tracking/Track/SPLIT_general_AssociatorByChi2");
	//	if (ifile ==0) dir= (TDirectory * ) file->Get("DQMData/Tracking/Track/general_AssociatorByChi2");
      }else if (type =="ByHits"){
      	dir = (TDirectory *) file->Get("DQMData/Tracking/Track/SPLIT_general_quickAssociatorByHits");
	if (ifile==0) dir = (TDirectory *) file->Get("DQMData/Tracking/Track/general_quickAssociatorByHits");
      }
      TH1F * h = (TH1F *) dir->Get(vHistoName[iname].c_str());
       if (DEBUG ) cout << __LINE__ <<" " <<  vHistoName[iname] << endl;
      TH1F *h1 = (TH1F*)  h->Clone();
      vHisto.push_back(h1);
    }
    histoMap.insert(std::pair < unsigned int, vector <TH1F*> >(ifile, vHisto) ); 
    
  }
  if (DEBUG) cout << __LINE__<< endl;
  outfile->cd();
  TLegend * leg ;
  if (type == "ByChi2") {
    leg = new TLegend (0.7,0.5, 0.95,0.4, "","brNDC");
  }else {
    leg = new TLegend (0.7,0.6, 0.95,0.5, "","brNDC");
  }    
  leg->SetFillColor(0);
  for (unsigned int iname=0; iname<vHistoName.size(); iname++){
    TCanvas *c = new TCanvas(); 
    char bufferE[100];
    char bufferM[100];
    char bufferR[100];
    char bufferN[100]; 
    stringstream name; 
    vector < TPaveStats *> vStat; 

    if (DEBUG) cout << __LINE__ << endl;
    for (unsigned int ifile = 0 ; ifile < filename.size(); ifile++){
      int icolor = ifile+1; 
      if (ifile >1) icolor+=1;
      if (DEBUG) cout << __LINE__<<  endl;
      name.str(""); 
      name << filename[ifile].substr(filename[ifile].find("_",filename[ifile].find("_")+1)+1,filename[ifile].find(".root")-filename[ifile].find("_",filename[ifile].find("_")+1)-1);
      sprintf(bufferN, name.str().c_str());
      if (iname == 0 ){
	if (ifile == 0 ) 	leg->AddEntry(histoMap[ifile][iname],"No splitting");
	if (ifile == 1) 	leg->AddEntry(histoMap[ifile][iname],"Cluster splitter");
	if (ifile == 2) 	leg->AddEntry(histoMap[ifile][iname],"Truth based splitter");
       }
      if (DEBUG) cout << __LINE__<<endl;
      if (ifile ==0){
	if (DEBUG ) cout << histoMap[ifile].size() << endl;
	histoMap[ifile][iname]->Draw(); 
	histoMap[ifile][iname]->GetYaxis()->SetRangeUser(0,1.05);
	c->SetLogx(0);	
	if ((string)histoMap[ifile][iname]->GetName() == "efficPt" || (string) histoMap[ifile][iname]->GetName() == "fakeratePt" ){
	  c->SetLogx(1);
	}
	histoMap[ifile][iname]->SetLineColor(icolor);
	TPaveStats *ptstatsDATA = new TPaveStats(0.8,0.7,.99,0.89,"brNDC");
	ptstatsDATA->SetName("stats");
	ptstatsDATA->SetBorderSize(2);
	ptstatsDATA->SetFillColor(kWhite);
	ptstatsDATA->SetTextAlign(12);
	ptstatsDATA->SetTextColor(icolor);
	sprintf(bufferE,"Entries = %1.0f",histoMap[ifile][iname]->GetEntries());
	sprintf(bufferM,"Mean = %1.2f",histoMap[ifile][iname]->GetMean());
	sprintf(bufferR,"RMS = %1.2f",histoMap[ifile][iname]->GetRMS());
	ptstatsDATA->AddText(bufferN);
	ptstatsDATA->AddText(bufferE);
	ptstatsDATA->AddText(bufferM);
	ptstatsDATA->AddText(bufferR);
	vStat.push_back(ptstatsDATA); 
	

// 	vector <TH1F *> vH = (histoMap[ifile]);
// 	vH[iname] ->GetName(); 
	if (DEBUG) cout << __LINE__<<endl;
      } else {
	histoMap[ifile][iname]->Draw("SAMES");
	histoMap[ifile][iname]->SetLineColor(icolor);
	if (DEBUG) cout << __LINE__<< endl;
	TPaveStats *ptstatsDATA = new TPaveStats(0.8,0.7-0.2*ifile,.99,0.89-0.2*ifile,"brNDC");
	ptstatsDATA->SetName("stats");
	ptstatsDATA->SetBorderSize(2);
	ptstatsDATA->SetFillColor(kWhite);
	ptstatsDATA->SetTextAlign(12);
	ptstatsDATA->SetTextColor(icolor);
	sprintf(bufferE,"Entries = %1.0f",histoMap[ifile][iname]->GetEntries());
	sprintf(bufferM,"Mean = %1.2f",histoMap[ifile][iname]->GetMean());
	sprintf(bufferR,"RMS = %1.2f",histoMap[ifile][iname]->GetRMS());
	ptstatsDATA->AddText(bufferN);
	ptstatsDATA->AddText(bufferE);
	ptstatsDATA->AddText(bufferM);
	ptstatsDATA->AddText(bufferR);
 	vStat.push_back(ptstatsDATA); 
	
	
      }
//       vStat[ifile] ->Draw();
      leg->Draw();
    }//file
    stringstream cname; 
    cname.str(""); 
    cname << vHistoName[iname];
    
    c->SetName(cname.str().c_str());
    c->Write();
    cname << ".png"; 

    c->SaveAs(cname.str().c_str());

  }//iname


//   for (map< unsigned int, vector <TH1F *> >::const_iterator it=vHistoName.begin(); it!=vHistoName.end(); it++){
//     for (unsigned int ih = 0; ih< it.second.size(); ih++){
//       (it.second)[ih]->Draw();
      
//     }
    
//   }
  outfile->Close(); 
  cout << "The End"<< endl;

}
