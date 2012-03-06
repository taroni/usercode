#include "Riostream.h" 
#include <memory>
#include <math.h>
#include <sstream>

#include <TBranch.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include "THashList.h"

using namespace std;

TH1F*openH1(TFile *f0)
{

 TDirectory *dir = (TDirectory *) f0->FindObjectAny("DQMData");
  dir ->cd();
  TDirectory *dir1 = (TDirectory *) dir->FindObjectAny("Run 100000");
  dir1 ->cd();
  TDirectory *dir2 = (TDirectory *) dir1->FindObjectAny("HLT_Summary");
  dir2 ->cd();
  TDirectory *dir3 = (TDirectory *) dir2->FindObjectAny("Run summary");
  dir3->cd();

  TH1F*h0= (TH1F *) dir3->Get("PassingBits_Summary_MinBiasRun");
  return h0;
}


ofstream rateFile, excltrigger, outlog;
ofstream singlerate,singlerate1, singlerate2, singlerate3,singlerate4, singlerate5, singlerate6,singlerate7;

void hltComparison(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  vector <string> vFileName;
  vector <TH1F *> vHisto;

  vector <double> evt;

  evt.push_back(794); //0
  evt.push_back(547) ; //1
  evt.push_back(763); //2
  evt.push_back(40); //3
  evt.push_back(553); //4
  evt.push_back(2075); //5
  evt.push_back(1859);  //6 
  evt.push_back(247);  //7
  evt.push_back(389); //8
  evt.push_back(409); //9
  evt.push_back(909); //10
  evt.push_back(557); //11
  evt.push_back(111); //12
  evt.push_back(149); //13
  evt.push_back(1561); //14
  evt.push_back(749); //15
  evt.push_back(794); //16

  string file0("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v6.1_HLT_V6/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
 
  string file1("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v6.1_HLT_V5/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file2("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v6.1_HLT_V3/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file3("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v6.1_HLT_V1/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file4("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v5.3_HLT_V2/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file5("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v5.3_HLT_V1/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file6("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v5.2_HLT_V7/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file7("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v5.2_HLT_V6/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file8("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v5.2_HLT_V5/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");

  string file9("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v5.2_HLT_V2/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");

  string file10("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v5.1_HLT_V3/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");

  string file11("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v4.3_HLT_V4/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file12("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v4.3_HLT_V3/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");

  string file13("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v4.2_HLT_V8/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file14("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v4.2_HLT_V7/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file15("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v4.2_HLT_V6/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  
  string file16("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/lucaroni/STREAM_HLT_120411/CMSSW_4_1_2/src/STREAM_VALIDATION/DiMuon/test/MinimumBias_Run2011A/_cdaq_physics_Run2011_5e32_v4.2_HLT_V5/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root");
  

  



  TFile * f0 = TFile::Open(file0.c_str());
  TFile * f1 = TFile::Open(file1.c_str());
  TFile * f2 = TFile::Open(file2.c_str());
  TFile * f3 = TFile::Open(file3.c_str());
  TFile * f4 = TFile::Open(file4.c_str());
  TFile * f5 = TFile::Open(file5.c_str());
  TFile * f6 = TFile::Open(file6.c_str());
  TFile * f7 = TFile::Open(file7.c_str());
  TFile * f8 = TFile::Open(file8.c_str());
  TFile * f9 = TFile::Open(file9.c_str());
  TFile * f10 = TFile::Open(file10.c_str());
  TFile * f11 = TFile::Open(file11.c_str());
  TFile * f12 = TFile::Open(file12.c_str());
  TFile * f13 = TFile::Open(file13.c_str());
  TFile * f14 = TFile::Open(file14.c_str());
  TFile * f15 = TFile::Open(file15.c_str());
  TFile * f16 = TFile::Open(file16.c_str());

  if (!f0 || !f1 || !f2 || !f3 || !f4 || !f5){
    cout << "at least one file not found" << endl;
    return;
  }

  cout << __LINE__ << endl;

 TH1F*h0=openH1(f0);
 TH1F*h1=openH1(f1);
 TH1F*h2=openH1(f2);
 TH1F*h3=openH1(f3);
 TH1F*h4=openH1(f4);
 TH1F*h5=openH1(f5);
 TH1F*h6=openH1(f6);
 TH1F*h7=openH1(f7);
 TH1F*h8=openH1(f8);
 TH1F*h9=openH1(f9);
 TH1F*h10=openH1(f10);
 TH1F*h11=openH1(f11);
 TH1F*h12=openH1(f12);
 TH1F*h13=openH1(f13);
 TH1F*h14=openH1(f14);
 TH1F*h15=openH1(f15);
 TH1F*h16=openH1(f16);

   vector <TH1F *> vhisto;

  vhisto.push_back(h0);
  vhisto.push_back(h1);
  vhisto.push_back(h2);
  vhisto.push_back(h3);
  vhisto.push_back(h4);
  vhisto.push_back(h5);
  vhisto.push_back(h6);
  vhisto.push_back(h7);
  vhisto.push_back(h8);
  vhisto.push_back(h9);
  vhisto.push_back(h10);
  vhisto.push_back(h11);
  vhisto.push_back(h12);
 vhisto.push_back(h13);
 vhisto.push_back(h14);
 vhisto.push_back(h15);
 vhisto.push_back(h16);
  cout << __LINE__ << endl;

  map < string, double > map0;
  map < string, double > map1;
  map < string, double > map2;
  map < string, double > map3;
  map < string, double > map4;
  map < string, double > map5;
  map < string, double > map6;
  map < string, double > map7;
  map < string, double > map8;
  map < string, double > map9;
  map < string, double > map10;
  map < string, double > map11;
  map < string, double > map12;
  map < string, double > map13;
  map < string, double > map14;
  map < string, double > map15;
  map < string, double > map16;


  vector < map < string,double > > vMap;

  vMap.push_back(map0);
  vMap.push_back(map1);
  vMap.push_back(map2);
  vMap.push_back(map3);
  vMap.push_back(map4);
  vMap.push_back(map5);
  vMap.push_back(map6);
  vMap.push_back(map7);
  vMap.push_back(map8);
  vMap.push_back(map9);
  vMap.push_back(map10);
  vMap.push_back(map11);
  vMap.push_back(map12);
  vMap.push_back(map13);
  vMap.push_back(map14);
  vMap.push_back(map15);
  vMap.push_back(map16);
 
  vector < ofstream > vOut;
  

  
  rateFile.open("averageRate.dat", ofstream::out);
  excltrigger.open("exclTrigger.dat", ofstream::out);
  singlerate.open("singlerate.dat", ofstream::out);
  outlog.open("hltComparison.log" , ofstream::out);
   
    cout << __LINE__ << endl;
    cout<<vhisto.size()<<endl;
 for (unsigned int iHisto=0; iHisto<vhisto.size(); iHisto++){
    int nlabels= 0 ;
    cout << "---------- " << vhisto[iHisto]->GetName()<<" ---------" << vhisto.size() <<endl;
    excltrigger << "--------------" << vhisto[iHisto]->GetName()<<" ---------" << endl;
    singlerate << "histo# " << iHisto <<  " --------------" << vhisto[iHisto]->GetName()<<" ---------" << endl;
    double avergRate=0;
    double npath =0;
        for (int iL = 0; iL <vhisto[iHisto]->GetXaxis()->GetNbins(); iL++){
      nlabels++; 
      if ((string) vhisto[iHisto]->GetXaxis()->GetBinLabel(iL)=="")continue;
      if ((string) vhisto[iHisto]->GetXaxis()->GetBinLabel(iL)=="DQM_FEDIntegrity"){
    	excltrigger << vhisto[iHisto]->GetXaxis()->GetBinLabel(iL) << endl;
    	continue;
     }
     if ((string) vhisto[iHisto]->GetXaxis()->GetBinLabel(iL)=="DQM_TriggerResults"){
    	excltrigger << vhisto[iHisto]->GetXaxis()->GetBinLabel(iL) << endl;
    	continue;
     }
     if ((string) vhisto[iHisto]->GetXaxis()->GetBinLabel(iL)=="HLTriggerFinalPath"){
    	excltrigger <<vhisto[iHisto]->GetXaxis()->GetBinLabel(iL) << endl;
    	continue;
     }
    if ((string) vhisto[iHisto]->GetXaxis()->GetBinLabel(iL)=="HLT_LogMonitor"){
    excltrigger << vhisto[iHisto]->GetXaxis()->GetBinLabel(iL) << endl;
    continue;
    }
    size_t found=1000;
    string str  = (vhisto[iHisto]->GetXaxis()->GetBinLabel(iL));
    found =str.find("_v");
    if (found<str.length()){ 
    string sbstr=str.substr (0,str.length()-3);
     // if (vhisto[iHisto]->GetBinContent(iL)>0)cout << __LINE__ <<  " "<< found << "  "<< str <<" " <<  sbstr << endl;
    	str=sbstr;
    	vhisto[iHisto]->GetXaxis()->SetBinLabel(iL, sbstr.c_str());
     }
      if (vhisto[iHisto]->GetBinContent(iL)<=0) continue;      
     singlerate << "histo# " << iHisto << " " << (string)vhisto[iHisto]->GetXaxis()->GetBinLabel(iL) << " " << (double)(vhisto[iHisto]->GetBinContent(iL))/evt[iHisto] << endl;
    vMap[iHisto].insert( pair < string, double > ((string)vhisto[iHisto]->GetXaxis()->GetBinLabel(iL),(double)(vhisto[iHisto]->GetBinContent(iL))/(evt[iHisto])));
    avergRate += (double)(vhisto[iHisto]->GetBinContent(iL))/evt[iHisto];
    npath++;
     }
    if (npath !=0) rateFile << vhisto[iHisto]->GetName() <<" average Rate " << avergRate << endl;
    }    
  cout << __LINE__ << endl;
  singlerate.close();
  //  int  maxvalue=0;
  //int  maxMap =0;
  //for (unsigned int imap=0; imap<vMap.size(); imap++){
  // if (vMap[imap].size()>=maxvalue) {
  //    maxvalue=vMap[imap].size();
  //   maxMap =imap;
  //   }
  // }
  
  vector <string> vLabel;
  for (unsigned int imap=0; imap<vMap.size(); imap++){
    for ( map < string, double>::const_iterator jt = vMap[imap].begin(); jt != vMap[imap].end(); jt++){
      if (imap ==0 ) {
	vLabel.push_back(jt->first);
      } else {
	int checkName =0;
	for (unsigned int ilabel = 0; ilabel<vLabel.size(); ilabel++){
	  if (vLabel[ilabel]== jt->first) checkName =1; 
	}
	if (checkName==0) vLabel.push_back(jt->first);
      }
    }
  }
  
  
  cout << __LINE__ << endl;
  int  maxvalue=vLabel.size();
  TH1F * hF0 = new TH1F ("_cdaq_physics_Run2011_5e32_v6.1_HLT_V6","_cdaq_physics_Run2011_5e32_v6.1_HLT_V6", maxvalue, 0, maxvalue);
  TH1F * hF1 = new TH1F ("_cdaq_physics_Run2011_5e32_v6.1_HLT_V5", "_cdaq_physics_Run2011_5e32_v6.1_HLT_V5", maxvalue, 0., maxvalue);
  TH1F * hF2 = new TH1F ("_cdaq_physics_Run2011_5e32_v6.1_HLT_V3","_cdaq_physics_Run2011_5e32_v6.1_HLT_V3", maxvalue, 0, maxvalue);
  TH1F * hF3 = new TH1F ("_cdaq_physics_Run2011_5e32_v6.1_HLT_V1","_cdaq_physics_Run2011_5e32_v6.1_HLT_V1", maxvalue, 0., maxvalue);
  TH1F * hF4 = new TH1F ("_cdaq_physics_Run2011_5e32_v5.3_HLT_V2", "_cdaq_physics_Run2011_5e32_v5.3_HLT_V2", maxvalue, 0., maxvalue);
  TH1F * hF5 = new TH1F ("_cdaq_physics_Run2011_5e32_v5.3_HLT_V1", "_cdaq_physics_Run2011_5e32_v5.3_HLT_V1", maxvalue, 0, maxvalue);
  TH1F * hF6 = new TH1F ("_cdaq_physics_Run2011_5e32_v5.2_HLT_V7", "_cdaq_physics_Run2011_5e32_v5.2_HLT_V7", maxvalue, 0, maxvalue);
  TH1F * hF7 = new TH1F ("_cdaq_physics_Run2011_5e32_v5.2_HLT_V6 ", "_cdaq_physics_Run2011_5e32_v5.2_HLT_V6 ", maxvalue, 0, maxvalue);
  TH1F * hF8 = new TH1F ("_cdaq_physics_Run2011_5e32_v5.2_HLT_V5", "_cdaq_physics_Run2011_5e32_v5.2_HLT_V5", maxvalue, 0, maxvalue);
  TH1F * hF9 = new TH1F ("_cdaq_physics_Run2011_5e32_v5.2_HLT_V2", "_cdaq_physics_Run2011_5e32_v5.2_HLT_V2", maxvalue, 0, maxvalue);
  TH1F * hF10 = new TH1F ("_cdaq_physics_Run2011_5e32_v5.1_HLT_V3", "_cdaq_physics_Run2011_5e32_v5.1_HLT_V3", maxvalue, 0, maxvalue);
  TH1F * hF11 = new TH1F ("_cdaq_physics_Run2011_5e32_v4.3_HLT_V4", "_cdaq_physics_Run2011_5e32_v4.3_HLT_V4", maxvalue, 0, maxvalue);
  TH1F * hF12 = new TH1F ("_cdaq_physics_Run2011_5e32_v4.3_HLT_V3 ", "_cdaq_physics_Run2011_5e32_v4.3_HLT_V3 ", maxvalue, 0, maxvalue);
  TH1F * hF13 = new TH1F ("_cdaq_physics_Run2011_5e32_v4.2_HLT_V8", "_cdaq_physics_Run2011_5e32_v4.2_HLT_V8", maxvalue, 0, maxvalue);
  TH1F * hF14 = new TH1F ("_cdaq_physics_Run2011_5e32_v4.2_HLT_V7", "_cdaq_physics_Run2011_5e32_v4.2_HLT_V7", maxvalue, 0, maxvalue);
  TH1F * hF15 = new TH1F ("_cdaq_physics_Run2011_5e32_v4.2_HLT_V6", "_cdaq_physics_Run2011_5e32_v4.2_HLT_V6", maxvalue, 0, maxvalue);
  TH1F * hF16 = new TH1F ("_cdaq_physics_Run2011_5e32_v4.2_HLT_V5", "_cdaq_physics_Run2011_5e32_v4.2_HLT_V5", maxvalue, 0, maxvalue);

 
  //TH1F * hF6 = new TH1F ("HLT139103", "HLT139103", maxvalue, 0., maxvalue);

  vector <TH1F*> vFhisto;
  cout << __LINE__ << endl;

  vFhisto.push_back(hF0);
  vFhisto.push_back(hF1);
  vFhisto.push_back(hF2);
  vFhisto.push_back(hF3);
  vFhisto.push_back(hF4);
  vFhisto.push_back(hF5);
  vFhisto.push_back(hF6);
  vFhisto.push_back(hF7);
  vFhisto.push_back(hF8);
  vFhisto.push_back(hF9);
  vFhisto.push_back(hF10);
  vFhisto.push_back(hF11);
  vFhisto.push_back(hF12);
  vFhisto.push_back(hF13);
  vFhisto.push_back(hF14);
  vFhisto.push_back(hF15);
  vFhisto.push_back(hF16);
 

  int bincount=0;
  for (unsigned int ilabel =0; ilabel<vLabel.size(); ilabel++){
    bincount++;
     // if(jt->second!=0)  vFhisto[maxMap]->Fill(bincount, jt->second);
    cout << __LINE__ << " " << bincount << " " << vLabel[ilabel] << endl;
    for (unsigned int imap=0; imap<vMap.size(); imap++){
      vFhisto[imap]->GetXaxis()->SetBinLabel(ilabel+1,vLabel[ilabel].c_str());
    }
    for (unsigned int imap=0; imap<vMap.size(); imap++){
      // if (imap==maxMap) continue;
      for ( map < string, double>::const_iterator it = vMap[imap].begin(); it != vMap[imap].end(); it++){
	outlog<< __LINE__ << " " << imap << " " << vLabel[ilabel] << "  " << it->first << " " << it->second<< endl;

	if (it->first!=vLabel[ilabel])continue;
	outlog<< __LINE__ << " " << imap << " " << vLabel[ilabel] << "  " << it->first << " " << it->second<< endl;
	// vFhisto[imap]->Fill(ilabel+1, it->second);
	// vFhisto[imap]->GetXaxis()->SetBinLabel(ilabel+1,it->first.c_str());
	vFhisto[imap]->Fill(it->first.c_str(), it->second);

      }// map iterator
    }//for map
  }//ilabel iterator

 
  TFile* outfile = new TFile("outfile.root","recreate");
  outfile->cd();
  for (unsigned int imap=0; imap<vMap.size(); imap++){
    vFhisto[imap]->Write();
  }


  excltrigger.close();
  outlog.close();

  cout << __LINE__<< endl;

  outfile -> Close();
  f0->Close();
  f1->Close();
  f2->Close();
  f3->Close();
  f4->Close();
  f5->Close();
  f6->Close();
  f7->Close();
  f8->Close();
  f9->Close();
  f10->Close();
  f11->Close();
f12->Close();
 f13->Close();
f14->Close();
f15->Close();
f16->Close();
}
