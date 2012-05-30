//to execute root -l hltComparisonMinBiasRun2012A.C+
#include "Riostream.h" 
#include <memory>
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>


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

#define DEBUG 1 

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

void hltComparisonMinBiasRun2012A(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  vector <string> vFileName;
  vector <TH1F *> vHisto;

  //Change the following array with the correct lumisections
  int ievt[] {670, 1096, 664, 984, 583, 298, 985, 1594, 1439, 5051, 69, 0, 0, 0, 1251, 1221, 1309, 424, 2557, 1075, 347, 443, 743, 302};
  vector <double> evt;

  ifstream inputfile("lista_key.txt",ifstream::in);
  if (inputfile.good()){
  }else {
    cout << "BAD file name. Does the file exist?"<< endl;
    return;
  }
 
  vector <TFile * > vFiles;
  vector <TH1F *> vhisto;
  vector <string> name; 
  vector < map < string,double > > vMap;
  int maxChars=256;
  char line[maxChars];
  int ijson =-1; 
  evt.clear();
  while (inputfile.good()){
    ijson++;
    inputfile.getline (line,256);	      
    cout << __LINE__ <<" " <<  line << endl;
    stringstream file; 
    file.str("");
    if (strlen(line) == 0 ) continue; 
    file <<  line << "/DQM_V0001_R000100000__Validation__MinBias__ALCARECO.root";
    TFile *f = TFile::Open(file.str().c_str());
//     stringstream filejson ; 
//     filejson << ijson<<".json"; 
    

    if (!f) continue;
    TH1F*h0=openH1(f);
    cout << __LINE__ << " file " << file.str() << " opened" << endl;
    if (h0->GetEntries() !=0){
      if (h0->GetBinContent(601)!=0) cout << "histo " << line << " has OVERFLOW " << endl;
      vFiles.push_back(f);
      vhisto.push_back(h0);
      name.push_back(line);
      map < string, double > map0;
      vMap.push_back(map0);
      evt.push_back(ievt[ijson]);  // an automatic determination of  lumisection need to be implemented
 //      cout << "histo with "  << h0->GetEntries() << " entries for "  << ievt[ijson] << " lumisection " << endl;
      if (ievt[ijson]==0) cout << "====> histo with zero lumisection" << endl;
    }
  }

  for (int i=0; i<evt.size(); i++)
    cout  << evt[i]<< " ";
  cout << endl;
  
  vector < ofstream > vOut;
  
  rateFile.open("averageRate.dat", ofstream::out);
  excltrigger.open("exclTrigger.dat", ofstream::out);
  singlerate.open("singlerate.dat", ofstream::out);
  outlog.open("hltComparison.log" , ofstream::out);
   
  for (unsigned int iHisto=0; iHisto<vhisto.size(); iHisto++){
    int nlabels= 0 ;
//     cout<< "---------- " << vhisto[iHisto]->GetName()<<" ---------" << vhisto.size() <<endl;
    excltrigger << "--------------" << vhisto[iHisto]->GetName()<<" ---------" << endl;
    singlerate << "histo# " << iHisto <<  " --------------" << vhisto[iHisto]->GetName()<<" ---------" << endl;
    double avergRate=0;
    double npath =0;
    for (int iL = 1; iL <=vhisto[iHisto]->GetXaxis()->GetNbins(); iL++){
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
       string sbstr=str.substr (0,found);
//        if (vhisto[iHisto]->GetBinContent(iL)>0)cout << __LINE__ <<  " "<< found << "  "<< str <<" " <<  sbstr << " " << vhisto[iHisto]->GetBinContent(iL)<< endl;
       str=sbstr;
       vhisto[iHisto]->GetXaxis()->SetBinLabel(iL, sbstr.c_str());
       if (vhisto[iHisto]->GetBinContent(iL)<=0) continue;      
       singlerate << "histo# " << iHisto << " " << iL << " "  << (string)vhisto[iHisto]->GetXaxis()->GetBinLabel(iL) << " " << (double)(vhisto[iHisto]->GetBinContent(iL))/evt[iHisto] 
		  << " " << vhisto[iHisto]->GetBinContent(iL) << " " <<evt[iHisto] <<   endl;
       vMap[iHisto].insert( pair < string, double > ((string)vhisto[iHisto]->GetXaxis()->GetBinLabel(iL),(double)(vhisto[iHisto]->GetBinContent(iL))/(evt[iHisto])));
       avergRate += (double)(vhisto[iHisto]->GetBinContent(iL))/evt[iHisto];
       npath++;
     }
    }
    if (npath !=0) rateFile << vhisto[iHisto]->GetName() <<" average Rate " << avergRate << endl;
  }    
//   cout << __LINE__ << endl;
  singlerate.close();
  
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
  
  vector <TH1F*> vFhisto;
//   cout << __LINE__ << endl;
 
  for (unsigned int iHisto=0; iHisto<vhisto.size(); iHisto++){
//     cout << __LINE__ << endl;
    int  maxvalue=vLabel.size();
    
    TH1F * hF0 = new TH1F (name[iHisto].c_str(),name[iHisto].c_str(), maxvalue, 0, maxvalue);
    vFhisto.push_back(hF0);

  }
  int bincount=0;
  for (unsigned int ilabel =0; ilabel<vLabel.size(); ilabel++){
    bincount++;
//     cout << __LINE__ << " " << bincount << " " << vLabel[ilabel] << endl;
    for (unsigned int imap=0; imap<vMap.size(); imap++){
      vFhisto[imap]->GetXaxis()->SetBinLabel(ilabel+1,vLabel[ilabel].c_str());
    }
    for (unsigned int imap=0; imap<vMap.size(); imap++){
      for ( map < string, double>::const_iterator it = vMap[imap].begin(); it != vMap[imap].end(); it++){
	outlog<< __LINE__ << " " << imap << " " << vLabel[ilabel] << "  " << it->first << " " << it->second<< endl;
	
	if (it->first!=vLabel[ilabel])continue;
	outlog<< __LINE__ << " " << imap << " " << vLabel[ilabel] << "  " << it->first << " " << it->second<< endl;
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

//   cout << __LINE__<< endl;

  outfile -> Close();
  for (unsigned int ifile = 0; ifile < vFiles.size(); ++ifile){
    vFiles[ifile]->Close();
  }

}
