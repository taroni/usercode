// TO RUN 
//root [0] .L moreActiveBits.C+ 
//root [1] moreActiveBits("outfile.root")
#include "Riostream.h" 
#include <memory>
#include <math.h>
#include <sstream>
#include <string>
#include <iostream>

#include <TBranch.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TKey.h>
#include "THashList.h"

using namespace std;
ofstream filedat;
ofstream fivefiledat;

void moreActiveBits(string filename) {
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  vector <TH1F*> vHisto;
   

  TFile * f = new TFile();
  f = TFile::Open(filename.c_str());
  if (!f){
    cout << "FILE NOT FOUND" << endl;
    return;
  }
  TIter hIt(f->GetListOfKeys());
  TKey *keyHisto;
  vector <string> hNames;
  while ((keyHisto=(TKey*)hIt())) {
    hNames.push_back(keyHisto->GetName());
  }

  string name;
  name = filename.substr(0,filename.length()-5);

  stringstream triggerfileName;
  triggerfileName.str("");
  triggerfileName << "trigger_" << name << ".dat";

  stringstream fivetriggerfileName;
  fivetriggerfileName.str("");
  fivetriggerfileName << "firsttriggers_" << name << ".dat";

  filedat.open(triggerfileName.str().c_str(), ofstream::out);
  fivefiledat.open(fivetriggerfileName.str().c_str(), ofstream::out);


  for(unsigned int iHisto=0;iHisto<hNames.size();iHisto++){
    TH1F *h = (TH1F*) f->FindObjectAny(hNames[iHisto].c_str());
    vHisto.push_back(h);

  }
  fivefiledat << "TRIGGER EXCLUDED FROM THE LIST: Alca_*, *_Random*, *_Physics"  << endl; //to check every time
  filedat << "TRIGGER EXCLUDED FROM THE LIST: Alca_*, *_Random*, *_Physics"  << endl; //to check every time
  fivefiledat << "from oldest (first) to most recent (last) trigger key" << endl;
  fivefiledat << endl;
  filedat << "from oldest (first) to most recent (last) trigger key" << endl;
  filedat <<  endl;
  for(unsigned int iHisto=0;iHisto<vHisto.size();iHisto++){
 
    fivefiledat << "#############################" << endl;
    fivefiledat << vHisto[iHisto]->GetName() << endl;
    fivefiledat << "-----------------------------" << endl;
    filedat << "#############################" << endl;
    filedat << vHisto[iHisto]->GetName() << endl;
    filedat << "-----------------------------" << endl;
    map < string, double > myMap;
    myMap.clear();
    cout <<iHisto << " " << vHisto[iHisto]->GetName()<< endl;
    for (unsigned int jBin=1; jBin<vHisto[iHisto]->GetNbinsX();jBin++){
      myMap.insert(pair < string,double > (vHisto[iHisto]->GetXaxis()->GetBinLabel(jBin),vHisto[iHisto]->GetBinContent(jBin) ));

    }//jBin    
    vector <pair  <  string,double  > >myTrigger;
    myTrigger.clear();
    for (map< string,double >::const_iterator jt= myMap.begin() ; jt!=myMap.end(); jt ++){
      double minValue=9999999.;
      double oldMaxValue=99999999;
      pair < string,double > currentTrigger;
      for (map< string,double >::const_iterator it= myMap.begin() ; it!=myMap.end(); it ++){
	bool triggerCheck = false ;
	if (myTrigger.size()!=0){
	  for (unsigned int itr = 0 ; itr< myTrigger.size(); itr++){
	    if (it->first == myTrigger[itr].first){
	      triggerCheck =true;
	    }
	  }//myIt
	}
	if (triggerCheck==true ) continue;
	if (it->second < minValue ){
	  minValue = it->second;
	  currentTrigger = make_pair (it->first , it->second);
	}
      }//it
      myTrigger.push_back(currentTrigger);
    }//jt
    fivefiledat<< "\t" <<"evt/LS" << endl;
    filedat<< "\t" <<"evt/LS" << endl;
    int mytrCount=0; 
    for (unsigned int itr = myTrigger.size()-1 ; itr != 0 ; itr--){
      if (myTrigger[itr].first.substr(0,5) == "AlCa_" ) continue;
      if (myTrigger[itr].first.substr(3,8) == "_Physics") continue;
      if (myTrigger[itr].first.substr(3,7) == "_Random")continue;
      filedat<< "\t" << myTrigger[itr].second << " \t " <<  myTrigger[itr].first << endl;
//       if (itr < myTrigger.size()-11)continue;
      if (mytrCount> 10 )continue;
      fivefiledat<< "\t" << myTrigger[itr].second << " \t " <<  myTrigger[itr].first << endl;
      mytrCount++;
    }
    filedat << endl;
    fivefiledat << endl;
  }//iHisto

  filedat.close();
  fivefiledat.close();


}
