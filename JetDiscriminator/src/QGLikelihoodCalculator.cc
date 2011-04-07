#include "TauAnalysis/FakeRate/interface/QGLikelihoodCalculator.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <vector>
#include "TMath.h"

using namespace std;

void getBins_int( int nBins_total, std::vector<double >& Lower, Double_t xmin, Double_t xmax, bool plotLog=true);




// constructor:

QGLikelihoodCalculator::QGLikelihoodCalculator( const std::string& fileName, int nPtBins ) {
  histoFile_ = TFile::Open(fileName.c_str());

  nPtBins_ = nPtBins;

}
QGLikelihoodCalculator::~QGLikelihoodCalculator( ) {
  histoFile_ ->Close();
//   cout << __PRETTY_FUNCTION__<< " " << __LINE__ << endl;
  histoFile_->~TFile();
}


// int QGLikelihoodCalculator::OpenFile(const std::string& fileName){
//   histoFile_ = TFile::Open(fileName.c_str());
//   if (!histoFile_ ) {
//     return 0;
//     // cout << "No File"  << endl;
//   }else {
//     return 1;
//   }
  
// }
float QGLikelihoodCalculator::computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD ) {

  return this->computeQGLikelihood( pt, nCharged, nNeutral, ptD, -1. );

}


float QGLikelihoodCalculator::computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD, float rmsCand ) {

// cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  " " ;
// cout << pt <<  " " <<nCharged << " " << nNeutral<< " " <<ptD;

// cout << endl;

  float ptMin = 0.;
  float ptMax = 0.;
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  const int nBinsPlusOne(nPtBins_+1);
  //Double_t ptBins[nBinsPlusOne];
  std::vector<double>  ptBins;
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
//fitTools::getBins_int( nBinsPlusOne, ptBins, 15., 1000. );
  getBins_int( nBinsPlusOne, ptBins, 15., 1000. );

// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " ptBins[nPtBins_ ] " <<  ptBins[nPtBins_ ] << endl;

  if( pt > ptBins[nPtBins_ ] ) {
   // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
   ptMin = ptBins[nPtBins_-1 ] ;
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  ptMax = ptBins[nPtBins_ ];
// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " "<<  ptMin<< " "<< ptMax << endl;
  } else {
    for( int iBin=0; iBin<nPtBins_; ++iBin ) {
// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " "<<ptBins[iBin]<< " "<<ptBins[iBin+1] << endl;
      if( pt>ptBins[iBin] && pt<=ptBins[iBin+1] ) {
        ptMin = ptBins[iBin];
        ptMax = ptBins[iBin+1];
// 	cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " "<<  ptMin<< " "<< ptMax << endl;
      } //if
    } //for
  } //else
  
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  if( ptMax==0. ) return -1.;


  char histoName[200];
  sprintf( histoName, "nCharged_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nCharged_gluon = (TH1F*)histoFile_->Get(histoName);
  sprintf( histoName, "nCharged_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nCharged_quark = (TH1F*)histoFile_->Get(histoName);

  sprintf( histoName, "nNeutral_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nNeutral_gluon = (nNeutral>0) ? (TH1F*)histoFile_->Get(histoName) : 0;
  sprintf( histoName, "nNeutral_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nNeutral_quark = (nNeutral>0) ? (TH1F*)histoFile_->Get(histoName) : 0;

  sprintf( histoName, "ptD_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_ptD_gluon = (ptD>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;
  sprintf( histoName, "ptD_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_ptD_quark = (ptD>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;

  sprintf( histoName, "rmsCand_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_rmsCand_gluon = (rmsCand>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;
  sprintf( histoName, "rmsCand_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_rmsCand_quark = (rmsCand>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;

  float gluonP = likelihoodProduct( nCharged, nNeutral, ptD, rmsCand, h1_nCharged_gluon, h1_nNeutral_gluon, h1_ptD_gluon, h1_rmsCand_gluon );
// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<h1_nNeutral_gluon->Integral() << endl;
  float quarkP = likelihoodProduct( nCharged, nNeutral, ptD, rmsCand, h1_nCharged_quark, h1_nNeutral_quark, h1_ptD_quark, h1_rmsCand_quark );

// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " << gluonP << " " << quarkP << endl;

  //float QGLikelihood = gluonP / (gluonP + quarkP );
  float QGLikelihood = quarkP / (gluonP + quarkP );

  return QGLikelihood;
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;

}


float QGLikelihoodCalculator::likelihoodProduct( float nCharged, float nNeutral, float ptD, float rmsCand, TH1F* h1_nCharged, TH1F* h1_nNeutral, TH1F* h1_ptD, TH1F* h1_rmsCand) {

  h1_nCharged->Scale(1./h1_nCharged->Integral("width"));
// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<h1_nCharged->Integral() << endl;
 if( h1_nNeutral!=0 ) h1_nNeutral->Scale(1./h1_nNeutral->Integral("width"));
// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<h1_nNeutral->Integral("width") << endl;
  if( h1_ptD!=0 )
    h1_ptD->Scale(1./h1_ptD->Integral("width"));
  if( h1_rmsCand!=0 )
    h1_rmsCand->Scale(1./h1_rmsCand->Integral("width"));

  float likeliProd =  h1_nCharged->GetBinContent(h1_nCharged->FindBin(nCharged));

  if( h1_nNeutral!=0 )
    likeliProd*=h1_nNeutral->GetBinContent(h1_nNeutral->FindBin(nNeutral));
// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<likeliProd << endl;
  if( h1_ptD!=0 )
    likeliProd*=h1_ptD->GetBinContent(h1_ptD->FindBin(ptD));
// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<likeliProd << endl;

  if( h1_rmsCand!=0 )
    likeliProd*=h1_rmsCand->GetBinContent(h1_rmsCand->FindBin(rmsCand));


// cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<likeliProd << endl;
  return likeliProd;

}


void getBins_int( int nBins_total, std::vector<double >&  Lower, Double_t xmin, Double_t xmax, bool plotLog) {
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;

  Double_t Lower_exact;
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  int nBins = nBins_total-1;
  // cout <<__PRETTY_FUNCTION__ << " " << __LINE__ <<" " << nBins<< endl;
  const double dx = (plotLog) ? pow((xmax / xmin), (1. / (double)nBins)) : ((xmax - xmin) / (double)nBins);
//   // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " "  << Lower->size() <<  endl;
//   if (Lower->size()==0){
   // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
    Lower.push_back(xmin);
   // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
//   }else{
//     Lower->at(0) = xmin;
//   // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
//   }
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  Lower_exact = Lower[0];
  for (int i = 1; i != nBins; ++i) {
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;

    if (plotLog) {
      Lower_exact *= dx;
//       Lower[i] = TMath::Ceil(Lower_exact);
      Lower.push_back(TMath::Ceil(Lower_exact));
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
    } else {
//       Lower[i] = TMath::Ceil(Lower[i-1] + dx);
      Lower.push_back(TMath::Ceil(Lower.at(i-1) + dx));
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
    }

  }
  // cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;

//   Lower[nBins] = xmax;
  Lower.push_back(xmax);

}

