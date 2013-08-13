#define dimensions_cxx
#include "dimensions.h"

void dimensions::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L dimensions.C
//      Root > dimensions t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   TFile * outfile = new TFile("outputSize.root", "RECREATE"); 
   

   Long64_t nentries = fChain->GetEntriesFast();
   map<int, TH1F*> hsize;
   map<int, TH1F*> hsizex;
   map<int, TH1F*> hsizey;
   stringstream name;
   
   for (int ilayer =1 ; ilayer < 4; ilayer++){
     for (int imodule =1; imodule < 9 ; imodule++){
       int histoID=10*ilayer + imodule; 
       name.str(""); 
       name << "clusterSize_layer"<<ilayer<< "_module"<<imodule <<"_bias150"; 
       TH1F * h = new TH1F (name.str().c_str(), name.str().c_str(), 100, 0, 100);
       hsize[histoID]=h;
       name.str(""); 
       name << "clusterSizeX_layer"<<ilayer<< "_module"<<imodule <<"_bias150"; 
       TH1F * hx = new TH1F (name.str().c_str(), name.str().c_str(), 100, 0, 100);
       hsizex[histoID]=hx;
       name.str(""); 
       name << "clusterSizeY_layer"<<ilayer<< "_module"<<imodule <<"_bias150"; 
       TH1F * hy = new TH1F (name.str().c_str(), name.str().c_str(), 100, 0, 100);
       hsizey[histoID]=hy;

     }
   }
   
   Long64_t nbytes = 0, nb = 0;
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
      //      std::cout << nb << endl;
      for (int ilayer =1 ; ilayer < 4; ilayer++){
	for (int imodule =1; imodule < 9 ; imodule++){
	  int histoID= 10* ilayer+ imodule; 
	  hsize[histoID] ->Fill (npix);
	  hsizex[histoID] ->Fill(clust_size_x);
	  hsizey[histoID] ->Fill(clust_size_y);
	  
	  
	}//imodule
      }//ilayer
   }
   
   TCanvas * canvas = new TCanvas(); 
   canvas->Draw(); 
   canvas->cd();
   canvas->SetLogy(1);
   outfile->cd(); 
   for (int ilayer =1 ; ilayer < 4; ilayer++){
     for (int imodule =1; imodule < 9 ; imodule++){
	  int histoID= 10* ilayer+ imodule; 

       hsize[histoID] ->Draw();
       name.str("");
       name << "clusterSize_layer"<<ilayer<< "_module"<<imodule<< ".png";
       canvas->SaveAs(name.str().c_str());
       
       hsizex[histoID] ->Draw();
       name.str("");
       name << "clusterSizeX_layer"<<ilayer<< "_module"<<imodule<< ".png";
       canvas->SaveAs(name.str().c_str());
       hsizey[histoID] ->Draw();
       name.str("");
       name << "clusterSizeY_layer"<<ilayer<< "_module"<<imodule<< ".png";
       canvas->SaveAs(name.str().c_str());
      
       hsize[histoID] ->Write();
       hsizex[histoID] ->Write();
       hsizey[histoID] ->Write();
     }//ilayer
   }//imodule   
   // if (Cut(ientry) < 0) continue;
   outfile->Close();
}
