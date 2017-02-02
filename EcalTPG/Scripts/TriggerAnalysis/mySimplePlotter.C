
{

  TCut formula = "rawTPData>0";
  TFile *_file0 = TFile::Open("ECALTPGtree_246926.root");
  TTree * tree = (TTree*) _file0->Get("EcalTPGAnalysis") ;
  TCanvas * c = new TCanvas(); 
  c->Draw();
  
  tree->Draw("trig_tower_sFGVB", formula);
  c->SaveAs("sfgvb.pdf") ; 
  c->SaveAs("sfgvb.png") ; 
  
  tree->Draw("trig_tower_sFGVB:ieta", formula , "colz") ;
  c->SaveAs("sfgvb_vs_eta.pdf") ; 
  c->SaveAs("sfgvb_vs_eta.png") ; 

  tree->Draw("ieta", formula&&"trig_tower_sFGVB==0") ;
  c->SaveAs("ieta_for_sfgvb.pdf") ; 
  c->SaveAs("ieta_for_sfgvb.png") ; 
 
  tree->Draw("sevlv:ieta", formula&&"sevlv!=0","colz") ;
  c->SaveAs("sevlv_vs_ieta.pdf") ; 
  c->SaveAs("sevlv_vs_ieta.png") ; 
  
  tree->Draw("sevlv:trig_tower_sFGVB", formula&&"sevlv!=0", "colz") ;
  c->SaveAs("sevlv_vs_trig_tower_sFGVB.pdf");
  c->SaveAs("sevlv_vs_trig_tower_sFGVB.png");
  
  tree->Draw("sevlv:trig_tower_sFGVB", formula&&"sevlv!=0 && fabs(ieta) < 18", "colz") ;
  c->SaveAs("sevlv_vs_trig_tower_sFGVB_barrel.pdf");
  c->SaveAs("sevlv_vs_trig_tower_sFGVB_barrel.png");
  
  tree->Draw("sevlv:trig_tower_sFGVB", formula&&"sevlv!=0 && fabs(ieta) >= 18", "colz") ;
  c->SaveAs("sevlv_vs_trig_tower_sFGVB_endcap.pdf");
  c->SaveAs("sevlv_vs_trig_tower_sFGVB_endcap.png");
  


  tree->Draw("ieta:iphi", formula&&"sevlv==1 && fabs(ieta)< 18 ", "colz") ;
  c->SaveAs("eta_vs_phi_sevlv1_barrel.pdf");
  c->SaveAs("eta_vs_phi_sevlv1_barrel.png");
  tree->Draw("ieta:iphi", formula&&"sevlv==2 && fabs(ieta)< 18 ", "colz") ;
  c->SaveAs("eta_vs_phi_sevlv2_barrel.pdf");
  c->SaveAs("eta_vs_phi_sevlv2_barrel.png");
  tree->Draw("ieta:iphi", formula&&"sevlv==3 && fabs(ieta)< 18 ", "colz") ;
  c->SaveAs("eta_vs_phi_sevlv3_barrel.pdf");
  c->SaveAs("eta_vs_phi_sevlv3_barrel.png");
  tree->Draw("ieta:iphi", formula&&"sevlv==4 && fabs(ieta)< 18 ", "colz") ;
  c->SaveAs("eta_vs_phi_sevlv4_barrel.pdf");
  c->SaveAs("eta_vs_phi_sevlv4_barrel.png");
  tree->Draw("ieta:iphi", formula&&"sevlv==5 && fabs(ieta)< 18 ", "colz") ;
  c->SaveAs("eta_vs_phi_sevlv5_barrel.pdf");
  c->SaveAs("eta_vs_phi_sevlv5_barrel.png");



}
