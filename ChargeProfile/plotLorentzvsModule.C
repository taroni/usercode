{
  gStyle->SetOptStat(0);
  // _file0->ls();
  TH1D * h300 = (TH1D * ) _file0->Get("tanLAvsModule_layer1_bias300");
  h300->Draw();
  h300->GetYaxis()->SetRangeUser(0.2,0.8);
  h300->Draw();
  TH1D * h100 = (TH1D * ) _file0->Get("tanLAvsModule_layer1_bias100");
  h300->SetMarkerStyle(20);
  h100->SetMarkerStyle(20);
  h100->SetMarkerColor(21);
  h100->SetMarkerColor(2);
  h100->SetLineColor(2);
  h100->Draw("SAME");
  TH1D * h70 = (TH1D * ) _file0->Get("tanLAvsModule_layer1_bias70");
  h70->SetMarkerStyle(20);
  h70->SetMarkerColor(4);
  h70->SetLineColor(4);
  h70->Draw("SAME");
  TLegend * leg = new TLegend(0.12, 0.75, 0.42, 0.85);
  leg->SetFillColor(0);
  leg->AddEntry(h100, "V_{b} = 100", "lpf");
  leg->AddEntry(h300, "V_{b} = 300", "lpf");
  leg->AddEntry(h70, "V_{b} = 70", "lpf");
  leg->Draw();


}
