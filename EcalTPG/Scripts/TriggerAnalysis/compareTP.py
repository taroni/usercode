import ROOT
import os
from array import array


from ROOT import TStyle, TLatex,  TH1F, TCanvas, TFile, TLegend


myStyle=TStyle()
myStyle.SetLabelSize(0.05,"xyz")
myStyle.SetTitleSize(0.06,"xyz")
myStyle.SetTitleOffset(1.0,"x")
myStyle.SetTitleOffset(1.2,"y")
myStyle.SetPalette(1)
myStyle.SetPadRightMargin(0.14)
myStyle.SetOptStat(0)


t1 = TLatex()
t1.SetNDC();
t1.SetTextAlign(26);
t1.SetTextSize(0.06);


##f1=TFile.Open("histoTPG_251643_eg12.root")
##f2=TFile.Open("histoTPG_000_eg12.root");
f1=TFile.Open("~/eos/cms/store/caf/user/taroni/TPG/histoTPG_254833_eg12.root")
#f2=TFile.Open("~/eos/cms/store/caf/user/taroni/TPG/newhistoTPG_251244_notranscorr_eg12.root")
f2=TFile.Open("~/eos/cms/store/caf/user/taroni/TPG/newhistoTPG_251244_withtranscorr_eg12.root")
f3=TFile.Open("etaphi_ee.root");

scaling=1.;

vec=[1.5,1.6,1.8,2.0,2.2,2.4,2.6,3.0]
eta1d_eep_1=TH1F("eta1d_eep_1","",7,array('d',vec));
eta1d_eep_2=TH1F("eta1d_eep_2","",7,array('d',vec));
eta1d_eep_rat=TH1F("eta1d_eep_rat","",7,array('d',vec));

eta1d_eem_1=TH1F("eta1d_eem_1","",7,array('d',vec));
eta1d_eem_2=TH1F("eta1d_eem_2","",7,array('d',vec));
eta1d_eem_rat=TH1F("eta1d_eem_rat","",7,array('d',vec));

eep_spect_1=f1.Get("TPEEPlus");
eem_spect_1=f1.Get("TPEEMinus");
eb_spect_1 =f1.Get("TPEB");
eep_spect_1.SetName("TPEEPlusData")
eem_spect_1.SetName("TPEEMinusData")
eb_spect_1.SetName("TPEBData")

eep_spect_1.Rebin(4); 
eem_spect_1.Rebin(4);
eb_spect_1.Rebin(4);


eep_spect_2=f1.Get("TPEmulEEPlus");
eem_spect_2=f1.Get("TPEmulEEMinus");
eb_spect_2 =f1.Get("TPEmulEB");
eep_spect_2.SetName("TPEmulEEPlusData")
eem_spect_2.SetName("TPEmulEEMinusData")
eb_spect_2.SetName("TPEmulEBData")

eep_spect_2.Rebin(4); 
eem_spect_2.Rebin(4);
eb_spect_2.Rebin(4);


eep_spect_mc_1=f2.Get("TPEEPlus");
eem_spect_mc_1=f2.Get("TPEEMinus");
eb_spect_mc_1 =f2.Get("TPEB");
eep_spect_mc_1.Rebin(4); 
eem_spect_mc_1.Rebin(4);
eb_spect_mc_1.Rebin(4);

eep_spect_mc_2=f2.Get("TPEmulEEPlus");
eem_spect_mc_2=f2.Get("TPEmulEEMinus");
eb_spect_mc_2 =f2.Get("TPEmulEB");
eep_spect_mc_2.Rebin(4); 
eem_spect_mc_2.Rebin(4);
eb_spect_mc_2.Rebin(4);


c1=TCanvas("c1","c1",10,10,1000,700);
c1.Draw()
c1.SetLogy(1)
c1.SetGridy(1)
c1.SetGridx(1)


eep_spect_1.SetTitle("EE+ TP Spectrum")
eep_spect_1.Scale(1./eep_spect_1.Integral())
eep_spect_1.Draw()
eep_spect_mc_1.Scale(1./eep_spect_mc_1.Integral())
eep_spect_mc_1.Draw("SAME")
eep_spect_1.SetLineColor(1)
eep_spect_1.SetLineWidth(2)
eep_spect_mc_1.SetLineWidth(2)
eep_spect_mc_1.SetLineColor(2)

leg = TLegend(0.6,0.85, 0.75, 0.7)
leg.SetFillColor(0)
#leg.AddEntry(eep_spect_1, "Run 251643", "l")
#leg.AddEntry(eep_spect_mc_1, "MC", "l")
leg.AddEntry(eep_spect_1, "Run 254833", "l")
leg.AddEntry(eep_spect_mc_1, "Run 251244", "l")

leg.Draw()

c1.Update()
c1.SaveAs("plots/tpSpectra_EEP.png")

ratio_eep_spect_1 = eep_spect_1.Clone();
ratio_eep_spect_1.Divide(eep_spect_mc_1);
ratio_eep_spect_1.SetTitle("Run 254833 / Run 251244 TP Spectrum ratio, EE+"); 
ratio_eep_spect_1.SetName("Run254833_Run251244_TP_Spectrum_EEP_ratio"); 
ratio_eep_spect_1.SetLineColor(1);
ratio_eep_spect_1.SetLineWidth(2);
ratio_eep_spect_1.Draw();
c1.SaveAs("plots/tpSpectraRatio_EEP.png")

c1.SetLogy(0)
eep_spect_1.Add(eep_spect_mc_1, -1)
#eep_spect_1.SetTitle("EE+ Data TP - MC TP Spectrum")
eep_spect_1.SetTitle("EE+ Run 254833 TP - Run 251244 TP Spectrum")
eep_spect_1.Draw()
c1.Update()
c1.SaveAs("plots/tpSpectraDiff_EEP.png")
eep_spect_1.SetTitle("EE+ TP Spectrum")
c1.SetLogy(1)


eep_spect_2.SetTitle("EE+ Emul TP Spectrum")
eep_spect_2.Scale(1./eep_spect_2.Integral())
eep_spect_2.Draw()
eep_spect_mc_2.Scale(1./eep_spect_mc_2.Integral())
eep_spect_mc_2.Draw("SAME")
eep_spect_2.SetLineColor(1)
eep_spect_2.SetLineWidth(2)
eep_spect_mc_2.SetLineWidth(2)
eep_spect_mc_2.SetLineColor(2)
leg.Draw()

c1.Update()
c1.SaveAs("plots/tpEmulSpectra_EEP.png")

ratio_eep_spect_2 = eep_spect_2.Clone();
ratio_eep_spect_2.Divide(eep_spect_mc_2);
ratio_eep_spect_2.SetTitle("Run 254833 / Run 251244 TP Emul Spectrum ratio, EE+"); 
ratio_eep_spect_2.SetName("Run254833_Run251244_TP_Spectrum_EEP_ratio"); 
ratio_eep_spect_2.Draw();
ratio_eep_spect_2.SetLineColor(1);
ratio_eep_spect_2.SetLineWidth(2);

c1.SaveAs("plots/tpEmulSpectraRatio_EEP.png")

c1.SetLogy(0)
eep_spect_2.Add(eep_spect_mc_2, -1)
#eep_spect_2.SetTitle("EE+ Emul Data TP - Emul MC TP Spectrum")
eep_spect_2.SetTitle("EE+ Emul Run 254833 TP - Emul Run 251244 TP Spectrum")
eep_spect_2.Draw()
c1.Update()
c1.SaveAs("plots/tpEmulSpectraDiff_EEP.png")
eep_spect_2.SetTitle("EE+ Emul TP Spectrum")
c1.SetLogy(1)


#### EE-

eem_spect_1.SetTitle("EE- TP Spectrum")
eem_spect_1.Scale(1./eem_spect_1.Integral())
eem_spect_1.Draw()
eem_spect_mc_1.Scale(1./eem_spect_mc_1.Integral())
eem_spect_mc_1.Draw("SAME")
eem_spect_1.SetLineColor(1)
eem_spect_1.SetLineWidth(2)
eem_spect_mc_1.SetLineWidth(2)
eem_spect_mc_1.SetLineColor(2)

leg.Draw()

c1.Update()
c1.SaveAs("plots/tpSpectra_EEM.png")
ratio_eem_spect_1 = eem_spect_1.Clone();
ratio_eem_spect_1.Divide(eem_spect_mc_1);
ratio_eem_spect_1.SetTitle("Run 254833 / Run 251244 TP Spectrum ratio, EE-"); 
ratio_eem_spect_1.SetName("Run254833_Run251244_TP_Spectrum_EEM_ratio"); 
ratio_eem_spect_1.Draw();
ratio_eem_spect_1.SetLineColor(1);
ratio_eem_spect_1.SetLineWidth(2);

c1.SaveAs("plots/tpSpectraRatio_EEM.png")


c1.SetLogy(0)
eem_spect_1.Add(eem_spect_mc_1, -1)
#eem_spect_1.SetTitle("EE- Data TP - MC TP Spectrum")
eem_spect_1.SetTitle("EE-  Run 254833 TP - Run 251244 TP Spectrum")
eem_spect_1.Draw()
c1.Update()
c1.SaveAs("plots/tpSpectraDiff_EEM.png")
eem_spect_1.SetTitle("EE- TP Spectrum")
c1.SetLogy(1)



eem_spect_2.SetTitle("EE- Emul TP Spectrum")
eem_spect_2.Scale(1./eem_spect_2.Integral())
eem_spect_2.Draw()
eem_spect_mc_2.Scale(1./eem_spect_mc_2.Integral())
eem_spect_mc_2.Draw("SAME")
eem_spect_2.SetLineColor(1)
eem_spect_2.SetLineWidth(2)
eem_spect_mc_2.SetLineWidth(2)
eem_spect_mc_2.SetLineColor(2)
leg.Draw()

c1.Update()
c1.SaveAs("plots/tpEmulSpectra_EEM.png")

ratio_eem_spect_2 = eem_spect_2.Clone();
ratio_eem_spect_2.Divide(eem_spect_mc_2);
ratio_eem_spect_2.SetTitle("Run 254833 / Run 251244 TP Emul Spectrum ratio, EE-"); 
ratio_eem_spect_2.SetName("Run254833_Run251244_TP_EmulSpectrum_EEM_ratio"); 
ratio_eem_spect_2.SetLineColor(1);
ratio_eem_spect_2.SetLineWidth(2);

ratio_eem_spect_2.Draw();
c1.SaveAs("plots/tpEmulSpectraRatio_EEM.png")

c1.SetLogy(0)
eem_spect_2.Add(eem_spect_mc_2, -1)
#eem_spect_2.SetTitle("EE- Emul Data TP - Emul MC TP Spectrum")
eem_spect_2.SetTitle("EE- Emul Run 254833 TP - Emul Run 251244 TP Spectrum")

eem_spect_2.Draw()
c1.Update()
c1.SaveAs("plots/tpEmulSpectraDiff_EEM.png")
eem_spect_2.SetTitle("EE- Emul TP Spectrum")
c1.SetLogy(1)



### EB
eb_spect_1.SetTitle("EB TP Spectrum")
eb_spect_1.Scale(1./eb_spect_1.Integral())
eb_spect_1.Draw()
eb_spect_mc_1.Scale(1./eb_spect_mc_1.Integral())
eb_spect_mc_1.Draw("SAME")
eb_spect_1.SetLineColor(1)
eb_spect_1.SetLineWidth(2)
eb_spect_mc_1.SetLineWidth(2)
eb_spect_mc_1.SetLineColor(2)

leg.Draw()

c1.Update()
c1.SaveAs("plots/tpSpectra_EB.png")

ratio_eb_spect_1 = eb_spect_1.Clone();
ratio_eb_spect_1.Divide(eb_spect_mc_1);
ratio_eb_spect_1.SetTitle("Run 254833 / Run 251244 TP Spectrum ratio, EB"); 
ratio_eb_spect_1.SetName("Run254833_Run251244_TP_Spectrum_EB_ratio"); 
ratio_eb_spect_1.Draw();
c1.SaveAs("plots/tpSpectraRatio_EB.png")

c1.SetLogy(0)
eb_spect_1.Add(eb_spect_mc_1, -1)
#eb_spect_1.SetTitle("EB Data TP - MC TP Spectrum")
eb_spect_1.SetTitle("EB Run 254833 TP - Run 251244 TP Spectrum")
eb_spect_1.Draw()
c1.Update()
c1.SaveAs("plots/tpSpectraDiff_EB.png")
eb_spect_1.SetTitle("EB TP Spectrum")
c1.SetLogy(1)



eb_spect_2.SetTitle("EB Emul TP Spectrum")
eb_spect_2.Scale(1./eb_spect_2.Integral())
eb_spect_2.Draw()
eb_spect_mc_2.Scale(1./eb_spect_mc_2.Integral())
eb_spect_mc_2.Draw("SAME")
eb_spect_2.SetLineColor(1)
eb_spect_2.SetLineWidth(2)
eb_spect_mc_2.SetLineWidth(2)
eb_spect_mc_2.SetLineColor(2)
leg.Draw()

c1.Update()
c1.SaveAs("plots/tpEmulSpectra_EB.png")
ratio_eb_spect_2 = eb_spect_2.Clone();
ratio_eb_spect_2.Divide(eb_spect_mc_2);
ratio_eb_spect_2.SetTitle("Run 254833 / Run 251244 TP Emul Spectrum ratio, EB"); 
ratio_eb_spect_2.SetName("Run254833_Run251244_TP_EmulSpectrum_EB_ratio"); 
ratio_eb_spect_2.Draw();
c1.SaveAs("plots/tpEmulSpectraRatio_EB.png")

c1.SetLogy(0)
eb_spect_2.Add(eb_spect_mc_2, -1)
eb_spect_2.SetTitle("EB Emul Data TP - Emul MC TP Spectrum")
eb_spect_2.Draw()
c1.Update()
c1.SaveAs("plots/tpEmulSpectraDiff_EB.png")
eb_spect_2.SetTitle("EB Emul TP Spectrum")
c1.SetLogy(1)
