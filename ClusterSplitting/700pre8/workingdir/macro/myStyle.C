#ifndef __mystyle_H__
#define __mystyle_H__

#include "TStyle.h"
#include "TPaveText.h"


TStyle * mystyle()
{

 TStyle * style_ = new TStyle("GoodStyle", "");

 style_ = new TStyle("DrawBaseStyle", "");
 style_->SetCanvasColor(0);
 style_->SetPadColor(0);
 style_->SetFrameFillColor(0);
 style_->SetStatColor(0);
 style_->SetOptStat(0);
 style_->SetTitleFillColor(0);
 style_->SetCanvasBorderMode(0);
 style_->SetPadBorderMode(0);
 style_->SetFrameBorderMode(0);
 style_->SetPadBottomMargin(0.1);
 style_->SetPadLeftMargin(0.1);
 style_->cd();

 // For the canvas:
   style_->SetCanvasBorderMode(0);
   style_->SetCanvasColor(kWhite);
   style_->SetCanvasDefH(600); //Height of canvas
   style_->SetCanvasDefW(600); //Width of canvas
   style_->SetCanvasDefX(0);   //POsition on screen
   style_->SetCanvasDefY(0);

 // For the Pad:
   style_->SetPadBorderMode(0);
   // style_->SetPadBorderSize(Width_t size = 1);
   style_->SetPadColor(kWhite);
   style_->SetPadGridX(false);
   style_->SetPadGridY(false);
   style_->SetGridColor(0);
   style_->SetGridStyle(3);
   style_->SetGridWidth(1);

 // For the frame:
   style_->SetFrameBorderMode(0);
   style_->SetFrameBorderSize(1);
   style_->SetFrameFillColor(0);
   style_->SetFrameFillStyle(0);
   style_->SetFrameLineColor(1);
   style_->SetFrameLineStyle(1);
   style_->SetFrameLineWidth(1);

//// For the histo:
//  // style_->SetHistFillColor(1);
//  // style_->SetHistFillStyle(0);
//  style_->SetHistLineColor(1);
//  style_->SetHistLineStyle(0);
//  style_->SetHistLineWidth(1);
//  // style_->SetLegoInnerR(Float_t rad = 0.5);
//  // style_->SetNumberContours(Int_t number = 20);

//  style_->SetEndErrorSize(2);
////  style_->SetErrorMarker(20);
//  style_->SetErrorX(0.);
//  
//  style_->SetMarkerStyle(20);

////For the fit/function:
//  style_->SetOptFit(1);
//  style_->SetFitFormat("5.4g");
//  style_->SetFuncColor(2);
//  style_->SetFuncStyle(1);
//  style_->SetFuncWidth(1);

////For the date:
//  style_->SetOptDate(0);
//  // style_->SetDateX(Float_t x = 0.01);
//  // style_->SetDateY(Float_t y = 0.01);

//// For the statistics box:
//  style_->SetOptFile(0);
//  style_->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
//  style_->SetStatColor(kWhite);
//  style_->SetStatFont(42);
//  style_->SetStatFontSize(0.025);
//  style_->SetStatTextColor(1);
//  style_->SetStatFormat("6.4g");
//  style_->SetStatBorderSize(1);
//  style_->SetStatH(0.1);
//  style_->SetStatW(0.15);
//  // style_->SetStatStyle(Style_t style = 1001);
//  // style_->SetStatX(Float_t x = 0);
//  // style_->SetStatY(Float_t y = 0);

 // Margins:
   style_->SetPadTopMargin(0.05);
   style_->SetPadBottomMargin(0.15);//0.13);
   style_->SetPadLeftMargin(0.15);//0.16);
   style_->SetPadRightMargin(0.05);//0.02);

 // For the Global title:

   style_->SetOptTitle(0);
   style_->SetTitleFont(42);
   style_->SetTitleColor(1);
   style_->SetTitleTextColor(1);
   style_->SetTitleFillColor(10);
   style_->SetTitleFontSize(0.05);
   // style_->SetTitleH(0); // Set the height of the title box
   // style_->SetTitleW(0); // Set the width of the title box
   // style_->SetTitleX(0); // Set the position of the title box
   // style_->SetTitleY(0.985); // Set the position of the title box
   // style_->SetTitleStyle(Style_t style = 1001);
   // style_->SetTitleBorderSize(2);

 // For the axis titles:

   style_->SetTitleColor(1, "XYZ");
   style_->SetTitleFont(42, "XYZ");
   style_->SetTitleSize(0.05, "XYZ");
   // style_->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
   // style_->SetTitleYSize(Float_t size = 0.02);
   style_->SetTitleXOffset(1.15);//0.9);
/*     style_->SetTitleYOffset(1.4); // => 1.15 if exponents */
   style_->SetTitleYOffset(1.6); // => 1.15 if exponents
   // style_->SetTitleYOffset(0.7);
   // style_->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

 // For the axis labels:

   style_->SetLabelColor(1, "XYZ");
   style_->SetLabelFont(42, "XYZ");
   style_->SetLabelOffset(0.007, "XYZ");
   style_->SetLabelSize(0.045, "XYZ");

 // For the axis:

   style_->SetAxisColor(1, "XYZ");
   style_->SetStripDecimals(kTRUE);
   style_->SetTickLength(0.03, "XYZ");
   style_->SetNdivisions(510, "XYZ");
   style_->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
   style_->SetPadTickY(1);

//// Change for log plots:
//  style_->SetOptLogx(0);
//  style_->SetOptLogy(0);
//  style_->SetOptLogz(0);

//// Postscript options:
//  style_->SetPaperSize(20.,20.);
//  // style_->SetLineScalePS(Float_t scale = 3);
//  // style_->SetLineStyleString(Int_t i, const char* text);
//  // style_->SetHeaderPS(const char* header);
//  // style_->SetTitlePS(const char* pstitle);

//  // style_->SetBarOffset(Float_t baroff = 0.5);
//  // style_->SetBarWidth(Float_t barwidth = 0.5);
//  // style_->SetPaintTextFormat(const char* format = "g");
//  // style_->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
//  // style_->SetTimeOffset(Double_t toffset);
//  // style_->SetHistMinimumZero(kTRUE);

//  // Additional settings for QCD-10-011
  style_->SetLegendBorderSize(0);


   return style_;

}



TStyle * rectStyle()
{
 TStyle * style_ = new TStyle("GoodStyle", "");



 style_ = new TStyle("DrawBaseStyle", "");
 style_->SetCanvasColor(0);
 style_->SetPadColor(0);
 style_->SetFrameFillColor(0);
 style_->SetStatColor(0);
 style_->SetOptStat(0);
 style_->SetTitleFillColor(0);
 style_->SetCanvasBorderMode(0);
 style_->SetPadBorderMode(0);
 style_->SetFrameBorderMode(0);
 style_->SetPadBottomMargin(0.12);
 style_->SetPadLeftMargin(0.12);
 style_->cd();

 // For the canvas:
 style_->SetCanvasBorderMode(0);
 style_->SetCanvasColor(kWhite);
 style_->SetCanvasDefX(0);   //POsition on screen
 style_->SetCanvasDefY(0);

 // For the Pad:
 style_->SetPadBorderMode(0);
 style_->SetPadColor(kWhite);
 style_->SetPadGridX(false);
 style_->SetPadGridY(false);
 style_->SetGridColor(0);
 style_->SetGridStyle(3);
 style_->SetGridWidth(1);

 // For the frame:
 style_->SetFrameBorderMode(0);
 style_->SetFrameBorderSize(1);
 style_->SetFrameFillColor(0);
 style_->SetFrameFillStyle(0);
 style_->SetFrameLineColor(1);
 style_->SetFrameLineStyle(1);
 style_->SetFrameLineWidth(1);


 // Margins:
 style_->SetPadTopMargin(0.05);
 style_->SetPadBottomMargin(0.15);//0.13);
 style_->SetPadLeftMargin(0.15);//0.16);
 style_->SetPadRightMargin(0.05);//0.02);

 // For the Global title:

 style_->SetOptTitle(0); /// TITLE DISAPPEARS
 style_->SetTitleFont(42);
 style_->SetTitleColor(1);
 style_->SetTitleTextColor(1);
 style_->SetTitleFillColor(10);
 style_->SetTitleFontSize(0.05);

 // For the axis titles:

 style_->SetTitleColor(1, "XYZ");
 style_->SetTitleFont(42, "XYZ");
 style_->SetTitleSize(0.05, "XYZ");
 style_->SetTitleXOffset(1.15);//0.9);
 style_->SetTitleYOffset(1.4); // => 1.15 if exponents

 // For the axis labels:

 style_->SetLabelColor(1, "XYZ");
 style_->SetLabelFont(42, "XYZ");
 style_->SetLabelOffset(0.007, "XYZ");
 style_->SetLabelSize(0.045, "XYZ");

 // For the axis:

 style_->SetAxisColor(1, "XYZ");
 style_->SetStripDecimals(kTRUE);
 style_->SetTickLength(0.03, "XYZ");
 style_->SetNdivisions(510, "XYZ");
 style_->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
 style_->SetPadTickY(1);

 return style_;

}







void setTitle(const std::string & lumi = "188")
{
//   TPaveText* cmslabel = new TPaveText( 0.145, 0.953, 0.6, 0.975, "brNDC" );
  TPaveText* cmsbox = new TPaveText( 0.13, 0.97, 0.6, .99999, "brNDC" );
  cmsbox->SetFillColor(0);
  cmsbox->SetBorderSize(0);
  cmsbox->Draw();
  TPaveText* cmslabel = new TPaveText( 0.13, 0.97, 0.6, .95, "brNDC" );
  // TPaveText* cmslabel = new TPaveText( 0.23, 0.953, 0.6, 0.975, "brNDC" );
 cmslabel->SetFillColor(kWhite);
 cmslabel->SetBorderSize(0);
 cmslabel->SetTextSize(0.035);
 //  if( legendQuadrant==0 ) cmslabel->SetTextAlign(11);
 cmslabel->SetTextSize(0.035);
 cmslabel->SetTextFont(62);

 std::string leftText;
 leftText = "";
 // leftText = "CMS Simulation - CMS Preliminary 2012";

 cmslabel->SetTextAlign(11); //center //11 align left
 // cmslabel->AddText(Form("%s, %s fb^{-1}", leftText.c_str(), lumi.c_str() ));
 cmslabel->AddText(Form("%s %s", leftText.c_str(), lumi.c_str() ));

 cmslabel->Draw();

 TPaveText* label_sqrt = new TPaveText(0.83,0.97,0.96,0.95, "brNDC");
 label_sqrt->SetFillColor(kWhite);
 label_sqrt->SetTextSize(0.035);
 label_sqrt->SetTextFont(42);
 label_sqrt->SetBorderSize(0);
 label_sqrt->SetTextAlign(31); // align right
 label_sqrt->AddText("#sqrt{s} = 8 TeV");

 label_sqrt->Draw();

}


#endif

