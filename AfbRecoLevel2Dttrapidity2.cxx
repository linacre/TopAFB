
/// Script for plotting Afb variables at the reco level
// S. Jindariani, J.Linacre, Y.Tu

#include <iostream>
using std::cout;
using std::endl;

#include "TRandom3.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TChain.h"
#include "TLegend.h"
#include "TColor.h"
#include "THStack.h"
#include "TCut.h"
#include "TText.h"
#include "TPaveText.h"
#include "RooUnfold/tdrstyle.C"

#include "RooUnfold/examples/AfbFinalUnfold.h"

std::string formatFloat(double x, const char* formatS) {
  std::string xS = Form(Form("%s", formatS),x);
  double xB = atof(xS.c_str());
  if (x>0 && xB==0){
    xS = Form(" %6.1g",x);
  }
  return xS;
}

//==============================================================================
// Global definitions
//==============================================================================

const Double_t _topScalingFactor=1.+(9824. - 10063.47)/9323.84;


void AfbRecoLevel2Dttrapidity2()
{

  setTDRStyle();
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cout.precision(3);

  TString observablename;
  TString xaxislabel;

  Float_t observable, ttrapidity2, absttrapidity2;
  Float_t observableMinus; 
  Double_t weight;
  Int_t Nsolns;

  int nVars = 9;

  for (Int_t iVar= 0; iVar < nVars; iVar++) {
    Initialize2DBinningttrapidity2(iVar);
    bool combineLepMinus = acceptanceName=="lepCosTheta" ? true : false;

    TH1D* hData= new TH1D ("Data", "Data",    nbins2D, xbins2D);
    TH1D* hBkg = new TH1D ("Background",  "Background",    nbins2D, xbins2D);
    TH1D* hTop = new TH1D ("Top",  "Top",    nbins2D, xbins2D);

    hData->Sumw2();
    hTop->Sumw2();
    hBkg->Sumw2();
    
    TChain *ch_data = new TChain("tree");
    ch_data->Add("data.root");
    ch_data->SetBranchAddress(observablename,    &observable);
    if( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
    ch_data->SetBranchAddress("weight",&weight);
    ch_data->SetBranchAddress("Nsolns",&Nsolns);
    ch_data->SetBranchAddress("ttRapidity2",&ttrapidity2);

    for (Int_t i= 0; i<ch_data->GetEntries(); i++) {
      ch_data->GetEntry(i);
      absttrapidity2=sqrt(ttrapidity2*ttrapidity2);
      if ( absttrapidity2 < 900 ) { // because sevents without top mass solution have ttRapidity2 = -999
        if(observablename=="lep_azimuthal_asymmetry2") observable = -cos(observable);
        fillUnderOverFlow(hData, sign(observable)*absttrapidity2, weight, Nsolns);    
        if (combineLepMinus) {
          fillUnderOverFlow(hData, sign(observableMinus)*absttrapidity2, weight, Nsolns);    
        }    
      }
    }

    TChain *ch_top = new TChain("tree");
    ch_top->Add("ttdil.root");
    ch_top->SetBranchAddress(observablename,    &observable);
    if( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
    ch_top->SetBranchAddress("weight",&weight);
    ch_top->SetBranchAddress("Nsolns",&Nsolns);
    ch_top->SetBranchAddress("ttRapidity2",&ttrapidity2);

    for (Int_t i= 0; i<ch_top->GetEntries(); i++) {
      ch_top->GetEntry(i);
      absttrapidity2=sqrt(ttrapidity2*ttrapidity2);
      if(observablename=="lep_azimuthal_asymmetry2") observable = -cos(observable);
      if ( absttrapidity2 < 900 ) { // because sevents without top mass solution have ttRapidity2 = -999
        fillUnderOverFlow(hTop, sign(observable)*absttrapidity2, weight, Nsolns);
        if (combineLepMinus) {
          fillUnderOverFlow(hTop, sign(observableMinus)*absttrapidity2, weight, Nsolns);
        }
      }
    }

    TChain *ch_bkg = new TChain("tree");
    ch_bkg->Add("ttotr.root");
    ch_bkg->Add("wjets.root");
    ch_bkg->Add("DYee.root");
    ch_bkg->Add("DYmm.root");
    ch_bkg->Add("DYtautau.root");
    ch_bkg->Add("tw.root");
    ch_bkg->Add("VV.root");
    ch_bkg->SetBranchAddress(observablename,    &observable);
    if( combineLepMinus ) ch_bkg->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
    ch_bkg->SetBranchAddress("weight",&weight);
    ch_bkg->SetBranchAddress("Nsolns",&Nsolns);
    ch_bkg->SetBranchAddress("ttRapidity2",&ttrapidity2);

    for (Int_t i= 0; i<ch_bkg->GetEntries(); i++) {
      ch_bkg->GetEntry(i);
      absttrapidity2=sqrt(ttrapidity2*ttrapidity2);
      if(observablename=="lep_azimuthal_asymmetry2") observable = -cos(observable);
      if ( absttrapidity2 < 900 ) { // because sevents without top mass solution have ttRapidity2 = -999
        fillUnderOverFlow(hBkg, sign(observable)*absttrapidity2, weight, Nsolns);
        if( combineLepMinus ) {
          fillUnderOverFlow(hBkg, sign(observableMinus)*absttrapidity2, weight, Nsolns);
        }
      }
    }

    hTop->Scale(_topScalingFactor);
    TH1D* hTop_bkgAdd= (TH1D*) hTop->Clone();
    hTop_bkgAdd->Add(hBkg);  

  //==================================================================
  // ============== Print the assymetry =============================
    cout<<"========= Variable:"<<observablename <<" ===================\n";
    Float_t Afb, AfbErr;
    const char* formatS = "%6.3f";

    GetAfb(hData, Afb, AfbErr);
    cout<<" Data: "<< formatFloat(Afb,formatS) <<" +/-  "<< formatFloat(AfbErr,formatS)<<"\n";
    GetAfb(hTop_bkgAdd, Afb, AfbErr);
    cout<<" True Top: "<<formatFloat(Afb,formatS) <<" +/-  "<< formatFloat(AfbErr,formatS)<<"\n";

    // hData->Scale(1./hData->Integral(),"width");
    // hTop->Scale(1./(hTop->Integral()+hBkg->Integral()),"width");
    // hBkg->Scale(1./(hTop->Integral()+hBkg->Integral()),"width");

    for (Int_t i= 1; i<=nbins2D; i++) {
      hData->SetBinContent(i, hData->GetBinContent(i)/hData->GetBinWidth(i) );
      hData->SetBinError(i, hData->GetBinError(i)/hData->GetBinWidth(i) );
      hTop->SetBinContent(i, hTop->GetBinContent(i)/hTop->GetBinWidth(i) );
      hTop->SetBinError(i, hTop->GetBinError(i)/hTop->GetBinWidth(i) );
      hBkg->SetBinContent(i, hBkg->GetBinContent(i)/hBkg->GetBinWidth(i) );
      hBkg->SetBinError(i, hBkg->GetBinError(i)/hBkg->GetBinWidth(i) );
    }

    TCanvas* c_test = new TCanvas("c_final","c_final",500,500);
    c_test->SetLeftMargin(0.2);

    hData->SetMinimum(0.0);
    hData->SetMaximum( 1.35* hData->GetMaximum());
    hData->GetXaxis()->SetTitle(yaxislabel + " #times sign(" + xaxislabel + ")");
    hData->GetYaxis()->SetTitle("Events/(" + yaxislabel + " #times sign(" + xaxislabel + ") )");
    hData->GetYaxis()->SetTitleOffset(1.6);
    hData->SetLineWidth(lineWidth);

    hTop->SetLineWidth(lineWidth);
    hTop->SetMarkerSize(0.0);
    hTop->SetLineColor(TColor::GetColorDark(kGreen));
    hTop->SetFillColor(TColor::GetColorDark(kGreen));
    hTop->SetFillStyle(3353);

    hBkg->SetLineWidth(lineWidth);
    hBkg->SetMarkerSize(0.0);
    hBkg->SetLineColor(kYellow);
    hBkg->SetFillColor(kYellow);

    THStack *hs = new THStack("hs","Stacked Top+BG");
    hs->Add(hBkg);
    hs->Add(hTop);
    hs->SetMinimum(0.0);
    hs->SetMaximum( 1.35* hs->GetMaximum());
    hs->Draw("hist");
    hs->GetXaxis()->SetTitle(yaxislabel + " #times sign(" + xaxislabel + ")");
    hs->GetYaxis()->SetTitle("Events/(" + yaxislabel + " #times sign(" + xaxislabel + ") )");
    hs->GetYaxis()->SetTitleOffset(1.6);

    hData->Draw("E same");

    //TLegend* leg1=new TLegend(0.6,0.62,0.9,0.838,NULL,"brNDC");
    TLegend* leg1=new TLegend(0.68, 0.75, 0.95, 0.93, NULL, "brNDC");
    leg1->SetFillStyle(0);
    leg1->SetEntrySeparation(100);  
    leg1->SetBorderSize(0);                                                                                 
    leg1->SetTextSize(0.03);
    leg1->AddEntry(hData, "Data");
    leg1->AddEntry(hTop,  "t#bar{t} (dileptonic)");                                                               
    leg1->AddEntry(hBkg,  "Background");                                                               
    leg1->Draw();                

    TPaveText *pt1 = new TPaveText(0.22, 0.88, 0.45, 0.91, "brNDC");
    pt1->SetName("pt1name");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);

    TText *blah;
    //blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at #sqrt{s}=7 TeV");
    blah = pt1->AddText("CMS, 5.0 fb^{-1} at #sqrt{s}=7 TeV");
    blah->SetTextSize(0.035);
    blah->SetTextAlign(11);

    pt1->Draw();
    c_test->SaveAs("2D_ttrapidity2_reco_"+acceptanceName+".pdf");
    c_test->SaveAs("2D_ttrapidity2_reco_"+acceptanceName+".C");

    delete hData;
    delete hBkg;
    delete hTop;
    delete hTop_bkgAdd;

  }

}

#ifndef __CINT__
int main () { AfbRecoLevel(); return 0; }  // Main program when run stand-alone
#endif
