
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
TString Region="";

void MttRecoLevel()
{

  setTDRStyle();
  tdrStyle->SetErrorX(0.0);
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cout.precision(3);

  TString observablename;
  TString xaxislabel;

  Float_t observable, observable_gen, ttmass, ttRapidity, tmass;
  Float_t observableMinus, observableMinus_gen; 
  Double_t weight;
  Int_t Nsolns;

  int nVars =8;

  for (Int_t iVar= 5; iVar < 6; iVar++) {
    Initialize1DBinning(iVar);
    xaxislabel = "M_{t#bar{t}}";
    bool combineLepMinus = acceptanceName=="lepCosTheta" ? true : false;

    //TH1D* hData= new TH1D ("Data", "Data",    12, 195, 1095);
    //TH1D* hBkg = new TH1D ("Background",  "Background",    12, 195, 1095);
    //TH1D* hTop = new TH1D ("Top",  "Top",    12, 195, 1095);

    TH1D* hData= new TH1D ("Data", "Data",    10, 345, 1095);
    TH1D* hBkg = new TH1D ("Background",  "Background",    10, 345, 1095);
    TH1D* hTop = new TH1D ("Top",  "Top",    10, 345, 1095);

    hData->Sumw2();
    hTop->Sumw2();
    hBkg->Sumw2();
    
    TChain *ch_data = new TChain("tree");
    ch_data->Add("data.root");
    ch_data->SetBranchAddress(observablename,    &observable);
    if( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
    ch_data->SetBranchAddress("weight",&weight);
    ch_data->SetBranchAddress("Nsolns",&Nsolns);
    ch_data->SetBranchAddress("tt_mass",&ttmass);
    ch_data->SetBranchAddress("ttRapidity",&ttRapidity);
    ch_data->SetBranchAddress("t_mass",&tmass);

    for (Int_t i= 0; i<ch_data->GetEntries(); i++) {
      ch_data->GetEntry(i);
      if ( (acceptanceName=="lepChargeAsym") || (acceptanceName=="lepAzimAsym") || (acceptanceName=="lepAzimAsym2") ) {
          // leptonic asymmetries don't need valid top mass solution
        fillUnderOverFlow(hData, ttmass, weight, Nsolns);    
      } else {
        if ( ttmass > 0 ) {
            // asymmetries with top properties are required to have a valid top mass solution
          fillUnderOverFlow(hData, ttmass, weight, Nsolns);    
        }
      }
      if (combineLepMinus) {
        // combine plus and minus
        if ( (acceptanceName=="lepChargeAsym") || (acceptanceName=="lepAzimAsym") || (acceptanceName=="lepAzimAsym2") ) {
            // leptonic asymmetries don't need valid top mass solution
          fillUnderOverFlow(hData, ttmass, weight, Nsolns);    
        } else {
          if ( ttmass > 0 ) {
              // asymmetries with top properties are required to have a valid top mass solution
            fillUnderOverFlow(hData, ttmass, weight, Nsolns);    
          }
        }
      }    
    }

    TChain *ch_top = new TChain("tree");
    ch_top->Add("ttdil.root");
    ch_top->SetBranchAddress(observablename,    &observable);
    if( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
    ch_top->SetBranchAddress("weight",&weight);
    ch_top->SetBranchAddress("Nsolns",&Nsolns);
    ch_top->SetBranchAddress("tt_mass",&ttmass);
    ch_top->SetBranchAddress("ttRapidity",&ttRapidity);
    ch_top->SetBranchAddress("t_mass",&tmass);

    for (Int_t i= 0; i<ch_top->GetEntries(); i++) {
      ch_top->GetEntry(i);
      if ( (acceptanceName=="lepChargeAsym") || (acceptanceName=="lepAzimAsym") || (acceptanceName=="lepAzimAsym2") ) {
          // leptonic asymmetries don't need valid top mass solution
        fillUnderOverFlow(hTop, ttmass, weight, Nsolns);    
      } else {
        if ( ttmass > 0 ) {
            // asymmetries with top properties are required to have a valid top mass solution
          fillUnderOverFlow(hTop, ttmass, weight, Nsolns);    
        }
      }
      if (combineLepMinus) {
        // combine plus and minus
        if ( (acceptanceName=="lepChargeAsym") || (acceptanceName=="lepAzimAsym") || (acceptanceName=="lepAzimAsym2") ) {
            // leptonic asymmetries don't need valid top mass solution
          fillUnderOverFlow(hTop, ttmass, weight, Nsolns);    
        } else {
          if ( ttmass > 0 ) {
              // asymmetries with top properties are required to have a valid top mass solution
            fillUnderOverFlow(hTop, ttmass, weight, Nsolns);    
          }
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
    ch_bkg->SetBranchAddress("tt_mass",&ttmass);
    ch_bkg->SetBranchAddress("ttRapidity",&ttRapidity);
    ch_bkg->SetBranchAddress("t_mass",&tmass);

    for (Int_t i= 0; i<ch_bkg->GetEntries(); i++) {
      ch_bkg->GetEntry(i);
      if ( (acceptanceName=="lepChargeAsym") || (acceptanceName=="lepAzimAsym") || (acceptanceName=="lepAzimAsym2") ) {
          // leptonic asymmetries don't need valid top mass solution
        fillUnderOverFlow(hBkg, ttmass, weight, Nsolns);    
      } else {
        if ( ttmass > 0 ) {
            // asymmetries with top properties are required to have a valid top mass solution
          fillUnderOverFlow(hBkg, ttmass, weight, Nsolns);    
        }
      }
      if (combineLepMinus) {
        // combine plus and minus
        if ( (acceptanceName=="lepChargeAsym") || (acceptanceName=="lepAzimAsym") || (acceptanceName=="lepAzimAsym2") ) {
            // leptonic asymmetries don't need valid top mass solution
          fillUnderOverFlow(hBkg, ttmass, weight, Nsolns);    
        } else {
          if ( ttmass > 0 ) {
              // asymmetries with top properties are required to have a valid top mass solution
            fillUnderOverFlow(hBkg, ttmass, weight, Nsolns);    
          }
        }
      }    
    }

    double topScalingFactor = (hData->Integral() - hBkg->Integral()) / hTop->Integral();
    //cout<<"topScalingFactor: "<<_topScalingFactor<<" "<<topScalingFactor<<" "<< _topScalingFactor - topScalingFactor <<endl;

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


    TCanvas* c_test = new TCanvas("c_final","c_final",500,500);
    c_test->SetLeftMargin(0.2);

    hData->SetMinimum(0.0);
    hData->SetMaximum( 1.1* hData->GetMaximum());
    hData->GetXaxis()->SetTitle(xaxislabel + " (GeV)" );
    hData->GetYaxis()->SetTitle("Events/(75 GeV)");
    hData->GetYaxis()->SetTitle("Events/(75 GeV)");
    hData->GetYaxis()->SetTitleOffset(1.6);
    hData->SetLineWidth(lineWidth);
    hData->SetMarkerSize(0.7);

    hTop->SetLineWidth(lineWidth);
    hTop->SetMarkerStyle(20);
    hTop->SetMarkerColor(TColor::GetColorDark(kRed));
    hTop->SetMarkerSize(0.0);
    hTop->SetLineColor(TColor::GetColorDark(kRed));
    //hTop->SetFillColor(10);
    hTop->SetFillStyle(0);

    hBkg->SetLineWidth(lineWidth);
    hBkg->SetMarkerStyle(20);
    hBkg->SetMarkerColor(TColor::GetColorDark(kBlue));
    hBkg->SetMarkerSize(0.0);
    hBkg->SetLineColor(TColor::GetColorDark(kBlue));
    hBkg->SetFillColor(TColor::GetColorDark(kBlue));
    hBkg->SetFillStyle(3353);

    THStack *hs = new THStack("hs","Stacked Top+BG");
    hs->Add(hBkg);
    hs->Add(hTop);
    hs->SetMinimum(0.0);
    hs->SetMaximum( 1.1* hs->GetMaximum());
    hs->Draw("hist");
    hs->GetXaxis()->SetTitle(xaxislabel + " (GeV)" );
    hs->GetYaxis()->SetTitle("Events/(75 GeV)");
    hs->GetYaxis()->SetTitle("Events/(75 GeV)");
    hs->GetYaxis()->SetTitleOffset(1.6);
    //hs->GetXaxis()->SetNdivisions(505,0);

    hData->Draw("E same");

    //TLegend* leg1=new TLegend(0.6,0.62,0.9,0.838,NULL,"brNDC");
    //TLegend* leg1=new TLegend(0.68, 0.75, 0.95, 0.93, NULL, "brNDC");
    TLegend *leg1 = new TLegend(0.66, 0.755, 0.95, 0.935, NULL, "brNDC");
    leg1->SetFillStyle(0);
    leg1->SetEntrySeparation(100);  
    leg1->SetBorderSize(0);                                                                                 
    leg1->SetTextSize(0.036);
    leg1->AddEntry(hData, "Data", "P");
    leg1->AddEntry(hTop,  "t#bar{t} (dileptonic)");                                                               
    leg1->AddEntry(hBkg,  "Background");                                                               
    leg1->Draw();                

    //TPaveText *pt1 = new TPaveText(0.22, 0.88, 0.45, 0.91, "brNDC");
    TPaveText *pt1 = new TPaveText(0.19, 0.94, 0.41, 0.98, "brNDC");
    pt1->SetName("pt1name");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);

    TText *blah;
    //blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
    blah = pt1->AddText("CMS, #sqrt{s} = 7 TeV, 5.0 fb^{-1}");
    blah->SetTextSize(0.036);
    blah->SetTextAlign(11);

    pt1->Draw();
    c_test->SaveAs("Mtt_reco.pdf");
    c_test->SaveAs("Mtt_reco.C");

    delete hData;
    delete hBkg;
    delete hTop;
    delete hTop_bkgAdd;

  }

}

#ifndef __CINT__
int main () { AfbRecoLevel(); return 0; }  // Main program when run stand-alone
#endif
