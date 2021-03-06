#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TLatex.h"
#include <fstream>
#include "tdrStyle.C"
#include "CommonFunctions.C"
#include <vector>

using namespace std;

TH1D* hnumerator;
TH1D* hdenominator;
TH1D* hacceptance;
TH1D* hacceptance_copy;
TH1D* hacceptance_statup;
TH1D* hacceptance_statdown;

TH2D* hnumerator2d_mtt;
TH2D* hdenominator2d_mtt;
TH2D* hnumerator2drebinned_mtt;
TH2D* hdenominator2drebinned_mtt;
TH2D* hacceptance2drebinned_mtt;

TH2D* hnumerator2d_ttpt;
TH2D* hdenominator2d_ttpt;
TH2D* hnumerator2drebinned_ttpt;
TH2D* hdenominator2drebinned_ttpt;
TH2D* hacceptance2drebinned_ttpt;

TH2D* hnumerator2d_ttrapidity2;
TH2D* hdenominator2d_ttrapidity2;
TH2D* hnumerator2drebinned_ttrapidity2;
TH2D* hdenominator2drebinned_ttrapidity2;
TH2D* hacceptance2drebinned_ttrapidity2;


void GetAfb(TH1D* h, Double_t &afb, Double_t  &afberr){
 
  Int_t nbins = h->GetNbinsX();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  //event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"width");
  //event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_plus_err,"width");
  event_total = event_plus + event_minus;
  
  //cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<endl;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
         + event_minus*event_minus*event_plus_err*event_plus_err)/
      (event_total*event_total*event_total*event_total));

}



void acceptanceplots(TString histname = "lepAzimAsym", bool drawnorm = false, TString FName1 = "results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", TString FName2 = "results/hist_noCuts.root"){
  setTDRStyle();


  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  
  hnumerator = (TH1D*)f_1->Get(Form("ttdil_h%sGen_allj_all", histname.Data())); 
  hnumerator2d_mtt = (TH2D*)f_1->Get(Form("ttdil_h%sGen2d_allj_all", histname.Data())); 
  hnumerator2d_ttpt = (TH2D*)f_1->Get(Form("ttdil_h%sttpTGen2d_allj_all", histname.Data())); 
  hnumerator2d_ttrapidity2 = (TH2D*)f_1->Get(Form("ttdil_h%sttRapidity2Gen2d_allj_all", histname.Data())); 

  std::cout << "Opening " << FName2.Data() << "\n";  
  TFile *f_2         = TFile::Open(FName2.Data());  
  hdenominator = (TH1D*)f_2->Get(Form("ttdil_h%sGen_allj_all", histname.Data()));
  hdenominator2d_mtt = (TH2D*)f_2->Get(Form("ttdil_h%sGen2d_allj_all", histname.Data())); 
  hdenominator2d_ttpt = (TH2D*)f_2->Get(Form("ttdil_h%sttpTGen2d_allj_all", histname.Data())); 
  hdenominator2d_ttrapidity2 = (TH2D*)f_2->Get(Form("ttdil_h%sttRapidity2Gen2d_allj_all", histname.Data())); 


  std::cout << "Opened " << Form("ttdil_h%sGen_allj_all", histname.Data()) << " and "<< Form("ttdil_h%sGen2d_allj_all", histname.Data()) <<"\n";

  Double_t pi = 3.141592653589793;
  Double_t bins_lepChargeAsym[] =  { -2., -0.8, -0.4, 0., 0.4, 0.8, 2.}; 
  Double_t bins_lepAzimAsym[] = {-1., -0.8, -0.4, 0., 0.4, 0.8, 1.}; 
  Double_t bins_lepAzimAsym2[] = {0., 4.*pi/20., 7.*pi/20., 10.*pi/20., 13.*pi/20., 16.*pi/20., pi}; 
  Double_t bins_topCosTheta[] = {-1., -0.7, -0.4, 0., 0.4, 0.7, 1.}; 
  Double_t bins_pseudorapiditydiff[] =  { -2., -1.0, -0.5, 0., 0.5, 1.0, 2.}; 
  Double_t bins_rapiditydiff[] =  { -2., -0.8, -0.3, 0., 0.3, 0.8, 2.}; 
  Double_t bins_rapiditydiffMarco[] =  { -2., -0.7, -0.3, 0., 0.3, 0.7, 2.}; 
  Double_t bins_lepCosTheta[] = {-1., -0.6, -0.3, 0., 0.3, 0.6, 1.}; 
  Double_t bins_topSpinCorr[] = {-1., -0.5, -0.2, 0., 0.2, 0.5, 1.}; 
  Double_t bins_lepCosOpeningAngle[] = {-1., -0.6, -0.3, 0., 0.3, 0.6, 1.}; 

  Double_t bins_lepChargeAsym_for2D[] =  { -2., 0., 2.};
  Double_t bins_lepAzimAsym_for2D[] = {-1., 0., 1.};
  Double_t bins_lepAzimAsym2_for2D[] = {0., pi/2., pi};
  Double_t bins_topCosTheta_for2D[] = {-1., 0., 1.};
  Double_t bins_pseudorapiditydiff_for2D[] =  { -2., 0., 2.};
  Double_t bins_rapiditydiff_for2D[] =  { -2., 0., 2.};
  Double_t bins_rapiditydiffMarco_for2D[] =  { -2., 0., 2.};
  Double_t bins_lepCosTheta_for2D[] = {-1., 0., 1.};
  Double_t bins_topSpinCorr_for2D[] = {-1., 0., 1.};
  Double_t bins_lepCosOpeningAngle_for2D[] = {-1., 0., 1.};

  Double_t binsmtt[] = {0., 410., 510., 1200.}; 
  Double_t binsttpt[] = {0., 24., 52., 300}; 
  Double_t binsttrapidity2[] = {0., 0.3, 0.7, 3.0}; 
  Double_t bins[7];
  Double_t binsfor2D[3];


  if(histname == "lepChargeAsym") memcpy(bins,bins_lepChargeAsym,7*8);
  if(histname == "lepAzimAsym") memcpy(bins,bins_lepAzimAsym,7*8);
  if(histname == "lepAzimAsym2") memcpy(bins,bins_lepAzimAsym2,7*8);
  if(histname == "topCosTheta") memcpy(bins,bins_topCosTheta,7*8);
  if(histname == "pseudorapiditydiff") memcpy(bins,bins_pseudorapiditydiff,7*8);
  if(histname == "rapiditydiff") memcpy(bins,bins_rapiditydiff,7*8);
  if(histname == "rapiditydiffMarco") memcpy(bins,bins_rapiditydiffMarco,7*8);
  if(histname == "lepCosTheta" || histname == "lepPlusCosTheta" || histname == "lepMinusCosTheta") memcpy(bins,bins_lepCosTheta,7*8);
  if(histname == "topSpinCorr") memcpy(bins,bins_topSpinCorr,7*8);
  if(histname == "lepCosOpeningAngle") memcpy(bins,bins_lepCosOpeningAngle,7*8);

  if(histname == "lepChargeAsym") memcpy(binsfor2D,bins_lepChargeAsym_for2D,3*8);
  if(histname == "lepAzimAsym") memcpy(binsfor2D,bins_lepAzimAsym_for2D,3*8);
  if(histname == "lepAzimAsym2") memcpy(binsfor2D,bins_lepAzimAsym2_for2D,3*8);
  if(histname == "topCosTheta") memcpy(binsfor2D,bins_topCosTheta_for2D,3*8);
  if(histname == "pseudorapiditydiff") memcpy(binsfor2D,bins_pseudorapiditydiff_for2D,3*8);
  if(histname == "rapiditydiff") memcpy(binsfor2D,bins_rapiditydiff_for2D,3*8);
  if(histname == "rapiditydiffMarco") memcpy(binsfor2D,bins_rapiditydiffMarco_for2D,3*8);
  if(histname == "lepCosTheta" || histname == "lepPlusCosTheta" || histname == "lepMinusCosTheta") memcpy(binsfor2D,bins_lepCosTheta_for2D,3*8);
  if(histname == "topSpinCorr") memcpy(binsfor2D,bins_topSpinCorr_for2D,3*8);
  if(histname == "lepCosOpeningAngle") memcpy(binsfor2D,bins_lepCosOpeningAngle_for2D,3*8);

  hnumerator = (TH1D*) hnumerator->Rebin(6,Form("numerator_%s", histname.Data()),bins);
  hdenominator = (TH1D*) hdenominator->Rebin(6,Form("denominator_%s", histname.Data()),bins);

  hnumerator2drebinned_mtt = new TH2D(Form("numerator_%s_mtt", histname.Data()),Form("numerator_%s_mtt", histname.Data()),2,binsfor2D,3, binsmtt);
  TAxis *xaxis = hnumerator2d_mtt->GetXaxis();
  TAxis *yaxis = hnumerator2d_mtt->GetYaxis();
  for (int j=1;j<=yaxis->GetNbins();j++) {
    for (int i=1;i<=xaxis->GetNbins();i++) {
      hnumerator2drebinned_mtt->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d_mtt->GetBinContent(i,j));
    }
  }
  hdenominator2drebinned_mtt = new TH2D(Form("denominator_%s_mtt", histname.Data()),Form("denominator_%s_mtt", histname.Data()),2,binsfor2D,3, binsmtt);
  TAxis *xaxisd = hdenominator2d_mtt->GetXaxis();
  TAxis *yaxisd = hdenominator2d_mtt->GetYaxis();
  for (int j=1;j<=yaxisd->GetNbins();j++) {
    for (int i=1;i<=xaxisd->GetNbins();i++) {
      hdenominator2drebinned_mtt->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d_mtt->GetBinContent(i,j));
    }
  }

  hnumerator2drebinned_ttpt = new TH2D(Form("numerator_%s_ttpt", histname.Data()),Form("numerator_%s_ttpt", histname.Data()),2,binsfor2D,3, binsttpt);
  xaxis = hnumerator2d_ttpt->GetXaxis();
  yaxis = hnumerator2d_ttpt->GetYaxis();
  for (int j=1;j<=yaxis->GetNbins();j++) {
    for (int i=1;i<=xaxis->GetNbins();i++) {
      hnumerator2drebinned_ttpt->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d_ttpt->GetBinContent(i,j));
    }
  }
  hdenominator2drebinned_ttpt = new TH2D(Form("denominator_%s_ttpt", histname.Data()),Form("denominator_%s_ttpt", histname.Data()),2,binsfor2D,3, binsttpt);
  xaxisd = hdenominator2d_ttpt->GetXaxis();
  yaxisd = hdenominator2d_ttpt->GetYaxis();
  for (int j=1;j<=yaxisd->GetNbins();j++) {
    for (int i=1;i<=xaxisd->GetNbins();i++) {
      hdenominator2drebinned_ttpt->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d_ttpt->GetBinContent(i,j));
    }
  }

  hnumerator2drebinned_ttrapidity2 = new TH2D(Form("numerator_%s_ttrapidity2", histname.Data()),Form("numerator_%s_ttrapidity2", histname.Data()),2,binsfor2D,3, binsttrapidity2);
  xaxis = hnumerator2d_ttrapidity2->GetXaxis();
  yaxis = hnumerator2d_ttrapidity2->GetYaxis();
  for (int j=1;j<=yaxis->GetNbins();j++) {
    for (int i=1;i<=xaxis->GetNbins();i++) {
      hnumerator2drebinned_ttrapidity2->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d_ttrapidity2->GetBinContent(i,j));
    }
  }
  hdenominator2drebinned_ttrapidity2 = new TH2D(Form("denominator_%s_ttrapidity2", histname.Data()),Form("denominator_%s_ttrapidity2", histname.Data()),2,binsfor2D,3, binsttrapidity2);
  xaxisd = hdenominator2d_ttrapidity2->GetXaxis();
  yaxisd = hdenominator2d_ttrapidity2->GetYaxis();
  for (int j=1;j<=yaxisd->GetNbins();j++) {
    for (int i=1;i<=xaxisd->GetNbins();i++) {
      hdenominator2drebinned_ttrapidity2->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d_ttrapidity2->GetBinContent(i,j));
    }
  }


  TString accepthistname = "accept_";
  accepthistname += histname;

  TFile *output = new TFile(Form("%s.root", accepthistname.Data()), "RECREATE");  

  hacceptance =  (TH1D*) hnumerator->Clone(accepthistname.Data());
  hacceptance->SetTitle(accepthistname.Data());
  hacceptance->Reset();
  hacceptance->Divide(hnumerator,hdenominator,1., 1.);

  hacceptance2drebinned_mtt =  (TH2D*) hnumerator2drebinned_mtt->Clone( Form("%s_mtt", accepthistname.Data()) );
  hacceptance2drebinned_mtt->Reset();
  hacceptance2drebinned_mtt->SetTitle(Form("%s_mtt", accepthistname.Data()));
  hacceptance2drebinned_mtt->Divide(hnumerator2drebinned_mtt,hdenominator2drebinned_mtt,1., 1.);

  hacceptance2drebinned_ttpt =  (TH2D*) hnumerator2drebinned_ttpt->Clone( Form("%s_ttpt", accepthistname.Data()) );
  hacceptance2drebinned_ttpt->Reset();
  hacceptance2drebinned_ttpt->SetTitle(Form("%s_ttpt", accepthistname.Data()));
  hacceptance2drebinned_ttpt->Divide(hnumerator2drebinned_ttpt,hdenominator2drebinned_ttpt,1., 1.);

  hacceptance2drebinned_ttrapidity2 =  (TH2D*) hnumerator2drebinned_ttrapidity2->Clone( Form("%s_ttrapidity2", accepthistname.Data()) );
  hacceptance2drebinned_ttrapidity2->Reset();
  hacceptance2drebinned_ttrapidity2->SetTitle(Form("%s_ttrapidity2", accepthistname.Data()));
  hacceptance2drebinned_ttrapidity2->Divide(hnumerator2drebinned_ttrapidity2,hdenominator2drebinned_ttrapidity2,1., 1.);

  hnumerator->SetLineColor(kBlue);
  hnumerator-> SetFillColor(0);
  hnumerator->SetMarkerColor(kBlue);
  hdenominator->SetLineColor(kRed);
  hdenominator->SetMarkerColor(kRed);
  hdenominator-> SetFillColor(0);
  hacceptance->SetLineColor(kBlack);
  hacceptance->SetMarkerColor(kBlack);
  hacceptance-> SetFillColor(0);

  gStyle->SetPaintTextFormat("6.4f");

  TCanvas *c1 = new TCanvas("c_acc", "c_acc", 500, 500); 
  c1->cd();

  hacceptance->SetMaximum(1.25*hacceptance->GetMaximum());
  if(hacceptance->GetMinimum() <0.15 *hacceptance->GetMaximum() ) hacceptance->SetMinimum(0.);  
  if(hacceptance->GetMinimum() > 0.) hacceptance->SetMinimum(0.75*hacceptance->GetMinimum() );  

  hacceptance->GetYaxis()->SetTitle("Acceptance #times Efficiency");
  hacceptance->GetYaxis()->SetDecimals(kTRUE);
  hacceptance->SetTitleSize(0.06, "XYZ");
  hacceptance->SetLabelSize(0.05, "XYZ");
  if(!histname.Contains("lepAzimAsym2") ) hacceptance->GetXaxis()->SetNdivisions(504,0);
  else hacceptance->GetXaxis()->SetNdivisions(506);
  hacceptance->GetYaxis()->SetTitleOffset(1.4);


  if(histname.Contains("lepChargeAsym") ) {
    hacceptance->GetXaxis()->SetTitle("#Delta|#eta_{l}|");
  }
  if(histname.Contains("lepAzimAsym2") ) {
    //hacceptance->GetXaxis()->SetTitle("#Delta#phi_{l+l-}");
    hacceptance->GetXaxis()->SetTitle("#Delta#phi_{l#lower[-0.4]{+}l#lower[-0.48]{-}}");
  }
  if(histname.Contains("lepCosTheta") ) {
    //hacceptance->GetXaxis()->SetTitle("cos(#theta_{l})");
    hacceptance->GetXaxis()->SetTitle("cos(^{}#theta_{l}#kern[-0.35]{*})");
  }
  if(histname.Contains("lepPlusCosTheta") ) {
    hacceptance->GetXaxis()->SetTitle("cos(#theta_{l+})");
  }
  if(histname.Contains("lepMinusCosTheta") ) {
    hacceptance->GetXaxis()->SetTitle("cos(#theta_{l-})");
  }
  if(histname.Contains("topSpinCorr") ) {
    //hacceptance->GetXaxis()->SetTitle("cos(#theta_{l+})cos(#theta_{l-})");
    hacceptance->GetXaxis()->SetTitle("cos(^{}#theta_{l#lower[-0.4]{+}}#kern[-1.38]{*}) cos(^{}#theta_{l#lower[-0.48]{-}}#kern[-1.0]{*})");
  }
  if(histname.Contains("lepCosOpeningAngle") ) {
    hacceptance->GetXaxis()->SetTitle("cos(#phi)");
  }
  if(histname.Contains("rapiditydiffMarco") ) {
    hacceptance->GetXaxis()->SetTitle("#Delta|y_{t}|");
  }

  double Asym1,Asym2,Asym3;
  double Asym1err,Asym2err,Asym3err;
  
  hacceptance_copy = new TH1D("hac","hac",6, bins);
  hacceptance_statup = new TH1D("hacup","hacup",6, bins);
  hacceptance_statdown = new TH1D("hacdown","hacdown",6, bins);
  TAxis *xaxis_ = hacceptance->GetXaxis();
  for (int i=1;i<=xaxis_->GetNbins();i++) {
    hacceptance_copy->Fill(xaxis_->GetBinCenter(i), hacceptance->GetBinContent(i));
    hacceptance_statup->Fill(xaxis_->GetBinCenter(i), 2.* hacceptance->GetBinError(i));
    hacceptance_statdown->Fill(xaxis_->GetBinCenter(i), hacceptance->GetBinContent(i) - hacceptance->GetBinError(i));
  }

  

  if(!drawnorm){

    THStack *hs = new THStack("hs_statband", "Stat band");
    hacceptance_statdown->SetLineColor(kWhite);
    hacceptance_statdown->SetFillColor(kWhite);
    hacceptance_statdown->SetFillStyle(0);
    hs->Add(hacceptance_statdown);
    tdrStyle->SetHatchesSpacing(0.6);
    hacceptance_statup->SetFillStyle(3335);
    hacceptance_statup->SetLineColor(kWhite);
    hacceptance_statup->SetFillColor(15);
    hs->Add(hacceptance_statup);

    hs->Draw();

    hs->SetMaximum(1.25*hs->GetMaximum()/1.05);
    if(hs->GetMinimum() <0.15 *hs->GetMaximum() ) hs->SetMinimum(0.);  
    else hs->SetMinimum(0.75*hs->GetMinimum() );  

    hs->GetXaxis()->SetTitle(hacceptance->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetDecimals(kTRUE);
    hs->GetYaxis()->SetTitle("Acceptance #times Efficiency");
    hs->GetYaxis()->SetDecimals(kTRUE);
    hs->GetXaxis()->SetTitleSize(0.06);
    hs->GetXaxis()->SetLabelSize(0.05);
    hs->GetYaxis()->SetTitleSize(0.06);
    hs->GetYaxis()->SetLabelSize(0.05);
    if(!histname.Contains("lepAzimAsym2") ) hs->GetXaxis()->SetNdivisions(504,0);
    else hs->GetXaxis()->SetNdivisions(506);
    hs->GetYaxis()->SetNdivisions(508);
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->GetXaxis()->SetTitleOffset(1.0);

    hacceptance->SetMarkerSize(1.5);
    //hacceptance->Draw("hist TEXT30E");
    //hacceptance->Draw("hist E");

    hacceptance->Draw("hist same");
  }
  else{

    hacceptance_copy->Scale(1./hacceptance_copy->Integral("width"));
    hnumerator->Scale(1./hnumerator->Integral(),"width");
    hdenominator->Scale(1./hdenominator->Integral(),"width");

    double max = 0.;
    if(hacceptance_copy->GetMaximum() > max) max = hacceptance_copy->GetMaximum();
    if(hnumerator->GetMaximum() > max) max = hnumerator->GetMaximum();
    if(hdenominator->GetMaximum() > max) max = hdenominator->GetMaximum();

    double min = 9999.;
    if(hacceptance_copy->GetMinimum() < min) min = hacceptance_copy->GetMinimum();
    if(hnumerator->GetMinimum() < min) min = hnumerator->GetMinimum();
    if(hdenominator->GetMinimum() < min) min = hdenominator->GetMinimum();
    min -= max*0.2;
    if(min < 0.2 * max ) min = 0.;

    max*=1.2;

    TH1 *frame=new TH1F("frame","",1000,bins[0],bins[6]);frame->SetMaximum(max);frame->SetMinimum(min);frame->Draw();
    frame->GetXaxis()->SetTitle( hacceptance->GetXaxis()->GetTitle() );
    frame->GetYaxis()->SetTitle( "Normalized to unit area");

    GetAfb(hacceptance_copy,Asym1, Asym1err);
    GetAfb(hnumerator,Asym2, Asym2err);
    GetAfb(hdenominator,Asym3, Asym3err);

    hacceptance_copy->Draw("histsame");
    hnumerator->Draw("histsame");
    hdenominator->Draw("histsame");
  }

  TLegend *leg = new TLegend(0.48,0.78,0.90,0.92);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.036);
  leg->SetFillStyle(0);
  if(drawnorm) leg->AddEntry(hacceptance, "acceptance","l");
  else { leg->AddEntry(hacceptance, "MC@NLO parton level","l"); leg->AddEntry(hacceptance_statup,    "Statistical uncertainty", "F"); }
  if(drawnorm) leg->AddEntry(hnumerator, "numerator","l");
  if(drawnorm) leg->AddEntry(hdenominator, "denominator","l");
  leg->Draw("same");

  //TPaveText *pt1 = new TPaveText(0.18, 0.86, 0.40, 0.92, "brNDC");
  //if(drawnorm) pt1 = new TPaveText(0.18, 0.77, 0.40, 0.92, "brNDC");
  TPaveText *pt1 = new TPaveText(0.155, 0.94, 0.41, 0.98, "brNDC");
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);

  TText *blah;
  //blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
  //blah = pt1->AddText("CMS, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
  blah = pt1->AddText("CMS Simulation, #sqrt{s} = 7 TeV");
  blah->SetTextSize(0.036);
  blah->SetTextAlign(11);  


  if(drawnorm) {
  blah = pt1->AddText("");

  TString Asym1_temp = formatFloat(Asym1,"%6.4f");
  Asym1_temp.ReplaceAll(" " , "" );
  Asym1_temp = TString("   Asym: ") +  Asym1_temp;
  blah = pt1->AddText(Asym1_temp.Data());
  blah->SetTextSize(0.036);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlack);  

  TString Asym2_temp = formatFloat(Asym2,"%6.4f");
  Asym2_temp.ReplaceAll(" " , "" );
  Asym2_temp = TString("   Asym: ") +  Asym2_temp;
  blah = pt1->AddText(Asym2_temp.Data());
  blah->SetTextSize(0.036);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlue);  

  TString Asym3_temp = formatFloat(Asym3,"%6.4f");
  Asym3_temp.ReplaceAll(" " , "" );
  Asym3_temp = TString("   Asym: ") +  Asym3_temp;
  blah = pt1->AddText(Asym3_temp.Data());
  blah->SetTextSize(0.036);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kRed);
  }
  

  pt1->Draw();

  c1->Print(Form("%s.pdf", accepthistname.Data()));

  hacceptance->Write();
  hdenominator->Write();
  hnumerator->Write();

  hdenominator2d_mtt->Write();
  hnumerator2d_mtt->Write();
  hacceptance2drebinned_mtt->Write();
  hdenominator2drebinned_mtt->Write();
  hnumerator2drebinned_mtt->Write();

  hdenominator2d_ttpt->Write();
  hnumerator2d_ttpt->Write();
  hacceptance2drebinned_ttpt->Write();
  hdenominator2drebinned_ttpt->Write();
  hnumerator2drebinned_ttpt->Write();

  hdenominator2d_ttrapidity2->Write();
  hnumerator2d_ttrapidity2->Write();
  hacceptance2drebinned_ttrapidity2->Write();
  hdenominator2drebinned_ttrapidity2->Write();
  hnumerator2drebinned_ttrapidity2->Write();

  f_1->Close();
  f_2->Close();
  output->Close();

}
