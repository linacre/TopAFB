#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include <fstream>
#include "tdrStyle.C"
#include "CommonFunctions.C"
#include <vector>

using namespace std;

TH1F* hnumerator;
TH1F* hdenominator;
TH1F* hacceptance;
TH2F* hnumerator2d;
TH2F* hdenominator2d;
//TH2F* hacceptance2d;
TH2F* hnumerator2drebinned;
TH2F* hdenominator2drebinned;
TH2F* hacceptance2drebinned;


void acceptanceplots(TString histname = "lepAzimAsym", bool drawnorm = false, TString FName1 = "results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", TString FName2 = "results/hist_noCuts.root"){
  setTDRStyle();

  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  
  hnumerator = (TH1F*)f_1->Get(Form("ttdil_h%sGen_allj_all", histname.Data())); 
  hnumerator2d = (TH2F*)f_1->Get(Form("ttdil_h%sGen2d_allj_all", histname.Data())); 
    
  std::cout << "Opening " << FName2.Data() << "\n";  
  TFile *f_2         = TFile::Open(FName2.Data());  
  hdenominator = (TH1F*)f_2->Get(Form("ttdil_h%sGen_allj_all", histname.Data()));
  hdenominator2d = (TH2F*)f_2->Get(Form("ttdil_h%sGen2d_allj_all", histname.Data())); 
  
  std::cout << "Opened " << Form("ttdil_h%sGen_allj_all", histname.Data()) << " and "<< Form("ttdil_h%sGen2d_allj_all", histname.Data()) <<"\n";
  
  Double_t bins1[] =  { -4., -2., -1., 0., 1., 2., 4.}; 
  Double_t bins2[] = {-1., -0.6, -0.3, 0., 0.3, 0.6, 1.}; 
  Double_t binsMtt[] = {0., 450., 550.,1200.}; 
  Double_t bins1forMtt[] = {-4., 0., 4.}; 
  Double_t bins2forMtt[] = {-1., 0., 1.}; 
  
  if(histname.Contains("lepChargeAsym") ||  histname.Contains("rapiditydiff")) {

  	hnumerator = (TH1F*) hnumerator->Rebin(6,Form("numerator_%s", histname.Data()),bins1);
  	hdenominator = (TH1F*) hdenominator->Rebin(6,Form("denominator_%s", histname.Data()),bins1);

  	hnumerator2drebinned = new TH2F(Form("numerator_%s_mtt", histname.Data()),Form("numerator_%s_mtt", histname.Data()),2,bins1forMtt,3, binsMtt);
  	TAxis *xaxis = hnumerator2d->GetXaxis();
  	TAxis *yaxis = hnumerator2d->GetYaxis();
  	for (int j=1;j<=yaxis->GetNbins();j++) {
  		for (int i=1;i<=xaxis->GetNbins();i++) {
  			hnumerator2drebinned->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d->GetBinContent(i,j));
  		}
  	}
  	
  	hdenominator2drebinned = new TH2F(Form("denominator_%s_mtt", histname.Data()),Form("denominator_%s_mtt", histname.Data()),2,bins1forMtt,3, binsMtt);
  	TAxis *xaxisd = hdenominator2d->GetXaxis();
  	TAxis *yaxisd = hdenominator2d->GetYaxis();
  	for (int j=1;j<=yaxisd->GetNbins();j++) {
  		for (int i=1;i<=xaxisd->GetNbins();i++) {
  			hdenominator2drebinned->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d->GetBinContent(i,j));
  		}
  	}
  	
  }
  
  else {
  	
  	hnumerator = (TH1F*) hnumerator->Rebin(6,Form("numerator_%s", histname.Data()),bins2);
  	hdenominator = (TH1F*) hdenominator->Rebin(6,Form("denominator_%s", histname.Data()),bins2);
  	
    	hnumerator2drebinned = new TH2F(Form("numerator_%s_mtt", histname.Data()),Form("numerator_%s_mtt", histname.Data()),2,bins2forMtt,3, binsMtt);
  	TAxis *xaxis = hnumerator2d->GetXaxis();
  	TAxis *yaxis = hnumerator2d->GetYaxis();
  	for (int j=1;j<=yaxis->GetNbins();j++) {
  		for (int i=1;i<=xaxis->GetNbins();i++) {
  			hnumerator2drebinned->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d->GetBinContent(i,j));
  		}
  	}
  	
  	hdenominator2drebinned = new TH2F(Form("denominator_%s_mtt", histname.Data()),Form("denominator_%s_mtt", histname.Data()),2,bins2forMtt,3, binsMtt);
  	TAxis *xaxisd = hdenominator2d->GetXaxis();
  	TAxis *yaxisd = hdenominator2d->GetYaxis();
  	for (int j=1;j<=yaxisd->GetNbins();j++) {
  		for (int i=1;i<=xaxisd->GetNbins();i++) {
  			hdenominator2drebinned->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d->GetBinContent(i,j));
  		}
  	}
    	
  }
  
  TString accepthistname = "accept_";
  accepthistname += histname;
  
  TFile *output = new TFile(Form("%s.root", accepthistname.Data()), "RECREATE");  
  
  hacceptance =  (TH1F*) hnumerator->Clone(accepthistname.Data());
  hacceptance->SetTitle(accepthistname.Data());
  hacceptance->Reset();
  hacceptance->Divide(hnumerator,hdenominator,1., 1.);
  
  //hacceptance2d =  (TH2F*) hnumerator2d->Clone( Form("%s_mtt_notrebinned", accepthistname.Data()) );
  //hacceptance2d->Reset();
  //hacceptance2d->SetTitle(Form("%s_mtt_notrebinned", accepthistname.Data()));
  //hacceptance2d->Divide(hnumerator2d,hdenominator2d,1., 1.);
  
  hacceptance2drebinned =  (TH2F*) hnumerator2drebinned->Clone( Form("%s_mtt", accepthistname.Data()) );
  hacceptance2drebinned->Reset();
  hacceptance2drebinned->SetTitle(Form("%s_mtt", accepthistname.Data()));
  hacceptance2drebinned->Divide(hnumerator2drebinned,hdenominator2drebinned,1., 1.);
  
  
    
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

  

  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  hacceptance->SetMaximum(1.25*hacceptance->GetMaximum());
  if(hacceptance->GetMinimum() <0.15 *hacceptance->GetMaximum() ) hacceptance->SetMinimum(0.);  
  if(hacceptance->GetMinimum() > 0.) hacceptance->SetMinimum(0.75*hacceptance->GetMinimum() );  

  if(!drawnorm){
  	hacceptance->Draw("hist TEXT00E");
  }
  else{
  	hacceptance->SetMaximum(0.45);
  	hacceptance->SetMinimum(0.); 
  	hacceptance->DrawNormalized("hist");
  	hnumerator->DrawNormalized("histsame");
  	hdenominator->DrawNormalized("histsame");
  }

  TLegend *leg = new TLegend(0.74,0.86,0.90,0.92);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.032);
  leg->SetFillStyle(0);
  leg->AddEntry(hacceptance, "acceptance","l");
  if(drawnorm) leg->AddEntry(hnumerator, "numerator","l");
  if(drawnorm) leg->AddEntry(hdenominator, "denominator","l");

  leg->Draw("same");


  TPaveText *pt1 = new TPaveText(0.18, 0.86, 0.40, 0.92, "brNDC");
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  
  TText *blah;
  blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
/*
  TString KS3_temp = formatFloat(KS3,"%6.2f");
  KS3_temp.ReplaceAll(" " , "" );
  KS3_temp = TString("   K-S: ") +  KS3_temp;
  blah = pt1->AddText(KS3_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kGreen+3);  

  TString KS1_temp = formatFloat(KS1,"%6.2f");
  KS1_temp.ReplaceAll(" " , "" );
  KS1_temp = TString("   K-S: ") +  KS1_temp;
  blah = pt1->AddText(KS1_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlue);  

  TString KS2_temp = formatFloat(KS2,"%6.2f");
  KS2_temp.ReplaceAll(" " , "" );
  KS2_temp = TString("   K-S: ") +  KS2_temp;
  blah = pt1->AddText(KS2_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kRed);  
*/  


  
  pt1->Draw();


  c1->Print(Form("%s.pdf", accepthistname.Data()));
  hacceptance->Write();
  hdenominator->Write();
  hnumerator->Write();
  
  //hacceptance2d->Write();
  hdenominator2d->Write();
  hnumerator2d->Write();

  hacceptance2drebinned->Write();
  hdenominator2drebinned->Write();
  hnumerator2drebinned->Write();

  f_1->Close();
  f_2->Close();
  output->Close();
  
}
