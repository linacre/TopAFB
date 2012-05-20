#include "TH1F.h"
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

TH1F* histo1;
TH1F* histo2;
TH1F* histo3;
TH1F* datahisto;



void acceptanceplots(TString histname = "lepAzimAsym", bool drawnorm = false, TString FName1 = "results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", TString FName2 = "results/hist_noCuts.root"){
  setTDRStyle();

  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  
  histo1 = (TH1F*)f_1->Get(Form("ttdil_h%sGen_allj_all", histname.Data())); 
    
  std::cout << "Opening " << FName2.Data() << "\n";  
  TFile *f_2         = TFile::Open(FName2.Data());  
  histo2 = (TH1F*)f_2->Get(Form("ttdil_h%sGen_allj_all", histname.Data()));    
  std::cout << "Opened " << Form("ttdil_h%sGen_allj_all", histname.Data()) << "\n";
  
  Double_t bins1[] =  { -4., -2., -1., 0., 1., 2., 4.}; 
  Double_t bins2[] = {-1., -0.6, -0.3, 0., 0.3, 0.6, 1.}; 

  if(histname.Contains("lepChargeAsym")) {
  	histo1 = (TH1F*) histo1->Rebin(6,Form("numerator_%s", histname.Data()),bins1);
  	histo2 = (TH1F*) histo2->Rebin(6,Form("denominator_%s", histname.Data()),bins1);
  }

  else {
  	histo1 = (TH1F*) histo1->Rebin(6,Form("numerator_%s", histname.Data()),bins2);
  	histo2 = (TH1F*) histo2->Rebin(6,Form("denominator_%s", histname.Data()),bins2);
  }
  
  TString accepthistname = "accept_";
  accepthistname += histname;
  
  TFile *output = new TFile(Form("%s.root", accepthistname.Data()), "RECREATE");  
  
  histo3 =  (TH1F*) histo1->Clone(accepthistname.Data());
  histo3->Reset();
  histo3->Divide(histo1,histo2,1., 1.);
    
  histo1->SetLineColor(kBlue);
  histo1-> SetFillColor(0);
  histo1->SetMarkerColor(kBlue);
  histo2->SetLineColor(kRed);
  histo2->SetMarkerColor(kRed);
  histo2-> SetFillColor(0);
  histo3->SetLineColor(kBlack);
  histo3->SetMarkerColor(kBlack);
  histo3-> SetFillColor(0);
  
  gStyle->SetPaintTextFormat("6.4f");

  

  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  histo3->SetMaximum(1.25*histo3->GetMaximum());
  if(histo3->GetMinimum() <0.15 *histo3->GetMaximum() ) histo3->SetMinimum(0.);  
  if(histo3->GetMinimum() > 0.) histo3->SetMinimum(0.75*histo3->GetMinimum() );  

  if(!drawnorm){
  	histo3->Draw("hist TEXT00E");
  }
  else{
  	histo3->SetMaximum(0.45);
  	histo3->SetMinimum(0.); 
  	histo3->DrawNormalized("hist");
  	histo1->DrawNormalized("histsame");
  	histo2->DrawNormalized("histsame");
  }

  TLegend *leg = new TLegend(0.74,0.86,0.90,0.92);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.032);
  leg->SetFillStyle(0);
  leg->AddEntry(histo3, "acceptance","l");
  if(drawnorm) leg->AddEntry(histo1, "numerator","l");
  if(drawnorm) leg->AddEntry(histo2, "denominator","l");

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
  histo3->Write();
  histo2->Write();
  histo1->Write();

  f_1->Close();
  f_2->Close();
  output->Close();
  
}
