
#include <iostream>
#include "AfbFinalUnfold.h"

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


#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"

#include "tdrstyle.C"

using std::cout;
using std::endl;


//==============================================================================
// Global definitions
//==============================================================================

const Double_t _topScalingFactor=8977.0/8191.0;


void AfbUnfoldExample()
{
  setTDRStyle();
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cout.precision(3);

 TRandom3* random = new TRandom3();                                                                                                        
 random->SetSeed(5);

 ofstream myfile;
 myfile.open ("summary_2Dunfolding.txt");
 cout.rdbuf(myfile.rdbuf());

 TString Region="";
 double xsection = 154.0;

 int kterm=4;
 int nVars =6;

 Float_t observable, observable_gen, weight, ttmass, ttmass_gen, ttRapidity, tmass;

  for (Int_t iVar= 0; iVar < nVars; iVar++) {

    Initialize2DBinning(iVar);

  TH1D* hData= new TH1D ("Data_BkgSub", "Data with background subtracted",    nbins2D, xbins2D);
  TH1D* hBkg = new TH1D ("Background",  "Background",    nbins2D, xbins2D);
  TH1D* hTop = new TH1D ("Top",  "Top",    nbins2D, xbins2D);
  TH1D* hTop_gen= new TH1D ("Top_Gen", "Top Gen", nbins2D, xbins2D); 
  TH1D* hData_unfolded= new TH1D ("Data_Unfold", "Data with background subtracted and unfolded", nbins2D, xbins2D);
 

  TH1D* hTrue= new TH1D ("true", "Truth",    nbins2D, xbins2D);
  TH1D* hMeas= new TH1D ("meas", "Measured", nbins2D, xbins2D);
  TH1D* hData_bkgSub;

  TMatrixD m_unfoldE(nbins2D,nbins2D);
  TMatrixD m_correctE(nbins2D,nbins2D);

  RooUnfoldResponse response (hMeas, hTrue);
  

  //  Now test with data and with BKG subtraction

  TChain *ch_bkg = new TChain("tree");
  TChain *ch_top = new TChain("tree");

  TChain *ch_data = new TChain("tree");
  
  TString path="../";
  
  ch_data->Add(path+"data.root");

  ch_top->Add(path+"ttdil_powheg.root");

  ch_bkg->Add(path+"ttotr.root");
  ch_bkg->Add(path+"wjets.root");
  ch_bkg->Add(path+"DYee.root");
  ch_bkg->Add(path+"DYmm.root");
  ch_bkg->Add(path+"DYtautau.root");
  ch_bkg->Add(path+"tw.root");
  ch_bkg->Add(path+"VV.root");

  ch_data->SetBranchAddress(observablename,    &observable);
  ch_data->SetBranchAddress("weight",&weight);
  ch_data->SetBranchAddress("tt_mass",&ttmass);
  ch_data->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_data->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_data->GetEntries(); i++) {
    ch_data->GetEntry(i);
    if (ttmass>0.0) hData->Fill(sign(observable)*ttmass,weight);
  }

  ch_bkg->SetBranchAddress(observablename,    &observable);
  ch_bkg->SetBranchAddress("weight",&weight);
  ch_bkg->SetBranchAddress("tt_mass",&ttmass);
  ch_bkg->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_bkg->SetBranchAddress("t_mass",&tmass);


  for (Int_t i= 0; i<ch_bkg->GetEntries(); i++) {
    ch_bkg->GetEntry(i);
    if (ttmass>0.0) hBkg->Fill(sign(observable)*ttmass,weight);
    }


  ch_top->SetBranchAddress(observablename,    &observable);
  ch_top->SetBranchAddress(observablename+"_gen",&observable_gen);
  ch_top->SetBranchAddress("weight",&weight);
  ch_top->SetBranchAddress("tt_mass",&ttmass);
  ch_top->SetBranchAddress("tt_mass_gen",&ttmass_gen);
  ch_top->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_top->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_top->GetEntries(); i++) {
    ch_top->GetEntry(i);
    if (ttmass>0.0) {
      response.Fill (sign(observable)*ttmass, sign(observable_gen)*ttmass_gen, weight);
      hTop->Fill(sign(observable)*ttmass,weight);
      hTop_gen->Fill(sign(observable_gen)*ttmass_gen,weight);
    }
  }
  
  hTop->Scale(_topScalingFactor);

  TCanvas* c_test = new TCanvas("c_final","c_final",500,500); 

  hData->SetLineWidth(lineWidth+2);
  hTop->SetLineWidth(lineWidth);

  hTop->SetLineColor(TColor::GetColorDark(kGreen));
  hTop->SetFillColor(TColor::GetColorDark(kGreen));

  hTop->SetFillStyle(3353);
  hTop_gen->SetFillStyle(3353);
  hBkg->SetLineColor(kYellow);
  hBkg->SetFillColor(kYellow);
  

  THStack *hs = new THStack("hs","Stacked Top+BG");

  hs->Add(hBkg);
  hs->Add(hTop);

  hs->SetMinimum(0.0);
  hs->SetMaximum( 2.0* hs->GetMaximum());
  hs->Draw("hist");
  hs->GetXaxis()->SetTitle("M_{tt}");
  hs->GetYaxis()->SetTitleOffset(1.3);
  hs->GetYaxis()->SetTitle("Events");

  hData->Draw("E same");

  TLegend* leg1=new TLegend(0.6,0.62,0.9,0.838,NULL,"brNDC");
  leg1->SetEntrySeparation(100);  
  leg1->SetFillColor(0);                                                                                                                  
  leg1->SetLineColor(0);                                                                                                                   
  leg1->SetBorderSize(0);                                                                                 
  leg1->SetTextSize(0.03);
  leg1->AddEntry(hData, "Data");
  leg1->AddEntry(hTop,  "t-tbar", "F");                                                               
  leg1->AddEntry(hBkg,  "Background", "F");                                                               
  leg1->Draw();                
  c_test->SaveAs("Mtt_"+observablename+Region+".pdf");


  hData_bkgSub= (TH1D*) hData->Clone();
  hData_bkgSub->Add(hBkg,-1.0);
 
  RooUnfoldSvd unfold (&response, hData_bkgSub, kterm); 
  unfold.Setup(&response,hData_bkgSub);                                                                                                         
  hData_unfolded = (TH1D*) unfold.Hreco();  
  m_unfoldE = unfold.Ereco(); 

  TFile *file = new TFile("../acceptance/powheg/accept_"+acceptanceName+".root");
  TH2D *acceptM_2d = (TH2D*) file->Get("accept_"+acceptanceName+"_mtt");

  TH1D* acceptM = new TH1D ("accept", "accept",    nbins2D, xbins2D);
  acceptM->SetBinContent(1,acceptM_2d->GetBinContent(1,3));
  acceptM->SetBinContent(2,acceptM_2d->GetBinContent(1,2));
  acceptM->SetBinContent(3,acceptM_2d->GetBinContent(1,1));

  acceptM->SetBinContent(4,acceptM_2d->GetBinContent(2,1));
  acceptM->SetBinContent(5,acceptM_2d->GetBinContent(2,2));
  acceptM->SetBinContent(6,acceptM_2d->GetBinContent(2,3));

  acceptM->Scale(1.0/acceptM->Integral());

  for (Int_t i= 1; i<=nbins2D; i++) {

    if (acceptM->GetBinContent(i)!=0) {
      hData_unfolded->SetBinContent(i, hData_unfolded->GetBinContent(i)*1.0/acceptM->GetBinContent(i));
      hData_unfolded->SetBinError  (i, hData_unfolded->GetBinError  (i)*1.0/acceptM->GetBinContent(i));
    }

    if (acceptM->GetBinContent(i)!=0) {
      hTop_gen->SetBinContent(i, hTop->GetBinContent(i)*1.0/acceptM->GetBinContent(i));
      hTop_gen->SetBinError  (i, hTop->GetBinError(i)  *1.0/acceptM->GetBinContent(i));
    }
  }
  

  double dataIntegral = hData_unfolded->Integral();
  hData_unfolded->Scale(xsection/hData_unfolded->Integral());
  hTop_gen->Scale(xsection/hTop_gen->Integral());
  
  for(int l=0;l<nbins2D;l++){
      for(int j=0;j<nbins2D;j++){
	double corr = 1.0 / ( acceptM->GetBinContent(l+1) * acceptM->GetBinContent(j+1) );
	corr = corr * pow(xsection / dataIntegral,2) ;
	m_correctE(l,j) = m_unfoldE(l,j)*corr;
      }
    }
  

  //==================================================================
  // ============== Print the assymetry =============================
  cout<<"========= Variable: "<<observablename <<"===================\n";
  Float_t Afb, AfbErr;

  GetAfb(hData, Afb, AfbErr);
  cout<<" Data: "<< Afb <<" +/-  "<< AfbErr<<"\n";

  GetAfb(hTop_gen, Afb, AfbErr);
  cout<<" True Top: "<< Afb <<" +/-  "<< AfbErr<<"\n";

  GetCorrectedAfb(hData_unfolded, m_correctE, Afb, AfbErr);
  cout<<" Unfolded: "<< Afb <<" +/-  "<< AfbErr<<"\n";

  //  TCanvas* c_test = new TCanvas("c_final","c_final",500,500); 
  hData_unfolded->GetXaxis()->SetTitle("M_{t#bar t}");
  hData_unfolded->GetYaxis()->SetTitle("d#sigma/dM_{t#bar t}");
  hData_unfolded->SetMinimum(0.0);
  hData_unfolded->SetMaximum( 2.0* hData_unfolded->GetMaximum());
  hData_unfolded->SetMarkerStyle(23);
  hData_unfolded->SetMarkerSize(2.0);
  hData_unfolded->Draw("E");
  hData_unfolded->SetLineWidth(lineWidth);
  hTop_gen->SetLineWidth(lineWidth);
  hTop_gen->SetLineColor(TColor::GetColorDark(kGreen));
  hTop_gen->SetFillColor(TColor::GetColorDark(kGreen));
  hTop_gen->Draw("hist same");
  hData_unfolded->Draw("E same");

  leg1=new TLegend(0.6,0.62,0.9,0.838,NULL,"brNDC");   
  leg1->SetEntrySeparation(100);                                                                                                          
  leg1->SetFillColor(0);                                                                                                                  
  leg1->SetLineColor(0);                                                                                                                   
  leg1->SetBorderSize(0);                                                                                                                  
  leg1->SetTextSize(0.03);                                                                              
  leg1->AddEntry(hData_unfolded, "( Data-BG ) Unfolded");                                                                                       
  leg1->AddEntry(hTop_gen,    "Powheg parton", "F");                                                               
  leg1->Draw();                

  c_test->SaveAs("finalplot_2D_unfolded_"+observablename+Region+".pdf");

  vector<double> afb_m;
  vector<double> afb_merr;
  GetAvsY(hData_unfolded, m_correctE, afb_m, afb_merr);  

  xbins2D[0]=200.0; xbins2D[1]=450; xbins2D[2]=550.0; xbins2D[3]=1500.0;
  TH1D* hAfbVsMtt = new TH1D ("AfbVsMtt",  "AfbVsMtt",  3, xbins2D);
  for (int nb=0; nb<3; nb++) {
      hAfbVsMtt->SetBinContent(nb+1,afb_m[nb]);
      hAfbVsMtt->SetBinError(nb+1,afb_merr[nb]);
    }
  
  //  GetAvsY(hTop_gen, m_unfoldE, afb_m, afb_merr);  
 
  //  TH1D* hTop_AfbVsMtt = new TH1D ("Top_AfbVsMtt",  "Top_AfbVsMtt",  3, xbins2D);
  //  for (int nb=0; nb<3; nb++) {
  //      hTop_AfbVsMtt->SetBinContent(nb+1,afb_m[nb]);
  //      hTop_AfbVsMtt->SetBinError(nb+1,afb_merr[nb]);
  //    }

  tdrStyle->SetErrorX(0.5);
  hAfbVsMtt->SetMinimum(-0.3);
  hAfbVsMtt->SetMaximum( 0.3);
  hAfbVsMtt->SetLineWidth( 2.0 );
  hAfbVsMtt->Draw("E");
  //  hTop_AfbVsMtt->SetLineColor(kGreen);
  hAfbVsMtt->GetYaxis()->SetTitle("A("+xaxislabel+")");
  hAfbVsMtt->GetYaxis()->SetTitleOffset(1.2);
  hAfbVsMtt->GetXaxis()->SetTitle("M_{t #bar t} GeV");
  //  hTop_AfbVsMtt->Draw("E same");
  c_test->SaveAs("AfbVsMtt_unfolded_"+observablename+Region+".pdf");

  }
  myfile.close();

}

#ifndef __CINT__
int main () { AfbUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
