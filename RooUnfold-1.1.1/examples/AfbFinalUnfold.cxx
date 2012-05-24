#include <iostream>
#include <fstream>
#include "AfbFinalUnfold.h"

#include "TROOT.h"
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


void AfbUnfoldExample()
{

  setTDRStyle();
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cout.precision(3);

  ofstream myfile;
  myfile.open ("summary_1Dunfolding.txt");
  cout.rdbuf(myfile.rdbuf());

  TRandom3* random = new TRandom3();                                                                                                        
  random->SetSeed(5);

  TString Region="";
 
  int kterm=3; 
  int nVars =8;

  Float_t observable, observable_gen, weight, ttmass, ttRapidity, tmass;

  for (Int_t iVar= 0; iVar < nVars; iVar++) {

    Initialize1DBinning(iVar);

  TH1D* hData= new TH1D ("Data_BkgSub", "Data with background subtracted",    nbins1D, xbins1D);
  TH1D* hBkg = new TH1D ("Background",  "Background",    nbins1D, xbins1D);
  TH1D* hTop = new TH1D ("Top",  "Top",    nbins1D, xbins1D);
  TH1D* hTop_gen = new TH1D ("Top Gen",  "Top Gen",    nbins1D, xbins1D);
  TH1D* hData_unfolded= new TH1D ("Data_Unfold", "Data with background subtracted and unfolded", nbins1D, xbins1D);

  TH1D* hTrue= new TH1D ("true", "Truth",    nbins1D, xbins1D);
  TH1D* hMeas= new TH1D ("meas", "Measured", nbins1D, xbins1D);
  TH1D* hData_bkgSub;

  hData->Sumw2();
  hTop->Sumw2();

  TMatrixD m_unfoldE (nbins1D,nbins1D);
  TMatrixD m_correctE(nbins1D,nbins1D);

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
    if ( (Region=="Signal") && (ttmass>450) )  
      hData->Fill(observable,weight);
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) 
      hData->Fill(observable,    weight);    
    if ( (Region=="") && (iVar<2) ) 
      hData->Fill(observable,    weight);   
  }

  ch_bkg->SetBranchAddress(observablename,    &observable);
  ch_bkg->SetBranchAddress("weight",&weight);
  ch_bkg->SetBranchAddress("tt_mass",&ttmass);
  ch_bkg->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_bkg->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_bkg->GetEntries(); i++) {
    ch_bkg->GetEntry(i);
    if ( (Region=="Signal") && (ttmass>450) )  
      hBkg->Fill(observable    ,weight);
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) 
      hBkg->Fill(observable,    weight);    
    if ( (Region=="") && (iVar<2) ) 
      hBkg->Fill(observable,    weight);   
  }

  ch_top->SetBranchAddress(observablename,    &observable);
  ch_top->SetBranchAddress(observablename+"_gen",&observable_gen);
  ch_top->SetBranchAddress("weight",&weight);
  ch_top->SetBranchAddress("tt_mass",&ttmass);
  ch_top->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_top->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_top->GetEntries(); i++) {
    ch_top->GetEntry(i);
    if ( (Region=="Signal") && (ttmass>450) ) {
      response.Fill (observable, observable_gen, weight);
      hTop->Fill(observable, weight);    
      hTop_gen->Fill(observable_gen, weight);    
    }
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) {
      response.Fill (observable, observable_gen, weight);
      hTop->Fill(observable, weight);    
      hTop_gen->Fill(observable_gen, weight);    
    }
    if ( (Region=="") && (iVar<2) ) {
      response.Fill (observable, observable_gen, weight);
      hTop->Fill(observable, weight);    
      hTop_gen->Fill(observable_gen, weight);    
    }
  }

  TChain *ch_top_mcnlo = new TChain("tree");
  ch_top_mcnlo->Add("../ttdil_mcnlo.root");

  ch_top_mcnlo->SetBranchAddress(observablename,    &observable);
  ch_top_mcnlo->SetBranchAddress(observablename+"_gen",&observable_gen);
  ch_top_mcnlo->SetBranchAddress("weight",&weight);
  ch_top_mcnlo->SetBranchAddress("tt_mass",&ttmass);
  ch_top_mcnlo->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_top_mcnlo->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_top_mcnlo->GetEntries(); i++) {
    ch_top_mcnlo->GetEntry(i);
    if ( (Region=="Signal") && (ttmass>450) ) {
      //      response.Fill (observable, observable_gen);
    }
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) {
      //      response.Fill (observable, observable_gen);
    }
    if ( (Region=="") && (iVar<2) ) {
      //      response.Fill (observable, observable_gen);
    }
  }


  hData_bkgSub= (TH1D*) hData->Clone();
  hData_bkgSub->Add(hBkg,-1.0);
 
  RooUnfoldSvd unfold (&response, hData_bkgSub, kterm); 
  unfold.Setup(&response,hData_bkgSub);                                                                                                         
  hData_unfolded = (TH1D*) unfold.Hreco();  
  m_unfoldE = unfold.Ereco(); 

  TFile *file = new TFile("../acceptance/powheg/accept_"+acceptanceName+".root");
  TH1D *acceptM = (TH1D*) file->Get("accept_"+acceptanceName);
  acceptM->Scale(1.0/acceptM->Integral());

  for (Int_t i= 1; i<=nbins1D; i++) {

    if (acceptM->GetBinContent(i)!=0) {
            hData_unfolded->SetBinContent(i, hData_unfolded->GetBinContent(i)*1.0/acceptM->GetBinContent(i));
            hData_unfolded->SetBinError  (i, hData_unfolded->GetBinError  (i)*1.0/acceptM->GetBinContent(i));
    }

    if (acceptM->GetBinContent(i)!=0) {
      hTop_gen->SetBinContent(i, hTop_gen->GetBinContent(i)*1.0/acceptM->GetBinContent(i));
      hTop_gen->SetBinError  (i, hTop_gen->GetBinError(i)  *1.0/acceptM->GetBinContent(i));
    }  
  } 

  double dataIntegral = hData_unfolded->Integral();
  hData_unfolded->Scale(xsection/hData_unfolded->Integral());
  hTop_gen->Scale(xsection/hTop_gen->Integral());
  
  for(int l=0;l<nbins1D;l++){
      for(int j=0;j<nbins1D;j++){
	double corr = 1.0 / ( acceptM->GetBinContent(l+1) * acceptM->GetBinContent(j+1) );
	corr = corr * pow(xsection / dataIntegral,2) ;
	m_correctE(l,j) = m_unfoldE(l,j) * corr;
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

  TCanvas* c_test = new TCanvas("c_final","c_final",500,500); 
  hData_unfolded->GetXaxis()->SetTitle(xaxislabel);
  hData_unfolded->GetYaxis()->SetTitle("d#sigma/d("+xaxislabel+")");
  hData_unfolded->SetMinimum(0.0);
  hData_unfolded->SetMaximum( 2.0* hData_unfolded->GetMaximum());
  hData_unfolded->SetMarkerStyle(23);
  hData_unfolded->SetMarkerSize(1.5);
  hData_unfolded->Draw("E");
  hData_unfolded->SetLineWidth(lineWidth);
  hTop_gen->SetLineWidth(lineWidth);
  hTop_gen->SetLineColor(TColor::GetColorDark(kGreen));
  hTop_gen->SetFillColor(TColor::GetColorDark(kGreen));
  hTop_gen->SetFillStyle(3353);
  hTop_gen->Draw("hist same");
  hData_unfolded->Draw("E same");

  TLegend* leg1=new TLegend(0.6,0.62,0.9,0.838,NULL,"brNDC");                                                                           
  leg1->SetEntrySeparation(100);                                                                                                          
  leg1->SetFillColor(0);                                                                                                                  
  leg1->SetLineColor(0);                                                                                                                   
  leg1->SetBorderSize(0);                                                                                                                  
  leg1->SetTextSize(0.03);                                                                                                                 
  leg1->AddEntry(hData_unfolded, "( Data - BG ) Unfolded");  
  leg1->AddEntry(hTop_gen,    "Powheg Parton", "F");                                                               
  leg1->Draw();                
  c_test->SaveAs("finalplot_unfolded_"+observablename+Region+".pdf");
  }

  myfile.close();
}

#ifndef __CINT__
int main () { AfbUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
