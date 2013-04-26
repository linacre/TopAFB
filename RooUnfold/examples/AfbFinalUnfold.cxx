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
#include "TPaveText.h"
#include "TLatex.h"

#include "TUnfold.h"

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"
#include "src/RooUnfoldTUnfold.h"

#include "tdrstyle.C"

using std::cout;
using std::endl;


//==============================================================================
// Global definitions
//==============================================================================

 // 0=SVD, 1=TUnfold via RooUnfold, 2=TUnfold
  int unfoldingType=0;

  TString Region="";
  Int_t kterm=3; 
  Double_t tau=1E-4;
  Int_t nVars =10;
  Int_t includeSys = 1;


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

  Float_t observable, observable_gen, ttmass, ttRapidity, tmass;
  Float_t observableMinus, observableMinus_gen; 
  Double_t weight;
  Int_t Nsolns;

  for (Int_t iVar= 0; iVar < nVars; iVar++) {

    Initialize1DBinning(iVar);
  bool combineLepMinus = iVar == 9 ? true : false;

  TH1D* hData= new TH1D ("Data_BkgSub", "Data with background subtracted",    nbins1D, xbins1D);
  TH1D* hBkg = new TH1D ("Background",  "Background",    nbins1D, xbins1D);
  TH1D* hTop = new TH1D ("Top",  "Top",    nbins1D, xbins1D);
  TH1D* hTop_gen = new TH1D ("Top Gen",  "Top Gen",    nbins1D, xbins1D);
  TH1D* hData_unfolded= new TH1D ("Data_Unfold", "Data with background subtracted and unfolded", nbins1D, xbins1D);
  
  double xbins1D_arccos[nbins1D+1];
  //xbins1D_arccos[0]=acos(-1.0); xbins1D_arccos[1]=acos(-0.6); xbins1D_arccos[2]=acos(-0.3); xbins1D_arccos[3]=acos(0.0); xbins1D_arccos[4]=acos(0.3); xbins1D_arccos[5]=acos(0.6); xbins1D_arccos[6]=acos(1.0);
  xbins1D_arccos[0]=acos(1.0); xbins1D_arccos[1]=acos(0.6); xbins1D_arccos[2]=acos(0.3); xbins1D_arccos[3]=acos(0.0); xbins1D_arccos[4]=acos(-0.3); xbins1D_arccos[5]=acos(-0.6); xbins1D_arccos[6]=acos(-1.0);
  TH1D* hData_unfolded_arccos= new TH1D ("Data_Unfold_arccos", "Data with background subtracted and unfolded arccos", nbins1D, xbins1D_arccos);
  TH1D* hTop_gen_arccos = new TH1D ("Top Gen arccos",  "Top Gen arccos",    nbins1D, xbins1D_arccos);

  TH1D* hTrue= new TH1D ("true", "Truth",    nbins1D, xbins1D);
  TH1D* hMeas= new TH1D ("meas", "Measured", nbins1D, xbins1D);

  TH2D* hTrue_vs_Meas= new TH2D ("true_vs_meas", "True vs Measured", nbins1D, xbins1D, nbins1D, xbins1D);

  TH1D* hData_bkgSub;

  hData->Sumw2();
  hBkg->Sumw2();
  hTop->Sumw2();
  hTop_gen->Sumw2();
  hData_unfolded->Sumw2();
  hTrue->Sumw2();
  hMeas->Sumw2();
  hTrue_vs_Meas->Sumw2();
  hData_unfolded_arccos->Sumw2();
  hTop_gen_arccos->Sumw2();
  

  TMatrixD m_unfoldE (nbins1D,nbins1D);
  TMatrixD m_correctE(nbins1D,nbins1D);


  //  Now test with data and with BKG subtraction

  TChain *ch_bkg = new TChain("tree");
  TChain *ch_top = new TChain("tree");

  TChain *ch_data = new TChain("tree");

  TString path="../";

  ch_data->Add(path+"data.root");

  ch_top->Add(path+"ttdil.root");

  ch_bkg->Add(path+"ttotr.root");
  ch_bkg->Add(path+"wjets.root");
  ch_bkg->Add(path+"DYee.root");
  ch_bkg->Add(path+"DYmm.root");
  ch_bkg->Add(path+"DYtautau.root");
  ch_bkg->Add(path+"tw.root");
  ch_bkg->Add(path+"VV.root");

  ch_data->SetBranchAddress(observablename,    &observable);
  if( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
  ch_data->SetBranchAddress("weight",&weight);
  ch_data->SetBranchAddress("Nsolns",&Nsolns);
  ch_data->SetBranchAddress("tt_mass",&ttmass);
  ch_data->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_data->SetBranchAddress("t_mass",&tmass);


  for (Int_t i= 0; i<ch_data->GetEntries(); i++) {
    ch_data->GetEntry(i);
    if ( (Region=="Signal") && (ttmass>450) )  
      fillUnderOverFlow(hData, observable, weight, Nsolns);
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) 
      fillUnderOverFlow(hData, observable, weight, Nsolns);    
    if ( (Region=="") && (iVar<2) ) 
      fillUnderOverFlow(hData, observable, weight, Nsolns);  

    if (combineLepMinus) {
	    if ( (Region=="Signal") && (ttmass>450) )  
	      fillUnderOverFlow(hData, observableMinus, weight, Nsolns);
	    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) 
	      fillUnderOverFlow(hData, observableMinus, weight, Nsolns);    
	    if ( (Region=="") && (iVar<2) ) 
	      fillUnderOverFlow(hData, observableMinus, weight, Nsolns);    
	} 
  }

  ch_bkg->SetBranchAddress(observablename,    &observable);
  if( combineLepMinus ) ch_bkg->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
  ch_bkg->SetBranchAddress("weight",&weight);
  ch_bkg->SetBranchAddress("Nsolns",&Nsolns);
  ch_bkg->SetBranchAddress("tt_mass",&ttmass);
  ch_bkg->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_bkg->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_bkg->GetEntries(); i++) {
    ch_bkg->GetEntry(i);
    if ( (Region=="Signal") && (ttmass>450) )  
      fillUnderOverFlow(hBkg, observable, weight, Nsolns);
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) 
      fillUnderOverFlow(hBkg, observable, weight, Nsolns);
    if ( (Region=="") && (iVar<2) ) 
      fillUnderOverFlow(hBkg, observable, weight, Nsolns);

    if (combineLepMinus) {
	    if ( (Region=="Signal") && (ttmass>450) )  
	      fillUnderOverFlow(hBkg, observableMinus, weight, Nsolns);
	    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) 
	      fillUnderOverFlow(hBkg, observableMinus, weight, Nsolns);
	    if ( (Region=="") && (iVar<2) ) 
	      fillUnderOverFlow(hBkg, observableMinus, weight, Nsolns);
	}
  }

  ch_top->SetBranchAddress(observablename,    &observable);
  ch_top->SetBranchAddress(observablename+"_gen",&observable_gen);
  if( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
  if( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms_gen",    &observableMinus_gen);
  ch_top->SetBranchAddress("weight",&weight);
  ch_top->SetBranchAddress("Nsolns",&Nsolns);
  ch_top->SetBranchAddress("tt_mass",&ttmass);
  ch_top->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_top->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_top->GetEntries(); i++) {
    ch_top->GetEntry(i);
    /* //for some reason it runs extremely slowly with this uncommented
    if ( (Region=="Signal") && (ttmass>450) ) {
      //response.Fill (observable, observable_gen, weight);
      fillUnderOverFlow(hTop, observable, weight, Nsolns);
      fillUnderOverFlow(hTop_gen, observable_gen, weight, Nsolns);
      fillUnderOverFlow(hMeas, observable, weight, Nsolns);
      fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
      fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
      if( combineLepMinus ) {
	      //response.Fill (observableMinus, observableMinus_gen, weight);
	      fillUnderOverFlow(hTop, observableMinus, weight, Nsolns);
	      fillUnderOverFlow(hTop_gen, observableMinus_gen, weight, Nsolns);
	      fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
	      fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
	      fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
	  }
    }
    */ 
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) {
      //response.Fill (observable, observable_gen, weight);
      fillUnderOverFlow(hTop, observable, weight, Nsolns);
      fillUnderOverFlow(hTop_gen, observable_gen, weight, Nsolns);
      fillUnderOverFlow(hMeas, observable, weight, Nsolns);
      fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
      fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
      if( combineLepMinus ) {
	      //response.Fill (observableMinus, observableMinus_gen, weight);
	      fillUnderOverFlow(hTop, observableMinus, weight, Nsolns);
	      fillUnderOverFlow(hTop_gen, observableMinus_gen, weight, Nsolns);
	      fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
	      fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
	      fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
	  }
    }
    if ( (Region=="") && (iVar<2) ) {
      //response.Fill (observable, observable_gen, weight);
      fillUnderOverFlow(hTop, observable, weight, Nsolns);
      fillUnderOverFlow(hTop_gen, observable_gen, weight, Nsolns);
      fillUnderOverFlow(hMeas, observable, weight, Nsolns);
      fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
      fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
      if( combineLepMinus ) {
	      //response.Fill (observableMinus, observableMinus_gen, weight);
	      fillUnderOverFlow(hTop, observableMinus, weight, Nsolns);
	      fillUnderOverFlow(hTop_gen, observableMinus_gen, weight, Nsolns);
	      fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
	      fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
	      fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
	  }
    }

    //if(i % 10000 == 0) cout<<i<<" "<<ch_top->GetEntries()<<endl;

  }

  RooUnfoldResponse response (hMeas, hTrue, hTrue_vs_Meas);
  //hTrue_vs_Meas = (TH2D*) response.Hresponse()->Clone(); 

  hData_bkgSub= (TH1D*) hData->Clone();
  hData_bkgSub->Add(hBkg,-1.0);
 
  if (unfoldingType==0)
    {
      RooUnfoldSvd unfold_svd (&response, hData_bkgSub, kterm); 
      unfold_svd.Setup(&response, hData_bkgSub);
      unfold_svd.IncludeSystematics(includeSys);
      hData_unfolded = (TH1D*) unfold_svd.Hreco();  
      m_unfoldE = unfold_svd.Ereco(); 
    }
  else 
    if (unfoldingType==1)
      {
	RooUnfoldTUnfold unfold_rooTUnfold (&response, hData_bkgSub, TUnfold::kRegModeCurvature); 
	unfold_rooTUnfold.Setup(&response, hData_bkgSub);
	unfold_rooTUnfold.FixTau(tau);
	unfold_rooTUnfold.IncludeSystematics(includeSys);
	hData_unfolded = (TH1D*) unfold_rooTUnfold.Hreco();  
	m_unfoldE = unfold_rooTUnfold.Ereco(); 
      }
  else
    if (unfoldingType==2)
      {
	TUnfold unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature);  
	unfold_TUnfold.SetInput(hTop);
	//Double_t biasScale=1.0;
	unfold_TUnfold.SetBias(hTop_gen);
	unfold_TUnfold.DoUnfold(tau);
	unfold_TUnfold.GetOutput(hData_unfolded);  

	
	TH2D* ematrix=unfold_TUnfold.GetEmatrix("ematrix","error matrix",0,0);
	for (Int_t cmi= 0; cmi<nbins1D; cmi++) {
	  for (Int_t cmj= 0; cmj<nbins1D; cmj++) {
	    m_unfoldE(cmi,cmj)= ematrix->GetBinContent(cmi+1,cmj+1);
	  }
	}    
      }
    else cout<<"Unfolding TYPE not Specified"<<"\n";


  if (unfoldingType==0){
    TCanvas* c_d = new TCanvas("c_d","c_d",500,500); 
    TH1D* dvec=unfold_svd.Impl()->GetD();
    dvec->Draw();
    c_d->SetLogy();
    c_d->SaveAs("D_"+observablename+Region+".pdf");
  }

  TCanvas* c_resp = new TCanvas("c_resp","c_resp");
  TH2D* hResp=(TH2D*) response.Hresponse();
  gStyle->SetPalette(1);
  hResp->GetXaxis()->SetTitle(xaxislabel+"_{gen}");
  hResp->GetYaxis()->SetTitle(xaxislabel);
  hResp->Draw("COLZ");
  c_resp->SaveAs("Response_"+observablename+Region+".pdf");
    
  TFile *file = new TFile("../acceptance/mcnlo/accept_"+acceptanceName+".root");
  TH1D *acceptM = (TH1D*) file->Get("accept_"+acceptanceName);
  acceptM->Scale(1.0/acceptM->Integral());

  TH1D *denominatorM = (TH1D*) file->Get("denominator_"+acceptanceName);

  for (Int_t i= 1; i<=nbins1D; i++) {

    if (acceptM->GetBinContent(i)!=0) {
            hData_unfolded->SetBinContent(i, hData_unfolded->GetBinContent(i)*1.0/acceptM->GetBinContent(i));
            hData_unfolded->SetBinError  (i, hData_unfolded->GetBinError  (i)*1.0/acceptM->GetBinContent(i));
            
            hData_unfolded_arccos->SetBinContent(nbins1D +1 - i, hData_unfolded->GetBinContent(i));
            hData_unfolded_arccos->SetBinError  (nbins1D +1 - i, hData_unfolded->GetBinError(i));
    }

    if (acceptM->GetBinContent(i)!=0) {
      hTop_gen->SetBinContent(i, hTop_gen->GetBinContent(i)*1.0/acceptM->GetBinContent(i));
      hTop_gen->SetBinError  (i, hTop_gen->GetBinError(i)  *1.0/acceptM->GetBinContent(i));
      
      hTop_gen_arccos->SetBinContent(nbins1D +1 - i, hTop_gen->GetBinContent(i) );
      hTop_gen_arccos->SetBinError  (nbins1D +1 - i, hTop_gen->GetBinError(i) );
    }  
  } 

  //scaling is now moved to after Afb is calculated
  //double dataIntegral = hData_unfolded->Integral();
  
  for(int l=0;l<nbins1D;l++){
      for(int j=0;j<nbins1D;j++){
	double corr = 1.0 / ( acceptM->GetBinContent(l+1) * acceptM->GetBinContent(j+1) );
	//corr = corr * pow(xsection / dataIntegral,2) ;
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
  
  GetAfb(denominatorM, Afb, AfbErr);
  cout<<" True Top from acceptance denominator: "<< Afb <<" +/-  "<< AfbErr<<"\n";

  GetCorrectedAfb(hData_unfolded, m_correctE, Afb, AfbErr);
  cout<<" Unfolded: "<< Afb <<" +/-  "<< AfbErr<<"\n";
  
  //GetCorrectedAfb(hData_unfolded_arccos, m_correctE, Afb, AfbErr);
  //cout<<" Unfolded (from arccos histo): "<< Afb <<" +/-  "<< AfbErr<<"\n";
  
  //GetAfbBinByBin(hData_unfolded);

  //GetAfb(hData_unfolded, Afb, AfbErr);
  //cout<<" Unfolded (ignoring correlation): "<< Afb <<" +/-  "<< AfbErr<<"\n";

  vector<double> afb_bins;
  vector<double> afb_bins;  
  GetCorrectedAfbBinByBin(hData_unfolded, m_correctE, afb_bins, afb_bins);


  //scale to total xsec with option "width",  so that differential xsec is plotted
  //hData_unfolded->Scale(xsection/hData_unfolded->Integral(),"width");
  //hTop_gen->Scale(xsection/hTop_gen->Integral(),"width");
  hData_unfolded->Scale(1./hData_unfolded->Integral(),"width");
  hTop_gen->Scale(1./hTop_gen->Integral(),"width");
  hData_unfolded_arccos->Scale(1./hData_unfolded_arccos->Integral(),"width");
  hTop_gen_arccos->Scale(1./hTop_gen_arccos->Integral(),"width");
  
  if(observablename=="lep_azimuthal_asymmetry") for(int i=1;i<nbins1D+1;i++){ cout<<i<<" bin = "<<hData_unfolded_arccos->GetBinContent(i)<<" +/- "<<hData_unfolded_arccos->GetBinError(i)<<endl; }
  else for(int i=1;i<nbins1D+1;i++){ cout<<i<<" bin = "<<hData_unfolded->GetBinContent(i)<<" +/- "<<hData_unfolded->GetBinError(i)<<endl; }


  TCanvas* c_test = new TCanvas("c_final","c_final",500,500); 
  if(observablename=="lep_azimuthal_asymmetry") {
  hData_unfolded_arccos->GetXaxis()->SetTitle("#Delta#phi_{l+l-}");
  hData_unfolded_arccos->GetYaxis()->SetTitle("1/#sigma d#sigma/d(#Delta#phi_{l+l-})");
  hData_unfolded_arccos->SetMinimum(0.0);
  hData_unfolded_arccos->SetMaximum( 2.0* hData_unfolded_arccos->GetMaximum());
  hData_unfolded_arccos->SetMarkerStyle(23);
  hData_unfolded_arccos->SetMarkerSize(1.5);
  hData_unfolded_arccos->Draw("E");
  hData_unfolded_arccos->SetLineWidth(lineWidth);
  hTop_gen_arccos->SetLineWidth(lineWidth);
  hTop_gen_arccos->SetLineColor(TColor::GetColorDark(kGreen));
  hTop_gen_arccos->SetFillColor(TColor::GetColorDark(kGreen));
  hTop_gen_arccos->SetFillStyle(3353);
  hTop_gen_arccos->Draw("hist same");
  hData_unfolded_arccos->Draw("E same");
  }
  else {
  hData_unfolded->GetXaxis()->SetTitle(xaxislabel);
  hData_unfolded->GetYaxis()->SetTitle("1/#sigma d#sigma/d("+xaxislabel+")");
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
  }

  TLegend* leg1=new TLegend(0.55,0.62,0.9,0.838,NULL,"brNDC");                                                                           
  leg1->SetEntrySeparation(100);                                                                                                          
  leg1->SetFillColor(0);                                                                                                                  
  leg1->SetLineColor(0);                                                                                                                   
  leg1->SetBorderSize(0);                    
  leg1->SetTextSize(0.03);
  if(observablename=="lep_azimuthal_asymmetry") {
  	leg1->AddEntry(hData_unfolded_arccos, "( Data - BG ) Unfolded");  
  	leg1->AddEntry(hTop_gen_arccos,    "SM parton level (mc@nlo)", "F"); 
  }
  else{
  	leg1->AddEntry(hData_unfolded, "( Data - BG ) Unfolded");  
  	leg1->AddEntry(hTop_gen,    "SM parton level (mc@nlo)", "F");                                                               
  }
  leg1->Draw();

  TPaveText *pt1 = new TPaveText(0.19, 0.85, 0.42, 0.89, "brNDC");
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);

  TText *blah;
  blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);
  pt1->Draw();

  c_test->SaveAs("finalplot_unfolded_"+observablename+Region+".pdf");
  }

  myfile.close();
}

#ifndef __CINT__
int main () { AfbUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
