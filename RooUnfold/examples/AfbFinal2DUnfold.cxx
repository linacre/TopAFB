
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

#include "src/RooUnfold.h"
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

 TRandom3* random = new TRandom3();                                                                                                        
 random->SetSeed(5);

 ofstream myfile;
 myfile.open ("summary_2Dunfolding.txt");
 cout.rdbuf(myfile.rdbuf());

 Float_t observable, observable_gen, ttmass, ttmass_gen, ttRapidity, tmass;
 Double_t weight;
 Int_t Nsolns;

  for (Int_t iVar= 0; iVar < nVars; iVar++) {


    Initialize2DBinning(iVar);

  TH1D* hData= new TH1D ("Data_BkgSub", "Data with background subtracted",    nbins2D, xbins2D);
  TH1D* hBkg = new TH1D ("Background",  "Background",    nbins2D, xbins2D);
  TH1D* hTop = new TH1D ("Top",  "Top",    nbins2D, xbins2D);
  TH1D* hTop_gen= new TH1D ("Top_Gen", "Top Gen", nbins2D, xbins2D); 
  TH1D* hData_unfolded= new TH1D ("Data_Unfold", "Data with background subtracted and unfolded", nbins2D, xbins2D);
 

  TH1D* hTrue= new TH1D ("true", "Truth",    nbins2D, xbins2D);
  TH1D* hMeas= new TH1D ("meas", "Measured", nbins2D, xbins2D);

  TH2D* hTrue_vs_Meas= new TH2D ("true_vs_meas", "True vs Measured", nbins2D, xbins2D, nbins2D, xbins2D);
  
  TH1D* hData_bkgSub;

  hData->Sumw2();
  hBkg->Sumw2();
  hTop->Sumw2();
  hTop_gen->Sumw2();
  hData_unfolded->Sumw2();
  hTrue->Sumw2();
  hMeas->Sumw2();      


  TMatrixD m_unfoldE(nbins2D,nbins2D);
  TMatrixD m_correctE(nbins2D,nbins2D);

  
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
  ch_data->SetBranchAddress("weight",&weight);
  ch_data->SetBranchAddress("Nsolns",&Nsolns);
  ch_data->SetBranchAddress("tt_mass",&ttmass);
  ch_data->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_data->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_data->GetEntries(); i++) {
    ch_data->GetEntry(i);
    if ( (Region=="Signal") && (ttmass>450) )  
      fillUnderOverFlow(hData, sign(observable)*ttmass, weight, Nsolns);
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) 
      fillUnderOverFlow(hData, sign(observable)*ttmass, weight, Nsolns);    
    if ( (Region=="") && (iVar<2) ) 
      fillUnderOverFlow(hData, sign(observable)*ttmass, weight, Nsolns);     
  }

  ch_bkg->SetBranchAddress(observablename,    &observable);
  ch_bkg->SetBranchAddress("weight",&weight);
  ch_bkg->SetBranchAddress("Nsolns",&Nsolns);
  ch_bkg->SetBranchAddress("tt_mass",&ttmass);
  ch_bkg->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_bkg->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_bkg->GetEntries(); i++) {
    ch_bkg->GetEntry(i);
    if ( (Region=="Signal") && (ttmass>450) )  
      fillUnderOverFlow(hBkg, sign(observable)*ttmass, weight, Nsolns);
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) 
      fillUnderOverFlow(hBkg, sign(observable)*ttmass, weight, Nsolns);
    if ( (Region=="") && (iVar<2) ) 
      fillUnderOverFlow(hBkg, sign(observable)*ttmass, weight, Nsolns);
  }

  ch_top->SetBranchAddress(observablename,    &observable);
  ch_top->SetBranchAddress(observablename+"_gen",&observable_gen);
  ch_top->SetBranchAddress("weight",&weight);
  ch_top->SetBranchAddress("Nsolns",&Nsolns);
  ch_top->SetBranchAddress("tt_mass",&ttmass);
  ch_top->SetBranchAddress("tt_mass_gen",&ttmass_gen);
  ch_top->SetBranchAddress("ttRapidity",&ttRapidity);
  ch_top->SetBranchAddress("t_mass",&tmass);

  for (Int_t i= 0; i<ch_top->GetEntries(); i++) {
    ch_top->GetEntry(i);
    if ( (Region=="") && (iVar>=2) && (ttmass>0) ) {
      //response.Fill (observable, observable_gen, weight);
      fillUnderOverFlow(hMeas, sign(observable)*ttmass, weight, Nsolns);
      fillUnderOverFlow(hTrue, sign(observable_gen)*ttmass_gen, weight, Nsolns);
      fillUnderOverFlow(hTrue_vs_Meas, sign(observable)*ttmass, sign(observable_gen)*ttmass_gen, weight, Nsolns);
    }
    if ( (Region=="") && (iVar<2) ) {
      //response.Fill (observable, observable_gen, weight);
      fillUnderOverFlow(hMeas, sign(observable)*ttmass, weight, Nsolns);
      fillUnderOverFlow(hTrue, sign(observable_gen)*ttmass_gen, weight, Nsolns);
      fillUnderOverFlow(hTrue_vs_Meas, sign(observable)*ttmass, sign(observable_gen)*ttmass_gen, weight, Nsolns);
    }
  }
  
  RooUnfoldResponse response (hMeas, hTrue, hTrue_vs_Meas);

  TCanvas* c_test = new TCanvas("c_final","c_final",500,500); 

  hData->SetLineWidth(lineWidth+2);
  hTop->SetLineWidth(lineWidth);

  hTrue->SetLineColor(TColor::GetColorDark(kGreen));
  hTrue->SetFillColor(TColor::GetColorDark(kGreen));

  hTop->SetFillStyle(3353);
  hTrue->SetFillStyle(3353);
  hBkg->SetLineColor(kYellow);
  hBkg->SetFillColor(kYellow);
  

  THStack *hs = new THStack("hs","Stacked Top+BG");

  hs->Add(hBkg);
  hs->Add(hTrue);

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
	unfold_TUnfold.SetInput(hMeas);
	//Double_t biasScale=1.0;
	unfold_TUnfold.SetBias(hTrue);
	unfold_TUnfold.DoUnfold(tau);
	unfold_TUnfold.GetOutput(hData_unfolded);  

	
	TH2D* ematrix=unfold_TUnfold.GetEmatrix("ematrix","error matrix",0,0);
	for (Int_t cmi= 0; cmi<nbins2D; cmi++) {
	  for (Int_t cmj= 0; cmj<nbins2D; cmj++) {
	    m_unfoldE(cmi,cmj)= ematrix->GetBinContent(cmi+1,cmj+1);
	  }
	}    
      }
    else cout<<"Unfolding TYPE not Specified"<<"\n";



  TFile *file = new TFile("../acceptance/mcnlo/accept_"+acceptanceName+".root");

  TH2D *acceptM_2d = (TH2D*) file->Get("accept_"+acceptanceName+"_mtt");

  TH1D* acceptM = new TH1D ("accept", "accept",    nbins2D, xbins2D);
  acceptM->SetBinContent(1,acceptM_2d->GetBinContent(1,3));
  acceptM->SetBinContent(2,acceptM_2d->GetBinContent(1,2));
  acceptM->SetBinContent(3,acceptM_2d->GetBinContent(1,1));

  acceptM->SetBinContent(4,acceptM_2d->GetBinContent(2,1));
  acceptM->SetBinContent(5,acceptM_2d->GetBinContent(2,2));
  acceptM->SetBinContent(6,acceptM_2d->GetBinContent(2,3));

  acceptM->Scale(1.0/acceptM->Integral());


  TH2D *denomM_2d = (TH2D*) file->Get("denominator_"+acceptanceName+"_mtt");
  TH1D* denomM = new TH1D ("denom", "denom",    nbins2D, xbins2D);
  
  denomM->SetBinContent(1,denomM_2d->GetBinContent(1,3));
  denomM->SetBinContent(2,denomM_2d->GetBinContent(1,2));
  denomM->SetBinContent(3,denomM_2d->GetBinContent(1,1));

  denomM->SetBinContent(4,denomM_2d->GetBinContent(2,1));
  denomM->SetBinContent(5,denomM_2d->GetBinContent(2,2));
  denomM->SetBinContent(6,denomM_2d->GetBinContent(2,3));  
  
  TH1D* denomM_0 = new TH1D ("denominator0", "denominator0",    2, -1500.,1500.);
  TH1D* denomM_1 = new TH1D ("denominator1", "denominator1",    2, -1500.,1500.);
  TH1D* denomM_2 = new TH1D ("denominator2", "denominator2",    2, -1500.,1500.);

  denomM_2->SetBinContent(1,denomM_2d->GetBinContent(1,3));
  denomM_1->SetBinContent(1,denomM_2d->GetBinContent(1,2));
  denomM_0->SetBinContent(1,denomM_2d->GetBinContent(1,1));

  denomM_0->SetBinContent(2,denomM_2d->GetBinContent(2,1));
  denomM_1->SetBinContent(2,denomM_2d->GetBinContent(2,2));
  denomM_2->SetBinContent(2,denomM_2d->GetBinContent(2,3));





  for (Int_t i= 1; i<=nbins2D; i++) {

    if (acceptM->GetBinContent(i)!=0) {
      hData_unfolded->SetBinContent(i, hData_unfolded->GetBinContent(i)*1.0/acceptM->GetBinContent(i));
      hData_unfolded->SetBinError  (i, hData_unfolded->GetBinError  (i)*1.0/acceptM->GetBinContent(i));
    }

    if (acceptM->GetBinContent(i)!=0) {
      hTrue->SetBinContent(i, hTrue->GetBinContent(i)*1.0/acceptM->GetBinContent(i));
      hTrue->SetBinError  (i, hTrue->GetBinError(i)  *1.0/acceptM->GetBinContent(i));
    }
  }
  

    for(int l=0;l<nbins2D;l++){
      for(int j=0;j<nbins2D;j++){
	double corr = 1.0 / ( acceptM->GetBinContent(l+1) * acceptM->GetBinContent(j+1) );
	//corr = corr * pow(xsection / dataIntegral,2) ;
	m_correctE(l,j) = m_unfoldE(l,j)*corr;
      }
    }
  

  //==================================================================
  // ============== Print the assymetry =============================
  cout<<"========= Variable: "<<observablename <<"===================\n";
  Float_t Afb, AfbErr;

  GetAfb(hData, Afb, AfbErr);
  cout<<" Data: "<< Afb <<" +/-  "<< AfbErr<<"\n";

  GetAfb(hTrue, Afb, AfbErr);
  cout<<" True Top: "<< Afb <<" +/-  "<< AfbErr<<"\n";
    
  GetAfb(denomM, Afb, AfbErr);
  cout<<" True Top from acceptance denominator: "<< Afb <<" +/-  "<< AfbErr<<"\n";
    
  Float_t AfbG[3];
    
  GetAfb(denomM_0, AfbG[0], AfbErr);
  cout<<" True Top 0 from acceptance denominator: "<< AfbG[0] <<" +/-  "<< AfbErr<<"\n";
  GetAfb(denomM_1, AfbG[1], AfbErr);
  cout<<" True Top 1 from acceptance denominator: "<< AfbG[1] <<" +/-  "<< AfbErr<<"\n";
  GetAfb(denomM_2, AfbG[2], AfbErr);
  cout<<" True Top 2 from acceptance denominator: "<< AfbG[2] <<" +/-  "<< AfbErr<<"\n";

  GetCorrectedAfb(hData_unfolded, m_correctE, Afb, AfbErr);
  cout<<" Unfolded: "<< Afb <<" +/-  "<< AfbErr<<"\n";

  //  hData_unfolded->Scale(1./hData_unfolded->Integral(),"width");
  //  hTrue->Scale(1./hTrue->Integral(),"width");
  
  TCanvas* c_test = new TCanvas("c_final","c_final",500,500); 
  hData_unfolded->GetXaxis()->SetTitle("M_{t#bar t}");
  hData_unfolded->GetYaxis()->SetTitle("d#sigma/dM_{t#bar t}");
  hData_unfolded->SetMinimum(0.0);
  hData_unfolded->SetMaximum( 2.0* hData_unfolded->GetMaximum());
  hData_unfolded->SetMarkerStyle(23);
  hData_unfolded->SetMarkerSize(2.0);
  hData_unfolded->Draw("E");
  hData_unfolded->SetLineWidth(lineWidth);
  hTrue->SetLineWidth(lineWidth);
  hTrue->SetLineColor(TColor::GetColorDark(kGreen));
  hTrue->SetFillColor(TColor::GetColorDark(kGreen));
  hTrue->Draw("hist same");
  hData_unfolded->Draw("E same");

  leg1=new TLegend(0.6,0.62,0.9,0.838,NULL,"brNDC");   
  leg1->SetEntrySeparation(100);                                                                                                          
  leg1->SetFillColor(0);                                                                                                                  
  leg1->SetLineColor(0);                                                                                                                   
  leg1->SetBorderSize(0);                                                                                                                  
  leg1->SetTextSize(0.03);                                                                              
  leg1->AddEntry(hData_unfolded, "( Data-BG ) Unfolded");                                                                                       
  leg1->AddEntry(hTrue,    "mc@nlo parton", "F");                                                               
  leg1->Draw();                

  c_test->SaveAs("Mtt_2D_unfolded_"+observablename+Region+".pdf");

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
 
   TH1D* hTop_AfbVsMtt = new TH1D ("Top_AfbVsMtt",  "Top_AfbVsMtt",  3, xbins2D);
   for (int nb=0; nb<3; nb++) {
       hTop_AfbVsMtt->SetBinContent(nb+1,AfbG[nb]);
       hTop_AfbVsMtt->SetBinError(nb+1,0);
   }
   
  tdrStyle->SetErrorX(0.5);
  hAfbVsMtt->SetMinimum(-0.3);
  hAfbVsMtt->SetMaximum( 0.3);
  hAfbVsMtt->SetLineWidth( 2.0 );
  hAfbVsMtt->Draw("E");
  hTop_AfbVsMtt->SetLineColor(kGreen);
  hTop_AfbVsMtt->SetMarkerColor(kGreen);
  hTop_AfbVsMtt->SetMarkerSize(0);
  hTop_AfbVsMtt->SetLineWidth( 2.0 );
  hAfbVsMtt->GetYaxis()->SetTitle(""+xaxislabel);
  hAfbVsMtt->GetYaxis()->SetTitleOffset(1.2);
  hAfbVsMtt->GetXaxis()->SetTitle("M_{t#bar t} (GeV/c^{2})");
  hTop_AfbVsMtt->Draw("E same");
  
  leg1=new TLegend(0.6,0.72,0.9,0.938,NULL,"brNDC");   
  leg1->SetEntrySeparation(100);                                                                                                          
  leg1->SetFillColor(0);                                                                                                                  
  leg1->SetLineColor(0);                                                                                                                   
  leg1->SetBorderSize(0);                                                                                                                  
  leg1->SetTextSize(0.03);                                                                              
  leg1->AddEntry(hAfbVsMtt, "data");                                                                                       
  leg1->AddEntry(hTop_AfbVsMtt,    "mc@nlo parton level");                                                               
  leg1->Draw();           
  
  
  c_test->SaveAs("AfbVsMtt_unfolded_"+observablename+Region+".pdf");

  }
  myfile.close();

}

#ifndef __CINT__
int main () { AfbUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
