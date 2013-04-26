
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
#include "tdrStyle.C"

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

const Double_t _ttMassCut=450.0;
//const Double_t _topScalingFactor=8988.5/8255.4;
const Double_t _topScalingFactor=9097.0/9344.0;

int lineWidth=3;

//==============================================================================
// Example Unfolding
//==============================================================================

void GetAfb(TH1D* h, Float_t &afb, Float_t  &afberr){
  Int_t nbins = h->GetNbinsX();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  //  afb = (event_plus-event_minus)/(event_plus+event_minus);
  //  afberr   = sqrt(4*(event_plus*event_minus)/(event_total*event_total*event_total));
  Double_t event_plus_err;
  Double_t event_minus_err;

  event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_total = event_plus + event_minus;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
		     + event_minus*event_minus*event_plus_err*event_plus_err)/
		  (event_total*event_total*event_total*event_total));

}

// This function is only needed for proper error estimation of the unfolded distributions, not used in this macro
void GetCorrectedAfb(TH1D* histogram, TMatrixD &covarianceM, Float_t &afb, Float_t  &afberr){

    //Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

    //Get histogram info
    int nbins = histogram->GetNbinsX();
    double n[16];
    for(int i=0;i<nbins;i++){
      n[i] = histogram->GetBinContent(i+1);
    }

    //Setup Alpha Vector
    double alpha[16], beta[16];
    for(int i=0;i<nbins;i++) if(i < nbins/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

    //Components of the error calculation
    double sum_n = 0.;
    double sum_alpha_n = 0.;
    for(int i=0;i<nbins;i++){
      sum_n += n[i];
      sum_alpha_n += alpha[i] * n[i];
    }

    double dfdn[16];
    for(int i=0;i<nbins;i++){
      dfdn[i] = ( alpha[i] * sum_n - sum_alpha_n ) / pow(sum_n,2);
    }

    //Error Calculation
    afberr = 0.;
    for(int i=0;i<nbins;i++){
      for(int j=0;j<nbins;j++){
	afberr += covarianceM(i,j) * dfdn[i] * dfdn[j];
      }
    }
    afberr = sqrt(afberr);

    //Calculate Afb
    afb = sum_alpha_n / sum_n;

    //    cout<<"AFB = "<<afb<<" "<<afberr<<endl;
}


void AfbRecoLevel(TString Region="")
{
  //setTDRStyle();
gStyle->SetOptFit();
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
cout.precision(2);
//gStyle->SetLabelSize(0.04,"xyz");
// gStyle->SetTitleSize(0.045,"xyz");
 gStyle->SetTitleColor(1, "XYZ");
 gStyle->SetTitleFont(42, "XYZ");
 gStyle->SetTitleSize(0.06, "X");
 gStyle->SetTitleSize(0.06, "Y");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(0.9);
  
  gStyle->SetTitleOffset( 0.1,"Y");
  
  TString observablename;
  TString xaxislabel;

  float xmin=-1.0;
  float xmax= 1.0;

  const int nbins=6;
  //const int nbins=40;
  double xbins[nbins+1];

  int nVars =9;

  Float_t observable, observable_gen, weight, ttmass, ttRapidity, tmass;

  for (Int_t iVar= 0; iVar < nVars; iVar++) {

  switch (iVar)
    {
      //   Lepton Charge Asymmetry
    case 0:
      {
      observablename="lep_charge_asymmetry";
      xaxislabel="|#eta_{l+}|-|#eta_{l-}|";
      xbins[0]=-3.0; xbins[1]=-2.0; xbins[2]=-1.0; xbins[3]=0.0; xbins[4]=1.0; xbins[5]=2.0; xbins[6]=3.0;
      xmin=-3.0;
      xmax= 3.0;
      break;
      }
  //   Lepton Azimuthal Asymmetry
    case 1:
      {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="cos(#Delta#phi_{l+l-})";
      xbins[0]=-1.0; xbins[1]=-0.6; xbins[2]=-0.3; xbins[3]=0.0; xbins[4]=0.3; xbins[5]=0.6; xbins[6]=1.0;
      xmin=-1.0;
      xmax= 1.0;
      break;
      }
  //   Top Charge Asymmetry
    case 2:
      {
      observablename="top_costheta_cms";
      xaxislabel="cos(#theta_{top})";
      xbins[0]=-1.0; xbins[1]=-0.6; xbins[2]=-0.3; xbins[3]=0.0; xbins[4]=0.3; xbins[5]=0.6; xbins[6]=1.0;
     
      xmin=-1.0;
      xmax= 1.0;
      break;
      }

  //   Top Polarization
    case 3:
      {
	observablename="lep_costheta_cms";
	//observablename="lepMinus_costheta_cms";
	xaxislabel="cos(#theta_{l}^{+})";
	//xaxislabel="cos(#theta_{l}^{-})";
      xbins[0]=-1.0; xbins[1]=-0.6; xbins[2]=-0.3; xbins[3]=0.0; xbins[4]=0.3; xbins[5]=0.6; xbins[6]=1.0;
      //xbins[0]=-1.0; xbins[1]=-0.95; xbins[2]=-0.90; xbins[3]=-0.85; xbins[4]=-0.80; xbins[5]=-0.75; xbins[6]=-0.70 ; xbins[7]=-0.65; xbins[8]=-0.60; xbins[9]=-0.55; xbins[10]=-0.50; xbins[11]=-0.45; xbins[12]=-0.40 ; xbins[13]=-0.35 ; xbins[14]=-0.30 ;xbins[15]=-0.25; xbins[16]=-0.20; xbins[17]=-0.15; xbins[18]=-0.10; xbins[19]=-0.05; xbins[20]=0.00; xbins[21]=0.05 ; xbins[22]=0.10; xbins[23]=0.15; xbins[24]=0.20; xbins[25]=0.25; xbins[26]=0.30; xbins[27]=0.35 ; xbins[28]=0.40 ; xbins[29]=0.45;  xbins[30]=0.50; xbins[31]=0.55 ; xbins[32]=0.60; xbins[33]=0.65; xbins[34]=0.70; xbins[35]=0.75; xbins[36]=0.80 ; xbins[37]=0.85; xbins[38]=0.90 ; xbins[39]=0.95; xbins[40]=1.00 ;
      xmin=-1.0;
      xmax= 1.0;
      break;
      }
  //   Top Spin Correlation
    case 4:
      {
      observablename="top_spin_correlation";
      xaxislabel="cos(#theta_{l+,n})cos(#theta_{l-,n})";
      xbins[0]=-1.0; xbins[1]=-0.6; xbins[2]=-0.3; xbins[3]=0.0; xbins[4]=0.3; xbins[5]=0.6; xbins[6]=1.0;
      xmin=-1.0;
      xmax= 1.0;
      break;
      }
 //   Top Asy I
    case 5:
      {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="|#eta_{top}|-|#eta_{tbar}|";
      xbins[0]=-6.0; xbins[1]=-2.0; xbins[2]=-1.0; xbins[3]=0.0; xbins[4]=1.0; xbins[5]=2.0; xbins[6]=6.0;
      xmin=-5.0;
      xmax= 5.0;
      break;
      }
      //   Top Asy I
    case 6:
      {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="(y_{top}-y_{tbar})(y_{top}+y_{tbar})";
      xbins[0]=-6.0; xbins[1]=-2.0; xbins[2]=-1.0; xbins[3]=0.0; xbins[4]=1.0; xbins[5]=2.0; xbins[6]=6.0;
      xmin=-6.0;
      xmax= 6.0;
      break;
      }
    case 7:
      {
	observablename="top_rapidtiydiff_Marco";
	xaxislabel="|y_{top}|-|y_{tbar}|";
	xbins[0]=-6.0; xbins[1]=-2.0; xbins[2]=-1.0; xbins[3]=0.0; xbins[4]=1.0; xbins[5]=2.0; xbins[6]=6.0;
	xmin=-6.0;
	xmax= 6.0;
	break;
      }
    case 8:
      {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="#Delta#phi_{l+l-}";
      xbins[0]=0.0; xbins[1]=0.52; xbins[2]=1.04; xbins[3]=1.56; xbins[4]=2.08; xbins[5]=2.6; xbins[6]=3.14;
      xmin=0.0;
      xmax=3.14;
      break;
      }
    default:
      {
      cout<<"Set the variable switch";
      }
    }



  TH1D* hData= new TH1D ("Data", "Data",    nbins, xbins);
  TH1D* hBkg = new TH1D ("Background",  "Background",    nbins, xbins);
  TH1D* hTop = new TH1D ("Top",  "Top",    nbins, xbins);

  //  Now test with data and with BKG subtraction

  hData->Sumw2();
  hTop->Sumw2();
  hBkg->Sumw2();

  const char* var=observablename;
  TCut baseline = "1";
  if  (Region=="Signal2")  
    baseline = "t_mass > 0 && tt_mass >450 && ((ttRapidity > 2.0)||(ttRapidity < -2.0))";
  if  (Region=="Signal")  
    baseline = "t_mass > 0 && tt_mass >450";
  if ( (Region=="") && (iVar>=2) ) 
    baseline = "t_mass > 0 && tt_mass > 0";
  if ( (Region=="") && (iVar<2) ) 
    baseline = "1";

  TChain *ch_data = new TChain("tree");
  ch_data->Add("data.root");
  ch_data->Draw(Form("%s >> %s", var, "Data"),     baseline*"weight");

  TChain *ch_top = new TChain("tree");
  ch_top->Add("ttdil.root");
  ch_top->Draw(Form("%s >> %s", var, "Top"), baseline*"weight");

  TChain *ch_bkg = new TChain("tree");
  ch_bkg->Add("ttotr.root");
  ch_bkg->Add("wjets.root");
  ch_bkg->Add("DYee.root");
  ch_bkg->Add("DYmm.root");
  ch_bkg->Add("DYtautau.root");
  ch_bkg->Add("tw.root");
  ch_bkg->Add("VV.root");
  ch_bkg->Draw(Form("%s >> %s", var, "Background"),       baseline*"weight");
  
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

  hData->SetLineWidth(lineWidth+2);
  hTop->SetLineWidth(lineWidth);

  hTop->SetLineColor(TColor::GetColorDark(kGreen));
  hTop->SetFillColor(TColor::GetColorDark(kGreen));

  hTop->SetFillStyle(3353);
  hBkg->SetLineColor(kYellow);
  hBkg->SetFillColor(kYellow);
  
  for (Int_t i= 1; i<=nbins; i++) {
     hData->SetBinContent(i, hData->GetBinContent(i)/hData->GetBinWidth(i) );
     hData->SetBinError(i, hData->GetBinError(i)/hData->GetBinWidth(i) );
     hTop->SetBinContent(i, hTop->GetBinContent(i)/hTop->GetBinWidth(i) );
     hTop->SetBinError(i, hTop->GetBinError(i)/hTop->GetBinWidth(i) );
     hBkg->SetBinContent(i, hBkg->GetBinContent(i)/hBkg->GetBinWidth(i) );
     hBkg->SetBinError(i, hBkg->GetBinError(i)/hBkg->GetBinWidth(i) );
  }
  
  // TH1D* hTotal_MC=(TH1D*)hTop->Clone();
//   hTotal_MC->Add(hBkg);
//   double KS3 = hTotal_MC->KolmogorovTest(hData,"UO");
//   std::cout <<"K-S powheg "<<KS3<< std::endl;

  THStack *hs = new THStack("hs","Stacked Top+BG");

  hs->Add(hBkg);
  hs->Add(hTop);

  hs->SetMinimum(0.0);
  hs->SetMaximum( 2.0* hs->GetMaximum());
  hs->Draw("hist");
  hs->GetXaxis()->SetTitle(xaxislabel);
  hs->GetYaxis()->SetTitleOffset(1.4);
  hs->GetYaxis()->SetTitle("Events/"+xaxislabel+"");

  hData->Draw("E same");

  TLegend* leg1=new TLegend(0.6,0.62,0.9,0.838,NULL,"brNDC");
  leg1->SetFillStyle(0);
  leg1->SetEntrySeparation(100);  
  leg1->SetBorderSize(0);                                                                                 
  leg1->SetTextSize(0.03);
  leg1->AddEntry(hData, "Data");
  leg1->AddEntry(hTop,  "t#bar{t} (dileptonic)");                                                               
  leg1->AddEntry(hBkg,  "Background");                                                               
  leg1->Draw();                

  TPaveText *pt1 = new TPaveText(0.19, 0.82, 0.42, 0.86, "brNDC");
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);

  TText *blah;
  blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at #sqrt{s}=7 TeV");
  blah->SetTextSize(0.035);
  blah->SetTextAlign(11);

  TPaveText *pt2 = new TPaveText(0.20, 0.70, 0.40, 0.75, "brNDC");
  pt2->SetName("pt2name");
  pt2->SetBorderSize(0);
  pt2->SetFillStyle(0);
  
 //  TText *blah2;

//   TString KS3_temp = formatFloat(KS3,"%6.2f");
//   KS3_temp.ReplaceAll(" " , "" );
//   KS3_temp = TString("   K-S: ") +  KS3_temp;
//   blah2 = pt2->AddText(KS3_temp.Data());
//   blah2->SetTextSize(0.032);
//   blah2->SetTextAlign(11);
//   blah2->SetTextColor(kGreen+3);

  pt1->Draw();
  // pt2->Draw();
  c_test->SaveAs("finalplot_"+observablename+Region+".pdf");

  delete hData;
  delete hBkg;
  delete hTop;
  delete hTop_bkgAdd;
  delete c_test;

  }

}

#ifndef __CINT__
int main () { AfbRecoLevel(); return 0; }  // Main program when run stand-alone
#endif
