#include <fstream>
#include "TH1D.h"
#include "TMatrixD.h"
#include "TMath.h"

double xsection = 154.0*0.06451;

int lineWidth=3;

TString observablename;
TString acceptanceName;
TString xaxislabel;
TString asymlabel;

float xmin=-1.0;
float xmax= 1.0;

const int nbins1D=6;
const int nbins2D=6;

double xbins1D[nbins1D+1];
double xbins2D[nbins2D+1];


Double_t stat_corr  [nbins1D]; //errors include syst error in the unfolding
Double_t stat_uncorr[nbins1D]; //errors do not include syst error in the unfolding
Double_t syst_corr  [nbins1D];

Float_t sign(Float_t t) 
{
  if( t >= 0.0 )
    return 1.0;
  else
    return -1.0;
}


void GetAfb(TH1D* h, Float_t &afb, Float_t  &afberr){

  Int_t nbins = h->GetNbinsX();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  //event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"");
  //event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_plus_err,"");
  event_total = event_plus + event_minus;

  //cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<endl;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
    + event_minus*event_minus*event_plus_err*event_plus_err)/
    (event_total*event_total*event_total*event_total));

}


void GetAfb_integratewidth(TH1D* h, Float_t &afb, Float_t  &afberr){

  Int_t nbins = h->GetNbinsX();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"width");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_plus_err,"width");
  event_total = event_plus + event_minus;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
    + event_minus*event_minus*event_plus_err*event_plus_err)/
    (event_total*event_total*event_total*event_total));

}


//void GetAfbBinByBin(TH1D* h, Float_t &afbbin, Float_t  &afberrbin){
void GetAfbBinByBin(TH1D* h){

  Int_t nbins = h->GetNbinsX();
  const int nbins2 = nbins/2 +1;
  Float_t event_minus[nbins2];
  Float_t event_plus[nbins2];
  Float_t event_total[nbins2];
  Double_t event_plus_err[nbins2];
  Double_t event_minus_err[nbins2];
  Double_t afbbin[nbins2];
  Double_t afberrbin[nbins2];

  Double_t event_plus_total = 0.;
  Double_t event_minus_total = 0.;
  Double_t event_total_total = 0.;

  for(int i=0;i<nbins2;i++){
    //event_minus[i]  = h-> IntegralAndError(i, i, event_plus_err[i],"");
    event_minus[i]  = h-> IntegralAndError(i, i, event_minus_err[i],"");
    event_minus_total += event_minus[i];
    //event_plus[i]   = h-> IntegralAndError(nbins+1-i, nbins+1-i, event_minus_err[i],"");
    event_plus[i]   = h-> IntegralAndError(nbins+1-i, nbins+1-i, event_plus_err[i],"");
    event_plus_total += event_plus[i];
    event_total[i] = event_plus[i] + event_minus[i];
    event_total_total += event_total[i];

    //cout<<event_minus[i]<<" "<<event_minus_err[i]<<" "<<event_plus[i]<<" "<<event_plus_err[i]<<" "<<event_total[i]<<endl;

    afbbin[i] = (event_plus[i]-event_minus[i])/(event_plus[i]+event_minus[i]);
    afberrbin[i]   = sqrt(4*(event_plus[i]*event_plus[i]*event_minus_err[i]*event_minus_err[i] 
      + event_minus[i]*event_minus[i]*event_plus_err[i]*event_plus_err[i])/
      (event_total[i]*event_total[i]*event_total[i]*event_total[i]));
    cout<<i<<" AFB = "<<afbbin[i]<<" +/- "<<afberrbin[i]<<endl;
  }
}


void GetAvsY(TH1* histogram, TMatrixD &covarianceM, vector<double> &myafb, vector<double> &myerr, ofstream& second_output_file){

  myafb.clear();
  myerr.clear();

    //Get info from histogram
  int nbins = histogram->GetNbinsX();
  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1);
  }

    //Output
  double afb[8], err[8];

    //Setup Some Needed Vectors
  double alpha[16], beta[16], dfdn[16];

    //Get Asymmetry for each Y bin
  for(int j=0;j<nbins/2;j++){

    int forBin = nbins/2 + j;
    int bacBin = nbins/2 - j - 1;

    for(int i=0;i<nbins;i++){
      if( i == forBin ){ 
        alpha[i] = 1; 
        beta[i] = 1;
      }else if( i == bacBin ){ 
        alpha[i] = -1; 
        beta[i] = 1;
      }else{
        alpha[i] = 0;
        beta[i] = 0;
      }
    }

    double sum = 0. , diff = 0.;
    for(int i=0;i<nbins;i++){
      sum += beta[i] * n[i];
      diff += alpha[i] * n[i];
    }

      //Calculate Everything
    if(sum > 0){ 

  //Error Calculation
      for(int i=0;i<nbins;i++){
        dfdn[i] = ( alpha[i] * sum - beta[i] * diff ) / pow(sum,2);
      }

      double afberr = 0.;
      for(int i=0;i<nbins;i++){
        for(int k=0;k<nbins;k++){
          afberr += covarianceM(i,k) * dfdn[i] * dfdn[k];
      //if(i==k) cout<<"DAH: "<<n[i]<<" "<<k<<" "<<covarianceM(i,k)<<endl;
        }
      }
      afberr = sqrt(afberr);

      err[j] = afberr;
      afb[j] = diff / sum; 

    }else{ 

      afb[j] = 0.; 
      err[j] = 0.;

    }
    myafb.push_back(afb[j]);
    myerr.push_back(err[j]);

    cout<<j<<" AFB = "<<afb[j]<<" +/- "<<err[j]<<endl;
    second_output_file << acceptanceName << " " << observablename << " AFB" << j << ": " << afb[j] << " +/- " << err[j] << endl;    	
  }
}

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




void GetCorrectedAfbBinByBin(TH1D* histogram, TMatrixD &covarianceM, vector<double> &myafb, vector<double> &myerr, ofstream& second_output_file){

    //Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

  myafb.clear();
  myerr.clear();

    //Get histogram info
  int nbins = histogram->GetNbinsX();
  const int nbins2 = nbins/2;

  Double_t afbbin[nbins2];
  Double_t afberrbin[nbins2];

  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1);
  }

    //Setup Alpha Vector
  double alpha[16], beta[16];
  for(int i=0;i<nbins;i++) if(i < nbins/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

    //Components of the error calculation
  double sum_n[nbins2];
  double sum_alpha_n[nbins2];
  double sum_n_total = 0.;
  double sum_alpha_n_total = 0.;


  for(int i=0;i<nbins2;i++){
    sum_n[i] = n[i] + n[nbins-1-i];
    sum_alpha_n[i] = alpha[i] * n[i] + alpha[nbins-1-i] * n[nbins-1-i];
    sum_n_total += sum_n[i];
    sum_alpha_n_total += sum_alpha_n[i];
  }

  double dfdn[16];
  for(int i=0;i<nbins;i++){
    int k = -999;
    if (i < nbins2) k = i;
    else k = nbins-1-i;
    dfdn[i] = ( alpha[i] * sum_n[k] - sum_alpha_n[k] ) / pow(sum_n[k],2);
  }

    //Error Calculation

  for(int k=0;k<nbins2;k++){
    afberrbin[k] = 0.;
    for(int i=0;i<nbins;i++){
      for(int j=0;j<nbins;j++){
        if( (i==k || i==nbins-1-k ) && (j==k || j==nbins-1-k ) ) {
          afberrbin[k] += covarianceM(i,j) * dfdn[i] * dfdn[j];
              //cout<<covarianceM(i,j)<<" "<<dfdn[i]<<" "<<dfdn[j]<<" "<<endl;
        }
      }
    }
    afberrbin[k] = sqrt(afberrbin[k]);
    afbbin[k] = sum_alpha_n[k] / sum_n[k];
    cout<<k<<" AFB = "<<afbbin[k]<<" +/- "<<afberrbin[k]<<endl;
    second_output_file << acceptanceName << " " << observablename << " AFB" << i << ": " << afbbin[k] << " +/- " << afberrbin[k] << endl;

    myafb.push_back(afbbin[k]);
    myerr.push_back(afberrbin[k]);
  }
  double afb = sum_alpha_n_total / sum_n_total;
    //cout<<"AFB = "<<afb<<endl;
}







void Initialize1DBinning(int iVar){


  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="|#eta_{l+}|-|#eta_{l-}|";
      acceptanceName="lepChargeAsym";
      asymlabel="A_{lepC}";
      xbins1D[0]=-2.0; xbins1D[1]=-0.8; xbins1D[2]=-0.4; xbins1D[3]=0.0; xbins1D[4]=0.4; xbins1D[5]=0.8; xbins1D[6]=2.0;
syst_corr[0] =  0.002272  ; stat_corr[0] =  0.005120  ; stat_uncorr[0] =  0.003483  ;
syst_corr[1] =  0.002421  ; stat_corr[1] =  0.006879  ; stat_uncorr[1] =  0.004675  ;
syst_corr[2] =  0.004388  ; stat_corr[2] =  0.007338  ; stat_uncorr[2] =  0.005284  ;
syst_corr[3] =  0.004508  ; stat_corr[3] =  0.007344  ; stat_uncorr[3] =  0.005309  ;
syst_corr[4] =  0.002648  ; stat_corr[4] =  0.006623  ; stat_uncorr[4] =  0.004689  ;
syst_corr[5] =  0.002842  ; stat_corr[5] =  0.004806  ; stat_uncorr[5] =  0.003518  ;
      xmin=-2.0;
      xmax= 2.0;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="#Delta#phi_{l+l-}";
      acceptanceName="lepAzimAsym2";
      asymlabel="A_{#Delta#phi}";
      Double_t pi = 3.141592653589793;
      xbins1D[0]=0.0; xbins1D[1]=4.*pi/20.; xbins1D[2]=7.*pi/20.; xbins1D[3]=10.*pi/20.; xbins1D[4]=13.*pi/20.; xbins1D[5]=16.*pi/20.; xbins1D[6]=pi;
syst_corr[0] =  0.005592  ; stat_corr[0] =  0.007795  ; stat_uncorr[0] =  0.005697  ;
syst_corr[1] =  0.004763  ; stat_corr[1] =  0.006603  ; stat_uncorr[1] =  0.004705  ;
syst_corr[2] =  0.002718  ; stat_corr[2] =  0.006543  ; stat_uncorr[2] =  0.004353  ;
syst_corr[3] =  0.002607  ; stat_corr[3] =  0.006642  ; stat_uncorr[3] =  0.004428  ;
syst_corr[4] =  0.005212  ; stat_corr[4] =  0.006648  ; stat_uncorr[4] =  0.004907  ;
syst_corr[5] =  0.005969  ; stat_corr[5] =  0.008155  ; stat_uncorr[5] =  0.006332  ;
      xmin=0.0;
      xmax=pi;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lep_costheta_cms";
      xaxislabel="cos(#theta_{l+})";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
syst_corr[0] =  0.026381  ; stat_corr[0] =  0.024436  ; stat_uncorr[0] =  0.017748  ;
syst_corr[1] =  0.020494  ; stat_corr[1] =  0.016636  ; stat_uncorr[1] =  0.012111  ;
syst_corr[2] =  0.012932  ; stat_corr[2] =  0.018346  ; stat_uncorr[2] =  0.012860  ;
syst_corr[3] =  0.009955  ; stat_corr[3] =  0.016185  ; stat_uncorr[3] =  0.011730  ;
syst_corr[4] =  0.019450  ; stat_corr[4] =  0.014609  ; stat_uncorr[4] =  0.010715  ;
syst_corr[5] =  0.028583  ; stat_corr[5] =  0.024646  ; stat_uncorr[5] =  0.016832  ;
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos(#theta_{l-})";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
syst_corr[0] =  0.023343  ; stat_corr[0] =  0.024218  ; stat_uncorr[0] =  0.017780  ;
syst_corr[1] =  0.008458  ; stat_corr[1] =  0.015829  ; stat_uncorr[1] =  0.012130  ;
syst_corr[2] =  0.005024  ; stat_corr[2] =  0.017465  ; stat_uncorr[2] =  0.012852  ;
syst_corr[3] =  0.009367  ; stat_corr[3] =  0.015869  ; stat_uncorr[3] =  0.011715  ;
syst_corr[4] =  0.012662  ; stat_corr[4] =  0.013843  ; stat_uncorr[4] =  0.010866  ;
syst_corr[5] =  0.014220  ; stat_corr[5] =  0.021370  ; stat_uncorr[5] =  0.017103  ;
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lep_costheta_cms";
      xaxislabel="cos(#theta_{l})";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
syst_corr[0] =  0.024398  ; stat_corr[0] =  0.016913  ; stat_uncorr[0] =  0.012560  ;
syst_corr[1] =  0.013633  ; stat_corr[1] =  0.012074  ; stat_uncorr[1] =  0.008569  ;
syst_corr[2] =  0.007663  ; stat_corr[2] =  0.012027  ; stat_uncorr[2] =  0.009090  ;
syst_corr[3] =  0.008886  ; stat_corr[3] =  0.010142  ; stat_uncorr[3] =  0.008288  ;
syst_corr[4] =  0.015734  ; stat_corr[4] =  0.010730  ; stat_uncorr[4] =  0.007630  ;
syst_corr[5] =  0.020589  ; stat_corr[5] =  0.017201  ; stat_uncorr[5] =  0.011998  ;
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos(#theta_{l+})cos(#theta_{l-})";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      xbins1D[0]=-1.0; xbins1D[1]=-0.5; xbins1D[2]=-0.2; xbins1D[3]=0.0; xbins1D[4]=0.2; xbins1D[5]=0.5; xbins1D[6]=1.0;
syst_corr[0] =  0.008637  ; stat_corr[0] =  0.014863  ; stat_uncorr[0] =  0.010240  ;
syst_corr[1] =  0.017095  ; stat_corr[1] =  0.030016  ; stat_uncorr[1] =  0.020079  ;
syst_corr[2] =  0.023907  ; stat_corr[2] =  0.043863  ; stat_uncorr[2] =  0.031697  ;
syst_corr[3] =  0.023001  ; stat_corr[3] =  0.042290  ; stat_uncorr[3] =  0.028571  ;
syst_corr[4] =  0.018808  ; stat_corr[4] =  0.028738  ; stat_uncorr[4] =  0.019560  ;
syst_corr[5] =  0.006292  ; stat_corr[5] =  0.011391  ; stat_uncorr[5] =  0.007670  ;
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="|y_{top}|-|y_{tbar}|";
      acceptanceName="rapiditydiffMarco";
      asymlabel="A_{C}";
      xbins1D[0]=-2.0; xbins1D[1]=-0.7; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.7; xbins1D[6]=2.0;
syst_corr[0] =  0.004265  ; stat_corr[0] =  0.006765  ; stat_uncorr[0] =  0.004693  ;
syst_corr[1] =  0.004710  ; stat_corr[1] =  0.011998  ; stat_uncorr[1] =  0.008558  ;
syst_corr[2] =  0.008643  ; stat_corr[2] =  0.014062  ; stat_uncorr[2] =  0.010686  ;
syst_corr[3] =  0.013091  ; stat_corr[3] =  0.014289  ; stat_uncorr[3] =  0.010699  ;
syst_corr[4] =  0.004617  ; stat_corr[4] =  0.011317  ; stat_uncorr[4] =  0.008462  ;
syst_corr[5] =  0.002602  ; stat_corr[5] =  0.006236  ; stat_uncorr[5] =  0.004700  ;
      xmin=-2.0;
      xmax= 2.0;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="cos(#theta_{top})";
      acceptanceName="topCosTheta";
      asymlabel="A_{FB}";
      xbins1D[0]=-1.0; xbins1D[1]=-0.7; xbins1D[2]=-0.4; xbins1D[3]=0.0; xbins1D[4]=0.4; xbins1D[5]=0.7; xbins1D[6]=1.0;
syst_corr[0] =  0.016123  ; stat_corr[0] =  0.033209  ; stat_uncorr[0] =  0.024855  ;
syst_corr[1] =  0.009248  ; stat_corr[1] =  0.013867  ; stat_uncorr[1] =  0.010505  ;
syst_corr[2] =  0.005435  ; stat_corr[2] =  0.010515  ; stat_uncorr[2] =  0.008204  ;
syst_corr[3] =  0.003678  ; stat_corr[3] =  0.010134  ; stat_uncorr[3] =  0.008252  ;
syst_corr[4] =  0.007670  ; stat_corr[4] =  0.012908  ; stat_uncorr[4] =  0.010373  ;
syst_corr[5] =  0.013827  ; stat_corr[5] =  0.032788  ; stat_uncorr[5] =  0.024581  ;
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 8:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="cos(#Delta#phi_{l+l-})";
      acceptanceName="lepAzimAsym";
      asymlabel="A_{#Delta#phi}";
      xbins1D[0]=-1.0; xbins1D[1]=-0.8; xbins1D[2]=-0.4; xbins1D[3]=0.0; xbins1D[4]=0.4; xbins1D[5]=0.8; xbins1D[6]=1.0;
      stat_corr[0] = 0.01; stat_corr [1] = 0.01;  stat_corr [2] = 0.01;  stat_corr [3] = 0.01; stat_corr [4] = 0.01;  stat_corr [5] = 0.01;  stat_corr [6] = 0.01; 
      stat_uncorr[0] = 0.00; stat_uncorr[1] = 0.00; stat_uncorr[2] = 0.00; stat_uncorr[3] = 0.00; stat_uncorr[4] = 0.00; stat_uncorr[5] = 0.00; stat_uncorr[6] = 0.00;
      syst_corr[0] = 0.02; syst_corr [1] = 0.02;  syst_corr [2] = 0.02;  syst_corr [3] = 0.02; syst_corr [4] = 0.02;  syst_corr [5] = 0.02;  syst_corr [6] = 0.02; 
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Asy I
    case 9:
    {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="|#eta_{top}|-|#eta_{tbar}|";
      acceptanceName="pseudorapiditydiff";
      asymlabel="A_{C}";
      xbins1D[0]=-4.0; xbins1D[1]=-1.0; xbins1D[2]=-0.5; xbins1D[3]=0.0; xbins1D[4]=0.5; xbins1D[5]=1.0; xbins1D[6]=4.0;
      stat_corr[0] = 0.01; stat_corr [1] = 0.01;  stat_corr [2] = 0.01;  stat_corr [3] = 0.01; stat_corr [4] = 0.01;  stat_corr [5] = 0.01;  stat_corr [6] = 0.01; 
      stat_uncorr[0] = 0.00; stat_uncorr[1] = 0.00; stat_uncorr[2] = 0.00; stat_uncorr[3] = 0.00; stat_uncorr[4] = 0.00; stat_uncorr[5] = 0.00; stat_uncorr[6] = 0.00;
      syst_corr[0] = 0.02; syst_corr [1] = 0.02;  syst_corr [2] = 0.02;  syst_corr [3] = 0.02; syst_corr [4] = 0.02;  syst_corr [5] = 0.02;  syst_corr [6] = 0.02; 
      xmin=-4.0;
      xmax= 4.0;
      break;
    }
    //   Top Asy II
    case 10:
    {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="(y_{top}-y_{tbar})(y_{top}+y_{tbar})";
      acceptanceName="rapiditydiff";
      asymlabel="A_{C}";
      xbins1D[0]=-4.0; xbins1D[1]=-0.8; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.8; xbins1D[6]=4.0;
      stat_corr[0] = 0.01; stat_corr [1] = 0.01;  stat_corr [2] = 0.01;  stat_corr [3] = 0.01; stat_corr [4] = 0.01;  stat_corr [5] = 0.01;  stat_corr [6] = 0.01; 
      stat_uncorr[0] = 0.00; stat_uncorr[1] = 0.00; stat_uncorr[2] = 0.00; stat_uncorr[3] = 0.00; stat_uncorr[4] = 0.00; stat_uncorr[5] = 0.00; stat_uncorr[6] = 0.00;
      syst_corr[0] = 0.02; syst_corr [1] = 0.02;  syst_corr [2] = 0.02;  syst_corr [3] = 0.02; syst_corr [4] = 0.02;  syst_corr [5] = 0.02;  syst_corr [6] = 0.02; 
      xmin=-4.0;
      xmax= 4.0;
      break;
    }
    default:
    {
      cout<<"Set the variable switch";
    }
  }
}

void Initialize2DBinning(int iVar){


  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepChargeAsym";
      asymlabel="A_{lepC}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.007318  ; stat_corr[0] =  0.008170  ; stat_uncorr[0] =  0.006160  ;
syst_corr[1] =  0.007694  ; stat_corr[1] =  0.017318  ; stat_uncorr[1] =  0.012322  ;
syst_corr[2] =  0.005939  ; stat_corr[2] =  0.024728  ; stat_uncorr[2] =  0.016786  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepAzimAsym2";
      asymlabel="A_{#Delta#phi}";
      Double_t pi = 3.141592653589793;
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.023651  ; stat_corr[0] =  0.009848  ; stat_uncorr[0] =  0.006565  ;
syst_corr[1] =  0.025452  ; stat_corr[1] =  0.019887  ; stat_uncorr[1] =  0.013103  ;
syst_corr[2] =  0.022379  ; stat_corr[2] =  0.024235  ; stat_uncorr[2] =  0.016475  ;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lep_costheta_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.014692  ; stat_corr[0] =  0.011057  ; stat_uncorr[0] =  0.007925  ;
syst_corr[1] =  0.048445  ; stat_corr[1] =  0.027972  ; stat_uncorr[1] =  0.019954  ;
syst_corr[2] =  0.068124  ; stat_corr[2] =  0.037687  ; stat_uncorr[2] =  0.026584  ;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.010814  ; stat_corr[0] =  0.010995  ; stat_uncorr[0] =  0.007913  ;
syst_corr[1] =  0.031719  ; stat_corr[1] =  0.027755  ; stat_uncorr[1] =  0.020015  ;
syst_corr[2] =  0.050730  ; stat_corr[2] =  0.037187  ; stat_uncorr[2] =  0.026815  ;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lep_costheta_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.012151  ; stat_corr[0] =  0.008145  ; stat_uncorr[0] =  0.005600  ;
syst_corr[1] =  0.039024  ; stat_corr[1] =  0.020343  ; stat_uncorr[1] =  0.014132  ;
syst_corr[2] =  0.058320  ; stat_corr[2] =  0.026816  ; stat_uncorr[2] =  0.018878  ;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="M_{t#bar t}";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.015401  ; stat_corr[0] =  0.014180  ; stat_uncorr[0] =  0.010028  ;
syst_corr[1] =  0.030914  ; stat_corr[1] =  0.037927  ; stat_uncorr[1] =  0.026826  ;
syst_corr[2] =  0.034871  ; stat_corr[2] =  0.050649  ; stat_uncorr[2] =  0.035982  ;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="M_{t#bar t}";
      acceptanceName="rapiditydiffMarco";
      asymlabel="A_{C}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.003017  ; stat_corr[0] =  0.009781  ; stat_uncorr[0] =  0.007485  ;
syst_corr[1] =  0.008211  ; stat_corr[1] =  0.024959  ; stat_uncorr[1] =  0.019329  ;
syst_corr[2] =  0.009798  ; stat_corr[2] =  0.033379  ; stat_uncorr[2] =  0.025751  ;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="topCosTheta";
      asymlabel="A_{FB}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.003138  ; stat_corr[0] =  0.010888  ; stat_uncorr[0] =  0.007558  ;
syst_corr[1] =  0.007280  ; stat_corr[1] =  0.027899  ; stat_uncorr[1] =  0.019511  ;
syst_corr[2] =  0.010199  ; stat_corr[2] =  0.036874  ; stat_uncorr[2] =  0.025993  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 8:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepAzimAsym";
      asymlabel="A_{#Delta#phi}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy I
    case 9:
    {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="pseudorapiditydiff";
      asymlabel="A_{C}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy II
    case 10:
    {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="rapiditydiff";
      asymlabel="A_{C}";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    default:
    {
      cout<<"Set the variable switch";
    }
  }
}

void Initialize2DBinningttpt(int iVar){


  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepChargeAsym";
      asymlabel="A_{lepC}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.005087  ; stat_corr[0] =  0.008183  ; stat_uncorr[0] =  0.005284  ;
syst_corr[1] =  0.005429  ; stat_corr[1] =  0.018370  ; stat_uncorr[1] =  0.011753  ;
syst_corr[2] =  0.011115  ; stat_corr[2] =  0.025259  ; stat_uncorr[2] =  0.016731  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepAzimAsym2";
      asymlabel="A_{#Delta#phi}";
      Double_t pi = 3.141592653589793;
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.007830  ; stat_corr[0] =  0.007155  ; stat_uncorr[0] =  0.005245  ;
syst_corr[1] =  0.014203  ; stat_corr[1] =  0.017064  ; stat_uncorr[1] =  0.011774  ;
syst_corr[2] =  0.028484  ; stat_corr[2] =  0.026298  ; stat_uncorr[2] =  0.016888  ;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lep_costheta_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.016340  ; stat_corr[0] =  0.009820  ; stat_uncorr[0] =  0.007676  ;
syst_corr[1] =  0.036816  ; stat_corr[1] =  0.027043  ; stat_uncorr[1] =  0.020872  ;
syst_corr[2] =  0.049544  ; stat_corr[2] =  0.039271  ; stat_uncorr[2] =  0.029880  ;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.013058  ; stat_corr[0] =  0.010809  ; stat_uncorr[0] =  0.007705  ;
syst_corr[1] =  0.023958  ; stat_corr[1] =  0.029437  ; stat_uncorr[1] =  0.020966  ;
syst_corr[2] =  0.032052  ; stat_corr[2] =  0.042381  ; stat_uncorr[2] =  0.030003  ;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lep_costheta_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.013735  ; stat_corr[0] =  0.007424  ; stat_uncorr[0] =  0.005438  ;
syst_corr[1] =  0.027114  ; stat_corr[1] =  0.020123  ; stat_uncorr[1] =  0.014792  ;
syst_corr[2] =  0.036131  ; stat_corr[2] =  0.028688  ; stat_uncorr[2] =  0.021172  ;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.007019  ; stat_corr[0] =  0.012939  ; stat_uncorr[0] =  0.009715  ;
syst_corr[1] =  0.020415  ; stat_corr[1] =  0.036131  ; stat_uncorr[1] =  0.027047  ;
syst_corr[2] =  0.025012  ; stat_corr[2] =  0.052149  ; stat_uncorr[2] =  0.038856  ;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="rapiditydiffMarco";
      asymlabel="A_{C}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.003071  ; stat_corr[0] =  0.010204  ; stat_uncorr[0] =  0.007877  ;
syst_corr[1] =  0.011181  ; stat_corr[1] =  0.028361  ; stat_uncorr[1] =  0.021528  ;
syst_corr[2] =  0.010222  ; stat_corr[2] =  0.041218  ; stat_uncorr[2] =  0.030818  ;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="topCosTheta";
      asymlabel="A_{FB}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.004387  ; stat_corr[0] =  0.010882  ; stat_uncorr[0] =  0.007908  ;
syst_corr[1] =  0.013393  ; stat_corr[1] =  0.029925  ; stat_uncorr[1] =  0.021610  ;
syst_corr[2] =  0.014671  ; stat_corr[2] =  0.043021  ; stat_uncorr[2] =  0.030935  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 8:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepAzimAsym";
      asymlabel="A_{#Delta#phi}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy I
    case 9:
    {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="pseudorapiditydiff";
      asymlabel="A_{C}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy II
    case 10:
    {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="rapiditydiff";
      asymlabel="A_{C}";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    default:
    {
      cout<<"Set the variable switch";
    }
  }
}

void Initialize2DBinningttrapidity2(int iVar){


  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepChargeAsym";
      asymlabel="A_{lepC}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.003463  ; stat_corr[0] =  0.007281  ; stat_uncorr[0] =  0.005302  ;
syst_corr[1] =  0.005698  ; stat_corr[1] =  0.015993  ; stat_uncorr[1] =  0.011634  ;
syst_corr[2] =  0.013340  ; stat_corr[2] =  0.021568  ; stat_uncorr[2] =  0.016385  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepAzimAsym2";
      asymlabel="A_{#Delta#phi}";
      Double_t pi = 3.141592653589793;
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.009670  ; stat_corr[0] =  0.006995  ; stat_uncorr[0] =  0.005269  ;
syst_corr[1] =  0.013223  ; stat_corr[1] =  0.014989  ; stat_uncorr[1] =  0.011671  ;
syst_corr[2] =  0.025492  ; stat_corr[2] =  0.020775  ; stat_uncorr[2] =  0.016472  ;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lep_costheta_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.015142  ; stat_corr[0] =  0.010353  ; stat_uncorr[0] =  0.007570  ;
syst_corr[1] =  0.034480  ; stat_corr[1] =  0.027519  ; stat_uncorr[1] =  0.020384  ;
syst_corr[2] =  0.048856  ; stat_corr[2] =  0.038596  ; stat_uncorr[2] =  0.028917  ;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.011975  ; stat_corr[0] =  0.010377  ; stat_uncorr[0] =  0.007601  ;
syst_corr[1] =  0.024868  ; stat_corr[1] =  0.028769  ; stat_uncorr[1] =  0.020488  ;
syst_corr[2] =  0.038707  ; stat_corr[2] =  0.041452  ; stat_uncorr[2] =  0.029096  ;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lep_costheta_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.012552  ; stat_corr[0] =  0.007371  ; stat_uncorr[0] =  0.005364  ;
syst_corr[1] =  0.028702  ; stat_corr[1] =  0.019493  ; stat_uncorr[1] =  0.014450  ;
syst_corr[2] =  0.043069  ; stat_corr[2] =  0.027454  ; stat_uncorr[2] =  0.020509  ;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.012554  ; stat_corr[0] =  0.011929  ; stat_uncorr[0] =  0.009414  ;
syst_corr[1] =  0.031478  ; stat_corr[1] =  0.033629  ; stat_uncorr[1] =  0.026276  ;
syst_corr[2] =  0.037265  ; stat_corr[2] =  0.048500  ; stat_uncorr[2] =  0.037520  ;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="rapiditydiffMarco";
      asymlabel="A_{C}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.007655  ; stat_corr[0] =  0.008960  ; stat_uncorr[0] =  0.006900  ;
syst_corr[1] =  0.010556  ; stat_corr[1] =  0.024765  ; stat_uncorr[1] =  0.018994  ;
syst_corr[2] =  0.013484  ; stat_corr[2] =  0.035883  ; stat_uncorr[2] =  0.026955  ;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="topCosTheta";
      asymlabel="A_{FB}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.007425  ; stat_corr[0] =  0.009803  ; stat_uncorr[0] =  0.006997  ;
syst_corr[1] =  0.011807  ; stat_corr[1] =  0.027521  ; stat_uncorr[1] =  0.019233  ;
syst_corr[2] =  0.014130  ; stat_corr[2] =  0.039478  ; stat_uncorr[2] =  0.027300  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 8:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepAzimAsym";
      asymlabel="A_{#Delta#phi}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy I
    case 9:
    {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="pseudorapiditydiff";
      asymlabel="A_{C}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy II
    case 10:
    {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="rapiditydiff";
      asymlabel="A_{C}";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    default:
    {
      cout<<"Set the variable switch";
    }
  }
}

void fillUnderOverFlow(TH1D *h1, float value, double weight, double Nsolns)
{
  double min = h1->GetXaxis()->GetXmin();
  double max = h1->GetXaxis()->GetXmax();

  if (value >= max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value <= min) value = h1->GetBinCenter(1);

  int bin_number = h1->FindBin(value);
  double orig_content = h1->GetBinContent(bin_number);
  double orig_error = h1->GetBinError(bin_number);

  //h1->Fill(value, weight);
  h1->SetBinContent( bin_number, orig_content+weight );
  h1->SetBinError( bin_number, sqrt( orig_error*orig_error + weight*weight*Nsolns ) );
}

//--------------------------------------------------------------------

void fillUnderOverFlow(TH2D *h2, float xvalue, float yvalue, double weight, double Nsolns)
{
  double maxx = h2->GetXaxis()->GetXmax();
  double minx = h2->GetXaxis()->GetXmin();
  double maxy = h2->GetYaxis()->GetXmax();
  double miny = h2->GetYaxis()->GetXmin();

  if (xvalue >= maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
  if (xvalue <= minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
  if (yvalue >= maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
  if (yvalue <= miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

  int bin_number = h2->FindBin(xvalue,yvalue);
  double orig_content = h2->GetBinContent(bin_number);
  double orig_error = h2->GetBinError(bin_number);

  //h2->Fill(xvalue, yvalue, weight);
  h2->SetBinContent( bin_number, orig_content+weight );
  h2->SetBinError( bin_number, sqrt( orig_error*orig_error + weight*weight*Nsolns ) );
}
