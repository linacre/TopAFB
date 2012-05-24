#include "TH1D.h"
#include "TMatrixD.h"
#include "TMath.h"

double xsection = 154.0;

int lineWidth=3;

TString observablename;
TString acceptanceName;
TString xaxislabel;

float xmin=-1.0;
float xmax= 1.0;

const int nbins1D=6;
const int nbins2D=6;

double xbins1D[nbins1D+1];
double xbins2D[nbins2D+1];


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

  event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_total = event_plus + event_minus;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
		     + event_minus*event_minus*event_plus_err*event_plus_err)/
		  (event_total*event_total*event_total*event_total));

}



void GetAvsY(TH1* histogram, TMatrixD &covarianceM, vector<double> &myafb, vector<double> &myerr){

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

void Initialize1DBinning(int iVar){


  switch (iVar)
    {
      //   Lepton Charge Asymmetry
    case 0:
      {
      observablename="lep_charge_asymmetry";
      xaxislabel="|#eta_{l+}|-|#eta_{l-}|";
      acceptanceName="lepChargeAsym";
      xbins1D[0]=-3.0; xbins1D[1]=-2.0; xbins1D[2]=-1.0; xbins1D[3]=0.0; xbins1D[4]=1.0; xbins1D[5]=2.0; xbins1D[6]=3.0;
      xmin=-3.0;
      xmax= 3.0;
      break;
      }
  //   Lepton Azimuthal Asymmetry
    case 1:
      {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="cos(#Delta#phi_{l+l-})";
      acceptanceName="lepAzimAsym";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      xmin=-1.0;
      xmax= 1.0;
      break;
      }
  //   Top Charge Asymmetry
    case 2:
      {
      observablename="top_costheta_cms";
      xaxislabel="cos(#theta_{top})";
      acceptanceName="topCosTheta";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      xmin=-1.0;
      xmax= 1.0;
      break;
      }

  //   Top Polarization
    case 3:
      {
      observablename="lep_costheta_cms";
      xaxislabel="cos(#theta_{l,n})";
      acceptanceName="lepPlusCosTheta";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      xmin=-1.0;
      xmax= 1.0;
      break;
      }
  //   Top Spin Correlation
    case 4:
      {
      observablename="top_spin_correlation";
      xaxislabel="cos(#theta_{l+,n})cos(#theta_{l-,n})";
      acceptanceName="topSpinCorr";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      xmin=-1.0;
      xmax= 1.0;
      break;
      }
 //   Top Asy I
    case 5:
      {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="|#eta_{top}|-|#eta_{tbar}|";
      acceptanceName="pseudorapiditydiff";
      xbins1D[0]=-6.0; xbins1D[1]=-2.0; xbins1D[2]=-1.0; xbins1D[3]=0.0; xbins1D[4]=1.0; xbins1D[5]=2.0; xbins1D[6]=6.0;
      xmin=-5.0;
      xmax= 5.0;
      break;
      }
      //   Top Asy II
    case 6:
      {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="(y_{top}-y_{tbar})(y_{top}+y_{tbar})";
      acceptanceName="rapiditydiff";
      xbins1D[0]=-6.0; xbins1D[1]=-2.0; xbins1D[2]=-1.0; xbins1D[3]=0.0; xbins1D[4]=1.0; xbins1D[5]=2.0; xbins1D[6]=6.0;
      xmin=-6.0;
      xmax= 6.0;
      break;
      }
      //   Top Asy III
    case 7:
      {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="|y_{top}|-|y_{tbar}|";
      acceptanceName="rapiditydiffMarco";
      xbins1D[0]=-3.0; xbins1D[1]=-2.0; xbins1D[2]=-1.0; xbins1D[3]=0.0; xbins1D[4]=1.0; xbins1D[5]=2.0; xbins1D[6]=3.0;
      xmin=-3.0;
      xmax= 3.0;
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
   //   Top Charge Asymmetry
    case 0:
      {
      observablename="top_costheta_cms";
      xaxislabel="A_{topFB}"; 
      acceptanceName="topCosTheta";
      xbins2D[0]=-1500.0; xbins2D[1]=-550.0; xbins2D[2]=-450.0; xbins2D[3]=0.0; xbins2D[4]=450; xbins2D[5]=550.0; xbins2D[6]=1500.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
      }

  //   Top Polarization
    case 1:
      {
      observablename="lep_costheta_cms";
      xaxislabel="P_{n}";
      acceptanceName="lepPlusCosTheta";
      xbins2D[0]=-1500.0; xbins2D[1]=-550.0; xbins2D[2]=-450.0; xbins2D[3]=0.0; xbins2D[4]=450; xbins2D[5]=550.0; xbins2D[6]=1500.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
      }
  //   Top Spin Correlation
    case 2:
      {
      observablename="top_spin_correlation";
      xaxislabel="A_{c1c2}";
      acceptanceName="topSpinCorr";
      xbins2D[0]=-1500.0; xbins2D[1]=-550.0; xbins2D[2]=-450.0; xbins2D[3]=0.0; xbins2D[4]=450; xbins2D[5]=550.0; xbins2D[6]=1500.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
      }
 //   Top Asy I
    case 3:
      {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="A_{I}";
      acceptanceName="pseudorapiditydiff";
      xbins2D[0]=-1500.0; xbins2D[1]=-550.0; xbins2D[2]=-450.0; xbins2D[3]=0.0; xbins2D[4]=450; xbins2D[5]=550.0; xbins2D[6]=1500.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
      }
      //   Top Asy I
    case 4:
      {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="A_{II}";
      acceptanceName="rapiditydiff";
      xbins2D[0]=-1500.0; xbins2D[1]=-550.0; xbins2D[2]=-450.0; xbins2D[3]=0.0; xbins2D[4]=450; xbins2D[5]=550.0; xbins2D[6]=1500.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
      }
      // Top Asy Marco
    case 5:
      {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="A_{III}";
      acceptanceName="rapiditydiffMarco";
      xbins2D[0]=-1500.0; xbins2D[1]=-550.0; xbins2D[2]=-450.0; xbins2D[3]=0.0; xbins2D[4]=450; xbins2D[5]=550.0; xbins2D[6]=1500.0;
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
