#include "TH1D.h"
#include "TMatrixD.h"
#include "TMath.h"

double xsection = 154.0*0.06451;

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




void GetCorrectedAfbBinByBin(TH1D* histogram, TMatrixD &covarianceM, vector<double> &myafb, vector<double> &myerr){

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
xaxislabel="cos(#theta^{+}_{l})";
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
        //   Top Polarization using negatively charged leptons
    case 8:
      {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos(#theta^{-}_{l})";
      acceptanceName="lepMinusCosTheta";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      xmin=-1.0;
      xmax= 1.0;
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
      xaxislabel="cos(#theta^{+}_{l})";
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
  //   Top Polarization using negatively charged leptons
    case 6:
      {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos(#theta^{-}_{l})";
      acceptanceName="lepMinusCosTheta";
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


void fillUnderOverFlow(TH1D *h1, float value, double weight, int Nsolns)
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
  h1->SetBinError( bin_number, sqrt( orig_error*orig_error + weight*weight*double(Nsolns) ) );
}

//--------------------------------------------------------------------

void fillUnderOverFlow(TH2D *h2, float xvalue, float yvalue, double weight, int Nsolns)
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
  h2->SetBinError( bin_number, sqrt( orig_error*orig_error + weight*weight*double(Nsolns) ) );
}
