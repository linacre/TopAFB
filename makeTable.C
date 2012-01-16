#include "CommonFunctions.C"
#include " TH1F.h"
void calculateA( TCut baseline, const char* var, bool isGen, float& event_total, float& event_plus, float& event_minus, float& event_corr,  float& asymm, float& asymm_err ){

  const char* temp = Form("(%s*%s_gen) >0", var, var);
  TCut corrSign = temp;
  cout <<temp <<endl;
  TCut plusSign  ;
  TCut minusSign ;
  if(!isGen){
    plusSign = Form("%s >0", var);
    minusSign = Form("%s <0", var);
 
  }
  else{
    plusSign = Form("%s_gen >0", var);
    minusSign = Form("%s_gen <0", var);
  }
 
  TH1F* hlepChargeA = new TH1F("lepChargeA", "lepChargeA", 4, -10, 10);
  TH1F* hlepChargeA_plus = new TH1F("lepChargeA_plus", "lepChargeA_plus", 4, -10, 10);
  TH1F* hlepChargeA_minus = new TH1F("lepChargeA_minus", "lepChargeA_minus", 4, -10, 10);
  TH1F* hlepChargeA_corr = new TH1F("lepChargeA_corr", "lepChargeA_corr", 4, -10, 10);
  tree->Draw(Form("%s >> %s", var, "lepChargeA"),       baseline);
  tree->Draw(Form("%s >> %s", var, "lepChargeA_corr"),  baseline && corrSign);
  
  if(!isGen){
    tree->Draw(Form("%s >> %s", var, "lepChargeA_plus") , baseline && plusSign);
    tree->Draw(Form("%s >> %s", var, "lepChargeA_minus"), baseline && minusSign);
  }
  else{
    tree->Draw(Form("%s >> %s", var, "lepChargeA_plus") , baseline && plusSign);
    tree->Draw(Form("%s >> %s", var, "lepChargeA_minus"), baseline && minusSign);
    
  }

  int nbins = lepChargeA->GetNbinsX();
 
  event_total= lepChargeA->Integral(0, nbins+1);
  event_plus  = lepChargeA_plus->Integral(0, nbins+1);
  event_minus = lepChargeA_minus->Integral(0, nbins+1);
  event_corr  = (lepChargeA_corr->Integral(0, nbins+1))/event_total;
  asymm = (event_plus-event_minus)/(event_plus+event_minus);
  asymm_err   = sqrt(4*(event_plus*event_minus)/(event_total*event_total*event_total));
  delete  hlepChargeA ;
  delete  hlepChargeA_plus;
  delete  hlepChargeA_minus;
  delete  hlepChargeA_corr;
  
}

void makeTable(bool isData=false){

  bool latex=true;
  std::string pmSign  = latex ? " $\\pm$ " : " &plusmn; ";
  std::string colSep  = latex ? " & " : " | ";
  std::string beginL  = latex ? ""   : "| ";
  std::string endL    = latex ? " \\\\ " : " | ";
  std::string mathSep = latex ? "$" : "";
  
  const char* formatS = "%6.3f";
  TCut baseline = "tt_mass <450 && ttRapidity >2.0";
  
  float lepChargeA    ;
  float lepChargeA_plus;
  float lepChargeA_minus;
  float lepChargeA_corr ;
  float lepChargeA_quant ;
  float lepChargeA_err  ;
  float lepChargeA_plus_gen ;
  float lepChargeA_minus_gen ;
  float lepChargeA_quant_gen ;
  float lepChargeA_err_gen   ;
  
  calculateA( baseline, "lep_charge_asymmetry", false, lepChargeA, lepChargeA_plus, lepChargeA_minus, lepChargeA_corr,  lepChargeA_quant, lepChargeA_err );
  calculateA( baseline, "lep_charge_asymmetry", true, lepChargeA, lepChargeA_plus_gen, lepChargeA_minus_gen, lepChargeA_corr,  lepChargeA_quant_gen, lepChargeA_err_gen );

  float topCosThetaA    ;
  float topCosThetaA_plus;
  float topCosThetaA_minus;
  float topCosThetaA_corr ;
  float topCosThetaA_quant ;
  float topCosThetaA_err  ;
  float topCosThetaA_plus_gen ;
  float topCosThetaA_minus_gen ;
  float topCosThetaA_quant_gen ;
  float topCosThetaA_err_gen   ;
  
  calculateA( baseline, "top_costheta_cms", false, topCosThetaA, topCosThetaA_plus, topCosThetaA_minus, topCosThetaA_corr,  topCosThetaA_quant, topCosThetaA_err );
  calculateA( baseline, "top_costheta_cms", true, topCosThetaA, topCosThetaA_plus_gen, topCosThetaA_minus_gen, topCosThetaA_corr,  topCosThetaA_quant_gen, topCosThetaA_err_gen );
  


  float lepCosThetaA    ;
  float lepCosThetaA_plus;
  float lepCosThetaA_minus;
  float lepCosThetaA_corr ;
  float lepCosThetaA_quant ;
  float lepCosThetaA_err  ;
  float lepCosThetaA_plus_gen ;
  float lepCosThetaA_minus_gen ;
  float lepCosThetaA_quant_gen ;
  float lepCosThetaA_err_gen   ;
  
  calculateA( baseline, "lep_costheta_cms", false, lepCosThetaA, lepCosThetaA_plus, lepCosThetaA_minus, lepCosThetaA_corr,  lepCosThetaA_quant, lepCosThetaA_err );
  calculateA( baseline, "lep_costheta_cms", true, lepCosThetaA, lepCosThetaA_plus_gen, lepCosThetaA_minus_gen, lepCosThetaA_corr,  lepCosThetaA_quant_gen, lepCosThetaA_err_gen );
  
  
  float topSpinCorrA    ;
  float topSpinCorrA_plus;
  float topSpinCorrA_minus;
  float topSpinCorrA_corr ;
  float topSpinCorrA_quant ;
  float topSpinCorrA_err  ;
  float topSpinCorrA_plus_gen ;
  float topSpinCorrA_minus_gen ;
  float topSpinCorrA_quant_gen ;
  float topSpinCorrA_err_gen   ;
  
  calculateA( baseline, "top_spin_correlation", false, topSpinCorrA, topSpinCorrA_plus, topSpinCorrA_minus, topSpinCorrA_corr,  topSpinCorrA_quant, topSpinCorrA_err );
  calculateA( baseline, "top_spin_correlation", true, topSpinCorrA, topSpinCorrA_plus_gen, topSpinCorrA_minus_gen, topSpinCorrA_corr,  topSpinCorrA_quant_gen, topSpinCorrA_err_gen );
  

  if(!isData){

  cout << "\\begin{tabular}{l| c  c  c  c c c}" << endl;
  cout << "\\hline" << endl;
  cout << beginL <<"Var "<< colSep << "     "  << colSep<< "Total Events"  << colSep  <<"Plus Sign"        << colSep <<"Minus Sign "     << colSep << "Prob of Corr Sign"  << colSep << "Asym" << endL <<endl;
  cout << beginL <<"Lep Charge Asy"<<colSep << " Gen  "     << colSep<< lepChargeA      << colSep  <<lepChargeA_plus_gen    << colSep <<lepChargeA_minus_gen  << colSep << "-"  << colSep <<formatFloat(lepChargeA_quant_gen, formatS) <<pmSign << formatFloat(lepChargeA_err_gen,formatS)<< endL <<endl;
  cout << beginL <<"Lep Charge Asy"<<colSep << " Reco "     << colSep << lepChargeA      << colSep  <<lepChargeA_plus    << colSep <<lepChargeA_minus  << colSep <<formatFloat(lepChargeA_corr,formatS) << colSep <<formatFloat(lepChargeA_quant,formatS) <<pmSign << formatFloat(lepChargeA_err, formatS) << endL <<endl;

  cout << beginL <<"Top AFB"<<colSep << " Gen  "     << colSep<< topCosThetaA      << colSep  <<topCosThetaA_plus_gen    << colSep <<topCosThetaA_minus_gen  << colSep << "-"  << colSep <<formatFloat(topCosThetaA_quant_gen, formatS) <<pmSign << formatFloat(topCosThetaA_err_gen,formatS)<< endL <<endl;
  cout << beginL <<"Top AFB"<<colSep << " Reco "     << colSep << topCosThetaA      << colSep  <<topCosThetaA_plus    << colSep <<topCosThetaA_minus  << colSep <<formatFloat(topCosThetaA_corr,formatS) << colSep <<formatFloat(topCosThetaA_quant,formatS) <<pmSign << formatFloat(topCosThetaA_err, formatS) << endL <<endl;

  cout << beginL <<"Top Polorization"<<colSep << " Gen  "     << colSep<< lepCosThetaA      << colSep  <<lepCosThetaA_plus_gen    << colSep <<lepCosThetaA_minus_gen  << colSep << "-"  << colSep <<formatFloat(lepCosThetaA_quant_gen, formatS) <<pmSign << formatFloat(lepCosThetaA_err_gen,formatS)<< endL <<endl;
  cout << beginL <<"Top Polorization"<<colSep << " Reco "     << colSep << lepCosThetaA      << colSep  <<lepCosThetaA_plus    << colSep <<lepCosThetaA_minus  << colSep <<formatFloat(lepCosThetaA_corr,formatS) << colSep <<formatFloat(lepCosThetaA_quant,formatS) <<pmSign << formatFloat(lepCosThetaA_err, formatS) << endL <<endl;

 cout << beginL <<"Top Spin Correlation"<<colSep << " Gen  "     << colSep<< topSpinCorrA      << colSep  <<topSpinCorrA_plus_gen    << colSep <<topSpinCorrA_minus_gen  << colSep << "-"  << colSep <<formatFloat(topSpinCorrA_quant_gen, formatS) <<pmSign << formatFloat(topSpinCorrA_err_gen,formatS)<< endL <<endl;
 cout << beginL <<"Top Spin Correlation"<<colSep << " Reco "     << colSep << topSpinCorrA      << colSep  <<topSpinCorrA_plus    << colSep <<topSpinCorrA_minus  << colSep <<formatFloat(topSpinCorrA_corr,formatS) << colSep <<formatFloat(topSpinCorrA_quant,formatS) <<pmSign << formatFloat(topSpinCorrA_err, formatS) << endL <<endl;



  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  }

  else{
    
  cout << "\\begin{tabular}{l| c  c  c  c c c}" << endl;
  cout << "\\hline" << endl;
  cout << beginL <<"Var "<< colSep<< "Total Events"  << colSep  <<"Plus Sign"        << colSep <<"Minus Sign "     << colSep << "Asym" << endL <<endl;
  cout << beginL <<"Lep Charge Asy" << colSep << lepChargeA      << colSep  <<lepChargeA_plus    << colSep <<lepChargeA_minus   << colSep <<formatFloat(lepChargeA_quant,formatS) <<pmSign << formatFloat(lepChargeA_err, formatS) << endL <<endl;

  cout << beginL <<"Top AFB"     << colSep << topCosThetaA      << colSep  <<topCosThetaA_plus    << colSep <<topCosThetaA_minus   << colSep <<formatFloat(topCosThetaA_quant,formatS) <<pmSign << formatFloat(topCosThetaA_err, formatS) << endL <<endl;


  cout << beginL <<"Top Polorization"   << colSep << lepCosThetaA      << colSep  <<lepCosThetaA_plus    << colSep <<lepCosThetaA_minus   << colSep <<formatFloat(lepCosThetaA_quant,formatS) <<pmSign << formatFloat(lepCosThetaA_err, formatS) << endL <<endl;


 cout << beginL <<"Top Spin Correlation"  << colSep << topSpinCorrA      << colSep  <<topSpinCorrA_plus    << colSep <<topSpinCorrA_minus   << colSep <<formatFloat(topSpinCorrA_quant,formatS) <<pmSign << formatFloat(topSpinCorrA_err, formatS) << endL <<endl;



  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
    
  }

}
