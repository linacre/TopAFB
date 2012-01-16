#include <iostream>
#include <vector>
#include "TStyle.h"
//#include "histtools.C"
#include "browseStacks.C"
#include "getMyHistosNames.C"
#include "TSystem.h"
#include "tdrStyle.C"
#include "DYest.C"
#include "CommonFunctions.C"
#include "EstimateSignalSpillage.C"

/*                                                                                                                                                                                                                                                            
  replacewFRest = replace the yields and errors for WJets, QCD with the FR estimate results                                                                                                                                                                   
  Does the Spillage subtraction as well                                                                                                                                                                                                                       
  replacewDYest = replace the yields and errors for DY with the data driven estimate                                                                                                                                                                          
  drawFullErrors = draw the errors in all their glory. Includes the systematic errors:                                                                                                                                                                        
  50% on the WJets estimate                                                                                                                                                                                                                                   
  100% on the QCD estimate                                                                                                                                                                                                                                    
  50% on DY estimate                                                                                                                                                                                                                                          
  50% on the backgrounds from MC                                                                                                                                                                                                                              
*/

//void makePSFile(const TString dataFName="results_data/hist_usePtGt2020_applyTriggers_hypDisamb_usejptJets_usetcMET_requireEcalEls_useOS_vetoHypMassLt12_require2BTag_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", const TString mcFName="results/hist_usePtGt2020_applyTriggers_hypDisamb_usejptJets_usetcMET_requireEcalEls_useOS_vetoHypMassLt12_require2BTag_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root",         float scaleMC=36.1, bool drawLogY = true, bool drawDiffs = false,

void makePSFile(const TString FName="results/hist_usePtGt2020_applyTriggers_hypDisamb_usepfMET_usepfJets_requireEcalEls_useOS_vetoHypMassLt12_require2BTag_sortJetCandidatesbyDR_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root",     
		int rebin=1, bool drawtprime = false,
		const TString pat="default", 
		float scaleMC=4.68, bool drawLogY = false, bool allhypOnly = false, bool drawDiffs = false,
                bool replacewFRest = false, bool replacewDYest = false,
                bool drawFullErrors = false, bool scalebTagMCbyDataDrivenEsts = false) {

  setTDRStyle();

  // hist::loadHist(FName.Data(),0,"*hnJet*");
  
  if(pat!="default") hist::loadHist(FName.Data(),0,pat.Data());
  else {
 //  hist::loadHist(FName.Data(),0,"*hnJet*");
//   hist::loadHist(FName.Data(),0,"*hnVtx*");
//   hist::loadHist(FName.Data(),0,"*hnBtagJet*");
//   hist::loadHist(FName.Data(),0,"*hmassltb_2j*");
//   hist::loadHist(FName.Data(),0,"*hmassllb_2j*");
//   hist::loadHist(FName.Data(),0,"*thefirstJetPt_2j*");
//   hist::loadHist(FName.Data(),0,"*thesecondJetPt_2j*");
//   hist::loadHist(FName.Data(),0,"*hjetPt_2j*");
//   hist::loadHist(FName.Data(),0,"*theSumJetPt_2j*");
//   hist::loadHist(FName.Data(),0,"*theSumBtagJetPt_2j*");
//   hist::loadHist(FName.Data(),0,"*hjetEta_2j*");
//   hist::loadHist(FName.Data(),0,"*htheleadinglepPt_2j*");
//   hist::loadHist(FName.Data(),0,"*hthesecondlepPt_2j*");
//   hist::loadHist(FName.Data(),0,"*theSumLepPt_2j*");
//   hist::loadHist(FName.Data(),0,"*hlepEta_2j*");
//   hist::loadHist(FName.Data(),0,"*hMET_2j*");
    hist::loadHist(FName.Data(),0,"*htopMass_2j*");
    hist::loadHist(FName.Data(),0,"*htopCosTheta_2j*");
    hist::loadHist(FName.Data(),0,"*hlepChargeAsym_2j*");
    hist::loadHist(FName.Data(),0,"*hlepCosTheta_2j*");
    hist::loadHist(FName.Data(),0,"*htopSpinCorr_2j*");
    hist::loadHist(FName.Data(),0,"*htopCosThetaGen_2j*");
    hist::loadHist(FName.Data(),0,"*hlepChargeAsymGen_2j*");
    hist::loadHist(FName.Data(),0,"*hlepCosThetaGen_2j*");
    hist::loadHist(FName.Data(),0,"*htopSpinCorrGen_2j*");
  }
  if(drawtprime && ! drawLogY) hist::scale("ttprime*", 10);

  //if(scalebTagMCbyDataDrivenEsts) {                                                                                                                                                                                                
  vector<TString> v_samples;
  vector<Color_t> v_colors;
  vector<TString> v_legEntries;
  vector<Style_t> v_styles;

  //the order you put the sample into the vector is the order                                                                                                                                                                                                 
  //that the histogram gets drawn. So if you put ttdil in first,                                                                                                                                                                                              
  //it will be drawn first, i.e, it will be drawn on the bottom                                                                                                                                                                                               
  //of the stack!                                                                                                                                                                                                                                             

  
  
  if(drawLogY) {
    v_samples.push_back("ttdil");
    v_colors.push_back(kRed+1);
    v_legEntries.push_back("t#bar{t}");
    v_styles.push_back(1001);
    
  }
  
   
  if(!drawLogY) {
  v_samples.push_back("wjets");
  v_colors.push_back(kGreen-3);
  if(!scalebTagMCbyDataDrivenEsts)
    v_legEntries.push_back("W#rightarrowl#nu");
  else
    v_legEntries.push_back("Non-W/Z prediction");
  v_styles.push_back(1001);    
  }


  //  v_samples.push_back("FR");
  // v_colors.push_back(kGreen-3);
  //v_legEntries.push_back("Non-W/Z prediction");
  //v_styles.push_back(1001);


  v_samples.push_back("VV");
  v_colors.push_back(kYellow-10);
  v_legEntries.push_back("VV");
  v_styles.push_back(1001);

  v_samples.push_back("tw");
  v_colors.push_back(kMagenta);
  v_legEntries.push_back("Single top");
  v_styles.push_back(1001);


  
  v_samples.push_back("DYtautau");
  v_colors.push_back(kAzure+8);
  v_legEntries.push_back("Z/#gamma*#rightarrow#tau^{+}#tau^{-}");
  v_styles.push_back(1001);


  std::vector<TString> v_prfxsToCombine;
  v_prfxsToCombine.clear();
  v_prfxsToCombine.push_back("DYee");
  v_prfxsToCombine.push_back("DYmm");
  hist::combineHists(v_prfxsToCombine, "DYeemm");


    
  v_samples.push_back("DYeemm");
  v_colors.push_back(kAzure - 2);
  //if(replacewDYest || scalebTagMCbyDataDrivenEsts)                                                                                                                                                                                                          
  //  v_legEntries.push_back("Z/#gamma*#rightarrowl^{+}l^{-} prediction");                                                                                                                                                                                    
  //else                                                                                                                                                                                                                                                      
  v_legEntries.push_back("Z/#gamma*#rightarrowl^{+}l^{-}");
  v_styles.push_back(1001);

  v_samples.push_back("ttotr");
  v_colors.push_back(kRed-1);
  v_legEntries.push_back("t#bar{t} (other)");
  v_styles.push_back(1001);
  
  
  if(!drawLogY) {
    v_samples.push_back("ttdil");
    v_colors.push_back(kRed+1);
    v_legEntries.push_back("t#bar{t} (dilepton)");
    v_styles.push_back(1001);
    
  }
    

  if(drawtprime) {
    v_samples.push_back("ttprime350");
    v_colors.push_back(kBlue-2);
    if(drawLogY) v_legEntries.push_back("t'#bar{t'} 350 GeV");
    else v_legEntries.push_back("10 #times t'#bar{t'} 350 GeV");
    v_styles.push_back(0);
  }

  
 
 
  //data should be the last thing put into the vectors                                                                                                                                                                                                        
  v_samples.push_back("data");
  v_colors.push_back(kBlack);
  v_legEntries.push_back("Data");
  v_styles.push_back(20);
 

  
  browseStacks(v_samples, v_colors, FName, v_legEntries, drawLogY, v_styles, drawFullErrors, drawDiffs, scaleMC, rebin, allhypOnly);

  hist::deleteHistos();

}
  
