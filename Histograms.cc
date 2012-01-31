#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <vector>
#include <string>
#include "TMath.h"
#include "TDirectory.h"
#include "topAFB_looper.h"
void topAFB_looper::bookHistos(const char *prefix, int nchannel, int nhists) {
  //  Book histograms...
  //  Naming Convention:
  //  Prefix comes from the sample and it is passed to the scanning function
  //  Suffix is "ee" "em" "em" "all" which depends on the final state
  //  For example: histogram named tt_hnJet_ee would be the Njet distribution
  //  for the ee final state in the ttbar sample.
  
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
    
  char jetbins[5][7]    = {"0", "1", "2", "3", "#geq 4"};
  char suffixall[4][4]  = {"ee", "mm", "em", "all"};
  char njetCh[4][5]     = {"0j", "1j", "2j", "allj"};
  

  
  //the plots not differ by number of jets
  for (int i=0; i<4; i++) {
    
    hnJet[i] = new TH1F(Form("%s_hnJet_%s",prefix,suffixall[i]),Form("%s_nJet_%s",prefix,suffixall[i]),10,0,10);	
    hnJet[i]->GetXaxis()->SetTitle("Number of jets");
    hnJet[i]->Sumw2();

    hnBtagJet[i] = new TH1F(Form("%s_hnBtagJet_%s",prefix,suffixall[i]),Form("%s_nBtagJet_%s",prefix,suffixall[i]),6,0,6);	
    hnBtagJet[i]->GetXaxis()->SetTitle("Number of b tagged jets");
    hnBtagJet[i]->Sumw2();

    hnVtx[i] = new TH1F(Form("%s_hnVtx_%s",prefix,suffixall[i]),Form("%s_nVtx_%s",prefix,suffixall[i]),18,0,18);	
    hnVtx[i]->GetXaxis()->SetTitle("Number of vertices");
    hnVtx[i]->Sumw2();

    for (int j = 0; j < 4; j++) {
      char suffix[7];
      sprintf(suffix, "%s_%s", njetCh[j], suffixall[i]);
      
      hlepChargeAsym_2d[i][j] = new TH2F(Form("%s_hlepChargeAsym2d_%s",prefix,suffix),Form("%s_lepChargeAsym2d_%s",prefix,suffix),4,-2,2, 4,-2,2);
      hlepChargeAsym_2d[i][j]->GetXaxis()->SetTitle("Charge_Asymmetry_lep_gen");
      hlepChargeAsym_2d[i][j]->GetYaxis()->SetTitle("Charge_Asymmetry_lep_rec");
      hlepChargeAsym_2d[i][j]->Sumw2();

      hlepAzimAsym_2d[i][j] = new TH2F(Form("%s_hlepAzimAsym2d_%s",prefix,suffix),Form("%s_lepAzimAsym2d_%s",prefix,suffix),4,-2,2, 4,-2,2);
      hlepAzimAsym_2d[i][j]->GetXaxis()->SetTitle("Azimuthal_Asymmetry_lep_gen");
      hlepAzimAsym_2d[i][j]->GetYaxis()->SetTitle("Azimuthal_Asymmetry_lep_rec");
      hlepAzimAsym_2d[i][j]->Sumw2();

      htopSpinCorr_2d[i][j] = new TH2F(Form("%s_htopSpinCorr2d_%s",prefix,suffix),Form("%s_topSpinCorr2d_%s",prefix,suffix),4,-2,2,4,-2,2);
      htopSpinCorr_2d[i][j]->GetXaxis()->SetTitle("Spin_Correlation_top_gen");
      htopSpinCorr_2d[i][j]->GetYaxis()->SetTitle("Spin_Correlation_top_rec");
      htopSpinCorr_2d[i][j]->Sumw2();
      
      htopCosTheta_2d[i][j] = new TH2F(Form("%s_htopCosTheta2d_%s",prefix,suffix),Form("%s_topCosTheta2d_%s",prefix,suffix),4,-2,2, 4,-2,2);
      htopCosTheta_2d[i][j]->GetXaxis()->SetTitle("Cos(theta)_top_gen");
      htopCosTheta_2d[i][j]->GetYaxis()->SetTitle("Cos(theta)_top_rec");
      htopCosTheta_2d[i][j]->Sumw2();

      hlepCosTheta_2d[i][j] = new TH2F(Form("%s_hlepCosTheta2d_%s",prefix,suffix),Form("%s_lepCosTheta2d_%s",prefix,suffix),4,-2,2, 4,-2,2);
      hlepCosTheta_2d[i][j]->GetXaxis()->SetTitle("Cos(theta)_lep_gen");
      hlepCosTheta_2d[i][j]->GetYaxis()->SetTitle("Cos(theta)_lep_rec");
      hlepCosTheta_2d[i][j]->Sumw2();
      
      
      hlepChargeAsym_gen[i][j] = new TH1F(Form("%s_hlepChargeAsymGen_%s",prefix,suffix),Form("%s_lepChargeAsymGen_%s",prefix,suffix),4,-2,2);
      hlepChargeAsym_gen[i][j]->GetXaxis()->SetTitle("Charge_Asymmetry_lep_gen");
      hlepChargeAsym_gen[i][j]->Sumw2();

      hlepAzimAsym_gen[i][j] = new TH1F(Form("%s_hlepAzimAsymGen_%s",prefix,suffix),Form("%s_lepAzimAsymGen_%s",prefix,suffix),4,-2,2);
      hlepAzimAsym_gen[i][j]->GetXaxis()->SetTitle("Azimuthal_Asymmetry_lep_gen");
      hlepAzimAsym_gen[i][j]->Sumw2();

      htopSpinCorr_gen[i][j] = new TH1F(Form("%s_htopSpinCorrGen_%s",prefix,suffix),Form("%s_topSpinCorrgen_%s",prefix,suffix),4,-2,2);
      htopSpinCorr_gen[i][j]->GetXaxis()->SetTitle("Spin_Correlation_top_gen");
      htopSpinCorr_gen[i][j]->Sumw2();
      
      htopCosTheta_gen[i][j] = new TH1F(Form("%s_htopCosThetaGen_%s",prefix,suffix),Form("%s_topCosThetaGen_%s",prefix,suffix),4,-2,2);
      htopCosTheta_gen[i][j]->GetXaxis()->SetTitle("Cos(theta)_top_gen");
      htopCosTheta_gen[i][j]->Sumw2();

      hlepCosTheta_gen[i][j] = new TH1F(Form("%s_hlepCosThetaGen_%s",prefix,suffix),Form("%s_lepCosThetaGen_%s",prefix,suffix),4,-2,2);
      hlepCosTheta_gen[i][j]->GetXaxis()->SetTitle("Cos(theta)_lep_gen");
      hlepCosTheta_gen[i][j]->Sumw2();
 
      hlepChargeAsym[i][j] = new TH1F(Form("%s_hlepChargeAsym_%s",prefix,suffix),Form("%s_lepChargeAsym_%s",prefix,suffix),4,-2,2);
      hlepChargeAsym[i][j]->GetXaxis()->SetTitle("Charge_Asymmetry_lep");
      hlepChargeAsym[i][j]->Sumw2();

      hlepAzimAsym[i][j] = new TH1F(Form("%s_hlepAzimAsym_%s",prefix,suffix),Form("%s_lepAzimAsym_%s",prefix,suffix),4,-2,2);
      hlepAzimAsym[i][j]->GetXaxis()->SetTitle("Azimuthal_Asymmetry_lep");
      hlepAzimAsym[i][j]->Sumw2();

      htopSpinCorr[i][j] = new TH1F(Form("%s_htopSpinCorr_%s",prefix,suffix),Form("%s_topSpinCorr_%s",prefix,suffix),4,-2,2);
      htopSpinCorr[i][j]->GetXaxis()->SetTitle("Spin_Correlation_top");
      htopSpinCorr[i][j]->Sumw2();
      
      htopCosTheta[i][j] = new TH1F(Form("%s_htopCosTheta_%s",prefix,suffix),Form("%s_topCosTheta_%s",prefix,suffix),4,-2,2);
      htopCosTheta[i][j]->GetXaxis()->SetTitle("Cos(theta)_top");
      htopCosTheta[i][j]->Sumw2();

      hlepCosTheta[i][j] = new TH1F(Form("%s_hlepCosTheta_%s",prefix,suffix),Form("%s_lepCosTheta_%s",prefix,suffix),4,-2,2);
      hlepCosTheta[i][j]->GetXaxis()->SetTitle("Cos(theta)_lep");
      hlepCosTheta[i][j]->Sumw2();
      
      httMass[i][j] = new TH1F(Form("%s_httMass_%s",prefix,suffix),Form("%s_ttMass_%s",prefix,suffix),200,200.,1200.);
      httMass[i][j]->GetXaxis()->SetTitle("TTBar Mass Estimate (GeV/c^{2})");
      httMass[i][j]->Sumw2();
      
      httMass_pull[i][j] = new TH1F(Form("%s_httMass_%s",prefix,suffix),Form("%s_ttMass_%s",prefix,suffix),200,-3,3);
      httMass_pull[i][j]->GetXaxis()->SetTitle("(Reco-Gen)/Gen");
      httMass_pull[i][j]->Sumw2();
      
      httMass_gen[i][j] = new TH1F(Form("%s_httMassGen_%s",prefix,suffix),Form("%s_ttMassGen_%s",prefix,suffix),200,200.,1200.);
      httMass_gen[i][j]->GetXaxis()->SetTitle("TTBar Mass Gen Estimate (GeV/c^{2})");
      httMass_gen[i][j]->Sumw2();

      httMass_2d[i][j] = new TH2F(Form("%s_httMass2d_%s",prefix,suffix),Form("%s_ttMass2d_%s",prefix,suffix),200,200.,1200., 200,200.,1200.);
      httMass_2d[i][j]->GetXaxis()->SetTitle("TTBar Mass Gen Estimate (GeV/c^{2})");
      httMass_2d[i][j]->GetYaxis()->SetTitle("TTBar Mass Reco Estimate (GeV/c^{2})");
      httMass_2d[i][j]->Sumw2();

      hllbbMass[i][j] = new TH1F(Form("%s_hllbbMass_%s",prefix,suffix),Form("%s_llbbMass_%s",prefix,suffix),200,200.,1200.);
      hllbbMass[i][j]->GetXaxis()->SetTitle("Mass(llbb) Estimate (GeV/c^{2})");
      hllbbMass[i][j]->Sumw2();
      

      httRapidity[i][j] = new TH1F(Form("%s_httRapidity_%s",prefix,suffix),Form("%s_ttRapidity_%s",prefix,suffix),100,-10,10);
      httRapidity[i][j]->GetXaxis()->SetTitle("TTBar Rapidity");
      httRapidity[i][j]->Sumw2();

      htopMass[i][j] = new TH1F(Form("%s_htopMass_%s",prefix,suffix),Form("%s_topMass_%s",prefix,suffix),100,0.,500.);
      htopMass[i][j]->GetXaxis()->SetTitle("Top Mass Estimate (GeV/c^{2})");
      htopMass[i][j]->Sumw2();
	
      hmassltb[i][j] =  new TH1F(Form("%s_hmassltb_%s",prefix,suffix),Form("%s_massltb_%s",prefix,suffix),75,0.,510.);
      hmassltb[i][j]->GetXaxis()->SetTitle("M_{l1b1} (GeV/c^{2})");
      hmassltb[i][j]->Sumw2();

      hmassllb[i][j] =  new TH1F(Form("%s_hmassllb_%s",prefix,suffix),Form("%s_massllb_%s",prefix,suffix),75,0.,510.);
      hmassllb[i][j]->GetXaxis()->SetTitle("M_{l2b2} (GeV/c^{2})");
      hmassllb[i][j]->Sumw2();
      
      hmassltb1Dmasscut[i][j] =  new TH1F(Form("%s_hmassltb1Dmasscut_%s",prefix,suffix),Form("%s_massltb1Dmasscut_%s",prefix,suffix),75,0.,510.);
      hmassltb1Dmasscut[i][j]->GetXaxis()->SetTitle("M_{l1b1} (GeV/c^{2}) for M_{l2b2} > 170 GeV/c^{2}");
      hmassltb1Dmasscut[i][j]->Sumw2();

      hmassllb1Dmasscut[i][j] =  new TH1F(Form("%s_hmassllb1Dmasscut_%s",prefix,suffix),Form("%s_massllb1Dmasscut_%s",prefix,suffix),75,0.,510.);
      hmassllb1Dmasscut[i][j]->GetXaxis()->SetTitle("M_{l2b2} (GeV/c^{2}) for M_{l1b1} > 170 GeV/c^{2}");
      hmassllb1Dmasscut[i][j]->Sumw2();

      htheSumJetPt[i][j] =  new TH1F(Form("%s_theSumJetPt_%s",prefix,suffix),Form("%s_theSumJetPt_%s",prefix,suffix),100,0.,600.);
      htheSumJetPt[i][j]->GetXaxis()->SetTitle("Sum of jet p_{T} (GeV/c)");
      htheSumJetPt[i][j]->Sumw2();

      htheSumBtagJetPt[i][j] =  new TH1F(Form("%s_theSumBtagJetPt_%s",prefix,suffix),Form("%s_theSumBtagJetPt_%s",prefix,suffix),100,0.,600.);
      htheSumBtagJetPt[i][j]->GetXaxis()->SetTitle("Sum of b tagged jet p_{T} (GeV/c)");
      htheSumBtagJetPt[i][j]->Sumw2();

      hthefirstJetPt[i][j] =  new TH1F(Form("%s_thefirstJetPt_%s",prefix,suffix),Form("%s_thefirstJetPt_%s",prefix,suffix),60,0.,360.);
      hthefirstJetPt[i][j]->GetXaxis()->SetTitle("First jet p_{T} (GeV/c)");
      hthefirstJetPt[i][j]->Sumw2();

      hthesecondJetPt[i][j] =  new TH1F(Form("%s_thesecondJetPt_%s",prefix,suffix),Form("%s_thesecondJetPt_%s",prefix,suffix),60,0.,360.);
      hthesecondJetPt[i][j]->GetXaxis()->SetTitle("Second jet p_{T} (GeV/c)");
      hthesecondJetPt[i][j]->Sumw2();
      
      htheleadinglepPt[i][j] = new TH1F(Form("%s_htheleadinglepPt_%s",prefix,suffix),Form("%s_theleadinglepPt_%s",prefix,suffix),60,0.,240.);
	  htheleadinglepPt[i][j]->GetXaxis()->SetTitle("Leading lepton p_{T} (GeV/c)");
      htheleadinglepPt[i][j]->Sumw2();

      hthesecondlepPt[i][j] = new TH1F(Form("%s_hthesecondlepPt_%s",prefix,suffix),Form("%s_thesecondlepPt_%s",prefix,suffix),60,0.,240.);
	  hthesecondlepPt[i][j]->GetXaxis()->SetTitle("Second lepton p_{T} (GeV/c)");
      hthesecondlepPt[i][j]->Sumw2();

      htheSumLepPt[i][j] =  new TH1F(Form("%s_theSumLepPt_%s",prefix,suffix),Form("%s_theSumLepPt_%s",prefix,suffix),120,0.,480.);
      htheSumLepPt[i][j]->GetXaxis()->SetTitle("Sum of lepton p_{T} (GeV/c)");
      htheSumLepPt[i][j]->Sumw2();

      hlepEta[i][j] = new TH1F(Form("%s_hlepEta_%s",prefix,suffix),Form("%s_lepEta_%s",prefix,suffix),60,-3.,3.);
      hlepEta[i][j]->GetXaxis()->SetTitle("Lepton #eta");
      hlepEta[i][j]->Sumw2();
      
      
      hjetPt[i][j] = new TH1F(Form("%s_hjetPt_%s",prefix,suffix),Form("%s_jetPt_%s",prefix,suffix),60,0.,360.);
      hjetPt[i][j]->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      hjetPt[i][j]->Sumw2();

      hjetEta[i][j] = new TH1F(Form("%s_hjetEta_%s",prefix,suffix),Form("%s_jetEta_%s",prefix,suffix),60,-3.,3.);
      hjetEta[i][j]->GetXaxis()->SetTitle("Jet #eta");
      hjetEta[i][j]->Sumw2();
      
     
      hMET[i][j] = new TH1F(Form("%s_hMET_%s",prefix,suffix),Form("%s_MET_%s",prefix,suffix),60,0.,360.);
      hMET[i][j]->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
      hMET[i][j]->Sumw2();

       
      hmasslb_2d[i][j] = new TH2F(Form("%s_hmasslb2d_%s",prefix,suffix),  Form("%s_masslb2d_%s" ,prefix,suffix),500,0,500,500,0,500);
      hmasslb_2d[i][j]->GetXaxis()->SetTitle("M_{l2b2}(GeV/c^{2})");
      hmasslb_2d[i][j]->GetYaxis()->SetTitle("M_{l1b1}(GeV/c^{2})");
      hmasslb_2d[i][j]->Sumw2();

      
      habcd_2d[i][j] = new TH2F(Form("%s_habcd2d_%s",prefix,suffix), Form("%s_habcd2d_%s" ,prefix,suffix),500,0,500,500,0,500);
      habcd_2d[i][j]->GetXaxis()->SetTitle("M_{l2b2}(GeV/c^{2})");
      habcd_2d[i][j]->GetYaxis()->SetTitle("M_{l1b1}(GeV/c^{2})");
      habcd_2d[i][j]->Sumw2();
	

      // generator level distributions

      httMassGluongenp[i][j] = new TH1F(Form("%s_httMassGluongenp_%s",prefix,suffix),Form("%s_ttMassGluongenp_%s",prefix,suffix),200, 0 , 1000);
      httMassGluongenp[i][j]->GetXaxis()->SetTitle("ttMassGluongenp ");
      httMassGluongenp[i][j]->Sumw2();
      
      httMassQuarkgenp[i][j] = new TH1F(Form("%s_httMassQuarkgenp_%s",prefix,suffix),Form("%s_ttMassQuarkgenp_%s",prefix,suffix),200, 0 , 1000);
      httMassQuarkgenp[i][j]->GetXaxis()->SetTitle("ttMassQuarkgenp ");
      httMassQuarkgenp[i][j]->Sumw2();

      httRapidityGluongenp[i][j] = new TH1F(Form("%s_httRapidityGluongenp_%s",prefix,suffix),Form("%s_ttRapidityGluongenp_%s",prefix,suffix),200, -10, 10);
      httRapidityGluongenp[i][j]->GetXaxis()->SetTitle("ttRapidityGluongenp ");
      httRapidityGluongenp[i][j]->Sumw2();
      
      httRapidityQuarkgenp[i][j] = new TH1F(Form("%s_httRapidityQuarkgenp_%s",prefix,suffix),Form("%s_ttRapidityQuarkgenp_%s",prefix,suffix),200,-10, 10);
      httRapidityQuarkgenp[i][j]->GetXaxis()->SetTitle("ttRapidityQuarkgenp ");
      httRapidityQuarkgenp[i][j]->Sumw2();

      hllbbRapidityGluongenp[i][j] = new TH1F(Form("%s_hllbbRapidityGluongenp_%s",prefix,suffix),Form("%s_llbbRapidityGluongenp_%s",prefix,suffix),200, -10, 10);
      hllbbRapidityGluongenp[i][j]->GetXaxis()->SetTitle("llbbRapidityGluongenp ");
      hllbbRapidityGluongenp[i][j]->Sumw2();
      
      hllbbRapidityQuarkgenp[i][j] = new TH1F(Form("%s_hllbbRapidityQuarkgenp_%s",prefix,suffix),Form("%s_llbbRapidityQuarkgenp_%s",prefix,suffix),200,-10, 10);
      hllbbRapidityQuarkgenp[i][j]->GetXaxis()->SetTitle("llbbRapidityQuarkgenp ");
      hllbbRapidityQuarkgenp[i][j]->Sumw2();


      htheSumBtagJetPtgenp[i][j] =  new TH1F(Form("%s_theSumBtagJetPtgenp_%s",prefix,suffix),Form("%s_theSumBtagJetPtgenp_%s",prefix,suffix),100,0.,500.);
      htheSumBtagJetPtgenp[i][j]->GetXaxis()->SetTitle("SumBtagJetPtgenp (GeV/c)");
      htheSumBtagJetPtgenp[i][j]->Sumw2();

      hthefirstBtagJetPtgenp[i][j] =  new TH1F(Form("%s_thefirstBtagJetPtgenp_%s",prefix,suffix),Form("%s_thefirstBtagJetPtgenp_%s",prefix,suffix),100,0.,500.);
      hthefirstBtagJetPtgenp[i][j]->GetXaxis()->SetTitle("firstBtagJetPtgenp (GeV/c)");
      hthefirstBtagJetPtgenp[i][j]->Sumw2();

      hthesecondBtagJetPtgenp[i][j] =  new TH1F(Form("%s_thesecondBtagJetPtgenp_%s",prefix,suffix),Form("%s_thesecondBtagJetPtgenp_%s",prefix,suffix),100,0.,500.);
      hthesecondBtagJetPtgenp[i][j]->GetXaxis()->SetTitle("secondBtagJetPtgenp (GeV/c)");
      hthesecondBtagJetPtgenp[i][j]->Sumw2();
      
      htheleadinglepPtgenp[i][j] = new TH1F(Form("%s_htheleadinglepPtgenp_%s",prefix,suffix),Form("%s_theleadinglepPtgenp_%s",prefix,suffix),50,0.,250.);
      htheleadinglepPtgenp[i][j]->GetXaxis()->SetTitle("leadingLeptongenp p_{T} (GeV/c)");
      htheleadinglepPtgenp[i][j]->Sumw2();

      hthesecondlepPtgenp[i][j] = new TH1F(Form("%s_hthesecondlepPtgenp_%s",prefix,suffix),Form("%s_thesecondlepPtgenp_%s",prefix,suffix),50,0.,250.);
      hthesecondlepPtgenp[i][j]->GetXaxis()->SetTitle("secondLeptongenp p_{T} (GeV/c)");
      hthesecondlepPtgenp[i][j]->Sumw2();

      htheSumLepPtgenp[i][j] =  new TH1F(Form("%s_theSumLepPtgenp_%s",prefix,suffix),Form("%s_theSumLepPtgenp_%s",prefix,suffix),100,0.,500.);
      htheSumLepPtgenp[i][j]->GetXaxis()->SetTitle("SumLepPtgenp (GeV/c)");
      htheSumLepPtgenp[i][j]->Sumw2();

      htheleadingNuPtgenp[i][j] = new TH1F(Form("%s_htheleadingNuPtgenp_%s",prefix,suffix),Form("%s_theleadingNuPtgenp_%s",prefix,suffix),50,0.,250.);
      htheleadingNuPtgenp[i][j]->GetXaxis()->SetTitle("leadingNugenp p_{T} (GeV/c)");
      htheleadingNuPtgenp[i][j]->Sumw2();

      hthesecondNuPtgenp[i][j] = new TH1F(Form("%s_hthesecondNuPtgenp_%s",prefix,suffix),Form("%s_thesecondNuPtgenp_%s",prefix,suffix),50,0.,250.);
      hthesecondNuPtgenp[i][j]->GetXaxis()->SetTitle("secondNugenp p_{T} (GeV/c)");
      hthesecondNuPtgenp[i][j]->Sumw2();

      hMETgenp[i][j] = new TH1F(Form("%s_hMETgenp_%s",prefix,suffix),Form("%s_METgenp_%s",prefix,suffix),50,0.,250.);
      hMETgenp[i][j]->GetXaxis()->SetTitle("METgenp (GeV)");
      hMETgenp[i][j]->Sumw2();

      htheSumLBPtgenp[i][j] =  new TH1F(Form("%s_theSumLBPtgenp_%s",prefix,suffix),Form("%s_theSumLBPtgenp_%s",prefix,suffix),100,0.,500.);
      htheSumLBPtgenp[i][j]->GetXaxis()->SetTitle("SumLBPtgenp (GeV/c)");
      htheSumLBPtgenp[i][j]->Sumw2();

      htheSumLBNPtgenp[i][j] =  new TH1F(Form("%s_theSumLBNPtgenp_%s",prefix,suffix),Form("%s_theSumLBNPtgenp_%s",prefix,suffix),100,0.,500.);
      htheSumLBNPtgenp[i][j]->GetXaxis()->SetTitle("SumLBNPtgenp (GeV/c)");
      htheSumLBNPtgenp[i][j]->Sumw2();


      hdRlbtruegenp[i][j] = new TH1F(Form("%s_hdRlbtruegenp_%s",prefix,suffix),Form("%s_dRlbtruegenp_%s",prefix,suffix),35,0.,3.5);
      hdRlbtruegenp[i][j]->GetXaxis()->SetTitle("dRlbtruegenp ");
      hdRlbtruegenp[i][j]->Sumw2();

      
      hdRlbfalsegenp[i][j] = new TH1F(Form("%s_hdRlbfalsegenp_%s",prefix,suffix),Form("%s_dRlbfalsegenp_%s",prefix,suffix),35,0.,3.5);
      hdRlbfalsegenp[i][j]->GetXaxis()->SetTitle("dRlbfalsegenp ");
      hdRlbfalsegenp[i][j]->Sumw2();

      hdRlbratiogenp[i][j] = new TH1F(Form("%s_hdRlbratiogenp_%s",prefix,suffix),Form("%s_dRlbratiogenp_%s",prefix,suffix),40,-2.,2);
      hdRlbratiogenp[i][j]->GetXaxis()->SetTitle("dRlbratiogenp ");
      hdRlbratiogenp[i][j]->Sumw2();
      
      htopptgenp[i][j] = new TH1F(Form("%s_htopptgenp_%s",prefix,suffix),Form("%s_topptgenp_%s",prefix,suffix),100, 0 , 500);
      htopptgenp[i][j]->GetXaxis()->SetTitle("topptgenp ");
      htopptgenp[i][j]->Sumw2();

      htopMassgenp[i][j] = new TH1F(Form("%s_htopMassgenp_%s",prefix,suffix),Form("%s_topMassgenp_%s",prefix,suffix),500, 0 , 500);
      htopMassgenp[i][j]->GetXaxis()->SetTitle("topMassgenp ");
      htopMassgenp[i][j]->Sumw2();
      
      htopptdrgenp_2d[i][j] = new TH2F(Form("%s_htopptdrgenp2d_%s",prefix,suffix), Form("%s_htopptdrgenp2d_%s" ,prefix,suffix),100,0,500,35,0,3.5);
      htopptdrgenp_2d[i][j]->GetXaxis()->SetTitle("top p_{T} (GeV/c)");
      htopptdrgenp_2d[i][j]->GetYaxis()->SetTitle("DR_lb");
      htopptdrgenp_2d[i][j]->Sumw2();
    
        
      hmasslbgenp_2d[i][j] = new TH2F(Form("%s_hmasslbgenp2d_%s",prefix,suffix), Form("%s_masslbgenp2d_%s" ,prefix,suffix),500,0,500,500,0,500);
      hmasslbgenp_2d[i][j]->GetXaxis()->SetTitle("Gen M_{l2b2}(GeV/c^{2})");
      hmasslbgenp_2d[i][j]->GetYaxis()->SetTitle("Gen M_{l1b1}(GeV/c^{2})");
      hmasslbgenp_2d[i][j]->Sumw2();

      hmasslbgenmatch1_2d[i][j] = new TH2F(Form("%s_hmasslbgenmatch12d_%s",prefix,suffix), Form("%s_masslbgenmatch12d_%s" ,prefix,suffix),500,0,500,500,0,500);
      hmasslbgenmatch1_2d[i][j]->GetXaxis()->SetTitle("Gen M_{l2b2}(GeV/c^{2})");
      hmasslbgenmatch1_2d[i][j]->GetYaxis()->SetTitle("Gen M_{l1b1}(GeV/c^{2})");
      hmasslbgenmatch1_2d[i][j]->Sumw2();
     
      hmasslbgenmatch_2d[i][j] = new TH2F(Form("%s_hmasslbgenmatch2d_%s",prefix,suffix), Form("%s_masslbgenmatch2d_%s" ,prefix,suffix),500,0,500,500,0,500);
      hmasslbgenmatch_2d[i][j]->GetXaxis()->SetTitle("Gen M_{l2b2}(GeV/c^{2})");
      hmasslbgenmatch_2d[i][j]->GetYaxis()->SetTitle("Gen M_{l1b1}(GeV/c^{2})");
      hmasslbgenmatch_2d[i][j]->Sumw2();
      
            
      //DYEst histos
      hdilMassWithMetDYEst[i][j] = new TH1F(Form("%s_hdilMassWithMetDYEst_%s",  prefix,suffix), "Di-lepton mass with MET for DY Estimation", 40, 0., 200.);
      hdilMassWithMetDYEst[i][j]->GetXaxis()->SetTitle("M_{ll}(GeV/c^{2}), with MET > 30 GeV cut");
      hdilMassWithMetDYEst[i][j]->Sumw2();

      hdilMassNoMetDYEst[i][j] = new TH1F(Form("%s_hdilMassNoMetDYEst_%s",  prefix,suffix), "Di-lepton mass without MET for DY Estimation", 40, 0., 200.);
      hdilMassNoMetDYEst[i][j]->GetXaxis()->SetTitle("M_{ll}(GeV/c^{2}), no MET cut");
      hdilMassNoMetDYEst[i][j]->Sumw2();

      hmetInDYEst[i][j] = new TH1F(Form("%s_hmetInDYEst_%s",  prefix,suffix), "MET in Z mass for DY Estimation", 40, 0., 200.);
      hmetInDYEst[i][j]->GetXaxis()->SetTitle("MET (GeV), inside Z mass window");
      hmetInDYEst[i][j]->Sumw2();

      hmetOutDYEst[i][j] = new TH1F(Form("%s_hmetOutDYEst_%s",  prefix,suffix), "MET outside Z mass for DY Estimation", 40, 0., 200.);
      hmetOutDYEst[i][j]->GetXaxis()->SetTitle("MET (GeV), outside Z mass window");
      hmetOutDYEst[i][j]->Sumw2();
    }
    
  }

 

}
