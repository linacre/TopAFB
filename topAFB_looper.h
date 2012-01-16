#ifndef topAFB_looper_h
#define topAFB_looper_h

#include <stdint.h>
#include "TFile.h"
#include "TTree.h"

#define NCHANNELS 1
#define NHISTS 1
#define JETPTCUT 30.0
// top mass                                                                                                                                                       
#include "../CORE/topmass/ttdilepsolve.cpp" 
#include "../CORE/topmass/getTopMassEstimate.icc" 
//class TChain;

class topAFB_looper
{
    public:
        topAFB_looper();
        ~topAFB_looper();
	enum FREnum   { e_qcd = 0, e_wjets };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain*, vector<TString> , string, 
			bool, float, float, bool, FREnum frmode = e_wjets);
	void bookHistos(const char* sample, int nchannels, int nhistsets);
	void FillHistograms(const unsigned int hypIdx, const vector<unsigned int> v_jets, const vector<unsigned int> v_jetsNoEtaCut,
					 const pair<float, float> p_met, const float weight, const string prefix) ;
	bool passbTagging(const unsigned int jet_idx, const string jetAlgo, const string bTagDiscriminator) ;
	double getFRWeight(const int hypIdx, SimpleFakeRate *mufr, SimpleFakeRate *elfr, FREnum frmode, bool isData);
	double getBFRWeight(const int hypIdx, 	vector<LorentzVector> & v_goodNonBtagJets_p4,vector<LorentzVector> & v_goodBtagJets_p4, bool isData);
	void fillUnderOverFlow(TH1F *h1, float value, float weight = 1.);
	void fillUnderOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
	//void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue);
	void fillOverFlow(TH1F *h1, float value, float weight = 1.);
	void fillOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
	void fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx);
	void fillHistos(TH2F *h2[4][4],float xvalue, float yvalue, float weight, int myType, int nJetsIdx);
	void fillHistos(TProfile *h2[4][4],float xvalue, float yvalue,  int myType, int nJetsIdx);
	
	int antimatch4vector(const LorentzVector &lvec, 
					    const vector<LorentzVector> &vec);
	double triggerEff(const int hypIdx);
	
    private:
        //
        // BABY NTUPLE VARIABLES
        //
        TFile *babyFile_;
        TTree *babyTree_;
	bool applyNoCuts;
	bool getVtxDistOnly;
	bool usePtGt2020;
	bool usePtGt2010;
	bool excludePtGt2020;
	bool applylepIDCuts;
	bool applyFOv1Cuts;
	bool applyFOv2Cuts;
	bool applyFOv3Cuts;
	bool applylepIsoCuts;
	bool applylepLooseIsoCuts;
	bool applyTriggers;
	bool vetoZmass;
	bool requireZmass;
	bool hypDisamb;
	bool useCorMET;
	bool usetcMET;
	bool usetcMET35X;
	bool usepfMET;
	bool vetoMET;
	bool vetoMET50;
	bool vetoProjectedMET;
	bool usejptJets;
	bool usecaloJets;
	bool usepfJets;
	bool veto1Jet;
	bool veto2Jets;
	bool requireEcalEls;
	bool useOS;
	bool useSS;
	bool applyAlignmentCorrection;
	bool vetoHypMassLt10;
	bool vetoHypMassLt12;
	bool scaleJESMETUp; 
	bool scaleJESMETDown; 
	bool estimateQCD;
	bool estimateWJets;
	bool requireBTag;
	bool require2BTag;
	bool sortJetCandidatesbyPt;
	bool sortJetCandidatesbyDR;
	bool applyLeptonJetInvMassCut450;
	bool applyTopSystEta;
	bool requireExact2BTag;
	bool applyHTCut;
	bool matchLeptonJetbyMaxDR;
	bool generalLeptonVeto;
	bool BTagAlgTCHE;
	bool createBabyNtuples;
	bool doBFR;
	bool applyMinMassLBCut;
	float globalJESRescale;
	FactorizedJetCorrector *jptL2L3Corr;
	FactorizedJetCorrector *pfL2L3Corr;
	FactorizedJetCorrector *jptL2L3ResidualCorr;
	FactorizedJetCorrector *pfL2L3ResidualCorr;
	ttdilepsolve *d_llsol;
	
        // event identification

        Int_t   run_;
        Int_t   ls_;
        Int_t   evt_;
	Int_t ndavtx_;
	float   t_mass_;
	float weight_;
	float massltb_;
	float massllb_;
	float dr_ltjet_gen_;
	float dr_lljet_gen_;
	float tt_mass_ ;
	float ttRapidity_ ;
	float lep_charge_asymmetry_ ;
	float top_spin_correlation_ ;
	float top_costheta_cms_    ;
	float lepPlus_costheta_cms_ ;
	float tt_mass_gen_ ;
	float ttRapidity_gen_ ;
	float lep_charge_asymmetry_gen_ ;
	float top_spin_correlation_gen_ ;
	float  top_costheta_cms_gen_    ;
	float lepPlus_costheta_cms_gen_ ;
	float massllbb_ ;
	float massllbb_gen_;
	//	float llbbRapidityQuark_gen_;
	//	float llbbRapidityGluon_gen_;
	TH1F* hnJet[4];                   // Njet distributions                                                                                                          
	TH1F* hnBtagJet[4];                   // NBTagjet distributions  
	TH1F* hnVtx[4];          
	//Top Mass Plots
	TH1F *httRapidity[4][4];     
	TH1F *httMass[4][4]; 
	TH1F *hllbbMass[4][4];
	TH1F *hllbbMass_gen[4][4];
	TH1F *httMass_pull[4][4]; 
	TH1F *httMass_gen[4][4];                                                                                                                                            
	TH1F *htopMass[4][4];
	TH1F *hmassllb[4][4];
	TH1F *hmassltb[4][4];
	TH1F *hmassllb1Dmasscut[4][4];
	TH1F *hmassltb1Dmasscut[4][4];
	TH1F *htheSumJetPt[4][4];
	TH1F *htheSumBtagJetPt[4][4];
	TH1F *hthefirstJetPt[4][4];
	TH1F *hthesecondJetPt[4][4];
	TH1F *htheleadinglepPt[4][4];
	TH1F *hthesecondlepPt[4][4];
	TH1F *hthesumlepPt[4][4];
        TH1F *hlepEta[4][4];
	TH1F *hjetPt[4][4];
	TH1F *hjetEta[4][4];
	TH1F *hMET[4][4];
	TH1F *htheSumLepPt[4][4];
	TH1F *htopCosTheta[4][4];
	TH1F *hlepCosTheta[4][4];
	TH1F *hlepChargeAsym[4][4];
	TH1F *htopSpinCorr[4][4];
	TH1F *htopCosTheta_gen[4][4];
	TH1F *hlepCosTheta_gen[4][4];
	TH1F *hlepChargeAsym_gen[4][4];
	TH1F *htopSpinCorr_gen[4][4];
	TH2F *htopCosTheta_2d[4][4];
	TH2F *hlepCosTheta_2d[4][4];
	TH2F *hlepChargeAsym_2d[4][4];
	TH2F *htopSpinCorr_2d[4][4];
	TH2F *httMass_2d[4][4];

	
	TH2F *hmasslb_2d[4][4];
	TH2F *habcd_2d[4][4];


	TH1F *hdRlbtruegenp[4][4];
	TH1F *hdRlbfalsegenp[4][4];
	TH1F *hdRlbratiogenp[4][4];
	TH1F *htopptgenp[4][4];
	
	TH1F *htopMassgenp[4][4];
	TH1F *htheSumBtagJetPtgenp[4][4];
	TH1F *hthefirstBtagJetPtgenp[4][4];
	TH1F *hthesecondBtagJetPtgenp[4][4];
	TH1F *htheleadinglepPtgenp[4][4];
	TH1F *hthesecondlepPtgenp[4][4];
	TH1F *htheSumLepPtgenp[4][4];
	TH1F *htheleadingNuPtgenp[4][4];
	TH1F *hthesecondNuPtgenp[4][4];
	TH1F *hMETgenp[4][4];
	TH1F *htheSumLBPtgenp[4][4];
	TH1F *htheSumLBNPtgenp[4][4];
	TH1F *httMassGluongenp[4][4];
	TH1F *httMassQuarkgenp[4][4];
	TH1F *httRapidityGluongenp[4][4];
	TH1F *httRapidityQuarkgenp[4][4];
	TH1F *hllbbRapidityGluongenp[4][4];
	TH1F *hllbbRapidityQuarkgenp[4][4];
	

	TH2F *hmasslbgenp_2d[4][4];
	TH2F *hmasslbgenmatch1_2d[4][4];
	TH2F *hmasslbgenmatch_2d[4][4];
	TH2F *htopptdrgenp_2d[4][4];

	// For the DY Estimation
	TH1F* hmetInDYEst[4][4];         // MET in Z window for the DY Estimation
	TH1F* hmetOutDYEst[4][4];        // MET outside Z window for the DY Estimation
	TH1F* hdilMassWithMetDYEst[4][4];// Dilepton mass with MET requirement for DY estimation
	TH1F* hdilMassNoMetDYEst[4][4];  // Dilepton mass without MET requirement for DY estimation

 private:
	Int_t   ngoodlep_;
        Int_t   ngoodel_;
        Int_t   ngoodmu_;



};



#endif
