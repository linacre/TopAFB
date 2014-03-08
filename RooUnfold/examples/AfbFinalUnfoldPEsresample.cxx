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
int unfoldingType = 2;

TString Region = "";
Int_t kterm = 3; //for SVD
Double_t tau = 0.005; //for TUnfold - this is a more reasonable default (1E-4 gives very little regularisation)

Int_t nVars = 8;
Int_t includeSys = 0;


void AfbFinalUnfoldPEsresample(Int_t iVar = 0, double scalettdil = 1., double scalettotr = 1., double scalewjets = 1., double scaleDY = 1., double scaletw = 1., double scaleVV = 1. )
{

    const int NPEs = 1000;
    Int_t lumiSF = 1;  //to simulate an integrated luminosity lumiSF times larger: multiplies number of events per PE and weights of ttdil events used for smearing matrix by lumiSF
    Int_t PEsize = 9084;
    PEsize*=lumiSF;
    bool split_sample = true;
    Double_t sample_split_factor = 2.;
    if (!split_sample) sample_split_factor = 1.;


    setTDRStyle();
    gStyle->SetOptFit();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cout.precision(3);

    Initialize1DBinning(iVar);
    TString summary_name = "summary_1Dunfolding_" + acceptanceName;

    if (!(scalettotr == 1. && scalewjets == 1. && scaleDY == 1. && scaletw == 1. && scaleVV == 1.))  summary_name = Form("summary_1Dunfolding_%i_%i_%i_%i_%i", int(10.*scalettotr + 0.5), int(10.*scalewjets + 0.5), int(10.*scaleDY + 0.5), int(10.*scaletw + 0.5), int(10.*scaleVV + 0.5));

    ofstream myfile;
    myfile.open (summary_name + ".txt");
    cout.rdbuf(myfile.rdbuf());

    // OGU 130516: add second output txt file with format easier to be pasted into google docs
    ofstream second_output_file;
    second_output_file.open(summary_name + "_formated.txt");

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);

    TRandom3 *random3_ = new TRandom3();


    const int nBkg = 7;
    TString path = "../";
    TString bkgroot[nBkg] = {"ttotr.root", "wjets.root", "DYee.root", "DYmm.root", "DYtautau.root", "tw.root", "VV.root"};
    double bkgSF[nBkg] = {scalettotr, scalewjets, scaleDY, scaleDY, scaleDY, scaletw, scaleVV};

    Float_t observable, observable_gen, ttmass, ttRapidity, tmass;
    Float_t observableMinus, observableMinus_gen;
    Double_t weight;
    Int_t evt;
    Int_t Nsolns;

    //for (Int_t iVar = 0; iVar < nVars; iVar++)
    //{
        //Initialize1DBinning(iVar);

        //use twice as many reco bins for TUnfold
        int recobinsmult = 2;
        if (unfoldingType == 0) recobinsmult = 1;
        const int nbins1Dreco = nbins1D * recobinsmult;
        double xbins1Dreco[nbins1Dreco + 1];

        for (int i = 0; i < nbins1Dreco + 1; ++i)
        {
            if (unfoldingType == 0) xbins1Dreco[i] = xbins1D[i];
            else xbins1Dreco[i] = (xbins1D[int(i / 2)] + xbins1D[int((i + 1) / 2)]) / 2.;  //reco-level bins from dividing gen-level bins in 2
            //else xbins1Dreco[i] = xmin + double(i) / nbins1Dreco * (xmax - xmin); //uniform reco-level bins
            //cout << xbins1Dreco[i] << endl;
        }

        bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;

        TH1D *hData = new TH1D ("Data_BkgSub", "Data with background subtracted",    nbins1Dreco, xbins1Dreco);
        TH1D *hBkg = new TH1D ("Background",  "Background",    nbins1Dreco, xbins1Dreco);
        TH1D *hData_unfolded = new TH1D ("Data_Unfold", "Data with background subtracted and unfolded", nbins1D, xbins1D);

        TH1D *hTrue = new TH1D ("true", "Truth",    nbins1D, xbins1D);
        TH1D *hTrue_PEs = new TH1D ("true_PEs", "Truth in PEs subsample", nbins1D, xbins1D);
        TH1D *hMeas = new TH1D ("meas", "Measured", nbins1Dreco, xbins1Dreco);

        TH2D *hTrue_vs_Meas = new TH2D ("true_vs_meas", "True vs Measured", nbins1Dreco, xbins1Dreco, nbins1D, xbins1D);

        TH1D *hData_bkgSub;

        TH1D *hPseudoData[NPEs];
        TH1D *hPseudoData_unfolded[NPEs];
        TMatrixD m_unfoldE_PEs[NPEs];
        TMatrixD m_correctE_PEs[NPEs];

        for (int iPE = 0; iPE < NPEs; ++iPE)
        {
            hPseudoData[iPE] = new TH1D ("", "", nbins1Dreco, xbins1Dreco);
            hPseudoData_unfolded[iPE] = new TH1D ("", "", nbins1D, xbins1D);
            m_unfoldE_PEs[iPE] = new TMatrixD(nbins1D, nbins1D);
            m_correctE_PEs[iPE] = new TMatrixD(nbins1D, nbins1D);
        }

        hData->Sumw2();
        hBkg->Sumw2();
        hData_unfolded->Sumw2();
        hTrue->Sumw2();
        hTrue_PEs->Sumw2();
        hMeas->Sumw2();
        hTrue_vs_Meas->Sumw2();

        TMatrixD m_unfoldE (nbins1D, nbins1D);
        TMatrixD m_correctE(nbins1D, nbins1D);
        TMatrixD m_unfoldcorr (nbins1D, nbins1D);


        //  Now test with data and with BKG subtraction

        TChain *ch_bkg[nBkg];
        TChain *ch_top = new TChain("tree");

        TChain *ch_data = new TChain("tree");


        ch_data->Add(path + "data.root");

        ch_top->Add(path + "ttdil.root");

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg] = new TChain("tree");
            ch_bkg[iBkg]->Add(path + bkgroot[iBkg]);
        }


        ch_data->SetBranchAddress(observablename,    &observable);
        if ( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
        ch_data->SetBranchAddress("weight", &weight);
        ch_data->SetBranchAddress("Nsolns", &Nsolns);
        ch_data->SetBranchAddress("tt_mass", &ttmass);
        ch_data->SetBranchAddress("ttRapidity", &ttRapidity);
        ch_data->SetBranchAddress("t_mass", &tmass);


        for (Int_t i = 0; i < ch_data->GetEntries(); i++)
        {
            ch_data->GetEntry(i);
            if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
            {
                // leptonic asymmetries don't need valid top mass solution
                fillUnderOverFlow(hData, observable, weight, Nsolns);
            }
            else
            {
                if ( ttmass > 0 )
                {
                    // asymmetries with top properties are required to have a valid top mass solution
                    fillUnderOverFlow(hData, observable, weight, Nsolns);
                }
            }
            if (combineLepMinus)
            {
                // combine plus and minus
                if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
                {
                    // leptonic asymmetries don't need valid top mass solution
                    fillUnderOverFlow(hData, observableMinus, weight, Nsolns);
                }
                else
                {
                    if ( ttmass > 0 )
                    {
                        // asymmetries with top properties are required to have a valid top mass solution
                        fillUnderOverFlow(hData, observableMinus, weight, Nsolns);
                    }
                }
            }
        }


        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {

            ch_bkg[iBkg]->SetBranchAddress(observablename,    &observable);
            if ( combineLepMinus ) ch_bkg[iBkg]->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
            ch_bkg[iBkg]->SetBranchAddress("weight", &weight);
            ch_bkg[iBkg]->SetBranchAddress("Nsolns", &Nsolns);
            ch_bkg[iBkg]->SetBranchAddress("tt_mass", &ttmass);
            ch_bkg[iBkg]->SetBranchAddress("ttRapidity", &ttRapidity);
            ch_bkg[iBkg]->SetBranchAddress("t_mass", &tmass);

            for (Int_t i = 0; i < ch_bkg[iBkg]->GetEntries(); i++)
            {
                ch_bkg[iBkg]->GetEntry(i);
                weight *= bkgSF[iBkg];
                if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
                {
                    // leptonic asymmetries don't need valid top mass solution
                    fillUnderOverFlow(hBkg, observable, weight, Nsolns);
                }
                else
                {
                    if ( ttmass > 0 )
                    {
                        // asymmetries with top properties are required to have a valid top mass solution
                        fillUnderOverFlow(hBkg, observable, weight, Nsolns);
                    }
                }
                if (combineLepMinus)
                {
                    if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
                    {
                        // leptonic asymmetries don't need valid top mass solution
                        fillUnderOverFlow(hBkg, observableMinus, weight, Nsolns);
                    }
                    else
                    {
                        if ( ttmass > 0 )
                        {
                            // asymmetries with top properties are required to have a valid top mass solution
                            fillUnderOverFlow(hBkg, observableMinus, weight, Nsolns);
                        }
                    }
                }
            }

        }

        ch_top->SetBranchAddress(observablename,    &observable);
        ch_top->SetBranchAddress(observablename + "_gen", &observable_gen);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms_gen",    &observableMinus_gen);
        ch_top->SetBranchAddress("weight", &weight);
        ch_top->SetBranchAddress("Nsolns", &Nsolns);
        ch_top->SetBranchAddress("tt_mass", &ttmass);
        ch_top->SetBranchAddress("ttRapidity", &ttRapidity);
        ch_top->SetBranchAddress("t_mass", &tmass);
        ch_top->SetBranchAddress("evt", &evt);

        Int_t NEvtsforPEs = 0;
        Int_t NEvtsforPEs_weighted = 0;
        Int_t prevevent = -999;
        for (Int_t i = 0; i < ch_top->GetEntries(); i++)
        {
            ch_top->GetEntry(i);
            if (evt != prevevent)
            {
                NEvtsforPEs++;
                if (weight > 0) NEvtsforPEs_weighted++;
                else NEvtsforPEs_weighted--;
            }
            prevevent = evt;
        }
        Int_t PEsize_unweighted = PEsize;
        Int_t PEsize_weighted = PEsize;
        double positive_event_SF = double(NEvtsforPEs_weighted) / double(NEvtsforPEs);
        if (NEvtsforPEs_weighted != NEvtsforPEs) PEsize_weighted = Int_t(0.5 + double(PEsize) / positive_event_SF );
        NEvtsforPEs /= sample_split_factor;

        const int Nevts = NEvtsforPEs;

        TH1D *hEvtSamplingMultiplicity = new TH1D ("EvtSamplingMultiplicity", "EvtSamplingMultiplicity", 5, 0, 5);
        double mean_event_duplication = double(PEsize_weighted * NPEs) / double(Nevts);
        TH1D *hEvtTotalMultiplicity = new TH1D ("hEvtTotalMultiplicity", "hEvtTotalMultiplicity", 2 * int(4.*sqrt(mean_event_duplication) + 0.5), int(mean_event_duplication + 0.5) - int(4.*sqrt(mean_event_duplication) + 0.5), int(mean_event_duplication + 0.5) + int(4.*sqrt(mean_event_duplication) + 0.5) );
        Int_t PENumEvts[NPEs] = {0};

        Int_t prevevt = -999;
        Int_t i_ev = -1;
        vector<Int_t> temp_event_multiplicity_vector;
        vector<Int_t> temp_iPEmapping_vector;

        for (Int_t i = 0; i < ch_top->GetEntries(); i++)
        {
            ch_top->GetEntry(i);
            weight *= scalettdil;

            if (evt != prevevt) i_ev++;


            if (!split_sample || i_ev < Nevts )
            {

                if (evt != prevevt)
                {
                    int randseed = i_ev + 1;
                    random3_->SetSeed(randseed);
                    temp_event_multiplicity_vector.clear();
                    temp_iPEmapping_vector.clear();
                    Int_t event_multiplicity_total = 0;

                    for (int iPE = 0; iPE < NPEs; ++iPE)
                    {
                        Int_t temp_event_multiplicity = random3_->Poisson( double(PEsize_weighted) / double(Nevts) );
                        hEvtSamplingMultiplicity->Fill(temp_event_multiplicity);
                        if (temp_event_multiplicity > 0)
                        {
                            temp_event_multiplicity_vector.push_back(temp_event_multiplicity);
                            temp_iPEmapping_vector.push_back(iPE);
                            if (weight > 0) PENumEvts[iPE] += temp_event_multiplicity;
                            else PENumEvts[iPE] -= temp_event_multiplicity;
                            event_multiplicity_total += temp_event_multiplicity;
                        }

                    }
                    hEvtTotalMultiplicity->Fill(event_multiplicity_total);

                }



                for (int inonzeroPE = 0; inonzeroPE < temp_event_multiplicity_vector.size(); ++inonzeroPE)
                {
                    Int_t PE_event_multiplicity = temp_event_multiplicity_vector[inonzeroPE];

                    if (PE_event_multiplicity > 0)
                    {

                        if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
                        {
                            fillUnderOverFlow(hPseudoData[ temp_iPEmapping_vector[ inonzeroPE ] ], observable, weight * double(PE_event_multiplicity), double(Nsolns) / double(PE_event_multiplicity));

                            if ( combineLepMinus )
                            {
                                fillUnderOverFlow(hPseudoData[ temp_iPEmapping_vector[ inonzeroPE ] ], observableMinus, weight * double(PE_event_multiplicity), double(Nsolns) / double(PE_event_multiplicity));
                            }
                        }
                        else
                        {
                            if ( ttmass > 0 )
                            {
                                fillUnderOverFlow(hPseudoData[ temp_iPEmapping_vector[ inonzeroPE ] ], observable, weight * double(PE_event_multiplicity), double(Nsolns) / double(PE_event_multiplicity));
                                if ( combineLepMinus )
                                {

                                    fillUnderOverFlow(hPseudoData[ temp_iPEmapping_vector[ inonzeroPE ] ], observableMinus, weight * double(PE_event_multiplicity), double(Nsolns) / double(PE_event_multiplicity));
                                }
                            }
                        }

                    }


                }

                if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
                {
                    fillUnderOverFlow(hTrue_PEs, observable_gen, weight, Nsolns);
                    if ( combineLepMinus )
                    {
                        fillUnderOverFlow(hTrue_PEs, observableMinus_gen, weight, Nsolns);
                    }
                }
                else
                {
                    if ( ttmass > 0 )
                    {
                        fillUnderOverFlow(hTrue_PEs, observable_gen, weight, Nsolns);
                        if ( combineLepMinus )
                        {
                            fillUnderOverFlow(hTrue_PEs, observableMinus_gen, weight, Nsolns);
                        }
                    }
                }

            }

            weight*=lumiSF;
            if(split_sample) weight /= 1. - 1./sample_split_factor;
            if (!split_sample || i_ev >= Nevts )
            {


                if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
                {
                    fillUnderOverFlow(hMeas, observable, weight, Nsolns);
                    fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
                    fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                    if ( combineLepMinus )
                    {
                        fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
                        fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
                        fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
                    }
                }
                else
                {
                    if ( ttmass > 0 )
                    {
                        fillUnderOverFlow(hMeas, observable, weight, Nsolns);
                        fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
                        fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                        if ( combineLepMinus )
                        {
                            fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
                            fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
                            fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
                        }
                    }
                }
            }

            prevevt = evt;
        }


        Int_t lastPE = NPEs;

        //cout << "last PE: " << lastPE << " " << Nevts << endl;

        RooUnfoldResponse response (hMeas, hTrue, hTrue_vs_Meas);

        hData_bkgSub = (TH1D *) hData->Clone();
        hData_bkgSub->Add(hBkg, -1.0);

        Double_t biasScale =  hData_bkgSub->Integral() / hMeas->Integral() ;
        hMeas->Scale(biasScale);

        if (unfoldingType == 0)
        {
            RooUnfoldSvd unfold_svd (&response, hData_bkgSub, kterm);
            unfold_svd.Setup(&response, hData_bkgSub);
            unfold_svd.IncludeSystematics(includeSys);
            hData_unfolded = (TH1D *) unfold_svd.Hreco();
            m_unfoldE = unfold_svd.Ereco();

            for (int iPE = 0; iPE < lastPE; ++iPE)
            {
                RooUnfoldSvd unfold_svd_PEs (&response, hPseudoData[iPE], kterm);
                unfold_svd_PEs.Setup(&response, hPseudoData[iPE]);
                unfold_svd_PEs.IncludeSystematics(includeSys);
                hPseudoData_unfolded[iPE] = (TH1D *) unfold_svd_PEs.Hreco();
                m_unfoldE_PEs[iPE] = unfold_svd_PEs.Ereco();
            }
        }
        else if (unfoldingType == 1)
        {
            RooUnfoldTUnfold unfold_rooTUnfold (&response, hData_bkgSub, TUnfold::kRegModeCurvature);
            unfold_rooTUnfold.Setup(&response, hData_bkgSub);
            unfold_rooTUnfold.FixTau(tau);
            unfold_rooTUnfold.IncludeSystematics(includeSys);
            hData_unfolded = (TH1D *) unfold_rooTUnfold.Hreco();
            m_unfoldE = unfold_rooTUnfold.Ereco();
        }
        else if (unfoldingType == 2)
        {

            TUnfold unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature);
            //TUnfold unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert);
            unfold_TUnfold.SetInput(hData_bkgSub);
            //unfold_TUnfold.SetBias(hTrue);  //doesn't make any difference, because if not set the bias distribution is automatically determined from hTrue_vs_Meas, which gives exactly hTrue

            /*if (doScanLCurve)
            {

                Int_t nScan = 3000;
                Int_t iBest;
                TSpline *logTauX, *logTauY;
                TGraph *lCurve;

                //restrict scan range to observed useful range of tau (5E-4 gives very weak regularisation, 2E-2 very strong)
                iBest = unfold_TUnfold.ScanLcurve(nScan, 5E-4, 2E-2, &lCurve, &logTauX, &logTauY);

                TCanvas *c_l = new TCanvas("c_l", "c_l", 1200, 600);
                c_l->Divide(3, 1);
                c_l->cd(1);
                logTauX->Draw("L");
                c_l->cd(2);
                logTauY->Draw("L");
                c_l->cd(3);
                lCurve->Draw("AC");
                c_l->SaveAs("Lcurve_" + acceptanceName + Region + ".pdf");


                std::cout << "tau=" << unfold_TUnfold.GetTau() << " iBest=" <<   iBest << "\n";

                tau = unfold_TUnfold.GetTau();

                //when an extreme value is chosen, normally it means scanLcurve didn't work very well, so reset tau to a reasonable value
                if (tau > 1E-2 || tau < 1E-3)
                {
                    tau = 0.005;
                    cout << "setting tau to 0.005" << endl;
                }

            }*/


            //biasScale = 0.0; //set biasScale to 0 when using kRegModeSize, or to compare with unfoldingType == 1
            //do the unfolding with calculated bias scale (N_data/N_MC), and tau from ScanLcurve if doScanLCurve=true. Note that the results will only be the same as unfoldingType == 1 with biasScale=0 and the same value of tau.
            cout << "bias scale for TUnfold: " << biasScale << endl;
            unfold_TUnfold.DoUnfold(tau, hData_bkgSub, biasScale);
            //unfold_TUnfold.DoUnfold(0.005,hData_bkgSub,biasScale);


            unfold_TUnfold.GetOutput(hData_unfolded);


            TH2D *ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
            TH2D *cmatrix = unfold_TUnfold.GetRhoIJ("cmatrix", "correlation matrix", 0, 0);
            for (Int_t cmi = 0; cmi < nbins1D; cmi++)
            {
                for (Int_t cmj = 0; cmj < nbins1D; cmj++)
                {
                    m_unfoldE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
                    m_unfoldcorr(cmi, cmj) = cmatrix->GetBinContent(cmi + 1, cmj + 1);
                }
            }

            for (int iPE = 0; iPE < lastPE; ++iPE)
            {
                TUnfold unfold_TUnfold_PEs (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature);
                unfold_TUnfold_PEs.SetInput(hPseudoData[iPE]);
                biasScale =  hPseudoData[iPE]->Integral() / hMeas->Integral() ;
                //cout<<"biasScale PE "<<iPE<<": "<<biasScale<<endl;
                hPseudoData[iPE]->Scale(1./biasScale);
                biasScale = 1.;
                unfold_TUnfold_PEs.DoUnfold(tau, hPseudoData[iPE], biasScale);
                unfold_TUnfold_PEs.GetOutput(hPseudoData_unfolded[iPE]);

                TH2D *ematrixPE = unfold_TUnfold_PEs.GetEmatrix("ematrix", "error matrix", 0, 0);
                //TH2D *cmatrixPE = unfold_TUnfold_PEs.GetRhoIJ("cmatrix", "correlation matrix", 0, 0);
                TMatrixD m_unfoldE_temp =  m_unfoldE_PEs[iPE];
                for (Int_t cmi = 0; cmi < nbins1D; cmi++)
                {
                    for (Int_t cmj = 0; cmj < nbins1D; cmj++)
                    {
                        m_unfoldE_temp(cmi, cmj) = ematrixPE->GetBinContent(cmi + 1, cmj + 1);
                    }
                }
                m_unfoldE_PEs[iPE] = m_unfoldE_temp;

            }


        }
        else cout << "Unfolding TYPE not Specified" << "\n";


        if (unfoldingType == 0)
        {
            TCanvas *c_d = new TCanvas("c_d", "c_d", 500, 500);
            TH1D *dvec = unfold_svd.Impl()->GetD();
            dvec->Draw();
            c_d->SetLogy();
            c_d->SaveAs("D_" + acceptanceName + Region + ".pdf");
        }

        TCanvas *c_resp = new TCanvas("c_resp", "c_resp");
        TH2D *hResp = (TH2D *) response.Hresponse();
        gStyle->SetPalette(1);
        hResp->GetXaxis()->SetTitle(xaxislabel);
        hResp->GetYaxis()->SetTitle(xaxislabel + "_{gen}");
        hResp->Draw("COLZ");
        c_resp->SaveAs("Response_" + acceptanceName + Region + ".pdf");

        TFile *file = new TFile("../acceptance/mcnlo/accept_" + acceptanceName + ".root");
        TH1D *acceptM = (TH1D *) file->Get("accept_" + acceptanceName);
        acceptM->Scale(1.0 / acceptM->Integral());

        TH1D *denominatorM = (TH1D *) file->Get("denominator_" + acceptanceName);

        for (Int_t i = 1; i <= nbins1D; i++)
        {

            if (acceptM->GetBinContent(i) != 0)
            {
                hData_unfolded->SetBinContent(i, hData_unfolded->GetBinContent(i) * 1.0 / acceptM->GetBinContent(i));
                hData_unfolded->SetBinError  (i, hData_unfolded->GetBinError  (i) * 1.0 / acceptM->GetBinContent(i));
                for (int iPE = 0; iPE < lastPE; ++iPE)
                {
                    hPseudoData_unfolded[iPE]->SetBinContent(i, hPseudoData_unfolded[iPE]->GetBinContent(i) * 1.0 / acceptM->GetBinContent(i));
                    hPseudoData_unfolded[iPE]->SetBinError  (i, hPseudoData_unfolded[iPE]->GetBinError  (i) * 1.0 / acceptM->GetBinContent(i));
                }
            }

            if (acceptM->GetBinContent(i) != 0)
            {
                hTrue->SetBinContent(i, hTrue->GetBinContent(i) * 1.0 / acceptM->GetBinContent(i));
                hTrue->SetBinError  (i, hTrue->GetBinError(i)  * 1.0 / acceptM->GetBinContent(i));
                hTrue_PEs->SetBinContent(i, hTrue_PEs->GetBinContent(i) * 1.0 / acceptM->GetBinContent(i));
                hTrue_PEs->SetBinError  (i, hTrue_PEs->GetBinError(i)  * 1.0 / acceptM->GetBinContent(i));
            }
        }

        for (int l = 0; l < nbins1D; l++)
        {
            for (int j = 0; j < nbins1D; j++)
            {
                double corr = 1.0 / ( acceptM->GetBinContent(l + 1) * acceptM->GetBinContent(j + 1) );
                //corr = corr * pow(xsection / dataIntegral,2) ;
                m_correctE(l, j) = m_unfoldE(l, j) * corr;
            }
        }


        for (int iPE = 0; iPE < lastPE; ++iPE)
        {
            TMatrixD m_unfoldE_temp =  m_unfoldE_PEs[iPE];
            TMatrixD m_correctE_temp =  m_correctE_PEs[iPE];
            for (int l = 0; l < nbins1D; l++)
            {
                for (int j = 0; j < nbins1D; j++)
                {
                    double corr = 1.0 / ( acceptM->GetBinContent(l + 1) * acceptM->GetBinContent(j + 1) );
                    m_correctE_temp(l, j) = m_unfoldE_temp(l, j) * corr;
                }
            }
            m_correctE_PEs[iPE] = m_correctE_temp;
        }


        //==================================================================
        // ============== Print the asymmetry =============================
        cout << "========= Variable: " << acceptanceName << "===================\n";

        Float_t Afb, AfbErr;

        GetAfb(hData, Afb, AfbErr);
        cout << " Data: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hTrue, Afb, AfbErr);
        cout << " True Top: " << Afb << " +/-  " << AfbErr << "\n";

        GetCorrectedAfb(hData_unfolded, m_correctE, Afb, AfbErr);
        cout << " Unfolded: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " Unfolded: " << Afb << " +/-  " << AfbErr << endl;

        GetAfb(denominatorM, Afb, AfbErr);
        cout << " True Top from acceptance denominator: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " True_Top_from_acceptance_denominator: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hTrue_PEs, Afb, AfbErr);
        Float_t Afb_true = Afb;
        cout << " True Top in PEs subsample: " << Afb << " +/-  " << AfbErr << "\n";

        TH1D *hPull = new TH1D ("pull", "pull", 100, -4, 4);
        int NevPEerr = int(sqrt(PEsize_weighted) + 0.5);
        int nBinsforNevPE = NevPEerr * ( int(100 / NevPEerr) > 0 ? int(100 / NevPEerr) : 1);
        TH1D *hNevPE = new TH1D ("NevPE", "NevPE", nBinsforNevPE, PEsize_unweighted - 4 * NevPEerr, PEsize_unweighted + 4 * NevPEerr );

        double meanAFB = 0.;
        double meanAFBerr = 0.;

        for (int iPE = 0; iPE < lastPE; ++iPE)
        {
            GetCorrectedAfb(hPseudoData_unfolded[iPE], m_correctE_PEs[iPE], Afb, AfbErr);
            hPull->Fill( (Afb - Afb_true) / AfbErr );
            hNevPE->Fill( PENumEvts[iPE] );
            meanAFB += Afb;
            meanAFBerr += AfbErr;
        }
        meanAFB /= double(lastPE);
        meanAFBerr /= double(lastPE);


        TH1D *hErr = new TH1D ("err", "err", 100, meanAFBerr * ( 1. - 4. / sqrt(PEsize_weighted) ), meanAFBerr * ( 1. + 4. / sqrt(PEsize_weighted) ) );
        TH1D *hAfb = new TH1D ("Afb", "Afb", 100, meanAFB - 4.*meanAFBerr, meanAFB + 4.*meanAFBerr);

        for (int iPE = 0; iPE < lastPE; ++iPE)
        {
            GetCorrectedAfb(hPseudoData_unfolded[iPE], m_correctE_PEs[iPE], Afb, AfbErr);
            hErr->Fill(AfbErr);
            hAfb->Fill(Afb);
        }

        hPull->GetXaxis()->SetTitle("A(" + xaxislabel + ") pull");
        hPull->GetYaxis()->SetTitle("PEs/bin");
        hNevPE->GetXaxis()->SetTitle("Number of events per PE");
        hNevPE->GetYaxis()->SetTitle("PEs/bin");
        hAfb->GetXaxis()->SetTitle("A(" + xaxislabel + ")");
        hAfb->GetYaxis()->SetTitle("PEs/bin");
        hErr->GetXaxis()->SetTitle("A(" + xaxislabel + ") uncertainty");
        hErr->GetYaxis()->SetTitle("PEs/bin");
        hEvtSamplingMultiplicity->GetXaxis()->SetTitle("PE event sampling multiplicity");
        hEvtSamplingMultiplicity->GetYaxis()->SetTitle("Number of occurances");
        hEvtTotalMultiplicity->GetXaxis()->SetTitle("Total event sampling multiplicity");
        hEvtTotalMultiplicity->GetYaxis()->SetTitle("Number of events");

        TCanvas *c_pull = new TCanvas("c_pull", "c_pull", 800, 800);
        gStyle->SetOptStat("e");
        gStyle->SetStatH(0.25);
        gStyle->SetStatW(0.30);
        gStyle->SetStatFormat("6.4g");
        gStyle->SetOptFit(1111);
        gStyle->SetFitFormat("6.4g");
        c_pull->Divide(2, 2);
        c_pull->cd(1);
        hPull->SetMaximum(1.3 * hPull->GetMaximum());
        hPull->Draw();
        hPull->Fit("gaus", "LEMV");
        c_pull->cd(2);
        hNevPE->SetMaximum(1.3 * hNevPE->GetMaximum());
        hNevPE->Draw();
        hNevPE->Fit("gaus", "LEMV");
        c_pull->cd(3);
        hAfb->SetMaximum(1.3 * hAfb->GetMaximum());
        hAfb->Draw();
        hAfb->Fit("gaus", "LEMV");
        cout<<"mean value: "<<hAfb->GetMean()<<" "<<hAfb->GetMeanError()<<endl;
        cout<<"mean RMS: "<<hAfb->GetRMS()<<" "<<hAfb->GetRMSError()<<endl;
        c_pull->cd(4);
        hErr->SetMaximum(1.3 * hErr->GetMaximum());
        hErr->Draw();
        hErr->Fit("gaus", "LEMV");
        cout<<"mean error: "<<hErr->GetMean()<<" "<<hErr->GetMeanError()<<endl;
        c_pull->SaveAs("Pull_" + acceptanceName + Region + ".pdf");
        c_pull->SaveAs("Pull_" + acceptanceName + Region + ".root");
        c_pull->SaveAs("Pull_" + acceptanceName + Region + ".C");


        TCanvas *c_Mult = new TCanvas("c_Mult", "c_Mult", 800, 400);
        c_Mult->Divide(2, 1);
        c_Mult->cd(1);
        c_Mult->SetLogy();
        hEvtSamplingMultiplicity->Draw();
        c_Mult->SetLogy();
        c_Mult->Update();
        c_Mult->cd(2);
        hEvtTotalMultiplicity->SetMaximum(1.3 * hEvtTotalMultiplicity->GetMaximum());
        hEvtTotalMultiplicity->Draw();
        hEvtTotalMultiplicity->Fit("gaus", "LEMV");
        c_Mult->SaveAs("SamplingMultiplicity_" + acceptanceName + Region + ".pdf");


        gStyle->SetOptStat(0);

        vector<double> afb_bins;
        vector<double> afb_bins_err;
        GetCorrectedAfbBinByBin(hData_unfolded, m_correctE, afb_bins, afb_bins_err, second_output_file);

        //scale to total xsec with option "width",  so that differential xsec is plotted
        hData_unfolded->Scale(1. / hData_unfolded->Integral(), "width");
        hTrue->Scale(1. / hTrue->Integral(), "width");


        for (int i = 1; i < nbins1D + 1; i++)
        {
            cout << i << " bin = " << hData_unfolded->GetBinContent(i) << " +/- " << hData_unfolded->GetBinError(i) << endl;
            second_output_file << acceptanceName << " " << observablename << " bin" << i << ": " << hData_unfolded->GetBinContent(i) << " +/- " << hData_unfolded->GetBinError(i) << endl;
        }


        TCanvas *c_test = new TCanvas("c_final", "c_final", 500, 500);

        hData_unfolded->GetXaxis()->SetTitle(xaxislabel);
        hData_unfolded->GetYaxis()->SetTitle("1/#sigma d#sigma/d(" + xaxislabel + ")");
        hData_unfolded->SetMinimum(0.0);
        hData_unfolded->SetMaximum( 2.0 * hData_unfolded->GetMaximum());
        hData_unfolded->SetMarkerStyle(23);
        hData_unfolded->SetMarkerSize(1.5);
        hData_unfolded->Draw("E");
        hData_unfolded->SetLineWidth(lineWidth);
        hTrue->SetLineWidth(lineWidth);
        hTrue->SetLineColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillStyle(3353);
        hTrue->Draw("hist same");
        hData_unfolded->Draw("E same");

        TLegend *leg1 = new TLegend(0.55, 0.62, 0.9, 0.838, NULL, "brNDC");
        leg1->SetEntrySeparation(100);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.03);

        leg1->AddEntry(hData_unfolded, "( Data - BG ) Unfolded");
        leg1->AddEntry(hTrue,    "SM parton level (mc@nlo)", "F");
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

        c_test->SaveAs("finalplot_unfolded_" + acceptanceName + Region + ".pdf");

        ch_data->Delete();

        ch_top->Delete();

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg]->Delete();
        }

    //} //end commented out iVar loop

    myfile.close();
    second_output_file.close();
}

#ifndef __CINT__
int main ()
{
    AfbUnfoldExample();    // Main program when run stand-alone
    return 0;
}
#endif
