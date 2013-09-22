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
int unfoldingType = 0;

TString Region = "";
Int_t kterm = 3; //for SVD
Double_t tau = 0.005; //for TUnfold - this is a more reasonable default (1E-4 gives very little regularisation)
bool doScanLCurve = false; //determine tau automatically when using unfoldingType=2, using scanLcurve (overrides value set above) - doesn't work very well
Int_t nVars = 8;
Int_t includeSys = 0;
bool checkErrors = true; //turn this on when making the final plots for the paper, to check the hard-coded systematics have been correctly entered
bool draw_truth_before_pT_reweighting = true; //turn this on when making the final plots for the paper (want to compare the data against the unweighted MC)
//bool drawTheory = true; //turn this on to show Bernreuther's predictions for AdeltaPhi and Ac1c2


void AfbUnfoldExample(double scalettdil = 1., double scalettotr = 1., double scalewjets = 1., double scaleDY = 1., double scaletw = 1., double scaleVV = 1. )
{

    setTDRStyle();
    gStyle->SetOptFit();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cout.precision(3);

    TString summary_name = "summary_1Dunfolding";

    if (!(scalettotr == 1. && scalewjets == 1. && scaleDY == 1. && scaletw == 1. && scaleVV == 1.))  summary_name = Form("summary_1Dunfolding_%i_%i_%i_%i_%i", int(10.*scalettotr + 0.5), int(10.*scalewjets + 0.5), int(10.*scaleDY + 0.5), int(10.*scaletw + 0.5), int(10.*scaleVV + 0.5));

    ofstream myfile;
    myfile.open (summary_name + ".txt");
    cout.rdbuf(myfile.rdbuf());

    // OGU 130516: add second output txt file with format easier to be pasted into google docs
    ofstream second_output_file;
    second_output_file.open(summary_name + "_formated.txt");

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);


    const int nBkg = 7;
    TString path = "../";
    TString bkgroot[nBkg] = {"ttotr.root", "wjets.root", "DYee.root", "DYmm.root", "DYtautau.root", "tw.root", "VV.root"};
    double bkgSF[nBkg] = {scalettotr, scalewjets, scaleDY, scaleDY, scaleDY, scaletw, scaleVV};

    Float_t observable, observable_gen, ttmass, ttRapidity, tmass;
    Float_t observableMinus, observableMinus_gen;
    Double_t weight;
    Int_t Nsolns;

    for (Int_t iVar = 0; iVar < nVars; iVar++)
    {

        Initialize1DBinning(iVar);

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
        TH1D *hMeas = new TH1D ("meas", "Measured", nbins1Dreco, xbins1Dreco);
        TH1D *denominatorM_nopTreweighting = new TH1D ("denominatorM_nopTreweighting", "denominatorM_nopTreweighting",    nbins1D, xbins1D);

        TH2D *hTrue_vs_Meas = new TH2D ("true_vs_meas", "True vs Measured", nbins1Dreco, xbins1Dreco, nbins1D, xbins1D);

        TH1D *hData_bkgSub;

        hData->Sumw2();
        hBkg->Sumw2();
        hData_unfolded->Sumw2();
        hTrue->Sumw2();
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
            // if ( (Region=="Signal") && (ttmass>450) ) {
            //   fillUnderOverFlow(hData, observable, weight, Nsolns);
            // } else if (Region=="") {
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
            // }
            if (combineLepMinus)
            {
                // combine plus and minus
                // if ( (Region=="Signal") && (ttmass>450) ) {
                //   fillUnderOverFlow(hData, observable, weight, Nsolns);
                // } else if (Region=="") {
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
                // }
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
                // if ( (Region=="Signal") && (ttmass>450) ) {
                //   fillUnderOverFlow(hBkg, observable, weight, Nsolns);
                // } else if (Region=="") {
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
                // }
                if (combineLepMinus)
                {
                    // combine plus and minus
                    // if ( (Region=="Signal") && (ttmass>450) ) {
                    //   fillUnderOverFlow(hBkg, observable, weight, Nsolns);
                    // } else if (Region=="") {
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
                    // }
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

        for (Int_t i = 0; i < ch_top->GetEntries(); i++)
        {
            ch_top->GetEntry(i);
            weight *= scalettdil;
            // if ( (Region=="Signal") && (ttmass>450) ) {
            //   fillUnderOverFlow(hMeas, observable, weight, Nsolns);
            //   fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
            //   fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
            //   if( combineLepMinus ) {
            //        //response.Fill (observableMinus, observableMinus_gen, weight);
            //        fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
            //        fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
            //        fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
            //      }
            // } else if (Region=="") {
            if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
            {
                //response.Fill (observable, observable_gen, weight);
                fillUnderOverFlow(hMeas, observable, weight, Nsolns);
                fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
                fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                if ( combineLepMinus )
                {
                    //response.Fill (observableMinus, observableMinus_gen, weight);
                    fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
                    fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
                    fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
                }
            }
            else
            {
                if ( ttmass > 0 )
                {
                    //response.Fill (observable, observable_gen, weight);
                    fillUnderOverFlow(hMeas, observable, weight, Nsolns);
                    fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
                    fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                    if ( combineLepMinus )
                    {
                        //response.Fill (observableMinus, observableMinus_gen, weight);
                        fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
                        fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
                        fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
                    }
                }
            }
            // }
            //if(i % 10000 == 0) cout<<i<<" "<<ch_top->GetEntries()<<endl;
        }

        RooUnfoldResponse response (hMeas, hTrue, hTrue_vs_Meas);
        //hTrue_vs_Meas = (TH2D*) response.Hresponse()->Clone();

        hData_bkgSub = (TH1D *) hData->Clone();
        hData_bkgSub->Add(hBkg, -1.0);

        Double_t biasScale =  hData_bkgSub->Integral() / hMeas->Integral() ;
        hMeas->Scale(biasScale);

        TCanvas *c_reco = new TCanvas("c_reco", "c_reco", 500, 500);

        hData->SetLineWidth(lineWidth + 2);

        hMeas->SetLineColor(TColor::GetColorDark(kGreen));
        hMeas->SetFillColor(TColor::GetColorDark(kGreen));
        hMeas->SetFillStyle(3353);

        hBkg->SetLineColor(kYellow);
        hBkg->SetFillColor(kYellow);


        THStack *hMC = new THStack("hMC", "Stacked Top+BG");

        hMC->Add(hBkg);
        hMC->Add(hMeas);

        hMC->SetMinimum(0.0);
        hMC->SetMaximum( 1.5 * hMC->GetMaximum());
        hMC->Draw("hist");
        hMC->GetXaxis()->SetTitle(xaxislabel);
        hMC->GetYaxis()->SetTitleOffset(1.3);
        hMC->GetYaxis()->SetTitle("Events/bin");

        hData->Draw("E same");

        TLegend *leg0 = new TLegend(0.58, 0.75, 0.9, 0.93, NULL, "brNDC");
        leg0->SetEntrySeparation(100);
        leg0->SetFillColor(0);
        leg0->SetLineColor(0);
        leg0->SetBorderSize(0);
        leg0->SetTextSize(0.03);
        leg0->SetFillStyle(0);
        leg0->AddEntry(hData, "Data");
        leg0->AddEntry(hMeas,  "MC@NLO reco level", "F");
        leg0->AddEntry(hBkg,  "Background", "F");
        leg0->Draw();
        c_reco->SaveAs("Reco_" + acceptanceName + Region + ".pdf");


        if (unfoldingType == 0)
        {
            RooUnfoldSvd unfold_svd (&response, hData_bkgSub, kterm);
            unfold_svd.Setup(&response, hData_bkgSub);
            unfold_svd.IncludeSystematics(includeSys);
            hData_unfolded = (TH1D *) unfold_svd.Hreco();
            m_unfoldE = unfold_svd.Ereco();
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

            if (doScanLCurve)
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

            }


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


        }
        else cout << "Unfolding TYPE not Specified" << "\n";

        //m_unfoldE.Print("f=%1.5g ");
        //m_unfoldcorr.Print("f=%1.5g ");

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
        c_resp->SaveAs("Response_" + acceptanceName + Region + ".C");
        c_resp->SaveAs("Response_" + acceptanceName + Region + ".root");

        TFile *file = new TFile("../acceptance/mcnlo/accept_" + acceptanceName + ".root");
        TH1D *acceptM = (TH1D *) file->Get("accept_" + acceptanceName);
        acceptM->Scale(1.0 / acceptM->Integral());

        TH1D *denominatorM = (TH1D *) file->Get("denominator_" + acceptanceName);


        TFile *file_nopTreweighting = new TFile("../acceptance/mcnlo_nopTreweighting/accept_" + acceptanceName + ".root");
        TH1D *denominatorM_nopTreweighting_raw = (TH1D *) file_nopTreweighting->Get("denominator_" + acceptanceName);




        for (Int_t i = 1; i <= nbins1D; i++)
        {

            if (acceptM->GetBinContent(i) != 0)
            {
                hData_unfolded->SetBinContent(i, hData_unfolded->GetBinContent(i) * 1.0 / acceptM->GetBinContent(i));
                hData_unfolded->SetBinError  (i, hData_unfolded->GetBinError  (i) * 1.0 / acceptM->GetBinContent(i));
            }

            if (acceptM->GetBinContent(i) != 0)
            {
                hTrue->SetBinContent(i, hTrue->GetBinContent(i) * 1.0 / acceptM->GetBinContent(i));
                hTrue->SetBinError  (i, hTrue->GetBinError(i)  * 1.0 / acceptM->GetBinContent(i));
            }

            denominatorM_nopTreweighting->SetBinContent(i, denominatorM_nopTreweighting_raw->GetBinContent(i));

        }

        denominatorM_nopTreweighting->Scale(1. / denominatorM_nopTreweighting->Integral(), "width");

        //scaling is now moved to after Afb is calculated
        //double dataIntegral = hData_unfolded->Integral();

        for (int l = 0; l < nbins1D; l++)
        {
            for (int j = 0; j < nbins1D; j++)
            {
                double corr = 1.0 / ( acceptM->GetBinContent(l + 1) * acceptM->GetBinContent(j + 1) );
                //corr = corr * pow(xsection / dataIntegral,2) ;
                m_correctE(l, j) = m_unfoldE(l, j) * corr;
            }
        }





        TProfile *theoryProfileCorr = new TProfile("thprofilecorrelated", "correlated data from theory file", nbins1D, xbins1D);
        TProfile *theoryProfileUnCorr = new TProfile("thprofileuncorrelated", "uncorrelated data from theory file", nbins1D, xbins1D);

        if (observablename == "lep_azimuthal_asymmetry2")
        {

            Float_t dphi, v1, v2, v3;
            Int_t ncols, nlines;
            nlines = 0;
            FILE *fp = fopen("theory/lhc7_mu1m_dphill_corr.dat", "r");
            while (1)
            {
                ncols = fscanf(fp, "%f %f %f %f", &dphi, &v1, &v2, &v3);
                if (ncols < 0) break;
                if (nlines < 5) printf("dphi=%8f, v=%8f\n", dphi, v3);
                theoryProfileCorr->Fill(dphi, v3);
                nlines++;
            }

            nlines = 0;
            fp = fopen("theory/lhc7_mu1m_dphill_uncorr.dat", "r");
            while (1)
            {
                ncols = fscanf(fp, "%f %f %f %f", &dphi, &v1, &v2, &v3);
                if (ncols < 0) break;
                if (nlines < 5) printf("dphi=%8f, v=%8f\n", dphi, v3);
                theoryProfileUnCorr->Fill(dphi, v3);
                nlines++;
            }
        }

        if (observablename == "top_spin_correlation")
        {

            Float_t c1c2, v1, v2, v3;
            Int_t ncols, nlines;
            nlines = 0;
            FILE *fp = fopen("theory/lhc7_mu1m_cos1cos2.dat", "r");
            while (1)
            {
                ncols = fscanf(fp, "%f %f %f %f", &c1c2, &v1, &v2, &v3);
                if (ncols < 0) break;
                if (nlines < 5) printf("c1c2=%8f, v=%8f\n", c1c2, v1);
                theoryProfileCorr->Fill(c1c2, v1);
                nlines++;
            }

            nlines = 0;
            fp = fopen("theory/lhc7_uncorr_mu1m_cos1cos2.dat", "r");
            while (1)
            {
                ncols = fscanf(fp, "%f %f", &c1c2, &v1);
                if (ncols < 0) break;
                if (nlines < 5) printf("c1c2=%8f, v=%8f\n", c1c2, v1);
                theoryProfileUnCorr->Fill(c1c2, v1);
                nlines++;
            }

        }



        //==================================================================
        // ============== Print the assymetry =============================
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

        //GetAfbBinByBin(hData_unfolded);

        //GetAfb(hData_unfolded, Afb, AfbErr);
        //cout<<" Unfolded (ignoring correlation): "<< Afb <<" +/-  "<< AfbErr<<"\n";


        if (observablename == "lep_azimuthal_asymmetry2" || observablename == "top_spin_correlation")
        {
            GetAfb_integratewidth( (TH1D *) theoryProfileCorr, Afb, AfbErr);
            cout << " Bernreuther correlated: " << Afb << " +/-  " << AfbErr << "\n";
            GetAfb_integratewidth( (TH1D *) theoryProfileUnCorr, Afb, AfbErr);
            cout << " Bernreuther uncorrelated: " << Afb << " +/-  " << AfbErr << "\n";
        }

        vector<double> afb_bins;
        vector<double> afb_bins_err;
        GetCorrectedAfbBinByBin(hData_unfolded, m_correctE, afb_bins, afb_bins_err, second_output_file);

        //scale to total xsec with option "width",  so that differential xsec is plotted
        //hData_unfolded->Scale(xsection/hData_unfolded->Integral(),"width");
        //hTrue->Scale(xsection/hTrue->Integral(),"width");


        hData_unfolded_clone = (TH1D *) hData_unfolded->Clone("hData_unfolded_clone");

        hData_unfolded->Scale(1. / hData_unfolded->Integral(), "width");
        hTrue->Scale(1. / hTrue->Integral(), "width");

        for (int i = 1; i < nbins1D + 1; i++)
        {
            cout << i << " bin = " << hData_unfolded->GetBinContent(i) << " +/- " << hData_unfolded->GetBinError(i) << endl;
            second_output_file << acceptanceName << " " << observablename << " bin" << i << ": " << hData_unfolded->GetBinContent(i) << " +/- " << hData_unfolded->GetBinError(i) << endl;
            //second_output_file << acceptanceName << " " << observablename << " truthbin" << i << ": " << hTrue->GetBinContent(i) << " +/- " << hTrue->GetBinError(i) << endl;
        }

        //calculate covariance matrix for normalised distribution
        for (int l = 0; l < nbins1D; l++)
        {
            for (int j = 0; j < nbins1D; j++)
            {
                m_correctE(l, j) = m_correctE(l, j) * (hData_unfolded->GetBinWidth(l + 1) * hData_unfolded->GetBinContent(l + 1) / hData_unfolded_clone->GetBinContent(l + 1)) * (hData_unfolded->GetBinWidth(j + 1) * hData_unfolded->GetBinContent(j + 1) / hData_unfolded_clone->GetBinContent(j + 1));
            }
        }

        m_correctE.Print("f=%1.5g ");

        //confirm covariance matrix for normalised distribution is correct by re-calculating Afb
        //GetCorrectedAfb_integratewidth(hData_unfolded, m_correctE, Afb, AfbErr);
        //cout << " Unfolded_after_scaling: " << Afb << " +/-  " << AfbErr << "\n";

        TH1D *hData_unfolded_minussyst;
        TH1D *hData_unfolded_plussyst;
        hData_unfolded_minussyst = (TH1D *) hData_unfolded->Clone("Data_unfolded_minussyst");
        hData_unfolded_plussyst = (TH1D *) hData_unfolded->Clone("Data_unfolded_plussyst");

        for (Int_t i = 1; i <= nbins1D; i++)
        {
            if (checkErrors)
            {
                if (includeSys)
                {
                    cout << "Difference between calculated and hard-coded stat errors: " << hData_unfolded->GetBinError(i) -  stat_corr[i - 1] << endl;
                }
                else
                {
                    cout << "Difference between calculated and hard-coded stat errors: " << hData_unfolded->GetBinError(i) -  stat_uncorr[i - 1] << endl;
                }
            }
            //hData_unfolded          ->SetBinError(i, stat_uncorr[i - 1]);  //running with includeSys = 0 means we can use the RooUnfold stat-only errors
            hData_unfolded_minussyst->SetBinContent(i, hData_unfolded->GetBinContent(i)
                                                    - sqrt(  pow(syst_corr[i - 1], 2)));  //hard-coded syst_corr now includes unfolding syst
            hData_unfolded_minussyst->SetBinError(i, 0);
            hData_unfolded_plussyst ->SetBinContent(i, 2 * sqrt( pow(syst_corr[i - 1], 2)));  //hard-coded syst_corr now includes unfolding syst
            hData_unfolded_plussyst ->SetBinError(i, 0);
        }

        THStack *hs = new THStack("hs_systband", "Systematic band");
        hData_unfolded_minussyst->SetLineColor(10);
        hData_unfolded_minussyst->SetFillColor(10);
        hData_unfolded_minussyst->SetFillStyle(0);
        hs->Add(hData_unfolded_minussyst);
        hData_unfolded_plussyst->SetFillStyle(3353);
        hData_unfolded_plussyst->SetLineColor(kWhite);
        hData_unfolded_plussyst->SetFillColor(15);
        hs->Add(hData_unfolded_plussyst);
        //hs->SetMinimum( 0 );
        hs->SetMinimum( hData_unfolded->GetMinimum() - ( 0.3 * hData_unfolded->GetMaximum() ) > 0.15 ? hData_unfolded->GetMinimum() - ( 0.3 * hData_unfolded->GetMaximum() ) : 0 );
        if (observablename == "lep_azimuthal_asymmetry2" || observablename == "top_spin_correlation") hs->SetMaximum(1.35 * hData_unfolded->GetMaximum());
        else hs->SetMaximum(1.3 * hData_unfolded->GetMaximum());


        TCanvas *c_test = new TCanvas("c_final", "c_final", 500, 500);

        theoryProfileCorr->SetLineColor(kBlue);
        theoryProfileCorr->SetLineWidth(2);
        theoryProfileCorr->SetMarkerStyle(1);

        theoryProfileUnCorr->SetLineColor(kBlue);
        theoryProfileUnCorr->SetLineWidth(2);
        theoryProfileUnCorr->SetLineStyle(2);
        theoryProfileUnCorr->SetMarkerStyle(1);

        hs->Draw();
        hs->GetXaxis()->SetTitle(xaxislabel);
        hs->GetYaxis()->SetTitle("1/#sigma d#sigma/d(" + xaxislabel + ")");
        //hData_unfolded->GetXaxis()->SetTitle(xaxislabel);
        //hData_unfolded->GetYaxis()->SetTitle("1/#sigma d#sigma/d("+xaxislabel+")");
        //hData_unfolded->SetMinimum(0.0);
        //hData_unfolded->SetMaximum( 2.0* hData_unfolded->GetMaximum());
        hData_unfolded->SetMarkerStyle(23);
        hData_unfolded->SetMarkerSize(1);
        hData_unfolded->SetFillStyle(0);
        hData_unfolded->Draw("E same");
        hData_unfolded->SetLineWidth(lineWidth);
        denominatorM_nopTreweighting->SetLineWidth(lineWidth);
        denominatorM_nopTreweighting->SetLineColor(TColor::GetColorDark(kRed));
        denominatorM_nopTreweighting->SetFillStyle(0);
        hTrue->SetLineWidth(lineWidth);
        hTrue->SetLineColor(TColor::GetColorDark(kRed));
        //hTrue->SetFillColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillStyle(0);
        if (!draw_truth_before_pT_reweighting) hTrue->Draw("hist same");
        else denominatorM_nopTreweighting->Draw("hist same");
        hData_unfolded->Draw("EP same");
        if (observablename == "lep_azimuthal_asymmetry2" || observablename == "top_spin_correlation")
        {
            theoryProfileUnCorr->Draw("hist same");
            theoryProfileCorr->Draw("hist same");
        }

        //TLegend* leg1=new TLegend(0.55,0.62,0.9,0.838,NULL,"brNDC");
        TLegend *leg1 = new TLegend(0.58, 0.75, 0.9, 0.93, NULL, "brNDC");
        leg1->SetEntrySeparation(0.1);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->SetTextSize(0.032);
        leg1->AddEntry(hData_unfolded, "( Data - BG ) unfolded");
        leg1->AddEntry(hData_unfolded_plussyst,    "Syst. uncertainty", "F");
        leg1->AddEntry(hTrue,    "MC@NLO parton level", "L");

        leg1->Draw();

        if (observablename == "lep_azimuthal_asymmetry2" || observablename == "top_spin_correlation")
        {
            TLegend *leg2 = new TLegend(0.18, 0.745, 0.45, 0.88, NULL, "brNDC");
            leg2->SetEntrySeparation(0.5);
            leg2->SetFillColor(0);
            leg2->SetLineColor(0);
            leg2->SetBorderSize(0);
            leg2->SetFillStyle(0);
            leg2->SetTextSize(0.032);
            leg2->AddEntry(theoryProfileCorr,  "#splitline{W.Bernreuther & Z.G.Si}{(SM, #mu=^{}m_{t})}", "L");
            leg2->AddEntry(theoryProfileUnCorr,  "#splitline{W.Bernreuther & Z.G.Si}{(uncorrelated, #mu=^{}m_{t})}", "L");
            leg2->Draw();
        }





        TPaveText *pt1 = new TPaveText(0.175, 0.885, 0.41, 0.91, "brNDC");
        pt1->SetName("pt1name");
        pt1->SetBorderSize(0);
        pt1->SetFillStyle(0);

        TText *blah;
        //blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
        blah = pt1->AddText("CMS, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
        blah->SetTextSize(0.032);
        blah->SetTextAlign(11);
        pt1->Draw();

        c_test->SaveAs("finalplot_unfolded_" + acceptanceName + Region + ".pdf");
        c_test->SaveAs("finalplot_unfolded_" + acceptanceName + Region + ".C");
        c_test->SaveAs("finalplot_unfolded_" + acceptanceName + Region + ".root");

        ch_data->Delete();

        ch_top->Delete();

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg]->Delete();
        }

    }

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
