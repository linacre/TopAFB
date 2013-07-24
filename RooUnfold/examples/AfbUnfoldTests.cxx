
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom3.h"
#include "TH1D.h"
#include "TUnfold.h"
#include "TMatrixD.h"

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"
#include "src/RooUnfoldTUnfold.h"

#include "AfbFinalUnfold.h"

#include "tdrstyle.C"


//==============================================================================
// Global definitions
//==============================================================================

// 0=SVD, 1=TUnfold via RooUnfold, 2=TUnfold
int unfoldingType = 0;

TString Region = "";


Int_t kterm = 3;
Double_t tau = 1E-4;
Int_t nPseudos = 10000;
Int_t includeSys = 0;

Int_t lineWidth = 5;


// "Pull" or "Linearity"
void AfbUnfoldTests(Int_t iVar = 0, TString TestType = "Pull")
{
#ifdef __CINT__
    gSystem->Load("libRooUnfold");
#endif

    setTDRStyle();
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit();
    gStyle->SetOptStat("emr");
    cout.precision(3);

    ofstream myfile;
    myfile.open ("summary_PEtest.txt");
    cout.rdbuf(myfile.rdbuf());

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);


    Initialize1DBinning(iVar);
    bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;

    TH1D *hTrue_before = new TH1D ("trueBeforeScaling", "Truth",    nbins1D, xbins1D);
    TH1D *hMeas_before = new TH1D ("measBeforeScaling", "Measured", nbins1D, xbins1D);

    TH1D *hTrue_after = new TH1D ("trueAfterScaling", "Truth",    nbins1D, xbins1D);
    TH1D *hMeas_after = new TH1D ("measAfterScaling", "Measured", nbins1D, xbins1D);

    TH1D *hSmeared = new TH1D ("smeared", "Smeared", nbins1D, xbins1D);
    TH1D *hUnfolded = new TH1D ("unfolded", "Unfolded", nbins1D, xbins1D);


    TH1D *AfbPull = new TH1D("h_afbpull", "Pulls for Afb", 50, -5, 5);

    TH2D *hTrue_vs_Meas = new TH2D ("true_vs_meas", "True vs Measured", nbins1D, xbins1D, nbins1D, xbins1D);


    hTrue_before->Sumw2();
    hMeas_before->Sumw2();
    hTrue_after->Sumw2();
    hMeas_after->Sumw2();
    hSmeared->Sumw2();
    hUnfolded->Sumw2();

    AfbPull->Sumw2();
    hTrue_vs_Meas->Sumw2();

    TMatrixD m_unfoldE(nbins1D, nbins1D);


    TH1F *h_pulls[nbins1D];
    TH1F *h_resd[nbins1D];
    for (int i = 0; i < nbins1D; i++)
    {
        TString name = "h_pull_";
        name += i;
        h_pulls[i] = new TH1F(name, name, 100, -5.0, 5.0);
        name = "h_resd_";
        name += i;
        h_resd[i] = new TH1F(name, name, 20, -1, 1);
    }


    TFile *file = new TFile("../ttdil.root");
    TTree *evtree = (TTree *) file->Get("tree");
    Int_t entries = (Int_t)evtree->GetEntries();
    cout << "RESPONSE: Number of Entries: " << entries << endl;

    Float_t observable, observable_gen, ttmass, ttRapidity, tmass;
    Float_t observableMinus, observableMinus_gen;
    Double_t weight;
    Int_t Nsolns;

    evtree->SetBranchAddress(observablename,    &observable);
    evtree->SetBranchAddress(observablename + "_gen", &observable_gen);
    if (observablename == "lep_azimuthal_asymmetry2") evtree->SetBranchAddress("lep_azimuthal_asymmetry_gen2", &observable_gen);
    if ( combineLepMinus ) evtree->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
    if ( combineLepMinus ) evtree->SetBranchAddress("lepMinus_costheta_cms_gen",    &observableMinus_gen);
    evtree->SetBranchAddress("weight", &weight);
    evtree->SetBranchAddress("Nsolns", &Nsolns);
    evtree->SetBranchAddress("tt_mass", &ttmass);
    evtree->SetBranchAddress("ttRapidity", &ttRapidity);
    evtree->SetBranchAddress("t_mass", &tmass);

    Float_t slope = 0.0;
    const int Nlin = 7;
    Float_t A_gen[Nlin], Aerr_gen[Nlin], A_unf[Nlin], Aerr_unf[Nlin], A_meas[Nlin], Aerr_meas[Nlin];
    Float_t A_pull[Nlin], A_pullwidth[Nlin], Aerr_pull[Nlin], Aerr_pullwidth[Nlin];

    for (int k = 0; k < Nlin; k++)
    {

        if ((TestType == "Pull") && (k == 1)) break;

        slope = -0.3 + 0.1 * k;

        cout << "slope =" << slope << "\n";

        hTrue_before->Reset();
        hMeas_before->Reset();
        hTrue_after->Reset();
        hMeas_after->Reset();
        hTrue_vs_Meas->Reset();
        AfbPull->Reset();


        if (observable < xmin) observable = xmin;
        if (observable > xmax) observable = xmax;

        if (observable_gen < xmin) observable_gen = xmin;
        if (observable_gen > xmax) observable_gen = xmax;

        if ( combineLepMinus )
        {

            if (observableMinus < xmin) observableMinus = xmin;
            if (observableMinus > xmax) observableMinus = xmax;

            if (observableMinus_gen < xmin) observableMinus_gen = xmin;
            if (observableMinus_gen > xmax) observableMinus_gen = xmax;
        }

        double asym_centre = (xmax + xmin) / 2.;

        for (Int_t i = 0; i < entries; i++)
        {
            evtree->GetEntry(i);
            double orig_weight = weight;
            //if(i % 10000 == 0) cout<<i<<" "<<ch_top->GetEntries()<<endl;

            if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
            {
                fillUnderOverFlow(hMeas_before, observable, weight, Nsolns);
                fillUnderOverFlow(hTrue_before, observable_gen, weight, Nsolns);
                fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                if (TestType == "Linearity") weight = weight * (1.0 + slope * (observable_gen - asym_centre) );
                fillUnderOverFlow(hMeas_after, observable, weight, Nsolns);
                fillUnderOverFlow(hTrue_after, observable_gen, weight, Nsolns);

            }
            else if ( ttmass > 0 )
            {
                fillUnderOverFlow(hMeas_before, observable, weight, Nsolns);
                fillUnderOverFlow(hTrue_before, observable_gen, weight, Nsolns);
                fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                if ( combineLepMinus )
                {
                    fillUnderOverFlow(hMeas_before, observableMinus, weight, Nsolns);
                    fillUnderOverFlow(hTrue_before, observableMinus_gen, weight, Nsolns);
                    fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
                }
                if (TestType == "Linearity") weight = weight * (1.0 + slope * (observable_gen - asym_centre) );
                fillUnderOverFlow(hMeas_after, observable, weight, Nsolns);
                fillUnderOverFlow(hTrue_after, observable_gen, weight, Nsolns);
                if ( combineLepMinus )
                {
                    if (TestType == "Linearity") weight = orig_weight * (1.0 + slope * (observableMinus_gen - asym_centre) );
                    fillUnderOverFlow(hMeas_after, observableMinus, weight, Nsolns);
                    fillUnderOverFlow(hTrue_after, observableMinus_gen, weight, Nsolns);
                }

            }
        }

        RooUnfoldResponse response (hMeas_before, hTrue_before, hTrue_vs_Meas);

        //scale to keep total yield constant
        hMeas_after->Scale( hMeas_before->Integral() / hMeas_after->Integral() );
        hTrue_after->Scale( hTrue_before->Integral() / hTrue_after->Integral() );


        TFile *accfile = new TFile("../acceptance/mcnlo/accept_" + acceptanceName + ".root");
        TH1D *acceptM = (TH1D *) accfile->Get("accept_" + acceptanceName);
        acceptM->Scale(1.0 / acceptM->Integral());

        for (Int_t accbin = 1; accbin <= nbins1D; accbin++)
        {

            if (acceptM->GetBinContent(accbin) != 0)
            {
                hTrue_after->SetBinContent(accbin, hTrue_after->GetBinContent(accbin) * 1.0 / acceptM->GetBinContent(accbin));
                hTrue_after->SetBinError  (accbin, hTrue_after->GetBinError  (accbin) * 1.0 / acceptM->GetBinContent(accbin));

            }
        }


        Float_t Afb, AfbErr;

        GetAfb(hTrue_after, Afb, AfbErr);
        A_gen[k] = Afb;
        Aerr_gen[k] = 0.0;
        cout << " True after re-weighting   : " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hMeas_after, Afb, AfbErr);
        A_meas[k] = Afb;
        Aerr_meas[k] = AfbErr;
        cout << " Measured after re-weighting   : " << Afb << " +/-  " << AfbErr << "\n";


        // Now do the pseudos

        Float_t trialAsym = 0.0, SumAsym = 0.0, SumErrAsym = 0.0;


        if (nPseudos > 1)
        {

            for (int i = 0; i < nPseudos; i++)
            {

                for (int j = 1; j < hMeas_after->GetNbinsX() + 1; j++)
                {
                    double fluct = random->Poisson(hMeas_after->GetBinContent(j));
                    //      double fluct = hMeas_after->GetBinContent(j);
                    hSmeared->SetBinError(j, sqrt(fluct));
                    hSmeared->SetBinContent(j, fluct);
                }


                if (unfoldingType == 0)
                {
                    RooUnfoldSvd unfold_svd (&response, hSmeared, kterm);
                    unfold_svd.Setup(&response, hSmeared);
                    //      unfold_svd.IncludeSystematics(includeSys);
                    hUnfolded = (TH1D *) unfold_svd.Hreco();
                    m_unfoldE = unfold_svd.Ereco();
                }
                else if (unfoldingType == 1)
                {
                    RooUnfoldTUnfold unfold_rooTUnfold (&response, hSmeared, TUnfold::kRegModeCurvature);
                    unfold_rooTUnfold.Setup(&response, hSmeared);
                    unfold_rooTUnfold.FixTau(tau);
                    //  unfold_rooTUnfold.IncludeSystematics(includeSys);
                    hUnfolded = (TH1D *) unfold_rooTUnfold.Hreco();
                    m_unfoldE = unfold_rooTUnfold.Ereco();
                }
                else if (unfoldingType == 2)
                {
                    TUnfold unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature);
                    //  Double_t biasScale=5.0;
                    //  unfold_TUnfold.SetBias(hTrue_before);
                    unfold_TUnfold.SetInput(hSmeared);
                    unfold_TUnfold.DoUnfold(tau);
                    unfold_TUnfold.GetOutput(hUnfolded);


                    TH2D *ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
                    for (Int_t cmi = 0; cmi < nbins1D; cmi++)
                    {
                        for (Int_t cmj = 0; cmj < nbins1D; cmj++)
                        {
                            m_unfoldE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
                        }
                    }
                }
                else cout << "Unfolding TYPE not Specified" << "\n";


                for (int l = 0; l < nbins1D; l++)
                {
                    for (int j = 0; j < nbins1D; j++)
                    {
                        double corr = 1.0 / ( acceptM->GetBinContent(l + 1) * acceptM->GetBinContent(j + 1) );
                        //corr = corr * pow(xsection / dataIntegral,2) ;
                        m_unfoldE(l, j) = m_unfoldE(l, j) * corr;
                    }
                }


                for (Int_t accbin = 1; accbin <= nbins1D; accbin++)
                {

                    if (acceptM->GetBinContent(accbin) != 0)
                    {
                        hUnfolded->SetBinContent(accbin, hUnfolded->GetBinContent(accbin) * 1.0 / acceptM->GetBinContent(accbin));
                        hUnfolded->SetBinError  (accbin, hUnfolded->GetBinError  (accbin) * 1.0 / acceptM->GetBinContent(accbin));

                    }
                }

                GetCorrectedAfb(hUnfolded, m_unfoldE, Afb, AfbErr);

                AfbPull -> Fill( (Afb - A_gen[k])  / AfbErr );

                SumAsym + = Afb;
                SumErrAsym + = AfbErr;


                for (int j = 0; j < nbins1D; j++)
                {
                    double pull = (hUnfolded->GetBinContent(j + 1) - hTrue_after->GetBinContent(j + 1)) / hUnfolded->GetBinError(j + 1);
                    h_pulls[j]->Fill(pull);
                    double resd = (hUnfolded->GetBinContent(j + 1) - hTrue_after->GetBinContent(j + 1)) / hTrue_after->GetBinContent(j + 1);
                    h_resd[j]->Fill(resd);
                }
            }

        }

        cout << "Average Asymmetry =" << SumAsym / nPseudos << " +/-  " << SumErrAsym / (nPseudos) << "\n";
        A_unf[k] = SumAsym / nPseudos;
        Aerr_unf[k] = SumErrAsym / nPseudos;

        A_pull[k] = AfbPull->GetMean();
        Aerr_pull[k] = AfbPull->GetMeanError();
        A_pullwidth[k] = AfbPull->GetRMS();
        Aerr_pullwidth[k] = AfbPull->GetRMSError();

    }

    TGraphErrors *Asym2D_TrueUnf = new TGraphErrors (Nlin, A_gen, A_unf, Aerr_gen, Aerr_unf);

    TGraphErrors *Asym2D_TrueMeas = new TGraphErrors (Nlin, A_gen, A_meas, Aerr_gen, Aerr_meas);

    TGraphErrors *Asym2D_PullWidth = new TGraphErrors (Nlin, A_gen, A_pullwidth, Aerr_gen, Aerr_pullwidth);

    TGraphErrors *Asym2D_Pull = new TGraphErrors (Nlin, A_gen, A_pull, Aerr_gen, Aerr_pull);

    if ((TestType == "Linearity"))
    {
        TCanvas *c_ttbar = new TCanvas("c_ttbar", "c_ttbar", 500, 500);
        Asym2D_TrueUnf->SetTitle(asymlabel);
        Asym2D_TrueUnf->SetMarkerStyle(23);
        Asym2D_TrueUnf->SetMarkerColor(kBlack);
        Asym2D_TrueUnf->SetMarkerSize(1.0);
        Asym2D_TrueUnf->GetXaxis()->SetTitle(asymlabel + " (true)");
        Asym2D_TrueUnf->GetYaxis()->SetTitle(asymlabel + " (unfolded)");
        Asym2D_TrueUnf->Draw("AP same");
        Asym2D_TrueUnf->Fit("pol1");
        c_ttbar->SaveAs(acceptanceName + "_LinearityCheck.pdf");
        c_ttbar->SaveAs(acceptanceName + "_LinearityCheck.C");

        TCanvas *c_Pull_lin = new TCanvas("c_Pull_lin", "c_Pull_lin", 500, 500);
        Asym2D_Pull->SetTitle(asymlabel);
        Asym2D_Pull->SetMarkerStyle(23);
        Asym2D_Pull->SetMarkerColor(kBlack);
        Asym2D_Pull->SetMarkerSize(1.0);
        Asym2D_Pull->GetXaxis()->SetTitle(asymlabel + " (true)");
        Asym2D_Pull->GetYaxis()->SetTitle(asymlabel + " pull");
        Asym2D_Pull->Draw("AP same");
        Asym2D_Pull->Fit("pol1");
        c_Pull_lin->SaveAs(acceptanceName + "_LinearityCheck_Pull.pdf");
        c_Pull_lin->SaveAs(acceptanceName + "_LinearityCheck_Pull.C");

        TCanvas *c_PullWidth_lin = new TCanvas("c_PullWidth_lin", "c_PullWidth_lin", 500, 500);
        Asym2D_PullWidth->SetTitle(asymlabel);
        Asym2D_PullWidth->SetMarkerStyle(23);
        Asym2D_PullWidth->SetMarkerColor(kBlack);
        Asym2D_PullWidth->SetMarkerSize(1.0);
        Asym2D_PullWidth->GetXaxis()->SetTitle(asymlabel + " (true)");
        Asym2D_PullWidth->GetYaxis()->SetTitle(asymlabel + " pull width");
        Asym2D_PullWidth->Draw("AP same");
        Asym2D_PullWidth->Fit("pol1");
        c_PullWidth_lin->SaveAs(acceptanceName + "_LinearityCheck_PullWidth.pdf");
        c_PullWidth_lin->SaveAs(acceptanceName + "_LinearityCheck_PullWidth.C");
    }
    else
    {
        TCanvas *c_pull = new TCanvas("c_pull", "c_pull", 500, 500);
        AfbPull->SetMarkerStyle(23);
        AfbPull->SetMarkerColor(kBlack);
        AfbPull->SetMarkerSize(1.0);
        AfbPull->GetXaxis()->SetTitle(asymlabel + " pull");
        AfbPull->GetYaxis()->SetTitle("Number of PEs / 0.2");
        AfbPull ->Fit("gaus");
        AfbPull ->Draw();
        c_pull->SaveAs(acceptanceName + "_Pull.pdf");
        c_pull->SaveAs(acceptanceName + "_Pull.C");


        TFile *plots = new TFile(acceptanceName + "_plots.root", "RECREATE");
        for (int i = 0; i < nbins1D; i++)
        {
            h_pulls[i] ->Write();
            h_resd[i] ->Write();
        }
        AfbPull ->Write();
    }

    myfile.close();

}

#ifndef __CINT__
int main ()
{
    AfbUnfoldLinearityTest();    // Main program when run stand-alone
    return 0;
}
#endif
