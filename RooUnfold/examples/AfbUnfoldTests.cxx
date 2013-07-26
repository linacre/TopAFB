
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom3.h"
#include "TH1D.h"
#include "TUnfold.h"
#include "TMatrixD.h"
#include <vector>

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"
#include "src/RooUnfoldTUnfold.h"

#include "AfbFinalUnfold.h"

#include "tdrstyle.C"
#include "../CommonFunctions.C"


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



//TestType: "Pull" or "Linearity"
//slopeOption: 0 = continuous reweighting, 1 = 6-binned reweighting, 2 = 2-binned reweighting (i.e. sign(var))
void AfbUnfoldTests(Int_t iVar = 0, TString TestType = "Pull", Int_t slopeOption = 0)
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
    //cout.rdbuf(myfile.rdbuf());
    ofstream second_output_file;
    second_output_file.open("summary_PEtest_formated.txt");

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);


    Initialize1DBinning(iVar);
    bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;
    const int nDiff = nbins1D / 2;

    TH1D *hEmpty = new TH1D ("Empty", "Empty",    nbins1D, xbins1D);

    TH1D *hTrue_before = new TH1D ("trueBeforeScaling", "Truth",    nbins1D, xbins1D);
    TH1D *hMeas_before = new TH1D ("measBeforeScaling", "Measured", nbins1D, xbins1D);

    TH1D *hTrue_after = new TH1D ("trueAfterScaling", "Truth",    nbins1D, xbins1D);
    TH1D *hMeas_after = new TH1D ("measAfterScaling", "Measured", nbins1D, xbins1D);

    TH1D *hSmeared = new TH1D ("smeared", "Smeared", nbins1D, xbins1D);
    TH1D *hUnfolded = new TH1D ("unfolded", "Unfolded", nbins1D, xbins1D);

    double pullMax = 5;
    int pullBins = 50;
    if (TestType == "Linearity") pullMax = 100;
    if (TestType == "Linearity") pullBins = 1000;

    TH1D *AfbPull[nDiff + 1];

    for (int iD = 0; iD < nDiff + 1; ++iD)
    {
        AfbPull[iD] = new TH1D("h_afbpull" + iD, "Pulls for Afb" + iD, pullBins, -pullMax, pullMax);
        AfbPull[iD]->Sumw2();
    }


    TH2D *hTrue_vs_Meas = new TH2D ("true_vs_meas", "True vs Measured", nbins1D, xbins1D, nbins1D, xbins1D);


    hTrue_before->Sumw2();
    hMeas_before->Sumw2();
    hTrue_after->Sumw2();
    hMeas_after->Sumw2();
    hSmeared->Sumw2();
    hUnfolded->Sumw2();

    hTrue_vs_Meas->Sumw2();

    TMatrixD m_unfoldE(nbins1D, nbins1D);


    TH1F *h_pulls[nbins1D];
    TH1F *h_resd[nbins1D];
    for (int i = 0; i < nbins1D; i++)
    {
        TString name = "h_pull_";
        name += i;
        h_pulls[i] = new TH1F(name, name, 50, -5.0, 5.0);
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
    //Float_t A_gen[Nlin], Aerr_gen[Nlin], A_unf[Nlin], Aerr_unf[Nlin], A_meas[Nlin], Aerr_meas[Nlin];
    //Float_t A_pull[Nlin], A_pullwidth[Nlin], Aerr_pull[Nlin], Aerr_pullwidth[Nlin];

    Float_t  A_meas[Nlin], Aerr_meas[Nlin];
    vector<vector<Float_t>> A_gen, Aerr_gen, A_unf, Aerr_unf;
    vector<vector<Float_t>> A_pull, A_pullwidth, Aerr_pull, Aerr_pullwidth;

    A_gen.clear();
    Aerr_gen.clear();
    A_unf.clear();
    Aerr_unf.clear();
    A_pull.clear();
    A_pullwidth.clear();
    Aerr_pull.clear();
    Aerr_pullwidth.clear();


    TH1D *hTrue_after_array[Nlin];
    TH1D *hMeas_after_array[Nlin];

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
        for (int iD = 0; iD < nDiff + 1; ++iD)
        {
            AfbPull[iD]->Reset();
        }


        double asym_centre = (xmax + xmin) / 2.;

        for (Int_t i = 0; i < entries; i++)
        {
            evtree->GetEntry(i);
            double orig_weight = weight;

            if (slopeOption == 1)
            {
                //fix the observable values to the bin centres so the acceptance function is unaffected by any reweighting
                observable =  hEmpty->GetBinCenter( hEmpty->FindBin( observable ) );
                observable_gen =  hEmpty->GetBinCenter( hEmpty->FindBin( observable_gen ) );
                if ( combineLepMinus )
                {
                    observableMinus =  hEmpty->GetBinCenter( hEmpty->FindBin( observableMinus ) );
                    observableMinus_gen =  hEmpty->GetBinCenter( hEmpty->FindBin( observableMinus_gen ) );
                }
            }

            //if(i % 10000 == 0) cout<<i<<" "<<ch_top->GetEntries()<<endl;

            if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
            {
                fillUnderOverFlow(hMeas_before, observable, weight, Nsolns);
                fillUnderOverFlow(hTrue_before, observable_gen, weight, Nsolns);
                fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                if (TestType == "Linearity" && slopeOption != 2) weight = weight * (1.0 + slope * (observable_gen - asym_centre) );
                if (TestType == "Linearity" && slopeOption == 2) weight = weight * (1.0 + 0.5 * slope * sign(observable_gen - asym_centre) );
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
                if (TestType == "Linearity" && slopeOption != 2) weight = weight * (1.0 + slope * (observable_gen - asym_centre) );
                if (TestType == "Linearity" && slopeOption == 2) weight = weight * (1.0 + 0.5 * slope * sign(observable_gen - asym_centre) );
                fillUnderOverFlow(hMeas_after, observable, weight, Nsolns);
                fillUnderOverFlow(hTrue_after, observable_gen, weight, Nsolns);
                if ( combineLepMinus )
                {
                    if (TestType == "Linearity" && slopeOption != 2) weight = orig_weight * (1.0 + slope * (observableMinus_gen - asym_centre) );
                    if (TestType == "Linearity" && slopeOption == 2) weight = orig_weight * (1.0 + 0.5 * slope * sign(observableMinus_gen - asym_centre) );
                    fillUnderOverFlow(hMeas_after, observableMinus, weight, Nsolns);
                    fillUnderOverFlow(hTrue_after, observableMinus_gen, weight, Nsolns);
                }

            }
        }

        RooUnfoldResponse response (hMeas_before, hTrue_before, hTrue_vs_Meas);

        //scale to keep total yield constant
        hMeas_after->Scale( hMeas_before->Integral() / hMeas_after->Integral() );
        hTrue_after->Scale( hTrue_before->Integral() / hTrue_after->Integral() );

        hTrue_after_array[k] = (TH1D *) hTrue_after->Clone();
        hMeas_after_array[k] = (TH1D *) hMeas_after->Clone();


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
        Float_t A_gen_k = Afb;
        Float_t Aerr_gen_k = 0.0;
        cout << " True after re-weighting   : " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hMeas_after, Afb, AfbErr);
        A_meas[k] = Afb;
        Aerr_meas[k] = AfbErr;
        cout << " Measured after re-weighting   : " << Afb << " +/-  " << AfbErr << "\n";


        // Now do the pseudos

        Float_t trialAsym = 0.0;
        vector <Float_t> SumAsym, SumErrAsym, SumTrueAsym, SumTrueErrAsym;
        SumAsym.clear();
        SumErrAsym.clear();
        SumTrueAsym.clear();
        SumTrueErrAsym.clear();

        for (int iD = 1; iD < nDiff + 1; ++iD)
        {
            SumAsym.push_back(0);
            SumErrAsym.push_back(0);
            SumTrueAsym.push_back(0);
            SumTrueErrAsym.push_back(0);
        }


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

                vector<double> afb2D;
                vector<double> afb2Derr;
                afb2D.clear();
                afb2Derr.clear();
                GetCorrectedAfbBinByBin(hUnfolded, m_unfoldE, afb2D, afb2Derr, second_output_file);

                vector<double> afbtrue2D;
                vector<double> afbtrue2Derr;
                afbtrue2D.clear();
                afbtrue2Derr.clear();
                GetCorrectedAfbBinByBin(hTrue_after, m_unfoldE, afbtrue2D, afbtrue2Derr, second_output_file);
                //true errors are much smaller (from denominator)
                afbtrue2Derr[0] = 0.0;
                afbtrue2Derr[1] = 0.0;
                afbtrue2Derr[2] = 0.0;

                AfbPull[0] -> Fill( (Afb - A_gen_k)  / AfbErr );
                SumAsym[0] + = Afb;
                SumErrAsym[0] + = AfbErr;
                SumTrueAsym[0] + = A_gen_k;
                SumTrueErrAsym[0] + = Aerr_gen_k;


                for (int iD = 1; iD < nDiff + 1; ++iD)
                {
                    AfbPull[iD]->Fill( (afb2D[iD - 1] - afbtrue2D[iD - 1])  / afb2Derr[iD - 1] );
                    SumAsym[iD] + = afb2D[iD - 1];
                    SumErrAsym[iD] + = afb2Derr[iD - 1];
                    SumTrueAsym[iD] + = afbtrue2D[iD - 1];
                    SumTrueErrAsym[iD] + = afbtrue2Derr[iD - 1];
                }



                for (int j = 0; j < nbins1D; j++)
                {
                    double pull = (hUnfolded->GetBinContent(j + 1) - hTrue_after->GetBinContent(j + 1)) / hUnfolded->GetBinError(j + 1);
                    h_pulls[j]->Fill(pull);
                    double resd = (hUnfolded->GetBinContent(j + 1) - hTrue_after->GetBinContent(j + 1)) / hTrue_after->GetBinContent(j + 1);
                    h_resd[j]->Fill(resd);
                }
            }

        }

        //cout << "Average Asymmetry =" << SumAsym / nPseudos << " +/-  " << SumErrAsym / (nPseudos) << "\n";

        vector<Float_t> temp_A_pull;
        vector<Float_t> temp_Aerr_pull;
        vector<Float_t> temp_A_pullwidth;
        vector<Float_t> temp_Aerr_pullwidth;
        temp_A_pull.clear();
        temp_Aerr_pull.clear();
        temp_A_pullwidth.clear();
        temp_Aerr_pullwidth.clear();

        for (int iD = 0; iD < nDiff + 1; ++iD)
        {
            temp_A_pull.push_back( AfbPull[iD]->GetMean() );
            temp_Aerr_pull.push_back( AfbPull[iD]->GetMeanError() );
            temp_A_pullwidth.push_back( AfbPull[iD]->GetRMS() );
            temp_Aerr_pullwidth.push_back( AfbPull[iD]->GetRMSError() );
            SumAsym[iD] /= nPseudos;
            SumErrAsym[iD] /= nPseudos;
            SumTrueAsym[iD] /= nPseudos;
            SumTrueErrAsym[iD] /= nPseudos;
        }

        A_unf.push_back( SumAsym );
        Aerr_unf.push_back( SumErrAsym );
        A_gen.push_back( SumTrueAsym );
        Aerr_gen.push_back( SumTrueErrAsym );

        A_pull.push_back( temp_A_pull );
        Aerr_pull.push_back( temp_Aerr_pull );
        A_pullwidth.push_back( temp_A_pullwidth );
        Aerr_pullwidth.push_back( temp_Aerr_pullwidth );

    }




    if ((TestType == "Linearity"))
    {
        //vector<vector<Float_t>> transposed_A_gen, transposed_Aerr_gen, transposed_A_unf, transposed_Aerr_unf;
        //vector<vector<Float_t>> transposed_A_pull, transposed_A_pullwidth, transposed_Aerr_pull, transposed_Aerr_pullwidth;

        //vector<Float_t> transposed_A_gen, transposed_Aerr_gen, transposed_A_unf, transposed_Aerr_unf;
        //vector<Float_t> transposed_A_pull, transposed_A_pullwidth, transposed_Aerr_pull, transposed_Aerr_pullwidth;
        //
        //transposed_A_gen.clear();
        //transposed_Aerr_gen.clear();
        //transposed_A_unf.clear();
        //transposed_Aerr_unf.clear();
        //transposed_A_pull.clear();
        //transposed_A_pullwidth.clear();
        //transposed_Aerr_pull.clear();
        //transposed_Aerr_pullwidth.clear();

        cout << " Got to line: " << __LINE__ << endl;

        //for (int iD = 0; iD < nDiff + 1; ++iD)
        //{

        /*
                vector<Float_t> tmp_A_gen, tmp_Aerr_gen, tmp_A_unf, tmp_Aerr_unf;
                vector<Float_t> tmp_A_pull, tmp_A_pullwidth, tmp_Aerr_pull, tmp_Aerr_pullwidth;

                tmp_A_gen.clear();
                tmp_Aerr_gen.clear();
                tmp_A_unf.clear();
                tmp_Aerr_unf.clear();
                tmp_A_pull.clear();
                tmp_A_pullwidth.clear();
                tmp_Aerr_pull.clear();
                tmp_Aerr_pullwidth.clear();
                */

        Float_t inclusive_A_gen[Nlin], inclusive_Aerr_gen[Nlin], inclusive_A_unf[Nlin], inclusive_Aerr_unf[Nlin];
        Float_t inclusive_A_pull[Nlin], inclusive_A_pullwidth[Nlin], inclusive_Aerr_pull[Nlin], inclusive_Aerr_pullwidth[Nlin];

        Float_t bin1_A_gen[Nlin], bin1_Aerr_gen[Nlin], bin1_A_unf[Nlin], bin1_Aerr_unf[Nlin];
        Float_t bin1_A_pull[Nlin], bin1_A_pullwidth[Nlin], bin1_Aerr_pull[Nlin], bin1_Aerr_pullwidth[Nlin];

        Float_t bin2_A_gen[Nlin], bin2_Aerr_gen[Nlin], bin2_A_unf[Nlin], bin2_Aerr_unf[Nlin];
        Float_t bin2_A_pull[Nlin], bin2_A_pullwidth[Nlin], bin2_Aerr_pull[Nlin], bin2_Aerr_pullwidth[Nlin];

        Float_t bin3_A_gen[Nlin], bin3_Aerr_gen[Nlin], bin3_A_unf[Nlin], bin3_Aerr_unf[Nlin];
        Float_t bin3_A_pull[Nlin], bin3_A_pullwidth[Nlin], bin3_Aerr_pull[Nlin], bin3_Aerr_pullwidth[Nlin];

        for (int k = 0; k < Nlin; k++)
        {

            inclusive_A_gen[k] = A_gen[k][0];
            inclusive_Aerr_gen[k] = Aerr_gen[k][0];
            inclusive_A_unf[k] = A_unf[k][0];
            inclusive_Aerr_unf[k] = Aerr_unf[k][0];
            inclusive_A_pull[k] = A_pull[k][0];
            inclusive_A_pullwidth[k] = A_pullwidth[k][0];
            inclusive_Aerr_pull[k] = Aerr_pull[k][0];
            inclusive_Aerr_pullwidth[k] = Aerr_pullwidth[k][0];

            bin1_A_gen[k] = A_gen[k][1];
            bin1_Aerr_gen[k] = Aerr_gen[k][1];
            bin1_A_unf[k] = A_unf[k][1];
            bin1_Aerr_unf[k] = Aerr_unf[k][1];
            bin1_A_pull[k] = A_pull[k][1];
            bin1_A_pullwidth[k] = A_pullwidth[k][1];
            bin1_Aerr_pull[k] = Aerr_pull[k][1];
            bin1_Aerr_pullwidth[k] = Aerr_pullwidth[k][1];

            bin2_A_gen[k] = A_gen[k][2];
            bin2_Aerr_gen[k] = Aerr_gen[k][2];
            bin2_A_unf[k] = A_unf[k][2];
            bin2_Aerr_unf[k] = Aerr_unf[k][2];
            bin2_A_pull[k] = A_pull[k][2];
            bin2_A_pullwidth[k] = A_pullwidth[k][2];
            bin2_Aerr_pull[k] = Aerr_pull[k][2];
            bin2_Aerr_pullwidth[k] = Aerr_pullwidth[k][2];

            bin3_A_gen[k] = A_gen[k][3];
            bin3_Aerr_gen[k] = Aerr_gen[k][3];
            bin3_A_unf[k] = A_unf[k][3];
            bin3_Aerr_unf[k] = Aerr_unf[k][3];
            bin3_A_pull[k] = A_pull[k][3];
            bin3_A_pullwidth[k] = A_pullwidth[k][3];
            bin3_Aerr_pull[k] = Aerr_pull[k][3];
            bin3_Aerr_pullwidth[k] = Aerr_pullwidth[k][3];


            /*
                        tmp_A_gen.push_back( A_gen[k][iD] );
                        tmp_Aerr_gen.push_back( Aerr_gen[k][iD] );
                        tmp_A_unf.push_back( A_unf[k][iD] );
                        tmp_Aerr_unf.push_back( Aerr_unf[k][iD] );
                        tmp_A_pull.push_back( A_pull[k][iD] );
                        tmp_A_pullwidth.push_back( A_pullwidth[k][iD] );
                        tmp_Aerr_pull.push_back( Aerr_pull[k][iD] );
                        tmp_Aerr_pullwidth.push_back( Aerr_pullwidth[k][iD] );
            */

        }

        //transposed_A_gen.push_back( tmp_A_gen );
        //transposed_Aerr_gen.push_back( tmp_Aerr_gen );
        //transposed_A_unf.push_back( tmp_A_unf );
        //transposed_Aerr_unf.push_back( tmp_Aerr_unf );
        //transposed_A_pull.push_back( tmp_A_pull );
        //transposed_A_pullwidth.push_back( tmp_A_pullwidth );
        //transposed_Aerr_pull.push_back( tmp_Aerr_pull );
        //transposed_Aerr_pullwidth.push_back( tmp_Aerr_pullwidth );

        //}



        TGraphErrors *Asym2D_TrueUnf = new TGraphErrors (Nlin, inclusive_A_gen, inclusive_A_unf, inclusive_Aerr_gen, inclusive_Aerr_unf);

        //TGraphErrors *Asym2D_TrueMeas = new TGraphErrors (Nlin, inclusive_A_gen, inclusive_A_meas, inclusive_Aerr_gen, inclusive_Aerr_meas);

        TGraphErrors *Asym2D_PullWidth = new TGraphErrors (Nlin, inclusive_A_gen, inclusive_A_pullwidth, inclusive_Aerr_gen, inclusive_Aerr_pullwidth);

        TGraphErrors *Asym2D_Pull = new TGraphErrors (Nlin, inclusive_A_gen, inclusive_A_pull, inclusive_Aerr_gen, inclusive_Aerr_pull);





        TGraphErrors *Asym2D_TrueUnfbin1 = new TGraphErrors (Nlin, bin1_A_gen, bin1_A_unf, bin1_Aerr_gen, bin1_Aerr_unf);
        TGraphErrors *Asym2D_PullWidthbin1 = new TGraphErrors (Nlin, bin1_A_gen, bin1_A_pullwidth, bin1_Aerr_gen, bin1_Aerr_pullwidth);
        TGraphErrors *Asym2D_Pullbin1 = new TGraphErrors (Nlin, bin1_A_gen, bin1_A_pull, bin1_Aerr_gen, bin1_Aerr_pull);

        TGraphErrors *Asym2D_TrueUnfbin2 = new TGraphErrors (Nlin, bin2_A_gen, bin2_A_unf, bin2_Aerr_gen, bin2_Aerr_unf);
        TGraphErrors *Asym2D_PullWidthbin2 = new TGraphErrors (Nlin, bin2_A_gen, bin2_A_pullwidth, bin2_Aerr_gen, bin2_Aerr_pullwidth);
        TGraphErrors *Asym2D_Pullbin2 = new TGraphErrors (Nlin, bin2_A_gen, bin2_A_pull, bin2_Aerr_gen, bin2_Aerr_pull);

        TGraphErrors *Asym2D_TrueUnfbin3 = new TGraphErrors (Nlin, bin3_A_gen, bin3_A_unf, bin3_Aerr_gen, bin3_Aerr_unf);
        TGraphErrors *Asym2D_PullWidthbin3 = new TGraphErrors (Nlin, bin3_A_gen, bin3_A_pullwidth, bin3_Aerr_gen, bin3_Aerr_pullwidth);
        TGraphErrors *Asym2D_Pullbin3 = new TGraphErrors (Nlin, bin3_A_gen, bin3_A_pull, bin3_Aerr_gen, bin3_Aerr_pull);



        TCanvas *c_ttbar = new TCanvas("c_ttbar", "c_ttbar", 500, 500);
        c_ttbar->Divide(2, 2);
        c_ttbar->cd(1);
        Asym2D_TrueUnf->SetTitle(asymlabel);
        Asym2D_TrueUnf->SetMarkerStyle(23);
        Asym2D_TrueUnf->SetMarkerColor(kBlack);
        Asym2D_TrueUnf->SetMarkerSize(0.6);
        Asym2D_TrueUnf->GetXaxis()->SetTitle(asymlabel + " inclusive (true)");
        Asym2D_TrueUnf->GetYaxis()->SetTitle(asymlabel + " inclusive (unfolded)");
        Asym2D_TrueUnf->Draw("AP same");
        Asym2D_TrueUnf->Fit("pol1");

        c_ttbar->cd(2);
        Asym2D_TrueUnfbin1->SetTitle(asymlabel);
        Asym2D_TrueUnfbin1->SetMarkerStyle(23);
        Asym2D_TrueUnfbin1->SetMarkerColor(kBlue);
        Asym2D_TrueUnfbin1->SetMarkerSize(0.6);
        Asym2D_TrueUnfbin1->GetXaxis()->SetTitle(asymlabel + " bin 1 (true)");
        Asym2D_TrueUnfbin1->GetYaxis()->SetTitle(asymlabel + " bin 1 (unfolded)");
        Asym2D_TrueUnfbin1->Draw("AP same");
        Asym2D_TrueUnfbin1->Fit("pol1");

        c_ttbar->cd(3);
        Asym2D_TrueUnfbin2->SetTitle(asymlabel);
        Asym2D_TrueUnfbin2->SetMarkerStyle(23);
        Asym2D_TrueUnfbin2->SetMarkerColor(kBlue);
        Asym2D_TrueUnfbin2->SetMarkerSize(0.6);
        Asym2D_TrueUnfbin2->GetXaxis()->SetTitle(asymlabel + " bin 2 (true)");
        Asym2D_TrueUnfbin2->GetYaxis()->SetTitle(asymlabel + " bin 2 (unfolded)");
        Asym2D_TrueUnfbin2->Draw("AP same");
        Asym2D_TrueUnfbin2->Fit("pol1");

        c_ttbar->cd(4);
        Asym2D_TrueUnfbin3->SetTitle(asymlabel);
        Asym2D_TrueUnfbin3->SetMarkerStyle(23);
        Asym2D_TrueUnfbin3->SetMarkerColor(kBlue);
        Asym2D_TrueUnfbin3->SetMarkerSize(0.6);
        Asym2D_TrueUnfbin3->GetXaxis()->SetTitle(asymlabel + " bin 3 (true)");
        Asym2D_TrueUnfbin3->GetYaxis()->SetTitle(asymlabel + " bin 3 (unfolded)");
        Asym2D_TrueUnfbin3->Draw("AP same");
        Asym2D_TrueUnfbin3->Fit("pol1");

        c_ttbar->SaveAs(acceptanceName + "_LinearityCheck.pdf");
        c_ttbar->SaveAs(acceptanceName + "_LinearityCheck.C");



        TCanvas *c_Pull_lin = new TCanvas("c_Pull_lin", "c_Pull_lin", 500, 500);
        c_Pull_lin->Divide(2, 2);
        c_Pull_lin->cd(1);
        Asym2D_Pull->SetTitle(asymlabel);
        Asym2D_Pull->SetMarkerStyle(23);
        Asym2D_Pull->SetMarkerColor(kBlack);
        Asym2D_Pull->SetMarkerSize(0.6);
        Asym2D_Pull->GetXaxis()->SetTitle(asymlabel + " inclusive (true)");
        Asym2D_Pull->GetYaxis()->SetTitle(asymlabel + " inclusive pull");
        Asym2D_Pull->Draw("AP same");
        Asym2D_Pull->Fit("pol1");

        c_Pull_lin->cd(2);
        Asym2D_Pullbin1->SetTitle(asymlabel);
        Asym2D_Pullbin1->SetMarkerStyle(23);
        Asym2D_Pullbin1->SetMarkerColor(kBlue);
        Asym2D_Pullbin1->SetMarkerSize(0.6);
        Asym2D_Pullbin1->GetXaxis()->SetTitle(asymlabel + " bin 1 (true)");
        Asym2D_Pullbin1->GetYaxis()->SetTitle(asymlabel + " bin 1 pull");
        Asym2D_Pullbin1->Draw("AP same");
        Asym2D_Pullbin1->Fit("pol1");

        c_Pull_lin->cd(3);
        Asym2D_Pullbin2->SetTitle(asymlabel);
        Asym2D_Pullbin2->SetMarkerStyle(23);
        Asym2D_Pullbin2->SetMarkerColor(kBlue);
        Asym2D_Pullbin2->SetMarkerSize(0.6);
        Asym2D_Pullbin2->GetXaxis()->SetTitle(asymlabel + " bin 2 (true)");
        Asym2D_Pullbin2->GetYaxis()->SetTitle(asymlabel + " bin 2 pull");
        Asym2D_Pullbin2->Draw("AP same");
        Asym2D_Pullbin2->Fit("pol1");

        c_Pull_lin->cd(4);
        Asym2D_Pullbin3->SetTitle(asymlabel);
        Asym2D_Pullbin3->SetMarkerStyle(23);
        Asym2D_Pullbin3->SetMarkerColor(kBlue);
        Asym2D_Pullbin3->SetMarkerSize(0.6);
        Asym2D_Pullbin3->GetXaxis()->SetTitle(asymlabel + " bin 3 (true)");
        Asym2D_Pullbin3->GetYaxis()->SetTitle(asymlabel + " bin 3 pull");
        Asym2D_Pullbin3->Draw("AP same");
        Asym2D_Pullbin3->Fit("pol1");

        c_Pull_lin->SaveAs(acceptanceName + "_LinearityCheck_Pull.pdf");
        c_Pull_lin->SaveAs(acceptanceName + "_LinearityCheck_Pull.C");




        TCanvas *c_PullWidth_lin = new TCanvas("c_PullWidth_lin", "c_PullWidth_lin", 500, 500);
        c_PullWidth_lin->Divide(2, 2);
        c_PullWidth_lin->cd(1);
        Asym2D_PullWidth->SetTitle(asymlabel);
        Asym2D_PullWidth->SetMarkerStyle(23);
        Asym2D_PullWidth->SetMarkerColor(kBlack);
        Asym2D_PullWidth->SetMarkerSize(0.6);
        Asym2D_PullWidth->GetXaxis()->SetTitle(asymlabel + " inclusive (true)");
        Asym2D_PullWidth->GetYaxis()->SetTitle(asymlabel + " inclusive pull width");
        Asym2D_PullWidth->Draw("AP same");
        Asym2D_PullWidth->Fit("pol1");

        c_PullWidth_lin->cd(2);
        Asym2D_PullWidthbin1->SetTitle(asymlabel);
        Asym2D_PullWidthbin1->SetMarkerStyle(23);
        Asym2D_PullWidthbin1->SetMarkerColor(kBlue);
        Asym2D_PullWidthbin1->SetMarkerSize(0.6);
        Asym2D_PullWidthbin1->GetXaxis()->SetTitle(asymlabel + " bin 1 (true)");
        Asym2D_PullWidthbin1->GetYaxis()->SetTitle(asymlabel + " bin 1 pull width");
        Asym2D_PullWidthbin1->Draw("AP same");
        Asym2D_PullWidthbin1->Fit("pol1");

        c_PullWidth_lin->cd(3);
        Asym2D_PullWidthbin2->SetTitle(asymlabel);
        Asym2D_PullWidthbin2->SetMarkerStyle(23);
        Asym2D_PullWidthbin2->SetMarkerColor(kBlue);
        Asym2D_PullWidthbin2->SetMarkerSize(0.6);
        Asym2D_PullWidthbin2->GetXaxis()->SetTitle(asymlabel + " bin 2 (true)");
        Asym2D_PullWidthbin2->GetYaxis()->SetTitle(asymlabel + " bin 2 pull width");
        Asym2D_PullWidthbin2->Draw("AP same");
        Asym2D_PullWidthbin2->Fit("pol1");

        c_PullWidth_lin->cd(4);
        Asym2D_PullWidthbin3->SetTitle(asymlabel);
        Asym2D_PullWidthbin3->SetMarkerStyle(23);
        Asym2D_PullWidthbin3->SetMarkerColor(kBlue);
        Asym2D_PullWidthbin3->SetMarkerSize(0.6);
        Asym2D_PullWidthbin3->GetXaxis()->SetTitle(asymlabel + " bin 3 (true)");
        Asym2D_PullWidthbin3->GetYaxis()->SetTitle(asymlabel + " bin 3 pull width");
        Asym2D_PullWidthbin3->Draw("AP same");
        Asym2D_PullWidthbin3->Fit("pol1");

        c_PullWidth_lin->SaveAs(acceptanceName + "_LinearityCheck_PullWidth.pdf");
        c_PullWidth_lin->SaveAs(acceptanceName + "_LinearityCheck_PullWidth.C");


        gStyle->SetOptStat(0);
        TCanvas *c_asymdist_lin = new TCanvas("c_asymdist_lin", "c_asymdist_lin", 1000, 500);
        c_asymdist_lin->Divide(4, 2);
        c_asymdist_lin->cd(1);
        hTrue_before->SetLineColor(TColor::GetColorDark(kRed));
        hTrue_before->SetLineWidth(1);
        hTrue_before->SetFillStyle(0);
        hTrue_before->GetXaxis()->SetTitle(yaxislabel + " #times sign(" + xaxislabel + ")");
        hTrue_before->GetYaxis()->SetTitle("Number of events");
        hTrue_before->Draw("hist");
        hMeas_before->SetLineColor(TColor::GetColorDark(kBlue));
        hMeas_before->SetLineWidth(1);
        hMeas_before->SetFillStyle(0);
        hMeas_before->GetXaxis()->SetTitle(yaxislabel + " #times sign(" + xaxislabel + ")");
        hMeas_before->GetYaxis()->SetTitle("Number of events");
        hMeas_before->Draw("hist same");

        TLegend *leg1 = new TLegend(0.60, 0.75, 0.9, 0.93, NULL, "brNDC");
        leg1->SetEntrySeparation(0.1);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->SetTextSize(0.03);
        leg1->AddEntry(hTrue_before,    "gen",  "L");
        leg1->AddEntry(hMeas_before,    "reco", "L");
        leg1->Draw();


        for (int k = 0; k < Nlin; ++k)
        {

            c_asymdist_lin->cd(k + 2);
            hTrue_after_array[k]->SetLineColor(TColor::GetColorDark(kRed));
            hTrue_after_array[k]->SetLineWidth(1);
            hTrue_after_array[k]->SetFillStyle(0);
            hTrue_after_array[k]->GetXaxis()->SetTitle(yaxislabel + " #times sign(" + xaxislabel + ")");
            hTrue_after_array[k]->GetYaxis()->SetTitle("Number of events");
            hTrue_after_array[k]->Draw("hist");
            hMeas_after_array[k]->SetLineColor(TColor::GetColorDark(kBlue));
            hMeas_after_array[k]->SetLineWidth(1);
            hMeas_after_array[k]->SetFillStyle(0);
            hMeas_after_array[k]->GetXaxis()->SetTitle(yaxislabel + " #times sign(" + xaxislabel + ")");
            hMeas_after_array[k]->GetYaxis()->SetTitle("Number of events");
            hMeas_after_array[k]->Draw("hist same");

            TPaveText *pt1 = new TPaveText(0.40, 0.76, 0.60, 0.91, "brNDC");
            pt1->SetName("pt1name");
            pt1->SetBorderSize(0);
            pt1->SetFillStyle(0);

            TText *blah;

            Float_t Afb, AfbErr;
            GetAfb(hTrue_after_array[k], Afb, AfbErr);

            TString Asym1_temp = formatFloat(Afb, "%6.2f");
            Asym1_temp.ReplaceAll(" " , "" );
            Asym1_temp = TString(" Asym (gen): ") +  Asym1_temp;
            blah = pt1->AddText(Asym1_temp.Data());
            blah->SetTextSize(0.03);
            blah->SetTextAlign(11);
            blah->SetTextColor(TColor::GetColorDark(kRed));

            GetAfb(hMeas_after_array[k], Afb, AfbErr);

            TString Asym2_temp = formatFloat(Afb, "%6.2f");
            Asym2_temp.ReplaceAll(" " , "" );
            Asym2_temp = TString(" Asym (reco): ") +  Asym2_temp;
            blah = pt1->AddText(Asym2_temp.Data());
            blah->SetTextSize(0.03);
            blah->SetTextAlign(11);
            blah->SetTextColor(TColor::GetColorDark(kBlue));

            pt1->Draw();

        }

        c_asymdist_lin->SaveAs(acceptanceName + "_LinearityCheck_AsymDists.pdf");
        c_asymdist_lin->SaveAs(acceptanceName + "_LinearityCheck_AsymDists.C");


    }
    else
    {


        TCanvas *c_pull = new TCanvas("c_pull", "c_pull", 500, 500);
        c_pull->Divide(2, 2);
        c_pull->cd(1);
        AfbPull[0]->SetMarkerStyle(23);
        AfbPull[0]->SetMarkerColor(kBlack);
        AfbPull[0]->SetMarkerSize(0.6);
        AfbPull[0]->GetXaxis()->SetTitle(asymlabel + " inclusive pull");
        AfbPull[0]->GetYaxis()->SetTitle("Number of PEs / 0.2");
        AfbPull[0] ->Fit("gaus");
        AfbPull[0] ->Draw();

        c_pull->cd(2);
        AfbPull[1]->SetMarkerStyle(23);
        AfbPull[1]->SetMarkerColor(kBlue);
        AfbPull[1]->SetMarkerSize(0.6);
        AfbPull[1]->GetXaxis()->SetTitle(asymlabel + " bin 1 pull");
        AfbPull[1]->GetYaxis()->SetTitle("Number of PEs / 0.2");
        AfbPull[1] ->Fit("gaus");
        AfbPull[1] ->Draw();

        c_pull->cd(3);
        AfbPull[2]->SetMarkerStyle(23);
        AfbPull[2]->SetMarkerColor(kBlue);
        AfbPull[2]->SetMarkerSize(0.6);
        AfbPull[2]->GetXaxis()->SetTitle(asymlabel + " bin 2 pull");
        AfbPull[2]->GetYaxis()->SetTitle("Number of PEs / 0.2");
        AfbPull[2] ->Fit("gaus");
        AfbPull[2] ->Draw();

        c_pull->cd(4);
        AfbPull[3]->SetMarkerStyle(23);
        AfbPull[3]->SetMarkerColor(kBlue);
        AfbPull[3]->SetMarkerSize(0.6);
        AfbPull[3]->GetXaxis()->SetTitle(asymlabel + " bin 3 pull");
        AfbPull[3]->GetYaxis()->SetTitle("Number of PEs / 0.2");
        AfbPull[3] ->Fit("gaus");
        AfbPull[3] ->Draw();

        c_pull->SaveAs(acceptanceName + "_Pull.pdf");
        c_pull->SaveAs(acceptanceName + "_Pull.C");



        TFile *plots = new TFile(acceptanceName + "_plots.root", "RECREATE");
        for (int i = 0; i < nbins2D; i++)
        {
            h_pulls[i] ->Write();
            h_resd[i] ->Write();
        }
        AfbPull[0] ->Write();
        AfbPull[1] ->Write();
        AfbPull[2] ->Write();
        AfbPull[3] ->Write();

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
