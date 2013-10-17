#include <iostream>
#include <fstream>
#include "AfbFinalUnfold.h"

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

#include "src/RooUnfold.h"
#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"

#include "tdrstyle.C"

using std::cout;
using std::endl;


//==============================================================================
// Global definitions
//==============================================================================

//const Double_t _topScalingFactor=1.+(9824. - 10070.94)/9344.25;  //needs to be changed from preselection ratio to ratio for events with a ttbar solution

// 0=SVD, 1=TUnfold via RooUnfold, 2=TUnfold
int unfoldingType = 2;

TString Region = "";
Int_t kterm = 3;  //note we used 4 here for ICHEP
Double_t tau = 3E-2;
Int_t nVars = 8;
Int_t includeSys = 0;
Int_t checkErrors = 1;
bool draw_truth_before_pT_reweighting = true;


void AfbUnfoldExample(TString Var2D = "mtt", double scalettdil = 1., double scalettotr = 1., double scalewjets = 1., double scaleDY = 1., double scaletw = 1., double scaleVV = 1.)
{
    setTDRStyle();
    gStyle->SetOptFit();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cout.precision(3);

    TString summary_name;
    if (Var2D == "mtt") summary_name = "summary_2Dunfolding";
    else if (Var2D == "ttrapidity2") summary_name = "summary_2Dunfolding_ttrapidity2";
    else if (Var2D == "ttpt") summary_name = "summary_2Dunfolding_ttpt";

    TString yaxisunit;
    if (Var2D == "mtt") yaxisunit = " (GeV/c^{2})";
    else if (Var2D == "ttrapidity2") yaxisunit = "";
    else if (Var2D == "ttpt") yaxisunit = " (GeV/c)";

    if (!(scalettotr == 1. && scalewjets == 1. && scaleDY == 1. && scaletw == 1. && scaleVV == 1.))  summary_name = summary_name + Form("_%i_%i_%i_%i_%i", int(10.*scalettotr + 0.5), int(10.*scalewjets + 0.5), int(10.*scaleDY + 0.5), int(10.*scaletw + 0.5), int(10.*scaleVV + 0.5));

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);

    ofstream myfile;
    myfile.open (summary_name + ".txt");
    //cout.rdbuf(myfile.rdbuf());

    // OGU 130516: add second output txt file with format easier to be pasted into google docs
    ofstream second_output_file;
    second_output_file.open(summary_name + "_formated.txt");

    const int nBkg = 7;
    TString path = "../";
    TString bkgroot[nBkg] = {"ttotr.root", "wjets.root", "DYee.root", "DYmm.root", "DYtautau.root", "tw.root", "VV.root"};
    double bkgSF[nBkg] = {scalettotr, scalewjets, scaleDY, scaleDY, scaleDY, scaletw, scaleVV};


    Float_t asymVar, asymVar_gen, tmass; //"asymVar" formerly known as "observable"
    Float_t asymVarMinus, asymVarMinus_gen;
    Float_t obs2D, obs2D_gen; //Note: change the name of obs2D to be something like "secondary variable"
    Double_t weight;
    Int_t Nsolns;

    for (Int_t iVar = 0; iVar < nVars; iVar++)
    {
	  ///////////////////////////////////////////////////////////////////////////////////////////
	  /////////////// 1. Set up all our histograms //////////////////////////////////////////////
	  
        //initialise 1D binning to get x values
        // Initialize1DBinning(iVar);
		// Float_t asym_centre = (xmax1D + xmin1D) / 2.;

        //initialise 2D binning
        if (Var2D == "mtt") Initialize2DBinning(iVar);
        else if (Var2D == "ttrapidity2") Initialize2DBinningttrapidity2(iVar);
        else if (Var2D == "ttpt") Initialize2DBinningttpt(iVar);
        bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;

		//Use proper 2D histograms, instead of the old 1D ones
        TH2D *hData = new TH2D ("Data_BkgSub", "Data with background subtracted",    nbinsx2D, xbins2D, nbinsy2D, ybins2D);
        TH2D *hBkg = new TH2D ("Background",  "Background",    nbinsx2D, xbins2D, nbinsy2D, ybins2D);
        TH2D *hData_unfolded = new TH2D ("Data_Unfold", "Data with background subtracted and unfolded", nbinsx2D, xbins2D, nbinsy2D, ybins2D);
        TH2D *hTrue = new TH2D ("true", "Truth",    nbinsx2D, xbins2D, nbinsy2D, ybins2D);
        TH2D *hMeas = new TH2D ("meas", "Measured", nbinsx2D, xbins2D, nbinsy2D, ybins2D);
		//for testing purposes
        TH2D *hData_bkgSub_rewrapped = new TH2D ("bkgsub", "bkgsub", nbinsx2D, xbins2D, nbinsy2D, ybins2D);
        //TH2D *hData_unwrapped_rewrapped = new TH2D ("lotsofwrap", "lotsofwrap", nbinsx2D, xbins2D, nbinsy2D, ybins2D);

		//Unwrapped histograms have n bins (where n = nx*ny), centered around the integers from 1 to n.
		TH1D *hData_unwrapped = new TH1D ("Data_BkgSub_Unwr", "Unwrapped data with background subtracted", nbinsunwrapped, 0.5, double(nbinsunwrapped)+0.5);
        TH1D *hBkg_unwrapped = new TH1D ("Background_Unwr",  "Background unwrapped",    nbinsunwrapped, 0.5, double(nbinsunwrapped)+0.5);
        TH1D *hData_unfolded_unwrapped = new TH1D ("Data_Unfold_Unwr", "Data with background subtracted and unfolded, unwrapped", nbinsunwrapped, 0.5, double(nbinsunwrapped)+0.5);
        TH1D *hTrue_unwrapped = new TH1D ("true_unwr", "Truth unwrapped",  nbinsunwrapped, 0.5, double(nbinsunwrapped)+0.5);
        TH1D *hMeas_unwrapped = new TH1D ("meas_unwr", "Measured unwrapped", nbinsunwrapped, 0.5, double(nbinsunwrapped)+0.5);

		//An n*n migration matrix, using the unwrapped binning on both axes
        TH2D *hTrue_vs_Meas = new TH2D ("true_vs_meas", "True vs Measured", nbinsunwrapped, 0.5, double(nbinsunwrapped)+0.5, nbinsunwrapped, 0.5, double(nbinsunwrapped)+0.5);

        TH1D *hData_bkgSub;

        hData->Sumw2();
        hBkg->Sumw2();
        hData_unfolded->Sumw2();
        hTrue->Sumw2();
        hMeas->Sumw2();
        hTrue_vs_Meas->Sumw2();

        hData_unwrapped->Sumw2();
        hBkg_unwrapped->Sumw2();
        hData_unfolded_unwrapped->Sumw2();
        hTrue_unwrapped->Sumw2();
        hMeas_unwrapped->Sumw2();


        TMatrixD m_unfoldE(nbinsunwrapped, nbinsunwrapped);
        TMatrixD m_correctE(nbinsunwrapped, nbinsunwrapped);

		/////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////// 2. Fill our histograms from the baby ntuples //////////////////

		//cout << "Filling histo for " << acceptanceName << "." << endl;
		//cout << "X-axis: " << xaxislabel << endl;
		//cout << "Y-axis: " << yaxislabel << endl;

		// Set up chains
        TChain *ch_bkg[nBkg];
        TChain *ch_top = new TChain("tree");
        TChain *ch_data = new TChain("tree");

        TString path = "../";

        ch_data->Add(path + "data.root");
        ch_top->Add(path + "ttdil.root");

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg] = new TChain("tree");
            ch_bkg[iBkg]->Add(path + bkgroot[iBkg]);
        }

		///// Load data from data chain, and fill hData //////////
        ch_data->SetBranchAddress(observablename,    &asymVar);
        if ( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &asymVarMinus);
        ch_data->SetBranchAddress("weight", &weight);
        ch_data->SetBranchAddress("Nsolns", &Nsolns);
        ch_data->SetBranchAddress("t_mass", &tmass);

        if (Var2D == "mtt")               ch_data->SetBranchAddress("tt_mass", &obs2D);
        else if (Var2D == "ttrapidity2")  ch_data->SetBranchAddress("ttRapidity2", &obs2D);
        else if (Var2D == "ttpt")         ch_data->SetBranchAddress("ttPt", &obs2D);


        for (Int_t i = 0; i < ch_data->GetEntries(); i++)
        {
            ch_data->GetEntry(i);
            obs2D = fabs(obs2D);
            if ( tmass > 0 )
            {
			  fillUnderOverFlow(hData, asymVar, obs2D, weight, Nsolns);
			  if (combineLepMinus) fillUnderOverFlow(hData, asymVarMinus, obs2D, weight, Nsolns);
            }
        }

		///// Load background MC from background chain, and fill h_bkg //////////
        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg]->SetBranchAddress(observablename,    &asymVar);
            if ( combineLepMinus ) ch_bkg[iBkg]->SetBranchAddress("lepMinus_costheta_cms",    &asymVarMinus);
            ch_bkg[iBkg]->SetBranchAddress("weight", &weight);
            ch_bkg[iBkg]->SetBranchAddress("Nsolns", &Nsolns);
            ch_bkg[iBkg]->SetBranchAddress("t_mass", &tmass);

            if (Var2D == "mtt")              ch_bkg[iBkg]->SetBranchAddress("tt_mass", &obs2D);
            else if (Var2D == "ttrapidity2") ch_bkg[iBkg]->SetBranchAddress("ttRapidity2", &obs2D);
            else if (Var2D == "ttpt")        ch_bkg[iBkg]->SetBranchAddress("ttPt", &obs2D);

            for (Int_t i = 0; i < ch_bkg[iBkg]->GetEntries(); i++)
            {
                ch_bkg[iBkg]->GetEntry(i);
                obs2D = fabs(obs2D);
                weight *= bkgSF[iBkg];
                if ( tmass > 0 )
                {
				  fillUnderOverFlow(hBkg, asymVar, obs2D, weight, Nsolns);
				  if (combineLepMinus)  fillUnderOverFlow(hBkg, asymVarMinus, obs2D, weight, Nsolns);
                }
            }
        }

		///// Load true top MC from top chain, and fill h_true and hTrue_vs_Meas ///////////
        ch_top->SetBranchAddress(observablename,    &asymVar);
        ch_top->SetBranchAddress(observablename + "_gen", &asymVar_gen);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms",    &asymVarMinus);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms_gen",    &asymVarMinus_gen);
        ch_top->SetBranchAddress("weight", &weight);
        ch_top->SetBranchAddress("Nsolns", &Nsolns);
        ch_top->SetBranchAddress("t_mass", &tmass);

        if (Var2D == "mtt")
        {
            ch_top->SetBranchAddress("tt_mass", &obs2D);
            ch_top->SetBranchAddress("tt_mass_gen", &obs2D_gen);
        }
        else if (Var2D == "ttrapidity2")
        {
            ch_top->SetBranchAddress("ttRapidity2", &obs2D);
            ch_top->SetBranchAddress("ttRapidity2_gen", &obs2D_gen);
        }
        else if (Var2D == "ttpt")
        {
            ch_top->SetBranchAddress("ttPt", &obs2D);
            ch_top->SetBranchAddress("ttPt_gen", &obs2D_gen);
        }

		int measbin = 0;
		int genbin = 0;

        for (Int_t i = 0; i < ch_top->GetEntries(); i++)
        {
            ch_top->GetEntry(i);
            obs2D = fabs(obs2D);
            obs2D_gen = fabs(obs2D_gen);
            weight *= scalettdil;
            if ( tmass > 0 )
            {
			  measbin = getUnwrappedBin(hMeas, asymVar, obs2D);
			  genbin  = getUnwrappedBin(hTrue, asymVar_gen, obs2D_gen);

			  fillUnderOverFlow(hMeas, asymVar, obs2D, weight, Nsolns);
			  fillUnderOverFlow(hTrue, asymVar_gen, obs2D_gen, weight, Nsolns);
			  fillUnderOverFlow(hTrue_vs_Meas, measbin, genbin, weight, Nsolns);
			  if ( combineLepMinus )
                {
				  measbin = getUnwrappedBin(hMeas, asymVarMinus, obs2D);
				  genbin  = getUnwrappedBin(hTrue, asymVarMinus_gen, obs2D_gen);

				  fillUnderOverFlow(hMeas, asymVarMinus, obs2D, weight, Nsolns);
				  fillUnderOverFlow(hTrue, asymVarMinus_gen, obs2D_gen, weight, Nsolns);
				  fillUnderOverFlow(hTrue_vs_Meas, measbin, genbin, weight, Nsolns);
                }
            }
        }

        RooUnfoldResponse response (hMeas, hTrue, hTrue_vs_Meas);

		///////////// Set up and draw some stuff ////////////////////////////////
		// Note: probably this section is broken by the switch to true 2D unfolding

		/*
        TCanvas *c_mtt = new TCanvas("c_mtt", "c_mtt", 500, 500);

        hData->SetLineWidth(lineWidth + 2);

        hTrue->SetLineColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillStyle(3353);

        hMeas->SetLineColor(TColor::GetColorDark(kGreen));
        hMeas->SetFillColor(TColor::GetColorDark(kGreen));
        hMeas->SetFillStyle(3353);

        hBkg->SetLineColor(kYellow);
        hBkg->SetFillColor(kYellow);


        THStack *hs = new THStack("hs", "Stacked Top+BG");  //does it even make sense to have a THStack anymore?

        hs->Add(hBkg);
        hs->Add(hMeas);

        hs->SetMinimum(0.0);
        hs->SetMaximum( 2.0 * hs->GetMaximum());
        hs->Draw("hist");
        hs->GetXaxis()->SetTitle(yaxislabel + yaxisunit);
        hs->GetYaxis()->SetTitleOffset(1.3);
        hs->GetYaxis()->SetTitle("Events");

        hData->Draw("E same");

        TLegend *leg1 = new TLegend(0.6, 0.62, 0.9, 0.838, NULL, "brNDC");
        leg1->SetEntrySeparation(100);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.03);
        leg1->SetFillStyle(0);
        leg1->AddEntry(hData, "Data");
        leg1->AddEntry(hMeas,  "MC@NLO reco level", "F");
        leg1->AddEntry(hBkg,  "Background", "F");
        leg1->Draw();
        c_mtt->SaveAs(Var2D + "_" + acceptanceName + Region + ".pdf");
		*/


		////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////// 3. Perform the unfolding procedure ////////////////////////////////////

		// First step: unwrap our TH2Ds into TH1Ds!
		unwrap2dhisto(hData, hData_unwrapped);
		unwrap2dhisto(hBkg,  hBkg_unwrapped);
		unwrap2dhisto(hTrue, hTrue_unwrapped);
		unwrap2dhisto(hMeas, hMeas_unwrapped);

        hData_bkgSub = (TH1D *) hData_unwrapped->Clone();
        hData_bkgSub->Add(hBkg_unwrapped, -1.0);

        // TCanvas *c_data = new TCanvas("c_data", "c_data");
		// hData_bkgSub->Draw();
        // c_data->SaveAs(Var2D + "_data_unwrapped_" + acceptanceName + Region + ".pdf");


		// Now let's actually do the unfolding.


        if (unfoldingType == 0)
	    {
		  //Dan made no changes here. Probably broken.
            RooUnfoldSvd unfold_svd (&response, hData_bkgSub, kterm);
            unfold_svd.Setup(&response, hData_bkgSub);
            unfold_svd.IncludeSystematics(includeSys);
            hData_unfolded = (TH1D *) unfold_svd.Hreco();
            m_unfoldE = unfold_svd.Ereco();
        }
        else if (unfoldingType == 1)
        {
		  //Dan made no changes here. Probably broken.
            RooUnfoldTUnfold unfold_rooTUnfold (&response, hData_bkgSub, TUnfold::kRegModeCurvature);
            unfold_rooTUnfold.Setup(&response, hData_bkgSub);
            unfold_rooTUnfold.FixTau(tau);
            unfold_rooTUnfold.IncludeSystematics(includeSys);
            hData_unfolded = (TH1D *) unfold_rooTUnfold.Hreco();
            m_unfoldE = unfold_rooTUnfold.Ereco();
        }
        else if (unfoldingType == 2)
        {
		  TUnfold unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone);  //need to set mode "None" here if regularizing by hand
            unfold_TUnfold.SetInput(hData_bkgSub);
            //Double_t biasScale=1.0;
            unfold_TUnfold.SetBias(hTrue_unwrapped);
			unfold_TUnfold.RegularizeBins2D(1,1,nbinsx2D,nbinsx2D,nbinsy2D,TUnfold::kRegModeCurvature);
            unfold_TUnfold.DoUnfold(tau,hData_bkgSub, 0.0);
            unfold_TUnfold.GetOutput(hData_unfolded_unwrapped);


            TH2D *ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
            for (Int_t cmi = 0; cmi < nbinsunwrapped; cmi++)
            {
                for (Int_t cmj = 0; cmj < nbinsunwrapped; cmj++)
                {
                    m_unfoldE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
                }
            }
        }
        else cout << "Unfolding TYPE not Specified" << "\n";


	   
		//Re-wrap 1D histogram into 2D
		rewrap1dhisto(hData_unfolded_unwrapped, hData_unfolded);
		rewrap1dhisto(hData_bkgSub, hData_bkgSub_rewrapped);
		// rewrap1dhisto(hData_unwrapped, hData_unwrapped_rewrapped);


		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////// 4. Output a bunch of histograms and tables //////////////////////////////////

        TCanvas *c_data = new TCanvas("c_data", "c_data");
        gStyle->SetPalette(1);
        hData_bkgSub_rewrapped->SetTitle("Data (background subtracted);"+xaxislabel+";"+yaxislabel);
        hData_bkgSub_rewrapped->Draw("COLZ");
        c_data->SetLogz();
        c_data->SaveAs("data_" + acceptanceName +"_" + Var2D + ".eps");
		hData_bkgSub->SetTitle("Unwrapped Data (background subtracted);Bin number;Entries per bin");
		hData_bkgSub->Draw();
        c_data->SaveAs("dataunwrapped_" + acceptanceName +"_" + Var2D + ".eps");

        hTrue->SetTitle("True MC;"+xaxislabel+";"+yaxislabel);
		hTrue->Draw("COLZ");
        c_data->SaveAs("true_" + acceptanceName +"_" + Var2D + ".eps");
		hTrue_unwrapped->SetTitle("Unwrapped Truth;Bin number;Entries per bin");
		hTrue_unwrapped->Draw();
        c_data->SaveAs("trueunwrapped_" + acceptanceName +"_" + Var2D + ".eps");

		hData_unfolded->SetTitle("Unfolded Data;"+xaxislabel+";"+yaxislabel);
		hData_unfolded->Draw("COLZ");
        c_data->SaveAs("unfolded_" + acceptanceName +"_" + Var2D + ".eps");


		TCanvas *c1 = new TCanvas("c1","c1",1200,400);
		c1->Divide(3,1);
		c1->cd(1);
		hData_bkgSub_rewrapped->Draw("COLZ");
		c1->cd(2);
		hData_unfolded->Draw("COLZ");
		c1->cd(3);
		hTrue->Draw("COLZ");
		c1->SaveAs("data_comparison_"+acceptanceName+"_"+Var2D+".pdf");

        if (unfoldingType == 0)
        {
            TCanvas *c_d = new TCanvas("c_d", "c_d", 500, 500);
            TH1D *dvec = unfold_svd.Impl()->GetD();
            dvec->Draw();
            c_d->SetLogy();
            c_d->SaveAs(Var2D + "_D_2D_" + acceptanceName + Region + ".pdf");
        }

        TCanvas *c_resp = new TCanvas("c_resp", "c_resp");
        TH2D *hResp = (TH2D *) response.Hresponse();
        gStyle->SetPalette(1);
		hResp->SetTitle("Migration matrix");
        hResp->GetXaxis()->SetTitle(yaxislabel + " and " + xaxislabel + ", unwrapped (reco)");
        hResp->GetYaxis()->SetTitle(yaxislabel + " and " + xaxislabel + ", unwrapped (gen)");
        hResp->Draw("COLZ");
        c_resp->SetLogz();
        c_resp->SaveAs(Var2D + "_Response_true2D_" + acceptanceName + Region + ".pdf");


        TFile *file = new TFile("../acceptance/mcnlo/accept_" + acceptanceName + ".root");

        TH2D *acceptM_2d = (TH2D *) file->Get("accept_" + acceptanceName + "_" + Var2D);

        TH2D *denomM_2d = (TH2D *) file->Get("denominator_" + acceptanceName + "_" + Var2D);



		/*
        TH1D *acceptM = new TH1D ("accept", "accept",    nbins2D, xbins2D);
        acceptM->SetBinContent(1, acceptM_2d->GetBinContent(1, 3));
        acceptM->SetBinContent(2, acceptM_2d->GetBinContent(1, 2));
        acceptM->SetBinContent(3, acceptM_2d->GetBinContent(1, 1));

        acceptM->SetBinContent(4, acceptM_2d->GetBinContent(2, 1));
        acceptM->SetBinContent(5, acceptM_2d->GetBinContent(2, 2));
        acceptM->SetBinContent(6, acceptM_2d->GetBinContent(2, 3));

        acceptM->Scale(1.0 / acceptM->Integral());
		


        TH1D *denomM = new TH1D ("denom", "denom",    nbins2D, xbins2D);

        denomM->SetBinContent(1, denomM_2d->GetBinContent(1, 3));
        denomM->SetBinContent(2, denomM_2d->GetBinContent(1, 2));
        denomM->SetBinContent(3, denomM_2d->GetBinContent(1, 1));

        denomM->SetBinContent(4, denomM_2d->GetBinContent(2, 1));
        denomM->SetBinContent(5, denomM_2d->GetBinContent(2, 2));
        denomM->SetBinContent(6, denomM_2d->GetBinContent(2, 3));

        TH1D *denomM_0 = new TH1D ("denominator0", "denominator0",    2, -1500., 1500.);
        TH1D *denomM_1 = new TH1D ("denominator1", "denominator1",    2, -1500., 1500.);
        TH1D *denomM_2 = new TH1D ("denominator2", "denominator2",    2, -1500., 1500.);

        denomM_2->SetBinContent(1, denomM_2d->GetBinContent(1, 3));
        denomM_1->SetBinContent(1, denomM_2d->GetBinContent(1, 2));
        denomM_0->SetBinContent(1, denomM_2d->GetBinContent(1, 1));

        denomM_0->SetBinContent(2, denomM_2d->GetBinContent(2, 1));
        denomM_1->SetBinContent(2, denomM_2d->GetBinContent(2, 2));
        denomM_2->SetBinContent(2, denomM_2d->GetBinContent(2, 3));
		



        TFile *file_nopTreweighting = new TFile("../acceptance/mcnlo_nopTreweighting/accept_" + acceptanceName + ".root");

        TH2D *denomM_nopTreweighting_2d = (TH2D *) file_nopTreweighting->Get("denominator_" + acceptanceName + "_" + Var2D);
        TH1D *denomM_nopTreweighting = new TH1D ("denomnopTreweighting", "denomnopTreweighting",    nbins2D, xbins2D);

        denomM_nopTreweighting->SetBinContent(1, denomM_nopTreweighting_2d->GetBinContent(1, 3));
        denomM_nopTreweighting->SetBinContent(2, denomM_nopTreweighting_2d->GetBinContent(1, 2));
        denomM_nopTreweighting->SetBinContent(3, denomM_nopTreweighting_2d->GetBinContent(1, 1));

        denomM_nopTreweighting->SetBinContent(4, denomM_nopTreweighting_2d->GetBinContent(2, 1));
        denomM_nopTreweighting->SetBinContent(5, denomM_nopTreweighting_2d->GetBinContent(2, 2));
        denomM_nopTreweighting->SetBinContent(6, denomM_nopTreweighting_2d->GetBinContent(2, 3));

        TH1D *denomM_nopTreweighting_0 = new TH1D ("denominator0nopTreweighting", "denominator0nopTreweighting",    2, -1500., 1500.);
        TH1D *denomM_nopTreweighting_1 = new TH1D ("denominator1nopTreweighting", "denominator1nopTreweighting",    2, -1500., 1500.);
        TH1D *denomM_nopTreweighting_2 = new TH1D ("denominator2nopTreweighting", "denominator2nopTreweighting",    2, -1500., 1500.);

        denomM_nopTreweighting_2->SetBinContent(1, denomM_nopTreweighting_2d->GetBinContent(1, 3));
        denomM_nopTreweighting_1->SetBinContent(1, denomM_nopTreweighting_2d->GetBinContent(1, 2));
        denomM_nopTreweighting_0->SetBinContent(1, denomM_nopTreweighting_2d->GetBinContent(1, 1));

        denomM_nopTreweighting_0->SetBinContent(2, denomM_nopTreweighting_2d->GetBinContent(2, 1));
        denomM_nopTreweighting_1->SetBinContent(2, denomM_nopTreweighting_2d->GetBinContent(2, 2));
        denomM_nopTreweighting_2->SetBinContent(2, denomM_nopTreweighting_2d->GetBinContent(2, 3));
		*/




        for (Int_t x = 1; x <= nbinsx2D; x++)
        {
		  for (Int_t y = 1; y<= nbinsy2D; y++)
		  {

            if (acceptM_2d->GetBinContent(x,y) != 0)
            {
			  hData_unfolded->SetBinContent(x,y, hData_unfolded->GetBinContent(x,y) * 1.0 / acceptM_2d->GetBinContent(x,y));
			  hData_unfolded->SetBinError  (x,y, hData_unfolded->GetBinError  (x,y) * 1.0 / acceptM_2d->GetBinContent(x,y));

			  hTrue->SetBinContent(x,y, hTrue->GetBinContent(x,y) * 1.0 / acceptM_2d->GetBinContent(x,y));
			  hTrue->SetBinError  (x,y, hTrue->GetBinError(x,y)  * 1.0 / acceptM_2d->GetBinContent(x,y));
            }
		  }
		}

		/*
		// Fix me at some point, please!
        for (int l = 0; l < nbins2D; l++)
        {
            for (int j = 0; j < nbins2D; j++)
            {
                double corr = 1.0 / ( acceptM->GetBinContent(l + 1) * acceptM->GetBinContent(j + 1) );
                //corr = corr * pow(xsection / dataIntegral,2) ;
                m_correctE(l, j) = m_unfoldE(l, j) * corr;
            }
        }
		*/

        //==================================================================
        // ============== Print the asymmetry =============================
		cout << "========= Variable: " << acceptanceName << "===================\n";
        Float_t Afb, AfbErr;

        GetAfb(hData, Afb, AfbErr);
        cout << " Data: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hTrue, Afb, AfbErr);
        cout << " True Top: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hData_unfolded, Afb, AfbErr); // Formerly GetCorrectedAfb
        cout << " Unfolded: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " Unfolded: " << Afb << " +/-  " << AfbErr << endl;

        GetAfb(denomM_2d, Afb, AfbErr);
        cout << " True Top from acceptance denominator: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " True_Top_from_acceptance_denominator: " << Afb << " +/-  " << AfbErr << "\n";

        vector<double> afb_m;
        vector<double> afb_merr;
        vector<double> afb_m_denom;
        vector<double> afb_merr_denom;
		cout << "From unfolded data:" << endl;
        GetAvsY2d(hData_unfolded, afb_m, afb_merr, second_output_file);
		cout << "From acceptance denominator:" << endl;
		GetAvsY2d(denomM_2d, afb_m_denom, afb_merr_denom, second_output_file);

		/*
        if (draw_truth_before_pT_reweighting)
        {
            GetAfb(denomM_nopTreweighting_0, AfbG[0], AfbErr);
            GetAfb(denomM_nopTreweighting_1, AfbG[1], AfbErr);
            GetAfb(denomM_nopTreweighting_2, AfbG[2], AfbErr);
        }
		*/

        TCanvas *c_afb = new TCanvas("c_afb", "c_afb", 500, 500);
        double ybinsForHisto[4] = {ybins2D[0], ybins2D[1], ybins2D[2], ybins2D[3]};
        if (Var2D == "mtt") ybinsForHisto[0] = 300.0;
        TH1D *hAfbVsMtt = new TH1D ("AfbVsMtt",  "AfbVsMtt",  3, ybinsForHisto);
        TH1D *hAfbVsMtt_statonly = new TH1D ("AfbVsMtt_statonly",  "AfbVsMtt_statonly",  3, ybinsForHisto);
        for (int nb = 0; nb < 3; nb++)
        {
            hAfbVsMtt->SetBinContent(nb + 1, afb_m[nb]);
            if (checkErrors)
            {
                if (includeSys)
                {
                    cout << "Difference between calculated and hard-coded stat errors: " << afb_merr[nb] - stat_corr[nb] << endl;
                }
                else
                {
                    cout << "Difference between calculated and hard-coded stat errors: " << afb_merr[nb] - stat_uncorr[nb] << endl;
                }
            }
            hAfbVsMtt->SetBinError(nb + 1,  sqrt( pow(afb_merr[nb], 2) + pow(syst_corr[nb], 2) ) );
            hAfbVsMtt_statonly->SetBinContent(nb + 1, afb_m[nb]);
            hAfbVsMtt_statonly->SetBinError(nb + 1, afb_merr[nb]);
        }

        //  GetAvsY(hTrue, m_unfoldE, afb_m, afb_merr);

        TH1D *hTop_AfbVsMtt = new TH1D ("Top_AfbVsMtt",  "Top_AfbVsMtt",  3, ybinsForHisto);
        for (int nb = 0; nb < 3; nb++)
        {
            hTop_AfbVsMtt->SetBinContent(nb + 1, afb_m_denom[nb]);
            hTop_AfbVsMtt->SetBinError(nb + 1, 0);
        }

        tdrStyle->SetErrorX(0.5);
        hAfbVsMtt->SetMinimum( hAfbVsMtt->GetMinimum() - 0.1 );
        hAfbVsMtt->SetMaximum( hAfbVsMtt->GetMaximum() + 0.1 );
        hAfbVsMtt->SetLineWidth( 2.0 );
        hAfbVsMtt->Draw("E");
        hAfbVsMtt_statonly->Draw("E1 same");
        hTop_AfbVsMtt->SetLineColor(TColor::GetColorDark(kRed));
        hTop_AfbVsMtt->SetMarkerColor(TColor::GetColorDark(kRed));
        hTop_AfbVsMtt->SetMarkerSize(0);
        hTop_AfbVsMtt->SetLineWidth( 2.0 );
        hAfbVsMtt->GetYaxis()->SetTitle(asymlabel);
        hAfbVsMtt->GetYaxis()->SetTitleOffset(1.2);
        hAfbVsMtt->GetXaxis()->SetTitle(yaxislabel + yaxisunit);
        hTop_AfbVsMtt->Draw("E same");

        leg1 = new TLegend(0.6, 0.72, 0.9, 0.938, NULL, "brNDC");
        leg1->SetEntrySeparation(100);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.03);
        leg1->SetFillStyle(0);
        leg1->AddEntry(hAfbVsMtt, "data");
        leg1->AddEntry(hTop_AfbVsMtt,    "MC@NLO parton level");
        leg1->Draw();

        TPaveText *pt1 = new TPaveText(0.18, 0.88, 0.41, 0.91, "brNDC");
        pt1->SetName("pt1name");
        pt1->SetBorderSize(0);
        pt1->SetFillStyle(0);

        TText *blah;
        //blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
        blah = pt1->AddText("CMS, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
        blah->SetTextSize(0.032);
        blah->SetTextAlign(11);
        pt1->Draw();

        c_afb->SaveAs("AfbVs" + Var2D + "_unfolded_" + acceptanceName + Region + ".pdf");
        c_afb->SaveAs("AfbVs" + Var2D + "_unfolded_" + acceptanceName + Region + ".root");
        c_afb->SaveAs("AfbVs" + Var2D + "_unfolded_" + acceptanceName + Region + ".C");

		/*
        TCanvas *c_mttu = new TCanvas("c_mttu", "c_mttu", 500, 500);

        hData_unfolded->Scale(1. / hData_unfolded->Integral(), "width");
        hTrue->Scale(1. / hTrue->Integral(), "width");

        hData_unfolded->GetXaxis()->SetTitle(yaxislabel + " #times sign(" + xaxislabel + ")");
        hData_unfolded->GetYaxis()->SetTitle("d#sigma/d(" + yaxislabel + " #times sign(" + xaxislabel + "))");
        hData_unfolded->SetMinimum(0.0);
        hData_unfolded->SetMaximum( 2.0 * hData_unfolded->GetMaximum());
        hData_unfolded->SetMarkerStyle(23);
        hData_unfolded->SetMarkerSize(2.0);
        hData_unfolded->Draw("E");
        hData_unfolded->SetLineWidth(lineWidth);
        hTrue->SetLineWidth(lineWidth);
        hTrue->SetLineColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillColor(TColor::GetColorDark(kGreen));
        hTrue->Draw("hist same");
        hData_unfolded->Draw("E same");

        leg1 = new TLegend(0.6, 0.62, 0.9, 0.838, NULL, "brNDC");
        leg1->SetEntrySeparation(100);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.03);
        leg1->SetFillStyle(0);
        leg1->AddEntry(hData_unfolded, "( Data-BG ) Unfolded");
        leg1->AddEntry(hTrue,    "MC@NLO parton level", "F");
        leg1->Draw();

        c_mttu->SaveAs(Var2D + "_2D_unfolded_" + acceptanceName + Region + ".pdf");

        TFile *output = new TFile(Form("DataMC_%s.root", Var2D.Data()), "UPDATE");

        TH1D *hDataMCratio  = (TH1D *) hData_unfolded->Clone("hDataMCratio" + Var2D + acceptanceName);
        hDataMCratio->SetTitle("hDataMCratio" + Var2D + acceptanceName);
        hDataMCratio->Reset();
        hDataMCratio->Divide(hData_unfolded, hTrue, 1., 1.);
        hDataMCratio->Write();
		
        output->Close();
		*/
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
