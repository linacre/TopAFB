{
gROOT->ProcessLine(".L drawHists.C+");
makePSFile("results/hist_usePtGt2020_hypDisamb_usepfMET_usepfJets_useOS_vetoHypMassLt12_requireBTag_sortJetCandidatesbyPt_generalLeptonVeto_createBabyNtuples_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root",5,0,"default",4.98,0,1,1);
}
