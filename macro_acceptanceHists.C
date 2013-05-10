{
gROOT->ProcessLine(".L acceptanceplots.C+");
acceptanceplots("lepAzimAsym2");
acceptanceplots("lepAzimAsym");
acceptanceplots("lepChargeAsym");
acceptanceplots("topCosTheta");
acceptanceplots("lepPlusCosTheta");
acceptanceplots("lepMinusCosTheta");
acceptanceplots("topSpinCorr");
acceptanceplots("rapiditydiff");
acceptanceplots("pseudorapiditydiff");
acceptanceplots("rapiditydiffMarco");
acceptanceplots("lepCosTheta");
}
