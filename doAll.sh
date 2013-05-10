#!/bin/bash

export execution_directory='/home/users/gutsche/software/TAS/AFB7TeV/CMSSW_5_3_2_patch4/src/Analysis_2012/TopAFB_lepton_energy_scale_down'

echo "Changing to $execution_directory"
cd $execution_directory

echo "Setting up environment"
source /code/osgcode/cmssoft/cmsset_default.sh CMSSW_1_6_12
export SCRAM_ARCH=slc5_amd64_gcc462
cmsenv

echo "executing doAll.C"
root -b -q doAll.C >& doAll.log

echo "extracting yields"
root -b -q macro_getYields.C >& getYields.log

echo "produce plots"
rootb -b -q macro_drawHists.C >& macro_drawHists.log

echo "executing doAll_nocuts.C"
root -b -q doAll_nocuts.C >& doAll_nocuts.log

echo "create acceptance histograms"
root -b -q macro_acceptanceHists.C >& macro_acceptanceHists.log
mkdir -p acceptance/mcnlo/
mv -f accept_*.* acceptance/mcnlo/

echo "run 1D unfolding"
cd RooUnfold
root -b -q macro_1DUnfolding.C >& macro_1DUnfolding.log

echo "run 2D unfolding"
root -b -q macro_2DUnfolding.C >& macro_2DUnfolding.log




