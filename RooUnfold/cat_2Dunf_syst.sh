#!/bin/bash
for file in $1/residuals_lepChargeAsym_mtt.txt $1/residuals_lepAzimAsym2_mtt.txt $1/residuals_lepPlusCosTheta_mtt.txt $1/residuals_lepMinusCosTheta_mtt.txt $1/residuals_lepCosTheta_mtt.txt $1/residuals_topSpinCorr_mtt.txt $1/residuals_rapiditydiffMarco_mtt.txt $1/residuals_topCosTheta_mtt.txt
do
	cat $file | grep inclusive | awk '{print $1"\t"$3"_"$4"\t"$7"\t"$9"\t"$11}';
	echo "";
	cat $file | grep bin | awk '{print $1"\t"$3"_"$4"\t"$7"\t"$9"\t"$11}';
	echo "";
	echo "";
	echo "";
done
echo ""
echo ""
echo ""
for file in $1/residuals_lepChargeAsym_ttpt.txt $1/residuals_lepAzimAsym2_ttpt.txt $1/residuals_lepPlusCosTheta_ttpt.txt $1/residuals_lepMinusCosTheta_ttpt.txt $1/residuals_lepCosTheta_ttpt.txt $1/residuals_topSpinCorr_ttpt.txt $1/residuals_rapiditydiffMarco_ttpt.txt $1/residuals_topCosTheta_ttpt.txt
do
	cat $file | grep inclusive | awk '{print $1"\t"$3"_"$4"\t"$7"\t"$9"\t"$11}';
	echo "";
	cat $file | grep bin | awk '{print $1"\t"$3"_"$4"\t"$7"\t"$9"\t"$11}';
	echo "";
	echo "";
	echo "";
done
echo ""
echo ""
echo ""
for file in $1/residuals_lepChargeAsym_ttrapidity2.txt $1/residuals_lepAzimAsym2_ttrapidity2.txt $1/residuals_lepPlusCosTheta_ttrapidity2.txt $1/residuals_lepMinusCosTheta_ttrapidity2.txt $1/residuals_lepCosTheta_ttrapidity2.txt $1/residuals_topSpinCorr_ttrapidity2.txt $1/residuals_rapiditydiffMarco_ttrapidity2.txt $1/residuals_topCosTheta_ttrapidity2.txt
do
	cat $file | grep inclusive | awk '{print $1"\t"$3"_"$4"\t"$7"\t"$9"\t"$11}';
	echo "";
	cat $file | grep bin | awk '{print $1"\t"$3"_"$4"\t"$7"\t"$9"\t"$11}';
	echo "";
	echo "";
	echo "";
done
echo ""
echo ""
echo ""
