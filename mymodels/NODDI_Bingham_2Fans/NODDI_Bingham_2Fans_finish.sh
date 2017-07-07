#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT

bindir=${FSLDEVDIR}/bin

Usage() {
    echo ""
    echo "Usage: NODDI_Bingham_2Fans_finish.sh <output_directory>"
    echo ""
    echo "expects to find all the estimatedParameters and nodif_brain_mask in subject directory"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

directory=$1
cd ${directory}

mv $directory/Param_0_samples.nii.gz $directory/fiso_samples.nii.gz
mv $directory/Param_1_samples.nii.gz $directory/fintra_samples.nii.gz
mv $directory/Param_2_samples.nii.gz $directory/ffib2_samples.nii.gz
mv $directory/Param_3_samples.nii.gz $directory/kappa1_samples.nii.gz
mv $directory/Param_4_samples.nii.gz $directory/beta1_samples.nii.gz
mv $directory/Param_5_samples.nii.gz $directory/th1_samples.nii.gz
mv $directory/Param_6_samples.nii.gz $directory/ph1_samples.nii.gz
mv $directory/Param_7_samples.nii.gz $directory/psi1_samples.nii.gz
mv $directory/Param_8_samples.nii.gz $directory/kappa2_samples.nii.gz
mv $directory/Param_9_samples.nii.gz $directory/beta2_samples.nii.gz
mv $directory/Param_10_samples.nii.gz $directory/th2_samples.nii.gz
mv $directory/Param_11_samples.nii.gz $directory/ph2_samples.nii.gz
mv $directory/Param_12_samples.nii.gz $directory/psi2_samples.nii.gz

${FSLDIR}/bin/fslmaths $directory/kappa1_samples.nii.gz -sub $directory/beta1_samples.nii.gz $directory/concentration_b1_samples.nii.gz 
${FSLDIR}/bin/fslmaths $directory/kappa2_samples.nii.gz -sub $directory/beta2_samples.nii.gz $directory/concentration_b2_samples.nii.gz 

Two_div_pi=0.636619772367581

$FSLDIR/bin/fslmaths $directory/fiso_samples.nii.gz -Tmean $directory/mean_fiso
$FSLDIR/bin/fslmaths $directory/fintra_samples.nii.gz -Tmean $directory/mean_fintra
$FSLDIR/bin/fslmaths $directory/ffib2_samples.nii.gz -Tmean $directory/mean_ffib2
$FSLDIR/bin/fslmaths $directory/kappa1_samples.nii.gz -Tmean $directory/mean_kappa1
$FSLDIR/bin/fslmaths $directory/beta1_samples.nii.gz -Tmean $directory/mean_beta1
$FSLDIR/bin/fslmaths $directory/concentration_b1_samples.nii.gz -Tmean $directory/mean_concentration_b1
$FSLDIR/bin/fslmaths $directory/kappa2_samples.nii.gz -Tmean $directory/mean_kappa2
$FSLDIR/bin/fslmaths $directory/beta2_samples.nii.gz -Tmean $directory/mean_beta2
$FSLDIR/bin/fslmaths $directory/concentration_b2_samples.nii.gz -Tmean $directory/mean_concentration_b2

$FSLDIR/bin/make_dyadic_vectors $directory/th1_samples $directory/ph1_samples $directory/nodif_brain_mask.nii.gz dyads1
$FSLDIR/bin/make_dyadic_vectors $directory/th2_samples $directory/ph2_samples $directory/nodif_brain_mask.nii.gz dyads2

${bindir}/utils/getFanningOrientation $directory/th1_samples $directory/ph1_samples $directory/psi1_samples $directory/nodif_brain_mask.nii.gz  $directory/fanningDir1
${bindir}/utils/getFanningOrientation $directory/th2_samples $directory/ph2_samples $directory/psi2_samples $directory/nodif_brain_mask.nii.gz  $directory/fanningDir2

${FSLDIR}/bin/fslmaths $directory/mean_concentration_b1 -recip -atan -mul $Two_div_pi $directory/ODIp1
${FSLDIR}/bin/fslmaths $directory/mean_kappa1 -recip -atan -mul $Two_div_pi $directory/ODIs1
${FSLDIR}/bin/fslmaths $directory/mean_kappa1 -recip $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/mean_concentration_b1 -recip $directory/ODItot_temp2
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -mul $directory/ODItot_temp2 $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -sqrt -abs -atan -mul $Two_div_pi $directory/ODItot1

${FSLDIR}/bin/fslmaths $directory/mean_concentration_b2 -recip -atan -mul $Two_div_pi $directory/ODIp2
${FSLDIR}/bin/fslmaths $directory/mean_kappa2 -recip -atan -mul $Two_div_pi $directory/ODIs2
${FSLDIR}/bin/fslmaths $directory/mean_kappa2 -recip $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/mean_concentration_b2 -recip $directory/ODItot_temp2
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -mul $directory/ODItot_temp2 $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -sqrt -abs -atan -mul $Two_div_pi $directory/ODItot2

rm $directory/ODItot_temp1.nii.gz
rm $directory/ODItot_temp2.nii.gz

#Dispersion Anisotropy Index
${FSLDIR}/bin/fslmaths $directory/mean_beta1 -div $directory/mean_concentration_b1 -atan -mul $Two_div_pi $directory/DA1
${FSLDIR}/bin/fslmaths $directory/mean_beta2 -div $directory/mean_concentration_b2 -atan -mul $Two_div_pi $directory/DA2

