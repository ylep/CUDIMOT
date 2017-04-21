#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT
#
# Script for initialise Psi angle and beta2kappa in NODDI-Bingham

Usage() {
    echo ""
    echo "Usage: initialise_Bingham <bindir> <DTI_output_dir> <mask_NIfTI_file><output_dir>"
    echo ""
    exit 1
}

[ "$4" = "" ] && Usage

bindir=$1
PathDTI=$2
mask=$3
outputdir=$4

${bindir}/initialise_Psi ${PathDTI}/dtifit_V1.nii.gz ${PathDTI}/dtifit_V2.nii.gz ${mask} ${outputdir}/initialPsi

#beta_to_kappa = 1 - (eigs(1)/eigs(2))^2;
${FSLDIR}/bin/fslmaths ${PathDTI}/dtifit_L3.nii.gz -div ${PathDTI}/dtifit_L2.nii.gz -sqr ${outputdir}/beta2kappa
${FSLDIR}/bin/fslmaths ${mask} -thr 0 -sub  ${outputdir}/beta2kappa  ${outputdir}/beta2kappa
