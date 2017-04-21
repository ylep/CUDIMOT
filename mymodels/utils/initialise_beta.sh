#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT
#
# Script for initialise beta Bingham parameter in NODDI-Bingham using kappa parameter and 
# beta2kappa relation values. 
# beta_to_kappa = 1 - (eigs(1)/eigs(2))^2;

Usage() {
    echo ""
    echo "Usage: initialise_beta <kappa_file> <beta2kappa_file> <output_file>"
    echo ""
    exit 1
}

[ "$3" = "" ] && Usage

kappa=$1
beta2kappa=$2
output=$3

${FSLDIR}/bin/fslmaths ${kappa} -mul ${beta2kappa} ${output}
