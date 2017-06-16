#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT
#
# Script for initialise k2 from k1. K2=k1*K2_to_k1

Usage() {
    echo ""
    echo "Usage: setK2.sh <k1_file> <k2_to_k1 file> <output_k2_file>"
    echo ""
    exit 1
}

[ "$3" = "" ] && Usage

k1=$1
k2_to_k1=$2
output=$3

${FSLDIR}/bin/fslmaths $k1 -mul $k2_to_k1 $output 
