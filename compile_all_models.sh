#!/bin/bash
set -e
models=`ls mymodels`
for m in $models; do
    echo
    echo ------------------------------
    echo -------- Compiling $m --------
    echo ------------------------------
    echo
    export modelname=$m
    make clean
    make install
done