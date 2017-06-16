#!/bin/bash
set -e
models=`ls mymodels -I utils`
for m in $models; do
    echo
    echo ------------------------------
    echo -------- Compiling $m --------
    echo ------------------------------
    echo
    export modelname=$m
    make cleanall
    make install
done