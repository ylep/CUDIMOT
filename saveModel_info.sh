#!/bin/bash
infoFile=bin/$modelname"@info"

echo > $infoFile
echo "---------------------------------------------" >> $infoFile
echo "---------------- NPARAMETERS  ---------------" >> $infoFile
echo "---------------------------------------------" >> $infoFile
cat mymodels/$modelname/modelparameters.h | grep define | grep -v INCLUDED >> $infoFile
cat mymodels/$modelname/modelparameters.cc | grep int >> $infoFile
echo "---------------------------------------------" >> $infoFile
echo >> $infoFile
echo >> $infoFile

echo "-------------------------------------------------------------" >> $infoFile
echo "--- Predicted Signal, Constraints & Derivatives Functions ---" >> $infoFile
echo "--- FILE: "$functionsFile" ---" >> $infoFile
echo "-------------------------------------------------------------" >> $infoFile
cat mymodels/$modelname/modelfunctions.h >> $infoFile
