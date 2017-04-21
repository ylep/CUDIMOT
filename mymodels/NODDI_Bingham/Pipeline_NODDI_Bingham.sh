#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT
#
# Pipeline for fitting NODDI-Bingham

premodel1=NODDI_Bingham_Init
premodel2=NODDI_Bingham_FitFrac
modelname=NODDI_Bingham
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FSLDIR}/lib
bindir=/home/fs0/moisesf/scratch/THESIS/microstructure_models/CUDIMOT/bin
modelsdir=/home/fs0/moisesf/scratch/THESIS/microstructure_models/CUDIMOT/cudimot/mymodels

Usage() {
    echo ""
    echo "Usage: Pipeline_NODDI_Bingham.sh <subject_directory> [options]"
    echo ""
    echo "expects to find data and nodif_brain_mask in subject directory"
    echo ""
    echo "<options>:"
    echo "-Q (name of the GPU(s) queue, default cuda.q (defined in environment variable: FSLGECUDAQ)"
    echo "-NJOBS (number of jobs to queue, the data is divided in NJOBS parts, usefull for a GPU cluster, default 4)"
    echo "--runMCMC (if you want to run MCMC)"
    echo "-b (burnin period, default 5000)"
    echo "-j (number of jumps, default 1250)"
    echo "-s (sample every, default 25)"
    echo ""
    exit 1
}

make_absolute(){
    dir=$1;
    if [ -d ${dir} ]; then
	OLDWD=`pwd`
	cd ${dir}
	dir_all=`pwd`
	cd $OLDWD
    else
	dir_all=${dir}
    fi
    echo ${dir_all}
}

[ "$1" = "" ] && Usage

subjdir=`make_absolute $1`
subjdir=`echo $subjdir | sed 's/\/$/$/g'`

echo "---------------------------------------------------------------------------------"
echo "------------------------------------ CUDIMOT ------------------------------------"
echo "----------------------------- MODEL: $modelname -----------------------------"
echo "---------------------------------------------------------------------------------"
echo subjectdir is $subjdir

start=`date +%s`

#parse option arguments
qsys=0
njobs=4
fudge=1
burnin=1000
njumps=1250
sampleevery=25
other=""
queue=""

shift
while [ ! -z "$1" ]
do
  case "$1" in
      -Q) queue="-q $2";shift;;
      -NJOBS) njobs=$2;shift;;
      -b) burnin=$2;shift;;
      -j) njumps=$2;shift;;
      -s) sampleevery=$2;shift;;
      *) other=$other" "$1;;
  esac
  shift
done

#Set options
opts="--fudge=$fudge --bi=$burnin --nj=$njumps --se=$sampleevery"
opts="$opts $other"

if [ $qsys -eq 0 ] && [ "x$SGE_ROOT" != "x" ]; then
	queue="-q $FSLGECUDAQ"
fi

#check that all required files exist

if [ ! -d $subjdir ]; then
	echo "subject directory $1 not found"
	exit 1
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}/data` -eq 0 ]; then
	echo "${subjdir}/data not found"
	exit 1
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif_brain_mask` -eq 0 ]; then
	echo "${subjdir}/nodif_brain_mask not found"
	exit 1
fi

if [ -e ${subjdir}.${modelname}/xfms/eye.mat ]; then
	echo "${subjdir} has already been processed: ${subjdir}.${modelname}." 
	echo "Delete or rename ${subjdir}.${modelname} before repeating the process."
	exit 1
fi

echo Making output directory structure

mkdir -p ${subjdir}.${modelname}/
mkdir -p ${subjdir}.${modelname}/diff_parts
mkdir -p ${subjdir}.${modelname}/logs
mkdir -p ${subjdir}.${modelname}/Dtifit
mkdir -p ${subjdir}.${modelname}/${premodel1}
mkdir -p ${subjdir}.${modelname}/${premodel1}/diff_parts
mkdir -p ${subjdir}.${modelname}/${premodel1}/logs
mkdir -p ${subjdir}.${modelname}/${premodel2}
mkdir -p ${subjdir}.${modelname}/${premodel2}/diff_parts
mkdir -p ${subjdir}.${modelname}/${premodel2}/logs
part=0

echo Copying files to output directory

${FSLDIR}/bin/imcp ${subjdir}/nodif_brain_mask ${subjdir}.${modelname}
if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif` = 1 ] ; then
    ${FSLDIR}/bin/fslmaths ${subjdir}/nodif -mas ${subjdir}/nodif_brain_mask ${subjdir}.${modelname}/nodif_brain
fi

# Specify Common Fixed Parameters
CFP_file=$subjdir.${modelname}/CFP
echo ${subjdir}/bvecs > $CFP_file
echo ${subjdir}/bvals >> $CFP_file

#Set more options
opts=$opts" --data=${subjdir}/data --maskfile=$subjdir.${modelname}/nodif_brain_mask --subjdir=$subjdir --forcedir --CFP=$CFP_file"

# Calculate S0 with the mean of the volumes with bval<50
bvals=`cat ${subjdir}/bvals`
mkdir -p ${subjdir}.${modelname}/temporal
pos=0
for i in $bvals; do 
    if [ $i -le 50 ]; then  
       	fslroi ${subjdir}/data  ${subjdir}.${modelname}/temporal/volume_$pos $pos 1    
    fi 
    pos=$(($pos + 1))
done
fslmerge -t ${subjdir}.${modelname}/temporal/S0s ${subjdir}.${modelname}/temporal/volume*
fslmaths ${subjdir}.${modelname}/temporal/S0s -Tmean ${subjdir}.${modelname}/S0
rm -rf ${subjdir}.${modelname}/temporal


##############################################################################
################################ First Dtifit  ###############################
##############################################################################
echo "Queue Dtifit"
PathDTI=${subjdir}.${modelname}/Dtifit
dtifit_command="${modelsdir}/utils/Run_dtifit.sh ${subjdir} ${subjdir}.${modelname} ${bindir}"
#SGE
dtifitProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathDTI/logs -N dtifit $dtifit_command`

#############################################################
##################### Grid Search Step ######################
#############################################################
echo "Queue Fitting process for model "${premodel1}
PathPreModel1=$subjdir.${modelname}/${premodel1}
PriorsFile=$modelsdir/${premodel1}/modelpriors

# Initialise Psi angle and beta2kappa
init_command="${modelsdir}/utils/initialise_Bingham.sh ${modelsdir}/utils ${PathDTI} ${subjdir}/nodif_brain_mask ${PathPreModel1}"
#SGE
initProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathPreModel1/logs -N ${premodel1}_initialisation -j $dtifitProcess $init_command`

# Fixed parameters
FixPFile=$PathPreModel1/FixP
echo ${PathPreModel1}/beta2kappa.nii.gz > $FixPFile
echo ${PathDTI}/dtifit_V1_th.nii.gz >> $FixPFile
echo ${PathDTI}/dtifit_V1_ph.nii.gz >> $FixPFile
echo ${PathPreModel1}/initialPsi.nii.gz >> $FixPFile
echo ${subjdir}.${modelname}/S0 >> $FixPFile

# Do GridSearch ... 
# Levenberg-Marquardt will find only local minimums. I need to check all these values
GridFile=$PathPreModel1/GridSearch
echo "search[0]=(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)" > $GridFile
echo "search[1]=(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)" >> $GridFile
echo "search[2]=(0.05,0.2,0.4,0.7,1,1.5,2,3,4,5,6,7,8,9,10,12,15,20,25,30,40,50,60)" >> $GridFile

# Split the dataset into parts
partsdir=$PathPreModel1/diff_parts
outputdir=$PathPreModel1
PreModel1Opts=$opts" --outputdir=$outputdir --partsdir=$partsdir --priors=$PriorsFile --FixP=$FixPFile --gridSearch=$GridFile --no_LevMar"

PreprocOpts=$PreModel1Opts" --idPart=0 --nParts=$njobs --logdir=$PathPreModel1/logs/preProcess"
preproc_command="$bindir/split_parts_${premodel1} $PreprocOpts"

#SGE
preProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathPreModel1/logs -N ${premodel1}_preproc -j $initProcess $preproc_command`

[ -f $PathPreModel1/commands.txt ] && rm $PathPreModel1/commands.txt
part=0
while [ $part -lt $njobs ]
do
    partzp=`$FSLDIR/bin/zeropad $part 4`
    
    Fitopts=$PreModel1Opts

    echo "$bindir/${premodel1} --idPart=$part --nParts=$njobs --logdir=$PathPreModel1/logs/${premodel1}_$partzp $Fitopts" >> $PathPreModel1/commands.txt
	    
    part=$(($part + 1))
done

#SGE
FitProcess=`${FSLDIR}/bin/fsl_sub $queue -N ${premodel1} -j $preProcess -t $PathPreModel1/commands.txt -l $PathPreModel1/logs`

PostprocOpts=$PreModel1Opts" --idPart=0 --nParts=$njobs --logdir=$PathPreModel1/logs/postProcess"
postproc_command="$bindir/merge_parts_${premodel1} $PostprocOpts"

#SGE
postProcess=`${FSLDIR}/bin/fsl_sub $queue -j $FitProcess -N ${premodel1}_postproc_gpu -l $PathPreModel1/logs $postproc_command`

#####################################################################
##################### Fit only Fractionalities ######################
#####################################################################
echo "Queue Fitting process for model "${premodel2}
PathPreModel2=$subjdir.${modelname}/${premodel2}
PriorsFile=$modelsdir/${premodel2}/modelpriors

# Initialise beta using beta2kappa from previous step
beta_command="${modelsdir}/utils/initialise_beta.sh ${PathPreModel1}/Param_2_samples ${PathPreModel1}/beta2kappa ${PathPreModel2}/beta"
#SGE
betaProcess=`${FSLDIR}/bin/fsl_sub $queue -j $postProcess -N ${premodel2}_beta_init -l $PathPreModel2/logs $beta_command`

# Create file to specify initialisation parameters (2 parameters: fiso,fintra)
InitializationFile=$PathPreModel2/InitializationParameters
echo $PathPreModel1/Param_0_samples > $InitializationFile
echo $PathPreModel1/Param_1_samples >> $InitializationFile

# Fixed parameters: kappa, beta, th, ph, psi, S0
FixPFile=$PathPreModel2/FixP
echo ${PathPreModel1}/Param_2_samples > $FixPFile
echo ${PathPreModel2}/beta >> $FixPFile
echo ${PathDTI}/dtifit_V1_th.nii.gz >> $FixPFile
echo ${PathDTI}/dtifit_V1_ph.nii.gz >> $FixPFile
echo ${PathPreModel1}/initialPsi.nii.gz >> $FixPFile
echo ${subjdir}.${modelname}/S0 >> $FixPFile


# Split the dataset into parts
partsdir=$PathPreModel2/diff_parts
outputdir=$PathPreModel2
PreModel2Opts=$opts" --outputdir=$outputdir --partsdir=$partsdir --init_params=$InitializationFile --priors=$PriorsFile --FixP=$FixPFile"

PreprocOpts=$PreModel2Opts" --idPart=0 --nParts=$njobs --logdir=$PathPreModel2/logs/preProcess"
preproc_command="$bindir/split_parts_${premodel2} $PreprocOpts"

#SGE
preProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathPreModel2/logs -N ${premodel2}_preproc -j $betaProcess $preproc_command`

[ -f $PathPreModel2/commands.txt ] && rm $PathPreModel2/commands.txt
part=0
while [ $part -lt $njobs ]
do
    partzp=`$FSLDIR/bin/zeropad $part 4`
    
    Fitopts=$PreModel2Opts

    echo "$bindir/${premodel2} --idPart=$part --nParts=$njobs --logdir=$PathPreModel2/logs/${premodel2}_$partzp $Fitopts" >> $PathPreModel2/commands.txt
	    
    part=$(($part + 1))
done

#SGE
FitProcess=`${FSLDIR}/bin/fsl_sub $queue -N ${premodel2} -j $preProcess -t $PathPreModel2/commands.txt -l $PathPreModel2/logs`

PostprocOpts=$PreModel2Opts" --idPart=0 --nParts=$njobs --logdir=$PathPreModel2/logs/postProcess"
postproc_command="$bindir/merge_parts_${premodel2} $PostprocOpts"

#SGE
postProcess=`${FSLDIR}/bin/fsl_sub $queue -j $FitProcess -N ${premodel2}_postproc_gpu -l $PathPreModel2/logs $postproc_command`

######################################################################################
######################### Fit all the parameters of the Model ########################
######################################################################################
echo "Queue Fitting process for model "${modelname}
PriorsFile=$modelsdir/${modelname}/modelpriors

# Create file to specify initialization parameters (7 parameters: fiso,fintra,kappa,beta,th,ph,psi)
InitializationFile=$subjdir.${modelname}/InitializationParameters
echo $PathPreModel2/Param_0_samples > $InitializationFile
echo $PathPreModel2/Param_1_samples >> $InitializationFile
echo ${PathPreModel1}/Param_2_samples >> $InitializationFile #kappa
echo ${PathPreModel2}/beta >> $InitializationFile #beta
echo ${PathDTI}/dtifit_V1_th.nii.gz >> $InitializationFile #th
echo ${PathDTI}/dtifit_V1_ph.nii.gz  >> $InitializationFile #ph
echo ${PathPreModel1}/initialPsi.nii.gz >> $InitializationFile

# S0 is a Fix parameter
FixPFile=$subjdir.${modelname}/FixP
echo  ${subjdir}.${modelname}/S0 > $FixPFile

partsdir=${subjdir}.${modelname}/diff_parts
outputdir=${subjdir}.${modelname}
ModelOpts=$opts" --outputdir=$outputdir --partsdir=$partsdir --init_params=$InitializationFile --FixP=$FixPFile --priors=$PriorsFile"

# Split the dataset into parts
PreprocOpts=$ModelOpts" --idPart=0 --nParts=$njobs --logdir=$subjdir.${modelname}/logs/preProcess"
preproc_command="$bindir/split_parts_${modelname} $PreprocOpts"

#SGE
preProcess=`${FSLDIR}/bin/fsl_sub $queue -j $postProcess -l ${subjdir}.${modelname}/logs -N ${modelname}_preproc $preproc_command`

[ -f ${subjdir}.${modelname}/commands.txt ] && rm ${subjdir}.${modelname}/commands.txt
part=0
while [ $part -lt $njobs ]
do
    partzp=`$FSLDIR/bin/zeropad $part 4`
    
    Fitopts=$ModelOpts
    
    echo "$bindir/${modelname} --idPart=$part --nParts=$njobs --logdir=$subjdir.${modelname}/logs/${modelname}_$partzp $Fitopts" >> ${subjdir}.${modelname}/commands.txt
    
    part=$(($part + 1))
done

#SGE
FitProcess=`${FSLDIR}/bin/fsl_sub $queue -N ${modelname} -j $preProcess -t ${subjdir}.${modelname}/commands.txt -l ${subjdir}.${modelname}/logs`

PostprocOpts=$ModelOpts" --idPart=0 --nParts=$njobs --logdir=$subjdir.${modelname}/logs/postProcess"
postproc_command="$bindir/merge_parts_${modelname} $PostprocOpts"

#SGE
postProcess=`${FSLDIR}/bin/fsl_sub $queue -j $FitProcess -N ${modelname}_postproc_gpu -l ${subjdir}.${modelname}/logs $postproc_command`

##########################################
### Calculate Dispersion Index & dyads ###
##########################################
finish_command="$modelsdir/${modelname}/${modelname}finish.sh ${subjdir}.${modelname}"
#SGE
finishProcess=`${FSLDIR}/bin/fsl_sub $queue -l ${subjdir}.${modelname}/logs -N ${modelname}_finish -j postProcess $finish_command`

endt=`date +%s`

runtime=$((endt-start))
echo Runtime $runtime
 
