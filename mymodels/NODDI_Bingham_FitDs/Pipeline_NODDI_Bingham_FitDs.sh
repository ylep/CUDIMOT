#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT
#
# Pipeline for fitting NODDI-Bingham

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
Usage() {
    echo ""
    echo "Usage: Pipeline_NODDI_Bingham_FitDs.sh <subject_directory> [options]"
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
    echo "--BIC (if you want to calculate BIC)"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

modelname=NODDI_Bingham_FitDs
step1=GS_FixFibreOrientation

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FSLDIR}/lib
BASEDIR=$(dirname "$0")
modelsdir=`make_absolute $BASEDIR/..`

if [ -z "$bindir" ];then
    echo "Please set the variable bindir with the path to your CUDIMOT binary files"
    exit 1
fi

subjdir=`make_absolute $1`
subjdir=`echo $subjdir | sed 's/\/$/$/g'`

echo "---------------------------------------------------------------------------------"
echo "------------------------------------ CUDIMOT ------------------------------------"
echo "----------------------------- MODEL: $modelname -----------------------------"
echo "---------------------------------------------------------------------------------"
echo subjectdir is $subjdir

start=`date +%s`

#parse option arguments
njobs=4
fudge=1
burnin=1000
njumps=1250
sampleevery=25
other=""
queue=""
lastStepModelOpts=""

shift
while [ ! -z "$1" ]
do
  case "$1" in
      -Q) queue="-q $2";shift;;
      -NJOBS) njobs=$2;shift;;
      -b) burnin=$2;shift;;
      -j) njumps=$2;shift;;
      -s) sampleevery=$2;shift;;
      --runMCMC) lastStepModelOpts=$lastStepModelOpts" --runMCMC";;
      --BIC) lastStepModelOpts=$lastStepModelOpts" --BIC";;
      *) other=$other" "$1;;
  esac
  shift
done

#Set options
opts="--fudge=$fudge --bi=$burnin --nj=$njumps --se=$sampleevery"
opts="$opts $other"

if [ "x$SGE_ROOT" != "x" ]; then
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
mkdir -p ${subjdir}.${modelname}/${step1}
mkdir -p ${subjdir}.${modelname}/${step1}/diff_parts
mkdir -p ${subjdir}.${modelname}/${step1}/logs

part=0

echo Copying files to output directory

${FSLDIR}/bin/imcp ${subjdir}/nodif_brain_mask ${subjdir}.${modelname}
if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif` = 1 ] ; then
    ${FSLDIR}/bin/fslmaths ${subjdir}/nodif -mas ${subjdir}/nodif_brain_mask ${subjdir}.${modelname}/nodif_brain
fi

# Specify Common Fixed Parameters
CFP_file=$subjdir.${modelname}/CFP
cp ${subjdir}/bvecs $subjdir.${modelname}
cp ${subjdir}/bvals $subjdir.${modelname}
echo  $subjdir.${modelname}/bvecs > $CFP_file
echo  $subjdir.${modelname}/bvals >> $CFP_file

# Specify Priors
PriorsFile=$modelsdir/$modelname/modelpriors

#Set more options
opts=$opts" --data=${subjdir}/data --maskfile=$subjdir.${modelname}/nodif_brain_mask --forcedir --CFP=$CFP_file"

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

# Specify Fixed parameters: S0
FixPFile=${subjdir}.${modelname}/FixP
echo ${subjdir}.${modelname}/S0 >> $FixPFile

##############################################################################
################################ First Dtifit  ###############################
##############################################################################
echo "Queue Dtifit"
PathDTI=${subjdir}.${modelname}/Dtifit
dtifit_command="${modelsdir}/utils/Run_dtifit.sh ${subjdir} ${subjdir}.${modelname} ${bindir}"
#SGE
dtifitProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathDTI/logs -N dtifit $dtifit_command`

## Model Parameters: fiso, fintra, kappa, beta, th, ph  psi ##
##############################################################
################ Grid Search + Fix th & ph ###################
##############################################################
echo "Queue $step1 process"
PathStep1=$subjdir.${modelname}/${step1}

# Initialise Psi angle and beta2kappa
init_command="${modelsdir}/utils/initialise_Bingham.sh ${modelsdir}/utils ${PathDTI} ${subjdir}/nodif_brain_mask ${PathStep1}"
#SGE
initProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathStep1/logs -N ${modelname}_initialisation -j $dtifitProcess $init_command`

# Create file to specify initialisation parameters
InitializationFile=$PathStep1/InitializationParameters
echo "" > $InitializationFile #fiso
echo "" >> $InitializationFile #fintra
echo "" >> $InitializationFile #kappa
echo "" >> $InitializationFile #beta
echo ${PathDTI}/dtifit_V1_th.nii.gz >> $InitializationFile #th
echo ${PathDTI}/dtifit_V1_ph.nii.gz >> $InitializationFile #ph
echo ${PathStep1}/initialPsi.nii.gz >> $InitializationFile #psi
echo "" >> $InitializationFile #Diso
echo "" >> $InitializationFile #Dparallel

# Do GridSearch (fiso,fintra,kappa,betta)
GridFile=$PathStep1/GridSearch
echo "search[0]=(0.01,0.1,0.2,0.3,0.8)" > $GridFile #fiso
echo "search[1]=(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)" >> $GridFile #fintra
echo "search[2]=(1,2,3,6,10,15,30,40,50)" >> $GridFile #kappa
echo "search[3]=(0.9,5,8,14,29,38)" >> $GridFile #beta

partsdir=$PathStep1/diff_parts
outputdir=$PathStep1
Step1Opts=$opts" --outputdir=$outputdir --partsdir=$partsdir --priors=$PriorsFile --FixP=$FixPFile --gridSearch=$GridFile --init_params=$InitializationFile --fixed=4,5,7,8"

postproc=`${modelsdir}/utils/jobs_wrapper.sh $PathStep1 $initProcess $modelname GS $njobs $Step1Opts`

######################################################################################
######################### Fit all the parameters of the Model ########################
######################################################################################
echo "Queue Fitting process"

# Create file to specify initialization parameters (5 parameters: fiso,fintra,kappa,th,ph)
InitializationFile=$subjdir.${modelname}/InitializationParameters
echo $PathStep1/Param_0_samples > $InitializationFile #fiso
echo $PathStep1/Param_1_samples >> $InitializationFile #fintra
echo $PathStep1/Param_2_samples >> $InitializationFile #kappa
echo $PathStep1/Param_3_samples >> $InitializationFile #beta
echo $PathStep1/Param_4_samples >> $InitializationFile #th
echo $PathStep1/Param_5_samples >> $InitializationFile #ph
echo $PathStep1/Param_6_samples >> $InitializationFile #psi
echo $PathStep1/Param_7_samples >> $InitializationFile #Diso
echo $PathStep1/Param_8_samples >> $InitializationFile #Dparallel

partsdir=${subjdir}.${modelname}/diff_parts
outputdir=${subjdir}.${modelname}
ModelOpts=$opts" --outputdir=$outputdir --partsdir=$partsdir --priors=$PriorsFile --FixP=$FixPFile --init_params=$InitializationFile $lastStepModelOpts --fixed=0,1"

postproc=`${modelsdir}/utils/jobs_wrapper.sh ${subjdir}.${modelname} $postproc $modelname FitProcess $njobs $ModelOpts`

#########################################
### Calculate Dispersion Index & dyads ###
##########################################
finish_command="${modelsdir}/${modelname}/${modelname}_finish.sh ${subjdir}.${modelname} ${modelsdir}"
#SGE
finishProcess=`${FSLDIR}/bin/fsl_sub $queue -l ${subjdir}.${modelname}/logs -N ${modelname}_finish -j $postproc $finish_command`

endt=`date +%s`
runtime=$((endt-start))
echo Runtime $runtime