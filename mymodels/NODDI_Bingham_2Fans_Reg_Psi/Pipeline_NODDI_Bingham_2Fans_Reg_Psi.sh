#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT
#
# Pipeline for fitting NODDI-Bingham 2 Fans with  - regularization 

bindir=${FSLDEVDIR}/bin

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
    echo "Usage: Pipeline_NODDI_Bingham_2Fans_Reg_Psi.sh <subject_directory> [options]"
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

modelname=NODDI_Bingham_2Fans_Reg_Psi
step1=GS_FixFibreOrientation

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FSLDIR}/lib

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
opts="--bi=$burnin --nj=$njumps --se=$sampleevery"
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
cp ${subjdir}/bvecs $subjdir.${modelname}
cp ${subjdir}/bvals $subjdir.${modelname}
cp ${subjdir}/bvecs_bvecs $subjdir.${modelname}
cp ${subjdir}/bvals_bvals $subjdir.${modelname}
cp ${subjdir}/LR_map $subjdir.${modelname}
CFP_file_Reg=$subjdir.${modelname}/CFP_Reg
echo  $subjdir.${modelname}/bvecs_bvecs > $CFP_file_Reg
echo  $subjdir.${modelname}/bvals_bvals >> $CFP_file_Reg
echo  $subjdir.${modelname}/LR_map >> $CFP_file_Reg
CFP_file=$subjdir.${modelname}/CFP
echo  $subjdir.${modelname}/bvecs > $CFP_file
echo  $subjdir.${modelname}/bvals >> $CFP_file

#Set more options
opts_Reg=$opts" --data=${subjdir}/data_dataSM --maskfile=$subjdir.${modelname}/nodif_brain_mask --forcedir --CFP=$CFP_file_Reg"
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
dtifit_command="${bindir}/Run_dtifit.sh ${subjdir} ${subjdir}.${modelname} ${bindir}"
#SGE
dtifitProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathDTI/logs -N dtifit $dtifit_command`

##############################################################
################ Grid Search + Fix th & ph ###################
############# Using normal NODDI_Bingham with 2 Fans #########
##############################################################
echo "Queue $step1 process"
PathStep1=$subjdir.${modelname}/${step1}

# Initialise Psi angle and beta2kappa
init_command="${bindir}/initialise_Bingham.sh ${bindir} ${PathDTI} ${subjdir}/nodif_brain_mask ${PathStep1}"
#SGE
initProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathStep1/logs -N ${modelname}_initialisation -j $dtifitProcess $init_command`

# Create file to specify initialisation parameters
InitializationFile=$PathStep1/InitializationParameters
echo "" > $InitializationFile #fiso
echo "" >> $InitializationFile #fintra
echo "" >> $InitializationFile #ffib2
echo "" >> $InitializationFile #kappa1
echo "" >> $InitializationFile #beta1
echo ${PathDTI}/dtifit_V1_th.nii.gz >> $InitializationFile #th
echo ${PathDTI}/dtifit_V1_ph.nii.gz >> $InitializationFile #ph
echo ${PathStep1}/initialPsi.nii.gz >> $InitializationFile #psi
echo "" >> $InitializationFile #kappa2
echo "" >> $InitializationFile #beta2
echo ${PathDTI}/dtifit_V2_th.nii.gz >> $InitializationFile #th2
echo ${PathDTI}/dtifit_V2_ph.nii.gz >> $InitializationFile #ph2
echo "" >> $InitializationFile #psi2

# Do GridSearch (fiso,fintra,kappa,betta)
GridFile=$PathStep1/GridSearch
echo "search[0]=(0.01,0.1,0.3,0.8)" > $GridFile #fiso
echo "search[1]=(0.3,0.4,0.5,0.6,0.8,0.9,1.0)" >> $GridFile #fintra
echo "search[2]=(0.1,0.2,0.3)" >> $GridFile #ffib2
echo "search[3]=(1,2,3,6,10,15,30,40,50)" >> $GridFile #kappa1
echo "search[4]=(0.9,5,8,14,29,38)" >> $GridFile #beta1
echo "search[8]=(1,2,3,6,10,15,30,40,50)" >> $GridFile #kappa2
echo "search[9]=(0.9,5,8,14,29,38)" >> $GridFile #beta2

partsdir=$PathStep1/diff_parts
outputdir=$PathStep1
Step1Opts=$opts" --outputdir=$outputdir --partsdir=$partsdir --FixP=$FixPFile --gridSearch=$GridFile --init_params=$InitializationFile --fixed=5,6,10,11"

postproc=`${bindir}/jobs_wrapper.sh $PathStep1 $initProcess NODDI_Bingham_2Fans GS $njobs $Step1Opts`

######################################################################################
######################### Fit all the parameters of the Model ########################
######################################################################################
echo "Queue Fitting process"

# Create file to specify initialization parameters (5 parameters: fiso,fintra,kappa,th,ph)
InitializationFile=$subjdir.${modelname}/InitializationParameters
echo $PathStep1/Param_0_samples > $InitializationFile #fiso
echo $PathStep1/Param_1_samples >> $InitializationFile #fintra
echo $PathStep1/Param_2_samples >> $InitializationFile #ffib2
echo $PathStep1/Param_3_samples >> $InitializationFile #kappa1
echo $PathStep1/Param_4_samples >> $InitializationFile #beta1
echo $PathStep1/Param_5_samples >> $InitializationFile #th1
echo $PathStep1/Param_6_samples >> $InitializationFile #ph1
echo $PathStep1/Param_7_samples >> $InitializationFile #psi1
echo $PathStep1/Param_8_samples >> $InitializationFile #kappa2
echo $PathStep1/Param_9_samples >> $InitializationFile #beta2
echo $PathStep1/Param_10_samples >> $InitializationFile #th2
echo $PathStep1/Param_11_samples >> $InitializationFile #ph2
echo $PathStep1/Param_12_samples >> $InitializationFile #psi2
echo $PathStep1/Param_0_samples >> $InitializationFile #LR_fiso
echo $PathStep1/Param_1_samples >> $InitializationFile #LR_fintra
echo $PathStep1/Param_2_samples >> $InitializationFile #LR_ffib2
echo $PathStep1/Param_3_samples >> $InitializationFile #LR_kappa1
echo $PathStep1/Param_4_samples >> $InitializationFile #LR_beta1
echo $PathStep1/Param_8_samples >> $InitializationFile #LR_kappa2
echo $PathStep1/Param_9_samples >> $InitializationFile #LR_beta2


partsdir=${subjdir}.${modelname}/diff_parts
outputdir=${subjdir}.${modelname}
ModelOpts=$opts_Reg" --outputdir=$outputdir --partsdir=$partsdir --FixP=$FixPFile --init_params=$InitializationFile $lastStepModelOpts"

postproc=`${bindir}/jobs_wrapper.sh ${subjdir}.${modelname} $postproc $modelname FitProcess $njobs $ModelOpts`

#########################################
### Calculate Dispersion Index & dyads ###
##########################################
finish_command="${bindir}/${modelname}_finish.sh ${subjdir}.${modelname}"
#SGE
finishProcess=`${FSLDIR}/bin/fsl_sub $queue -l ${subjdir}.${modelname}/logs -N ${modelname}_finish -j $postproc $finish_command`

endt=`date +%s`
runtime=$((endt-start))
echo Runtime $runtime