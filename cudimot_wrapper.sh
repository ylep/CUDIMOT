#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT

modelname=dtimodel
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FSLDIR}/lib
bindir=/home/fs0/moisesf/scratch/THESIS/microstructure_models/CUDIMOT/bin

Usage() {
    echo ""
    echo "Usage: cuditmot <subject_directory> <file_Common_Fixed_Parameters> [options]"
    echo ""
    echo "expects to find data and nodif_brain_mask in subject directory"
    echo ""
    echo "<options>:"
    #echo "-QSYS (Queue System, 0 use fsl_sub: FMRIB, 1 TORQUE (default): WashU)"
    echo "-Q (name of the GPU(s) queue, default cuda.q (defined in environment variable: FSLGECUDAQ)"
    #echo "-Q (name of the GPU(s) queue, default cuda.q for QSYS=0 and no queue for QSYS=1)"	
    echo "-NJOBS (number of jobs to queue, the data is divided in NJOBS parts, usefull for a GPU cluster, default 4)"
    echo "-w (ARD weight, more weight means less secondary fibres per voxel, default 1)"
    echo "-b (burnin period, default 5000)"
    echo "-j (number of jumps, default 1250)"
    echo "-s (sample every, default 25)"
    #echo "-cfp (File with alist of names and files for Common Fixed Parameters in this model)"
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

[ "$2" = "" ] && Usage

subjdir=`make_absolute $1`
subjdir=`echo $subjdir | sed 's/\/$/$/g'`
common_fixed_params=$2

echo "---------------------------------------------"
echo "------------------ CUDIMOT ------------------"
echo "---------------------------------------------"
echo subjectdir is $subjdir

start=`date +%s`

#parse option arguments
qsys=0
njobs=4
fudge=1
burnin=5000
njumps=1250
sampleevery=25
other=""
queue=""

shift;shift
while [ ! -z "$1" ]
do
  case "$1" in
      -QSYS) qsys=$2;shift;;
      -Q) queue="-q $2";shift;;
      -NJOBS) njobs=$2;shift;;
      -w) fudge=$2;shift;;
      -b) burnin=$2;shift;;
      -j) njumps=$2;shift;;
      -s) sampleevery=$2;shift;;
      *) other=$other" "$1;;
  esac
  shift
done

#Set options
opts="--fudge=$fudge --bi=$burnin --nj=$njumps --se=$sampleevery --CFP=$common_fixed_params"
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

if [ -e ${subjdir}.cudimot/xfms/eye.mat ]; then
	echo "${subjdir} has already been processed: ${subjdir}.cudimot." 
	echo "Delete or rename ${subjdir}.cudimot before repeating the process."
	exit 1
fi

echo Making cudimot directory structure

mkdir -p ${subjdir}.cudimot/
mkdir -p ${subjdir}.cudimot/diff_parts
mkdir -p ${subjdir}.cudimot/logs
part=0

#mkdir -p ${subjdir}.cudimot/logs/logs_gpu
partsdir=${subjdir}.cudimot/diff_parts

echo Copying files to cudimot directory

${FSLDIR}/bin/imcp ${subjdir}/nodif_brain_mask ${subjdir}.cudimot
if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif` = 1 ] ; then
    ${FSLDIR}/bin/fslmaths ${subjdir}/nodif -mas ${subjdir}/nodif_brain_mask ${subjdir}.cudimot/nodif_brain
fi

#Set more default options
opts=$opts" --data=${subjdir}/data --maskfile=$subjdir.cudimot/nodif_brain_mask --subjdir=$subjdir --partsdir=$partsdir --outputdir=$subjdir.cudimot --forcedir"

# Split the dataset in parts
echo Pre-processing stage
	PreprocOpts=$opts" --idPart=0 --nParts=$njobs --logdir=$subjdir.cudimot/logs/preProcess"
	preproc_command="$bindir/split_parts $PreprocOpts"

	if [ $qsys -eq 0 ]; then
		#SGE
		preProcess=`${FSLDIR}/bin/fsl_sub -T 40 -l ${subjdir}.cudimot/logs -N cudimot_preproc $preproc_command`
	else
		#TORQUE
		echo $prepoc_command > ${subjdir}.cudimot/temp
		torque_command="qsub -V $queue -l nodes=1:ppn=1:gpus=1,walltime=00:40:00 -N cudimot_preproc -o ${subjdir}.cudimot/logs -e ${subjdir}.cudimot/logs"
		preProcess=`exec $torque_command ${subjdir}.cudimot/temp | awk '{print $1}' | awk -F. '{print $1}'`
		rm ${subjdir}.cudimot/temp
		sleep 10
	fi



echo Queuing Fitting model processing stage

	[ -f ${subjdir}.cudimot/commands.txt ] && rm ${subjdir}.cudimot/commands.txt

	part=0
	while [ $part -lt $njobs ]
	do
	    	partzp=`$FSLDIR/bin/zeropad $part 4`
	    
		Fitopts=$opts

		#${FSLDIR}/bin/cudimot
		echo "$bindir/${modelname} --idPart=$part --nParts=$njobs --logdir=$subjdir.cudimot/logs/cudimot_$partzp $Fitopts" >> ${subjdir}.cudimot/commands.txt
	    
	    	part=$(($part + 1))
	done

	if [ $qsys -eq 0 ]; then
		#SGE
		FitProcess=`${FSLDIR}/bin/fsl_sub $queue -N cudimot -j $preProcess -t ${subjdir}.cudimot/commands.txt -l ${subjdir}.cudimot/logs `
	else
		#TORQUE
		taskfile=${subjdir}.cudimot/commands.txt
		echo "command=\`cat "$taskfile" | head -\$PBS_ARRAYID | tail -1\` ; exec \$command" > ${subjdir}.cudimot/temp
		tasks=`wc -l $taskfile | awk '{print $1}'`
		sge_tasks="-t 1-$tasks"
		#PBS -t x-y: x and y are the array bounds
		torque_command="qsub -V $queue -l nodes=1:ppn=1:gpus=1,walltime=3:00:00,pmem=16gb -N cudimot -o ${subjdir}.cudimot/logs -e ${subjdir}.cudimot/logs -W depend=afterok:$preProcess $sge_tasks"
		FitProcess=`exec $torque_command ${subjdir}.cudimot/temp | awk '{print $1}' | awk -F. '{print $1}'`
		rm ${subjdir}.cudimot/temp
		sleep 10
	fi


echo Queuing Post-processing stage
# Needs the parent directory where all the output parts are stored $subjdir.cudimot
PostprocOpts=$opts" --idPart=0 --nParts=$njobs --logdir=$subjdir.cudimot/logs/postProcess"

#${FSLDIR}/bin/
postproc_command="$bindir/merge_parts $PostprocOpts"

if [ $qsys -eq 0 ]; then
	#SGE
	postProcess=`${FSLDIR}/bin/fsl_sub $queue -T 120 -j $FitProcess -N cudimot_postproc_gpu -l ${subjdir}.cudimot/logs $postproc_command`
else
	#TORQUE
	echo $post_command > ${subjdir}.cudimot/temp
	torque_command="qsub -V $queue -l nodes=1:ppn=1:gpus=1,walltime=00:40:00 -N cudimot_postproc_gpu -o ${subjdir}.cudimot/logs -e ${subjdir}.cudimot/logs -W depend=afterokarray:$FitProcess"
	postProcess=`exec $torque_command ${subjdir}.cudimot/temp | awk '{print $1}' | awk -F. '{print $1}'`
        rm ${subjdir}.cudimot/temp
	sleep 10
fi

endt=`date +%s`

runtime=$((endt-start))
echo Runtime $runtime

