# This is how I managed to compile CUDIMOT under Ubuntu 22.04
# -- yann.leprince@cea.fr 2026-03-18

# Preparatory step: install the CUDA Toolkit version 10.2 to
# /volatile/opt/cuda-10.2. This is necessary to get a working libcudart and
# libdevcudart (they are not included in FSL). We still have to use nvcc from
# the system installation of the cuda toolkit (nvidia-cuda-toolkit version
# 11.5.1-1ubuntu1) because cuda-10.2/bin/nvcc does not support GCC > 8.
wget https://developer.download.nvidia.com/compute/cuda/10.2/Prod/local_installers/cuda_10.2.89_440.33.01_linux.run
sh cuda_10.2.89_440.33.01_linux.run --silent --toolkit --toolkitpath=/volatile/opt/cuda-10.2 --override --librarypath=/volatile/opt/cuda-10.2-libs

# Set up FSL envionment (same as fsl_init)
export FSLDIR=/drf/local/fsl-6.0.7.6
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
. "${FSLDIR}/etc/fslconf/fsl.sh"

# Compile
mkdir -p fsldevdir
make modelname=NODDI_Watson FSLDEVDIR=$(pwd)/fsldevdir FSLCONFDIR=${FSLDIR}/config CUDA=/volatile/opt/cuda-10.2 NVCC=/usr/bin/nvcc
