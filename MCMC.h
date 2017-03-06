#ifndef CUDIMOT_MCMC_H_INCLUDED
#define CUDIMOT_MCMC_H_INCLUDED

/**
 *
 * \class MVMV
 *
 * \brief A class for fitting a model to some dataset on a GPU using MCMC algorithm
 *
 * \author Moises Hernandez-Fernandez - FMRIB Image Analysis Group
 *
 * \date March 2017
 *
 *
 *  Copyright (C) 2005 University of Oxford 
 */

/* CCOPYRIGHT */

#include <curand_kernel.h>
#include "compileOptions/grid.h"

namespace Cudimot{
  template <typename T>
  class MCMC{

  private:

    /** 
     * The number of voxels to fit in a part
     */
    int nvoxFit_part;

    /** 
     * Total number of jumps in MCMC to be discarded at the beginning (default is 5000)
     */
    int nburnin; 
    
    /** 
     * Number of jumps in MCMC after burnin period (default is 1250)
     */
    int njumps;
    
    /** 
     * Number of samples to record per parameter and voxel
     */
    int nsamples;

    /** 
     * Number of jumps between sample recording
     */
    int sampleevery;

    /** 
     * Number of jumps between updates of the proposal standard deviations
     */
    int updateproposalevery;

    /** 
     * If activated, use Rician noise modelling
     */
    bool RicianNoise;

    /**
     * State of several random number generators on the GPU
     */
    curandState* randStates;
    

  public:

    /**
     * Constructor
     * @param nvoxFitpart Number of voxel of the data to fit
     */
    MCMC(int nvoxFitpart);

    /**
     * Run MCMC algorithm on the GPU
     * @param nvox Number of voxels in the dataset
     * @param nmeas Number of measurements in the dataset
     * @param CFP_size Size (without measurements) of the common fixed parameters of the model
     * @param meas Measurements of all the voxels of the dataset (on GPU)
     * @param params Value of the parameters to estimate of the model for all the voxels of the dataset (on GPU)
     * @param CFP Common (to all the voxels) fixed parameters of the model (on GPU). CFP_size*nmeas
     * @param samples Samples of the parameters estimation will be stored here (on the GPU)
     */
    void run( int nvox, int nmeas, int CMP_size,
	      T* meas, T* params, T* CFP,
	      T* samples);
  };
}

#endif
