#ifndef CUDIMOT_LEVENBERG_MARQUARDT_H_INCLUDED
#define CUDIMOT_LEVENBERG_MARQUARDT_H_INCLUDED

/**
 *
 * \class Levenberg_Marquardt
 *
 * \brief A class for fitting a model to some dataset on a GPU using Levenberg-Marquardt algorithm
 *
 * \author Moises Hernandez-Fernandez - FMRIB Image Analysis Group
 *
 * \date March 2017
 *
 *
 *  Copyright (C) 2005 University of Oxford 
 */

/* CCOPYRIGHT */

#include "compileOptions/grid.h"

namespace Cudimot{
  template <typename T>
  class Levenberg_Marquardt{

  private:
    /**
     * Maximum number of iterations
     */
    int max_iterations;
    
    /**
     * If activeted, use Marquardt contribution
     */
    bool Marquardt;
 
  public:
    
    /**
     * Constructor
     */
    Levenberg_Marquardt();

    /**
     * Run Levenberg-Marquard on the GPU
     * @param nvox Number of voxels in the dataset
     * @param nmeas Number of measurements in the dataset
     * @param CFP_size Size (without measurements) of the common fixed parameters of the model
     * @param meas Measurements of all the voxels of the dataset (on GPU)
     * @param params Value of the parameters to estimate of the model for all the voxels of the dataset (on GPU)
     * @param CFP Common (to all the voxels) fixed parameters of the model (on GPU). CFP_size*nmeas
     */
    void run( int nvox, int nmeas,int CFP_size,
	      T* meas, T* params,
	      T* CFP);
  };
}

#endif
