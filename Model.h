#ifndef CUDIMOT_MODEL_H_INCLUDED
#define CUDIMOT_MODEL_H_INCLUDED

/**
 *
 * \class Model
 *
 * \brief A class for managing the Model specified by the Designer
 *
 * This class contains the information (parameters) of the diffusion MRI model specified by a designer.
 *
 * \author Moises Hernandez-Fernandez - FMRIB Image Analysis Group
 *
 * \date March 2017
 *
 *
 *  Copyright (C) 2005 University of Oxford 
 */

/* CCOPYRIGHT */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "utils.h"
#include "cudimotoptions.h"
#include "mymodels/mymodel.h"

namespace Cudimot{
  
  template <typename T>
  class Parameters;
  
  template <typename T>
  class Model{
    friend class Parameters<T>;
    
  private:
    /**
     * Number of parameters in the model (to estimate)
     */
    int nparams;

    /**
     * Vector with a default values for initializing the model parameters
     */
    vector<T> params_init;

    /**
     * If activated, the Designer has provided default values for initializing the model parameters
     */
    bool provided_params_init;
    
    ////////////////////////////////////////////////////////////////////////
    /// FIXED PARAMETERS different for each voxel, nifti volumes (ex.T1) ///
    ////////////////////////////////////////////////////////////////////////
 
    /**
     * Number of Fixed Parameters in the model (different for each voxel)
     */
    int nFP; 

    
    /////////////////////////////////////////////////////////////////
    /// COMMON (to all voxels) FIXED PARAMETERS (bvals,bvecs,...) ///
    /////////////////////////////////////////////////////////////////
    
    /**
     * Number of Common Fixed Parameters in the model (common to all voxels)
     */
    int nCFP;
    
    /**
     * Total size of Common Fixed Parameters (without counting measurements)
     */
    int CFP_Tsize; 

    /**
     * Vector with the size of each Common Fixed Parameter (without counting measurements)
     */
    vector<int> CFP_sizes;
    /////////////////////////////////////////////////////////////////
    
    ///////////////////////
    /// Priors for MCMC ///
    ///////////////////////

    /**
     * This method reads the model configuration file and sets the class attributes. File name must be given by the Designer.
     */
    void Modelparser();

    /**
     * This method writes in a file the model information (for the final User)
     */
    void WriteModelSpecification(string file);
    
  public:
    /**
     * Constructor
     */
    Model();

    /**
     * Destructor
     */
    ~Model();

    /** 
     * @return The number of parameters in the model (to estimate)
     */
    int getNparams();

    /**
     * @return Indicates if the Designer has provided default values for initializing the model parameters
     */
    bool initProvided();

    /**
     * @param A number to identify a parameter of the model 
     * @return The default value for initializing one of the model parameters
     */
    T getParam_init(int id_param);
  };
}

#endif
