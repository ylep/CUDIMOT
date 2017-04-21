/* getPredictedSignal.cu

   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
   
   Copyright (C) 2005 University of Oxford */

/* CCOPYRIGHT */

// Method to calculate the model predicted signal on the GPU of a group of voxels given the model parameters.

#include "getPredictedSignal.h"
#include "functions_gpu.h"
#include "modelparameters.h"
#include "modelfunctions.h"

namespace Cudimot{

#define THREADS_VOXEL 32 // Multiple of 32: Threads collaborating to compute a voxel. Do not change this, otherwise Synchronization will be needed
  
  template <typename T>
  __global__ void getPredictedSignal_kernel(
					    int nvox, // nvoxels
					    int nmeas, // nmeasurements
					    int CFP_Tsize, //size*M-measurements
					    int FixP_Tsize, // fixed params: size*N-voxels
					    T* parameters, // model parameters 
					    T* CFP_global, // common fixed model parameters
					    T* FixP, // fixed model parameters
					    T* PredictedSignal)
  {
    // 1 block of threads process several voxels
    // Each warp processes 1 voxel
    int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
    int idVOX_inBlock =  threadIdx.x/THREADS_VOXEL;
    int idSubVOX= threadIdx.x%THREADS_VOXEL;
    bool leader = (idSubVOX==0);  // Some steps are performed by only one thread of the warp
    
    ////////// DYNAMIC SHARED MEMORY ///////////
    extern __shared__ double shared[];	      	//Size:
    
    T* CFP = (T*)shared;		 	//nmeas*CMP_Tsize
    T* params = (T*) &CFP[nmeas*CFP_Tsize]; 	//NPARAMS*VOXELS_BLOCK
    ////////////////////////////////////////////
    
    /// Copy common fixed model parameters to Shared Memory ///
    if(threadIdx.x==0){ // only one thread of the whole block. Common to all voxels
      for(int i=0;i<nmeas*CFP_Tsize;i++){
	CFP[i]=CFP_global[i];
      }
    }
    __syncthreads();
    ///////////////////////////////////////////////////////////
    
    ///////// each voxel/warp of the block points to its data///////////
    params = &params[idVOX_inBlock*NPARAMS];
    PredictedSignal = &PredictedSignal[idVOX*nmeas]; //Global
    FixP = &FixP[idVOX*FixP_Tsize]; // Global memory
    ////////////////////////////////////////////////////////////////////

    /// Ititialise shared values of each voxel: only the leader///
    if(leader){ 
      #pragma unroll
      for(int par=0;par<NPARAMS;par++){
	params[par]=parameters[idVOX*NPARAMS+par];
      }
    }
    __syncthreads();
    //sync?   

        
    int idMeasurement=idSubVOX;
    int nmeas2compute = nmeas/THREADS_VOXEL;
    if (idSubVOX<(nmeas%THREADS_VOXEL)) nmeas2compute++;
    
    for(int iter=0;iter<nmeas2compute;iter++){
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred=Predicted_Signal(NPARAMS,params,myCFP,FixP);
      PredictedSignal[idMeasurement]=pred;
      idMeasurement+=THREADS_VOXEL;
    }
  }
  
  
  template <typename T>
  getPredictedSignal<T>::getPredictedSignal(){}
  
  template <typename T>
  void getPredictedSignal<T>::run(
				   int nvox, int nmeas,
				   int CFP_size, int FixP_size,
				   T* params, T* CFP, T* FixP,
				   T* PredictedSignal) 
  {
    
    long int amount_shared_mem = 0;
    amount_shared_mem += (nmeas*CFP_size)*sizeof(T); // CFP
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Parameters
    
    cout << "Shared Memory used in PredictedSignal kernel: " << amount_shared_mem << endl;
    
    int threads_block = VOXELS_BLOCK * THREADS_VOXEL;
    int nblocks=(nvox/VOXELS_BLOCK);
    if(nvox%VOXELS_BLOCK) nblocks++;
    
    getPredictedSignal_kernel<T><<<nblocks,threads_block,amount_shared_mem>>>(nvox,nmeas,CFP_size,FixP_size,params,CFP,FixP,PredictedSignal);
    sync_check("getPredictedSignal Kernel");
  }
  
  // Explicit Instantiations of the template
  template class getPredictedSignal<float>;
  template class getPredictedSignal<double>;
}
