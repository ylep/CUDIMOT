/* Levenberg_Marquardt.cu

   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
   
   Copyright (C) 2005 University of Oxford */

/* CCOPYRIGHT */

#include "Levenberg_Marquardt.h"
#include "utils.h"
#include "cudimotoptions.h"
#include "mymodels/mymodel.h"
#include MYMODEL_FUNCTIONS

namespace Cudimot{

#define CFTOL 1.0e-8
#define LTOL 1.0e20
#define EPS_gpu 2.0e-16 //Losely based on NRinC 20.1
  
  __device__ inline bool zero_cf_diff_conv(double* cfo,double* cfn){
    return(2.0*fabs(*cfo-*cfn) <= CFTOL*(fabs(*cfo)+fabs(*cfn)+EPS_gpu));
  }
  
  template <typename T>
  __device__ void Cost_Function(int idSubVOX,
				int nmeas,
				int CFP_Tsize,
				T* measurements,
				T* parameters,
				T* fixed_params,
				T* CFP,
				double* result)
  {
    int idMeasurement=idSubVOX;
    T accumulated_error=0;

    int nmeas2compute = nmeas/THREADS_VOXEL;
    if (idSubVOX<(nmeas%THREADS_VOXEL)) nmeas2compute++;

    for(int iter=0;iter<nmeas2compute;iter++){
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred_error=Predicted_Signal(NPARAMS,parameters,fixed_params,myCFP);
      
      pred_error=pred_error-measurements[idMeasurement];
      accumulated_error+=pred_error*pred_error;
      idMeasurement+=THREADS_VOXEL;
    }
    
    for (int offset=THREADS_VOXEL/2; offset>0; offset/=2){
      accumulated_error+= __shfl_down(accumulated_error,offset);
    }
    if(idSubVOX==0){
      *result=accumulated_error;
    }
  }
  
  template <typename T>
  __device__ void Calculate_Gradient(int idSubVOX,
				     int nmeas,
				     int CFP_Tsize,
				     T* measurements,
				     T* parameters,
				     T* fixed_params,
				     T* CFP,
				     T* Gradient)
  {
    int idMeasurement=idSubVOX;
    T myderivatives[NPARAMS];
    
    int max_iters = nmeas/THREADS_VOXEL;
    if(nmeas%THREADS_VOXEL) max_iters++;
    
    if(idSubVOX==0){
      // #pragma unroll NPARAMS ... nvcc does not like this
      // warning: extra characters in the unroll pragma (expected a single positive integer), ignoring pragma for this loop
#pragma unroll
      for(int p=0;p<NPARAMS;p++){
	Gradient[p]=0;
      }
    }
    
    for(int iter=0;iter<max_iters;iter++){
      for(int p=0;p<NPARAMS;p++){
     	myderivatives[p]=0; // Maybe idMeasurement > nmeas
      }
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred_error=0;

      if(idMeasurement<nmeas){
	pred_error=Predicted_Signal(NPARAMS,parameters,fixed_params,myCFP);
	pred_error=pred_error-measurements[idMeasurement];
	
	Partial_Derivatives(NPARAMS,parameters,fixed_params,myCFP,myderivatives);
      }
      
      // #pragma unroll NPARAMS ... nvcc does not like this
#pragma unroll 
      for(int p=0;p<NPARAMS;p++){
	myderivatives[p]=2.0*pred_error*myderivatives[p];
	for (int offset=THREADS_VOXEL/2; offset>0; offset/=2){
	  myderivatives[p]+= __shfl_down(myderivatives[p],offset);
	}
      }
      if(idSubVOX==0){
	// #pragma unroll NPARAMS ... nvcc does not like this
#pragma unroll
	for(int p=0;p<NPARAMS;p++){
	  Gradient[p]+=myderivatives[p];
	}
      }
      idMeasurement+=THREADS_VOXEL;
    }  
  }

  template <typename T>
  __device__ void Calculate_Hessian(int idSubVOX,
				    int nmeas,
				    int CFP_Tsize,
				    T* measurements,
				    T* parameters,
				    T* fixed_params,
				    T* CFP,
				    T* Hessian)
  {
    int idMeasurement=idSubVOX;
    T myderivatives[NPARAMS];
    
    int max_dir = nmeas/THREADS_VOXEL;
    if(nmeas%THREADS_VOXEL) max_dir++;
    
    if(idSubVOX==0){
      for(int p=0;p<NPARAMS;p++){
	for(int p2=0;p2<NPARAMS;p2++){
	  Hessian[p*NPARAMS+p2]=0;
	}
      }
    }

    for(int iter=0;iter<max_dir;iter++){

      T* myCFP = &CFP[idMeasurement*CFP_Tsize];

      for(int p=0;p<NPARAMS;p++){
	myderivatives[p]=0;
      }
      T pred_error=0;

      if(idMeasurement<nmeas){
	pred_error=Predicted_Signal(NPARAMS,parameters,fixed_params,myCFP);
	pred_error=pred_error-measurements[idMeasurement];
	Partial_Derivatives(NPARAMS,parameters,fixed_params,myCFP, myderivatives);
      }
      
// #pragma unroll NPARAMS ... nvcc does not like this
#pragma unroll
      for(int p=0;p<NPARAMS;p++){
	for(int p2=0;p2<NPARAMS;p2++){
	  T element = 2.0 * myderivatives[p] * myderivatives[p2];
	  for (int offset=THREADS_VOXEL/2; offset>0; offset/=2){
	    element+= __shfl_down(element,offset);
	  }
	  if(idSubVOX==0){
	    Hessian[p*NPARAMS+p2]+=element;
	  }
	}
      }
      idMeasurement+=THREADS_VOXEL;
    }  
  }
  
  template <typename T>
  __device__ void LUsolver(int idSubVOX,
			   T* Hessian,
			   T* Gradient,
			   T* Solution){
    
    // If NPARAMS > 32 the current version of the method will fail !!
    // Need to generalise
    
    T col_elems[NPARAMS];
    T pivot;
    
    // Initialise Matrix. Each thread contains a column of the Hessian and one thread the Gradient column:   Matrix = [Hessian | Gradient]
    if (idSubVOX<NPARAMS){
      for(int p=0;p<NPARAMS;p++){
	col_elems[p] = Hessian[p*NPARAMS+idSubVOX];
      }
    }else if(idSubVOX==NPARAMS){
      for(int p=0;p<NPARAMS;p++){
	col_elems[p] = Gradient[p];
      }
    }else{
      for(int p=0;p<NPARAMS;p++){
	col_elems[p] = 0;
      }
    }
    
    // Solve in two steps: 
    // Forward step: Zero's under diagonal
    // Backward step: Zero's above diagonal

    // Forward step
    for (int col=0; col<NPARAMS; col++){
      // Divide row by diagonal element (1 in the diagonal)
      // Cannot have negative numbers in the diagonal -r/-r = +1
      pivot = col_elems[col];
      pivot = __shfl(pivot,col);
      col_elems[col] = col_elems[col]/pivot; 
      
      // Eliminate all terms under diagonal element (1)
      // Pivot is the element to make zero, Pivot-Pivot*1 = 0
      // This_row = This_row - Pivot * row_of_diagonal_element
      for (int row=col+1; row<NPARAMS; row++){
	pivot  = col_elems[row];
	pivot  = __shfl(pivot,col);
	col_elems[row] -= pivot*col_elems[col];
      }
    }

    // Backward step
    for (int col=NPARAMS-1; col>0; col--) {
      // Eliminate all terms above diagonal element
      for (int row=0; row<col; row++) {
	pivot  = col_elems[row];
	pivot  = __shfl(pivot,col);
	col_elems[row] -= pivot*col_elems[col];
      }
    }

    if(idSubVOX==NPARAMS){
      for(int p=0;p<NPARAMS;p++){
	Solution[p] = col_elems[p];
      }
    }
    
  }
  
  template <typename T, bool MARQUARDT>
  __global__ void levenberg_kernel(
				   int nvox, // nvoxels
				   int nmeas, // nmeasurements
				   int CFP_Tsize, //size*M-measurements
				   T* meas, // measurements
				   T* parameters, // model parameters 
				   T* CFP, // common fixed model parameters
				   int nmax_iters)
{
  // 1 block of threads process several voxels
  // Each warp processes 1 voxel
  int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
  int idVOX_inBlock =  threadIdx.x/THREADS_VOXEL;
  int idSubVOX= threadIdx.x%THREADS_VOXEL;
  bool leader = (idSubVOX==0);  // Some steps are performed by only one thread of the warp
    
  ////////// DYNAMIC SHARED MEMORY ///////////
  extern __shared__ double shared[];				//Size:
  double* pcf = (double*) shared;				//VOXELS_BLOCK 
  double* ncf = (double*) &pcf[VOXELS_BLOCK];			//VOXELS_BLOCK
  double* lambda = (double*) &ncf[VOXELS_BLOCK];		//VOXELS_BLOCK
  double* olambda = (double*) &lambda[VOXELS_BLOCK];		//VOXELS_BLOCK
  
  T* S_CFP = (T*) &olambda[VOXELS_BLOCK];			//nmeas*CMP_Tsize
  T* params = (T*) &S_CFP[nmeas*CFP_Tsize]; 			//NPARAMS*VOXELS_BLOCK
  T* Gradient = (T*) &params[NPARAMS*VOXELS_BLOCK];		//NPARAMS*VOXELS_BLOCK
  T* Hessian = (T*) &Gradient[NPARAMS*VOXELS_BLOCK];		//NPARAMS*NPARAMS*VOXELS_BLOCK
  T* step = (T*) &Hessian[NPARAMS*NPARAMS*VOXELS_BLOCK];	//NPARAMS*VOXELS_BLOCK
  
  int* success = (int*) &step[NPARAMS*VOXELS_BLOCK];		//VOXELS_BLOCK
  int* end = (int*) &success[VOXELS_BLOCK];			//VOXELS_BLOCK
  ////////////////////////////////////////////
  
  /// Copy common fixed model parameters to Shared Memory ///
  if(threadIdx.x==0){ // only one thread of the whole block. Common to all voxels
    for(int i=0;i<nmeas*CFP_Tsize;i++){
      S_CFP[i]=CFP[i];
    }
  }
  ///////////////////////////////////////////////////////////

  ///////// each voxel/warp of the block points to its data///////////
  meas = &meas[idVOX*nmeas]; //Global memory
  pcf = &pcf[idVOX_inBlock];
  ncf = &ncf[idVOX_inBlock];
  lambda = &lambda[idVOX_inBlock];
  olambda = &olambda[idVOX_inBlock];
  params = &params[idVOX_inBlock*NPARAMS];
  Gradient = &Gradient[idVOX_inBlock*NPARAMS];
  Hessian = &Hessian[idVOX_inBlock*NPARAMS*NPARAMS];
  step = &step[idVOX_inBlock*NPARAMS];
  success = &success[idVOX_inBlock];
  end = &end[idVOX_inBlock];

  int iter=0;
    
  /// Ititialise shared values of each voxel: only the leader///
  if(leader){ 
    *end=false;
    *success=true;
    *lambda=0.1;
    *olambda= 0.0;    
    *ncf=0.0;
    for(int i=0;i<NPARAMS;i++){
      params[i]=parameters[idVOX*NPARAMS+i];
    }
  }
  __syncthreads();
  ///////////////////////////////////////////

  T* fixed_params=new T[1]; // TODO: add feature
  Cost_Function(idSubVOX,nmeas,CFP_Tsize,meas,params,fixed_params,S_CFP,pcf);
  
  while (!( (*success) && iter++>=nmax_iters)){ 
    //if success we don't increase niter (first condition is true)
    //function cost has been decreased, we have advanced.
    if(*success){
      Calculate_Gradient(idSubVOX,nmeas,CFP_Tsize,meas,params,fixed_params,S_CFP,Gradient);
      Calculate_Hessian(idSubVOX,nmeas,CFP_Tsize,meas,params,fixed_params,S_CFP,Hessian);
    }
    
    if(leader){
      for (int i=0; i<NPARAMS; i++){
	if(MARQUARDT)
	  Hessian[(i*NPARAMS)+i]=((1+(*lambda))/(1+(*olambda)))*Hessian[i*NPARAMS+i];	//Levenberg-Marquardt
	else
	  Hessian[(i*NPARAMS)+i]+=(*lambda)-(*olambda);	//Levenberg
      }
    }
    
    LUsolver(idSubVOX,Hessian,Gradient,step);

    if(leader){
      for(int i=0;i<NPARAMS;i++){
	step[i]=params[i]-step[i];
      }
    }
    __syncthreads();
    
    Cost_Function(idSubVOX,nmeas,CFP_Tsize,meas,step,fixed_params,S_CFP,ncf);

    if(leader){
      if ( *success = ((*ncf) < (*pcf))){ 
	*olambda = 0.0;
	for(int i=0;i<NPARAMS;i++){
	  params[i]=step[i];
	}
	*lambda=(*lambda)/10.0;
	
	if (zero_cf_diff_conv(pcf,ncf)){
	  *end=true;
	}
	*pcf=*ncf;
      }else{
	*olambda=*lambda;
	*lambda=(*lambda)*10.0;
	if(*lambda > LTOL){ 
	 *end=true;
	}
      }
    }	
    __syncthreads();
    if(*end) break;		
  }
  if(leader){
    for(int i=0;i<NPARAMS;i++){
      parameters[idVOX*NPARAMS+i]=params[i];
    }
  }
  __syncthreads();
}
  
  
  template <typename T>
  Levenberg_Marquardt<T>::Levenberg_Marquardt(){
    cudimotOptions& opts = cudimotOptions::getInstance();
    max_iterations=opts.iterLevMar.value();
    Marquardt=opts.useMarquardt.value();
  }
  
  template <typename T>
  void Levenberg_Marquardt<T>::run(
				   int nvox, int nmeas,
				   int CFP_size,
				   T* meas,
				   T* params,
				   T* CFP) 
  {
    
    long int amount_shared_mem = 0;
    amount_shared_mem += 4*VOXELS_BLOCK*sizeof(double); // Levenberg parameters
    amount_shared_mem += (nmeas*CFP_size)*sizeof(T); // CFP
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Parameters
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Gradient
    amount_shared_mem += (NPARAMS*NPARAMS*VOXELS_BLOCK)*sizeof(T); // Hessian
    amount_shared_mem += 2*VOXELS_BLOCK*sizeof(int); // Levenberg parameters
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // step
    
    cout << "Shared Memory used in Levenberg-Marquardt kernel: " << amount_shared_mem << endl;
    
    int threads_block = VOXELS_BLOCK * THREADS_VOXEL;
    int nblocks=(nvox/VOXELS_BLOCK);
    if(nvox%VOXELS_BLOCK) nblocks++;
    
    if(Marquardt){
      levenberg_kernel<T,true><<<nblocks,threads_block,amount_shared_mem>>>(nvox,nmeas,CFP_size,meas,params,CFP,max_iterations);
    }else{
      levenberg_kernel<T,false><<<nblocks,threads_block,amount_shared_mem>>>(nvox,nmeas,CFP_size,meas,params,CFP,max_iterations);
    }
    sync_check("Levenberg_Marquardt Kernel");
  }
  
  // Explicit Instantiations of the template
  template class Levenberg_Marquardt<float>;
  template class Levenberg_Marquardt<double>;
}
