/* MCMC.cu

   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
   
   Copyright (C) 2005 University of Oxford */

/* CCOPYRIGHT */

// Markoc Chain Monte Carlo method

#include <curand.h>
#include <curand_kernel.h>
#include "MCMC.h"
#include "utils.h"
#include "cudimotoptions.h"
#include "mymodels/mymodel.h"
#include MYMODEL_FUNCTIONS

namespace Cudimot{

#define maxfloat 1e10f
  
  __constant__ float LowerLimits [NPARAMS];
  __constant__ float UpperLimits [NPARAMS];
  
  // SHFL double precision
  /* __device__ inline double __shfl_down(double var, unsigned int srcLane, int width=32) {
    int2 a = *reinterpret_cast<int2*>(&var);
    a.x = __shfl_down(a.x, srcLane, width);
    a.y = __shfl_down(a.y, srcLane, width);
    return *reinterpret_cast<double*>(&a);
    }*/
  
  __device__ inline float log_gpu(float x){return logf(x);}
  __device__ inline double log_gpu(double x){return log(x);}
  __device__ inline float exp_gpu(float x){return expf(x);}
  __device__ inline double exp_gpu(double x){return exp(x);}
  
  // Returns the natural log of the 0th order modified Bessel function of first kind for an argument x
  // Follows the exponential implementation of the Bessel function in Numerical Recipes, Ch. 6
  __device__ inline float logIo(const float x){
    float y,b;
    b=fabsf(x);
    if (b<3.75f){
      float a=x/3.75f;
      a*=a;
      //Bessel function evaluation
      y=1.0f+a*(3.5156229f+a*(3.0899424f+a*(1.2067492f+a*(0.2659732f+a*(0.0360768f+a*0.0045813f)))));
      y=logf(y);
    }else{
      float a=3.75f/b; 
      //Logarithm of Bessel function
      y=b+logf((0.39894228f+a*(0.01328592f+a*(0.00225319f+a*(-0.00157565f+a*(0.00916281f+a*(-0.02057706f+a*(0.02635537f+a*(-0.01647633f+a*0.00392377f))))))))/sqrt(b));
    }
    return y;
  }
  
  __device__ inline double logIo(const double x){
    double y,b;
    b=fabs(x);
    if (b<3.75){
      double a=x/3.75;
      a*=a;
      //Bessel function evaluation
      y=1.0+a*(3.5156229+a*(3.0899424+a*(1.2067492+a*(0.2659732+a*(0.0360768+a*0.0045813)))));
      y=log(y);
    }
    else{
      double a=3.75/b; 
      //Logarithm of Bessel function
      y=b+log((0.39894228+a*(0.01328592+a*(0.00225319+a*(-0.00157565+a*(0.00916281+a*(-0.02057706+a*(0.02635537+a*(-0.01647633+a*0.00392377))))))))/sqrt(b));
    }
    return y;
  }


  template <typename T>
  __device__ inline void Propose(int par, T* params, T* old, T* propSD, curandState* localrandState){
    *old=params[par];
    params[par] = params[par] + curand_normal(localrandState)*propSD[par];
  }

  template <typename T>
  __device__ inline int Check_limits(int idpar, T* params){
    //if(params[idpar]<LowerLimits[idpar]) return 0;
    //if(params[idpar]>UpperLimits[idpar]) return 0;
    return 1;
  }
  
  
  template <typename T>
  __device__ inline void Check_prior(){
     // check lower limit

     // check upper limit

     // Gaussian (mean,sd)
     // if (param<=0.0f || m_f[fibre]>=1.0f)
    
     // Gamma (alpha,beta)
     
     //(1.0-alpha)* log(param) + beta*param;

     // ARD - Automatic Relevance Determination

     // Beta (alpha,beta)
    
    // sin()

  }
  
  
  template <typename T, bool RICIAN_NOISE>
  __device__ void  Compute_Likelihood(int idSubVOX,
				      int nmeas,
				      int CFP_Tsize,
				      T* measurements,
				      T* parameters,
				      T* tau,
				      T* fixed_params,
				      T* CFP,
				      T* likelihood)
  {
    int idMeasurement=idSubVOX;
    int nmeas2compute = nmeas/THREADS_VOXEL;
    if (idSubVOX<(nmeas%THREADS_VOXEL)) nmeas2compute++;
    
    T accumulated_error=0;
    for(int dir=0;dir<nmeas2compute;dir++){
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred_error=Predicted_Signal(NPARAMS,parameters,fixed_params,myCFP);
      if(RICIAN_NOISE){
	T meas = measurements[idMeasurement];
	pred_error=log_gpu(meas)+(-0.5*(*tau)*(meas*meas+pred_error*pred_error)+logIo((*tau)*pred_error*meas));
	accumulated_error+=pred_error;
      }else{
	pred_error=pred_error-measurements[idMeasurement];
	accumulated_error+=pred_error*pred_error;
      }
      idMeasurement+=THREADS_VOXEL;
    }
   
    for (int offset=THREADS_VOXEL/2; offset>0; offset/=2){
      accumulated_error+= __shfl_down(accumulated_error,offset);
    }
    
    if(idSubVOX==0){
      if(RICIAN_NOISE){
	*likelihood = -nmeas*log_gpu(*tau)-accumulated_error;
      }else{
	*likelihood = (nmeas/2.0f)*log_gpu(accumulated_error/2.0f);
      }
    }
  }
  
  template <typename T>
  __device__ inline int Compute_test_energy(T* new_energy, T* old_energy, T* prior, T* likelihood, curandState* localrandState){
    (*old_energy) = (*new_energy);
    (*new_energy) = (*prior)+ (*likelihood);
    
    T tmp=exp_gpu((*old_energy)-(*new_energy));
    
    return (tmp>curand_uniform(localrandState));
  }
  

  template <typename T, bool RICIAN_NOISE, bool RECORDING>
  __global__ void mcmc_kernel(
			      curandState* randstate, // to generate random numbers
			      int nvox, // num voxels
			      int nmeas, // num measurements
			      int CFP_Tsize, //size*M-measurements
			      int niters,
			      int nsamples, // num samples per parameter
			      int sampleevery, // record a sample every x iterations
			      int updateproposalevery, // update SD proposals every x iters
			      T* meas, // measurements
			      T* parameters, // model parameters 
			      T* CFP, // common fixed model parameters
			      T* samples) // to record parameters samples
  {
    // 1 block of threads process several voxels
    // Each warp processes 1 voxel
    int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
    int idVOX_inBlock =  threadIdx.x/THREADS_VOXEL;
    int idSubVOX= threadIdx.x%THREADS_VOXEL;
    bool leader = (idSubVOX==0);  // Some steps are performed by only one thread of the warp
    
    ////////// DYNAMIC SHARED MEMORY ///////////
    extern __shared__ double shared[];		     			// Size: 
    curandState* localrandState = (curandState*)shared;		

    T* S_CFP = (T*) &localrandState[VOXELS_BLOCK]; 	       	        // nmeas*CFP_Tsize
    T* params = (T*) &S_CFP[nmeas*CFP_Tsize]; 				// NPARAMS*VOXELS_BLOCK 
    T* priors = (T*) &params[NPARAMS*VOXELS_BLOCK]; 			// NPARAMS*VOXELS_BLOCK
    T* propSD = (T*) &priors[NPARAMS*VOXELS_BLOCK]; 			// NPARAMS*VOXELS_BLOCK
    
    T* likelihood = (T*) &propSD[NPARAMS*VOXELS_BLOCK]; 		// VOXELS_BLOCK 
    T* Tprior = (T*) &likelihood[VOXELS_BLOCK]; 			// VOXELS_BLOCK 
    T* energy = (T*) &Tprior[VOXELS_BLOCK];				// VOXELS_BLOCK 
    T* tau = (T*) &energy[VOXELS_BLOCK];			       	// VOXELS_BLOCK 

    T* old_param = (T*) &tau[VOXELS_BLOCK];		 	      	// VOXELS_BLOCK 
    T* old_energy =  (T*) &old_param[VOXELS_BLOCK];		       	// VOXELS_BLOCK 

    int* naccepted = (int*) &old_energy[VOXELS_BLOCK];			//NPARAMS*VOXELS_BLOCK
    int* nrejected = (int*) &naccepted[NPARAMS*VOXELS_BLOCK];		//NPARAMS*VOXELS_BLOCK
    ////////////////////////////////////////////
    
    /// Copy common fixed model parameters to Shared Memory ///
    if(threadIdx.x==0){ // only one thread of the whole block. Common to all voxels
      for(int i=0;i<nmeas*CFP_Tsize;i++){
	S_CFP[i]=CFP[i];
      }
    }
    ///////////////////////////////////////////////////////////
    
    ///////// each voxel/warp of the block points to its data///////////
    meas = &meas[idVOX*nmeas]; // Global memory
    if(RECORDING){
      samples = &samples[idVOX*NPARAMS*nsamples]; //Global memory
    }
    localrandState = (curandState*)&localrandState[idVOX_inBlock];
    params = &params[idVOX_inBlock*NPARAMS];
    priors = &priors[idVOX_inBlock*NPARAMS];
    propSD = &propSD[idVOX_inBlock*NPARAMS];
    likelihood = &likelihood[idVOX_inBlock];
    Tprior = &Tprior[idVOX_inBlock];
    energy = &energy[idVOX_inBlock];
    tau = &tau[idVOX_inBlock];
    old_param = &old_param[idVOX_inBlock];
    old_energy = &old_energy[idVOX_inBlock];
    naccepted = &naccepted[idVOX_inBlock*NPARAMS];
    nrejected = &nrejected[idVOX_inBlock*NPARAMS];
    
    /// Ititialise shared values of each voxel: only the leader///
    if(leader){ 
      *localrandState = randstate[idVOX];
      for(int par=0;par<NPARAMS;par++){
      	params[par]=parameters[idVOX*NPARAMS+par];
	naccepted[par]=0;
	nrejected[par]=0;
	priors[par]=0;
	propSD[par]=params[par]/10.0; // add by the user
      }
      *tau=0;
      *Tprior=0;
    }
    __syncthreads();
    ///////////////////////////////////////////
    
    T* fixed_params=new T[1]; // TODO: add feature
    
    Compute_Likelihood<T,RICIAN_NOISE>(idSubVOX,nmeas,CFP_Tsize,meas,params,tau,fixed_params,S_CFP,likelihood);
    if(leader){
      *energy=(*Tprior)+(*likelihood);
    }
    
    for (int iter=0; iter<niters; iter++){
      for (int par=0; par<NPARAMS; par++){
	int criteria=0;
	if(leader){
	  Propose(par,params,old_param,propSD,localrandState);
	  criteria=Check_limits(par,params);
	}
	criteria = __shfl(criteria,0);
	
	if(criteria){
	  Compute_Likelihood<T,RICIAN_NOISE>(idSubVOX,nmeas,CFP_Tsize,meas,params,tau,fixed_params,S_CFP,likelihood);
	  if(leader){
	    criteria=Compute_test_energy(energy,old_energy,Tprior,likelihood,localrandState);
	    if(criteria){
	      naccepted[par]++;
	    }else{
	      nrejected[par]++;
	      params[par]=(*old_param);
	      *energy=*old_energy;
	    }
	  }
	  criteria = __shfl(criteria,0); // better than __syncthreads(); ?
	}else{
	  nrejected[par]++;
	  params[par]=(*old_param);
	}
      }
      
      if(RECORDING){
	if((!(iter%sampleevery))&&(leader)){
	  int nsamp=iter/sampleevery;
	  for (int par=0; par<NPARAMS; par++){
	    samples[par*nsamples+nsamp]=params[par];
	  }
	}
      }

      if(!RECORDING){  // deactivated when not recording ?
	if((iter>0)&&(!(iter%updateproposalevery))&&(leader)){
	  for (int par=0; par<NPARAMS; par++){
	    propSD[par]=sqrt((naccepted[par]+1.0)/(nrejected[par]+1.0)); // not too big ?
	    propSD[par]=min(propSD[par],maxfloat);
	  }
	}
      }
    }
    if(leader){
      randstate[idVOX]=*localrandState; // save state, otherwise random numbers will be repeated (start at the same point)
    }
  }
  
  __global__ void setup_randoms_kernel(curandState* randstate, double seed){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    curand_init(seed,id,0,&randstate[id]);
  }
  
  template <typename T>
  MCMC<T>::MCMC(int nvoxFitpart){
    cudimotOptions& opts = cudimotOptions::getInstance();
    nvoxFit_part=nvoxFitpart;
    nburnin=opts.nburn.value();
    njumps=opts.njumps.value();
    nsamples=opts.njumps.value()/opts.sampleevery.value();
    sampleevery=opts.sampleevery.value();
    updateproposalevery=opts.sampleevery.value();
    RicianNoise=opts.rician.value();
    
    // Initialise Randoms
    int blocks_Rand = nvoxFit_part/256;
    if(nvoxFit_part%256) blocks_Rand++;
    cudaMalloc((void**)&randStates, blocks_Rand*256*sizeof(curandState));
    dim3 Dim_Grid_Rand(blocks_Rand,1);
    dim3 Dim_Block_Rand(256,1);
    srand(opts.seed.value());  //randoms seed
    setup_randoms_kernel<<<Dim_Grid_Rand,Dim_Block_Rand>>>(randStates,rand());
    sync_check("Setup_Randoms_kernel");
  }
  
  template <typename T>
  void MCMC<T>::run(
		    int nvox, int nmeas,
		    int CFP_size,
		    T* meas,
		    T* params,
		    T* CFP,
		    T* samples) 
  {
    long int amount_shared_mem = 0;
    amount_shared_mem += VOXELS_BLOCK*sizeof(curandState); // curandState
    amount_shared_mem += (nmeas*CFP_size)*sizeof(T); // CFP
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Parameters
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Priors
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // PropSD
    amount_shared_mem += (4*VOXELS_BLOCK)*sizeof(T); // Likelihod,TPrior,Energy, Tau
    amount_shared_mem += (2*VOXELS_BLOCK)*sizeof(T); // old_param, old_energy
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(int); // naccepted
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(int); // nrejected

    cout << "Shared Memory used in MCMC kernel: " << amount_shared_mem << endl;
    
    int threads_block = VOXELS_BLOCK * THREADS_VOXEL;
    int nblocks=(nvox/VOXELS_BLOCK);
    if(nvox%VOXELS_BLOCK) nblocks++;
    
    // Burn-In
    if(RicianNoise){
      mcmc_kernel<T,true,false><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nvox,nmeas,CFP_size,nburnin,nsamples,sampleevery,updateproposalevery,meas,params,CFP,samples);
    }else{
      mcmc_kernel<T,false,false><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nvox,nmeas,CFP_size,nburnin,nsamples,sampleevery,updateproposalevery,meas,params,CFP,samples);
    }
    sync_check("MCMC Kernel: burnin step");
    
    // Recordig
    if(RicianNoise){
      mcmc_kernel<T,true,true><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nvox,nmeas,CFP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,CFP,samples);
    }else{
      mcmc_kernel<T,false,true><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nvox,nmeas,CFP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,CFP,samples);
    }
    sync_check("MCMC Kernel: recording step"); 
    
    printf("--------------------- MCMC completed ------------\n");
  }

  // Explicit Instantiations of the template
  template class MCMC<float>;
  template class MCMC<double>;
}
