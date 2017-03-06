/*  cudimot.cc

    Moises Hernandez-Fernandez - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */


#include <sys/time.h>
#include "compileOptions/type.h"
#include "cudimotoptions.h"
#include "init_gpu.h"
#include "dMRI_Data.h"
#include "Model.h"
#include "Parameters.h"
#include "Levenberg_Marquardt.h"
#include "MCMC.h"

using namespace Cudimot;

double timeval_diff(struct timeval *a, struct timeval *b){
  return (double)(a->tv_sec +(double)a->tv_usec/1000000) - (double)(b->tv_sec +(double)b->tv_usec/1000000);
}

int main(int argc, char *argv[]){
  
  struct timeval t1,t2;
  double time;
  gettimeofday(&t1,NULL); 
  
  // Setup logging:
  Log& logger = LogSingleton::getInstance();
  cudimotOptions& opts = cudimotOptions::getInstance();
  opts.parse_command_line(argc,argv,logger);
  srand(opts.seed.value());  //randoms seed
  
  init_gpu();
  
  // Encapsulate dMRI data
  dMRI_Data<MyType> data;
  
  Model<MyType> model;
  Parameters<MyType> params(model,data);
  Levenberg_Marquardt<MyType> methodLM; // add options
  MCMC<MyType> methodMCMC(data.getNvoxFit_part());
  // Data and parameters are divided into parts => process each part
  for(int part=0;part<data.getNparts();part++){
    int part_size=0;
    MyType* meas=data.getMeasPart(part,part_size);
    
    MyType* parameters_part = params.getParametersPart(part);
    
    methodLM.run(part_size,data.getNmeas(),
		 params.getTsize_CFP(),
		 meas,parameters_part,
		 params.getCFP());
    
    if(!opts.runMCMC.value()){
      params.copyParamsPartGPU2Host(part);
    }else{
      methodMCMC.run(part_size,data.getNmeas(),
		     params.getTsize_CFP(),
		     meas,parameters_part,
		     params.getCFP(),
		     params.getSamples());
      
      params.copySamplesPartGPU2Host(part); 
    }
  }
  if(!opts.runMCMC.value()){
    params.copyParams2Samples();
  }
  params.writeSamples();
  
  gettimeofday(&t2,NULL);
  time=timeval_diff(&t2,&t1); 
  cout << endl << "Part processed in: " << time << " seconds" << endl;	
  
}

