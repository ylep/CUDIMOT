/* Parameters.cu
   
   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
   
   Copyright (C) 2005 University of Oxford */

/* CCOPYRIGHT */

#include "Parameters.h"

namespace Cudimot{
  
  template <typename T>
  Parameters<T>::Parameters(Model<T> model,dMRI_Data<T> dMRI_data):
    nparams(model.nparams),
    nFP(model.nFP), nCFP(model.nCFP),CFP_Tsize(model.CFP_Tsize),
    nvox(dMRI_data.nvox), nmeas(dMRI_data.nmeas),nparts(dMRI_data.nparts),
    size_part(dMRI_data.size_part),size_last_part(dMRI_data.size_last_part),
    nvoxFit_part(dMRI_data.nvoxFit_part)
  {
    
    Log& logger = LogSingleton::getInstance();
    cudimotOptions& opts = cudimotOptions::getInstance();
    
    //////////////////////////////////////////////////////
    /// Initialise parameters
    /// The user can provide nifti files for some parameters
    //////////////////////////////////////////////////////
    params_host=new T[nvox*nparams];
    
    if (opts.init_params.set()){
      // If volumes provided for initialize parameters
      for(int idParam=0;idParam<nparams;idParam++){
	// Read binary file with values for this part, logfile
	string name_file;
	name_file.append(opts.partsdir.value());
	name_file.append("/part_");
	name_file.append(num2str(opts.idPart.value()));
	name_file.append("/Param_");
	name_file.append(num2str(idParam));
	
	ifstream in;
	long nbytes;
	int nvox_file,nmeas_file;
	in.open(name_file.data(), ios::in | ios::binary);
	in.read((char*)&nvox_file, 4);
	in.read((char*)&nmeas_file, 4);
	in.read((char*)&nbytes, sizeof(long));
	
	if(nvox!=nvox_file || nmeas_file!=1){
	  cerr << "CUDIMOT Error: The amount of data in the input file " <<  name_file << " for initializing the parameters is not correct" << endl;
	  exit(-1);
	}
	
	Matrix Parameters_init;
	Parameters_init.ReSize(1,nvox);
	in.read((char*)&Parameters_init(1,1),nbytes);
	in.close();
	
	for(int i=0;i<nvox;i++){
	  params_host[i*nparams+idParam]=Parameters_init(1,i+1);
	}
      }
      
    }else{
      // If NOT volumes provided for initialize parameters, try default values
      for(int idParam=0;idParam<nparams;idParam++){
	for(int i=0;i<nvox;i++){
	  if(model.initProvided()){
	    params_host[i*nparams+idParam]=model.getParam_init(idParam);
	  }else{
	    params_host[i*nparams+idParam]=0;
	  }
	}
      }
    }
    
    //////////////////////////////////////////////////////
    /// Read Common Fixed Parameters (kxM, M:measurements)
    /// Provided by the user
    /// If not provided use the defaut values (given by the model-designer)
    //////////////////////////////////////////////////////
    CFP_host= new T[nmeas*CFP_Tsize];
    
    if (opts.CFP.set()){
      int nCFP_set = 0; //to check that all the values are set
      int cumulative_size=0;
      string filename(opts.CFP.value());
      ifstream file(filename.data());
      if (file.is_open()){
	string line;
	while(getline(file,line)){
       	  if (!line.empty()){ // if is empty, what should I do?
	    // Read file with values
	    Matrix values;
	    values=read_ascii_matrix(line);
	    int param_size = model.CFP_sizes[nCFP_set];
	    if(values.Nrows()!=param_size){
	      cerr << "CUDIMOT Error: Common Fixed Parameter number " << nCFP_set << " and file " << line << " do not match the dimensions specified for this model" << endl;
	      exit(-1);
	    }
	    if(values.Ncols()!=nmeas){
	      cerr << "CUDIMOT Error:  Common Fixed Parameter number " << nCFP_set << " and file " << line << " do not match the number of measurements: " << nmeas << endl;
	      exit(-1);
	    }
	    
	    for(int i=0;i<param_size;i++){
	      for (int j=0;j<nmeas;j++){
		CFP_host[cumulative_size+i+j*CFP_Tsize]=values(i+1,j+1);
	      }
	    }
	    nCFP_set++;
	    cumulative_size+=param_size;
	    
	  }
	} //end lines
	
	if(nCFP_set!=nCFP){
	  cerr << "CUDIMOT Error: The number of common fixed parameter provided in file: " << filename.data() << " is not correct for this model. The number of common fixed parameter must be " << nCFP << endl;
	  exit(-1);
	}
      }else{
	cerr << "CUDIMOT Error: Unable to open Common Fixed Parameters file: " << filename.data() << endl; 
	exit(-1);
      }
    
    }else{
      // No CFP provided ?
      
    }
    //////////////////////////////////////////////////////
    
        
    /// Allocate GPU memory
    cudaMalloc((void**)&params_gpu,nvoxFit_part*nparams*sizeof(T));
    cudaMalloc((void**)&CFP_gpu,CFP_Tsize*nmeas*sizeof(T));
    sync_check("Allocating Model Parameters on GPU\n");
    
    // Copy Parameters from host to GPU
    cudaMemcpy(CFP_gpu,CFP_host,CFP_Tsize*nmeas*sizeof(T),cudaMemcpyHostToDevice);
    sync_check("Copying Common-Fixed Model Parameters to GPU\n");
    
    /// If MCMC: allocate memory in host and GPU for samples;
    if(opts.runMCMC.value()){
      nsamples=(opts.njumps.value()/opts.sampleevery.value());   
    }else{
      nsamples=1;
    }
    samples_host = new T[nsamples*nparams*nvox];
    cudaMalloc((void**)&samples_gpu,nvoxFit_part*nparams*nsamples*sizeof(T));
  }
  
  template <typename T>
  Parameters<T>::~Parameters(){}
  
  template <typename T>
  T* Parameters<T>::getParametersPart(int part){
    if(part>=nparts){
      cerr << "CUDIMOT Error: Trying to get an incorrect part of the Parameters: " << part << ". There are only " << nparts << " parts and index starts at 0." << endl;
      exit(-1);
    }
    
    int initial_pos=part*size_part*nparams;
    // Copy from host to GPU
    cudaMemcpy(params_gpu,&params_host[initial_pos],nvoxFit_part*nparams*sizeof(T),cudaMemcpyHostToDevice);
    sync_check("Copying Model Parameters to GPU\n");

    return params_gpu;
  }
 
  template <typename T>
  T* Parameters<T>::getSamples(){
    return samples_gpu;
  }

  template <typename T>
  int Parameters<T>::getTsize_CFP(){
    return CFP_Tsize;
  }
  
  template <typename T>
  T* Parameters<T>::getCFP(){
    return CFP_gpu;
  }
  
  template <typename T>
  void Parameters<T>::copyParamsPartGPU2Host(int part){

    if(part>=nparts){
      cerr << "CUDIMOT Error: Trying to store an incorrect part of Parameters: " << part << ". There are only " << nparts << " parts and index starts at 0." << endl;
      exit(-1);
    }
    
    int size=size_part; // this ignores the extra voxels added
    if(part==(nparts-1)){
      size=size_last_part; // this ignores the extra voxels added
    }
    int initial_pos=part*size_part*nparams;
    
    cudaMemcpy(&params_host[initial_pos],params_gpu,size*nparams*sizeof(T),cudaMemcpyDeviceToHost);
     sync_check("Copying Model Parameters from GPU\n");
     
     samples_host=params_host; 
  }
  
  template <typename T>
  void Parameters<T>::copyParams2Samples(){
    samples_host=params_host;
  }
  
  
  template <typename T>
  void Parameters<T>::copySamplesPartGPU2Host(int part){
    if(part>=nparts){
      cerr << "CUDIMOT Error: Trying to store an incorrect part of Samples: " << part << ". There are only " << nparts << " parts and index starts at 0." << endl;
      exit(-1);
    }
    
    int size=size_part; // this ignores the extra voxels added
    if(part==(nparts-1)){
      size=size_last_part; // this ignores the extra voxels added
    }
    int initial_pos=part*size_part*nparams*nsamples;
    
    cudaMemcpy(&samples_host[initial_pos],samples_gpu,size*nparams*nsamples*sizeof(T),cudaMemcpyDeviceToHost);
    sync_check("Copying Samples from GPU\n");
  }
  
  template <typename T>
  void Parameters<T>::writeSamples(){
    Log& logger = LogSingleton::getInstance();
    cudimotOptions& opts = cudimotOptions::getInstance();

    // Create a Matrix [nsamples X nvoxels] for each parameter
    vector<Matrix> samples;
    samples.resize(nparams);
    
    for(int par=0;par<nparams;par++){
      samples[par].ReSize(nsamples,nvox);
      samples[par]=0;
    }
    // Copy samples to each Matrix
    for(int vox=0;vox<nvox;vox++){  
      for(int par=0;par<nparams;par++){
	for(int sam=0;sam<nsamples;sam++){
	  samples[par](sam+1,vox+1)=samples_host[vox*nparams*nsamples+par*nsamples+sam];
	}
      }
    }
    
    // Write to file
    for(int par=0;par<nparams;par++){
      string file_name;
      file_name.append(opts.partsdir.value());
      file_name.append("/part_");
      file_name.append(num2str(opts.idPart.value()));
      file_name.append("/Param_"+num2str(par)+"_samples");
      ofstream out;
      out.open(file_name.data(), ios::out | ios::binary);
      out.write((char*)&nvox,4); // number of voxels
      out.write((char*)&nsamples,4); // number of measurements
      long size=nvox*nsamples*sizeof(Real); //need Real here (NEWMAT Object!)
      out.write((char*)&size,sizeof(long)); // number of bytes
      out.write((char*)&samples[par](1,1),size);
      out.close();
    }
  }
  
  // Explicit Instantiations of the template
  template class Parameters<float>;
  template class Parameters<double>;
}
