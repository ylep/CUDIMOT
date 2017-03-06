/*  split_parts.cc

    Moises Hernandez-Fernandez - FMRIB Image Analysis Group
    
    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include <sys/stat.h>
#include "newimage/newimageall.h"
#include "compileOptions/type.h"
#include "cudimotoptions.h"
#include "Model.h"

void save_part(Matrix data, string path, string name, int idpart){
  int nvox = data.Ncols();
  int nmeas = data.Nrows();
 
  string file_name;
  file_name = path+num2str(idpart)+"/"+name;
  
  ofstream out;
  out.open(file_name.data(), ios::out | ios::binary);
  out.write((char*)&nvox,4); // number of voxels
  out.write((char*)&nmeas,4); // number of measurements
  long size=nvox*nmeas*sizeof(Real);
  out.write((char*)&size,sizeof(long)); // number of bytes
  out.write((char*)&data(1,1),size); // information
  out.close(); 
}

using namespace Cudimot;

int main(int argc, char *argv[]){
  Log& logger = LogSingleton::getInstance();
  cudimotOptions& opts = cudimotOptions::getInstance();
  opts.parse_command_line(argc,argv,logger);
  
  NEWIMAGE::volume4D<MyType> data;
  NEWIMAGE::volume<MyType> mask;
  read_volume4D(data,opts.datafile.value());
  read_volume(mask,opts.maskfile.value());

  Matrix dataM;
  dataM=data.matrix(mask);
  
  int nmeas = dataM.Nrows();
  if(nmeas<=0){
    cerr << "CUDIMOT Error: The number of diffusion-weighted measurements must be greater than 0 in the input volume\n" << endl;
    exit (EXIT_FAILURE);
  }
  
  int nvoxels=dataM.Ncols();
  if(nvoxels<=0){
    cerr << "CUDIMOT Error: The number of voxels must be greater than 0" << endl;
    exit (EXIT_FAILURE);
  }

  // Create directories for the different parts
  for(int i=0;i<(opts.nParts.value());i++){
    string dirpath=opts.partsdir.value()+"/part_"+num2str(i);
    
    mkdir(dirpath.data(),0700);
  }
  
  int size_part=nvoxels/opts.nParts.value();
  
  Matrix data_part;
  string out_path;
  string out_name("data");
  out_path.append(opts.partsdir.value());
  out_path.append("/part_");
  for(int i=0;i<(opts.nParts.value()-1);i++){
    data_part = dataM.SubMatrix(1,nmeas,i*size_part+1,(i+1)*size_part);
    save_part(data_part,out_path,out_name,i);
  }

  // last part
  data_part = dataM.SubMatrix(1,nmeas,(opts.nParts.value()-1)*size_part+1,nvoxels);
  save_part(data_part,out_path,out_name,(opts.nParts.value()-1));

  //////////////////////////////////////////////////////
  /// Initialization of parameters
  /// The user can provide nifti files for some parameters
  /// Divide into different parts
  //////////////////////////////////////////////////////
  Model<MyType> model;
  int nparams=model.getNparams();
  
  if (opts.init_params.set()){
    string filename(opts.init_params.value());
    ifstream file(filename.data());
    if (file.is_open()){
      string line;
      int id_param=-1;
      while(getline(file,line)){
	id_param++;
	if (!line.empty()){
	  // Read volume with values fot this parameter
	  string name_file(line);
	  NEWIMAGE::volume4D<MyType> param_vals;
	  read_volume4D(param_vals,name_file);

	  Matrix paramM;
	  paramM=param_vals.matrix(mask);
	  
	  int size4dim = paramM.Nrows();
	  if(size4dim!=1){
	    cerr << "CUDIMOT Error: The volume used for initilize the parameters: " << name_file << " must be a 3D volume \n" << endl;
	    exit (EXIT_FAILURE);
	  }
	  
	  if(nvoxels!=paramM.Ncols()){
	    cerr << "CUDIMOT Error: The number of voxels in the data and the volume used for initilizing the parameters: " << name_file << " does not match\n" << endl;
	    exit (EXIT_FAILURE);
	  }
	  
	  Matrix param_part;
	  string out_name_p;
	  out_name_p.append("Param_");
	  out_name_p.append(num2str(id_param));
	  
	  for(int i=0;i<(opts.nParts.value()-1);i++){
	    param_part = paramM.SubMatrix(1,1,i*size_part+1,(i+1)*size_part);
	    save_part(param_part,out_path,out_name_p,i);
	  }
	  // last part
	  param_part = paramM.SubMatrix(1,1,(opts.nParts.value()-1)*size_part+1,nvoxels);
	  save_part(param_part,out_path,out_name_p,(opts.nParts.value()-1));
	  
	}else{
	  cout << "empty" << endl;
	  // Empty line, initialise this parameter with default value if provided or zeros otherwise
	  Matrix param_part;
	  string out_name_p;
	  out_name_p.append("Param_");
	  out_name_p.append(num2str(id_param));
	  
	  param_part.ReSize(1,size_part);
	  cout << "PARTS " <<  opts.nParts.value() << endl;
	  for(int i=0;i<(opts.nParts.value()-1);i++){
	    if(model.initProvided()){
	      param_part = model.getParam_init(id_param);
	    }else{
	      param_part = 0;
	    }
	    save_part(param_part,out_path,out_name_p,i);
	  }
	  // last part
	  int size_last_part=nvoxels-((opts.nParts.value()-1)*size_part);
	  param_part.ReSize(1,size_last_part); 
	  if(model.initProvided()){
	    param_part = model.getParam_init(id_param);
	  }else{
	    param_part = 0;
	  }
	  save_part(param_part,out_path,out_name_p,(opts.nParts.value()-1));
	}
	
      } //end lines
      
      if(id_param!=nparams){
	cerr << "CUDIMOT Error: The number of volumes (lines) provided for initializing the parameters in: " << filename.data() << " does not match the number of parameters of this model: " << nparams << ". If a parameter does not need initialization, its line can be empty.\n" << endl;
      }
      
    }else{
      cerr << "CUDIMOT Error: Unable to open Initialization Parameter file: " << filename.data() << endl; 
      exit(-1);
    }
  }else{
    // Not initialization file provided. Initialise with default parameter values if provided or zeros otherwise. But then not file division is needed.

  }
}
