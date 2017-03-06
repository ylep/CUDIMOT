/*  merge_parts.cc

    Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include <sys/stat.h>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "compileOptions/type.h"
#include "cudimotoptions.h"
#include "Model.h"


using namespace Cudimot;


void join_Parts(NEWIMAGE::volume<MyType> mask, string directory_in, string name_in, string name_out, int nsamples, int nParts, float max, float min){
    
  Matrix result(nsamples,0);
  Matrix part;

  for(int i=0;i<nParts;i++){
    
    std::string file_name;
    file_name.assign(directory_in);
    file_name += num2str(i);
    file_name += "/"; 
    file_name += name_in; 
    
    ifstream in;
    long nbytes_file;
    int nvox_file,nsamples_file;
    in.open(file_name.data(), ios::in | ios::binary);
    in.read((char*)&nvox_file, 4);
    in.read((char*)&nsamples_file, 4);
    in.read((char*)&nbytes_file, sizeof(long));
    if(nvox_file<=0 || nsamples_file<=0 || nsamples_file!=nsamples ){
      cerr << "CUDIMOT Error: The amount of data in the intermediate output file: " << file_name.data() << " is not correct." << endl;
      exit(-1);
    }
    part.ReSize(nsamples_file,nvox_file);
    in.read((char*)&part(1,1), nbytes_file);
    in.close();
    result |= part;
  }
  
  NEWIMAGE::volume4D<MyType> tmp;
  tmp.setmatrix(result,mask);
  if(max==-10) max=tmp.max();
  if(min==-10) min=tmp.min(); 
  tmp.setDisplayMaximumMinimum(max,min);
  save_volume4D(tmp,name_out);
}

//////////////////////////////////////////////////////////
//       MERGE THE OUTPUTS FILES OF CUDIMOT
//////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  // Setup logging:
  Log& logger = LogSingleton::getInstance();
  cudimotOptions& opts = cudimotOptions::getInstance();
  opts.parse_command_line(argc,argv,logger);
  
  NEWIMAGE::volume<MyType> mask;
  read_volume(mask,opts.maskfile.value());
  
  int nsamples=0;
  if(!opts.runMCMC.value()){
    nsamples=1; // LM
  }else{
    nsamples=opts.njumps.value()/opts.sampleevery.value();
  }
  if(nsamples<=0){
    cerr << "CUDIMOT Error: The number of samples must be greater than 0" << endl;
    exit (EXIT_FAILURE);
  }

  Model<MyType> model;
  int nparams = model.getNparams();

  string path_in;
  path_in.append(opts.partsdir.value());
  path_in.append("/part_");

  string path_out;
  path_out.append(opts.outputdir.value());

  for(int par=0;par<nparams;par++){
    string file_name = "Param_" + num2str(par) + "_samples";
    std::string output_file=path_out+"/"+file_name;
            
    join_Parts(mask,path_in,file_name,output_file,nsamples,opts.nParts.value(),-10,-10);
  }
  return 0;
}

