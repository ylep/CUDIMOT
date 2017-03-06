/* Model.cu
   
   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
   
   Copyright (C) 2005 University of Oxford */

/* CCOPYRIGHT */

// This class contains the information of the Model

#include "Model.h"

namespace Cudimot{
  
  // This method reads the list of parameters from a file. File name must be given by the Designer.
  template <typename T>
  void Model<T>::Modelparser(){
    string input_file(MYMODEL_PARAMETERS);
    ifstream file_parameters(input_file.data());
    string mark_start="NEW_MODEL";
    string mark_comment="//";
    string mark_NPARAMS="NPARAMS=";
    string Pinit_mark="P_init["; // mark for defining initialization values
    string CFP_size_mark="CFP_size[";
    string mark_end="]";
    string delimiter=",";
    string line;
    
    nparams=0;
    provided_params_init=false;
    nCFP=0;
    CFP_Tsize=0;
    
    // File reading steps 0:Ignore, 1:Reading
    int reader_state=0;
    if (file_parameters.is_open()){
      while(getline(file_parameters,line)){
	
       	if(!line.compare(mark_start)) reader_state=1; // Starting MARK
	
	if(reader_state){ 
	  // Reading useful information
	  int pos0, pos1;
	  
	  if((pos0=line.find(mark_comment))>=0){
	    // If line commented out using "//", ignore the line
	    continue;
	  }
	  
	  //if((pos0=line.find(mark_NPARAMS))>=0){ 
	  // Read Number of Parameters
	  //pos0=pos0+mark_NPARAMS.length();
	  //line.erase(0,pos0);
	  //nparams=atoi(line.data());
	  //}
	  nparams=NPARAMS;
	  
	  if((pos0=line.find(Pinit_mark))>=0){ 
	    // Read Parameter Initializarion
	    provided_params_init=true;
	    pos0=pos0+Pinit_mark.length();
	    pos1=line.rfind(mark_end);
	    int len = pos1-pos0;
	    string vals=line.substr(pos0,len);
	    string token;
	    while ((pos0 = vals.find(delimiter))>=0){
	      params_init.push_back(strtod(vals.substr(0,pos0).data(),NULL));
	      vals.erase(0,pos0+delimiter.length());
	    }
	    params_init.push_back(strtod(vals.data(),NULL));
	  }
	  
    	  if((pos0=line.find(CFP_size_mark))>=0){ 
	    // Read common fixed param Sizes
	    pos0=pos0+CFP_size_mark.length();
	    pos1=line.rfind(mark_end);
	    int len = pos1-pos0;
	    string vals=line.substr(pos0,len);
	    string token;
	    while ((pos0 = vals.find(delimiter))>=0){
	      int size=atoi(vals.substr(0,pos0).data());
	      CFP_sizes.push_back(size);
	      CFP_Tsize+=size;
	      nCFP++;
	      vals.erase(0,pos0+delimiter.length());
	    }
	    int size=atoi(vals.data());
	    CFP_sizes.push_back(size);
	    CFP_Tsize+=size;
	    nCFP++;
	  }
	}
      }
      file_parameters.close();
      
      // All the infomation needed ?
      // CHECK
      if(!nparams){
	string input_file(MYMODEL_PARAMETERS);
	cerr << "CUDIMOT Error: Number of parameters must be greater than 0. Please check your parameter specification: " << input_file.data() << endl;
	exit(-1);
      }

    }else{
      cerr << "CUDIMOT Error: Unable to open Parameters file configuration: " << input_file.data() << endl; 
      exit(-1);
    }
  }

  template <typename T>
  Model<T>::Model(){
    cudimotOptions& opts = cudimotOptions::getInstance();
    /// Read text file with parameters information
    Modelparser();
  }
  
  template <typename T>
  Model<T>::~Model(){}

  template <typename T>
  int Model<T>:: getNparams(){
    return nparams;
  }
  
  template <typename T>
  bool Model<T>:: initProvided(){
    return provided_params_init;
  }
  
  template <typename T>
  T Model<T>:: getParam_init(int id_param){
    return params_init[id_param];
  }
  
  // Explicit Instantiations of the template
  template class Model<float>;
  template class Model<double>;
}
