/*  cudimotOptions.h
    
    Moises Hernandez-Fernandez FMRIB Image Analysis Group
    
    Copyright (C) 1999-2010 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(cudimotOptions_h)
#define cudimotOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"

using namespace Utilities;

namespace Cudimot {
  
  class cudimotOptions {
  public:
    static cudimotOptions& getInstance();
    ~cudimotOptions() { delete gopt; }
    
    Option<bool> verbose;
    Option<bool> help;
    Option<string> logdir;
    Option<bool> forcedir;
    Option<string> datafile;
    Option<string> maskfile;
    Option<string> subjdir;
    Option<string> partsdir;
    Option<string> outputdir;
    Option<int> idPart;
    Option<int> nParts;
    Option<float> fudge;
    Option<int> njumps;
    Option<int> nburn;
    Option<int> nburn_noard;
    Option<int> sampleevery;
    Option<int> updateproposalevery;
    Option<int> seed;
    Option<bool> runLevMar;
    Option<bool> useMarquardt;
    Option<int> iterLevMar;
    Option<bool> runMCMC;
    Option<bool> no_ard;
    Option<bool> all_ard;
    Option<bool> rician;
    Option<string> CFP;
    Option<string> init_params;
    
    void parse_command_line(int argc, char** argv,  Log& logger);
  
  private:
    cudimotOptions();  
    const cudimotOptions& operator=(cudimotOptions&);
    cudimotOptions(cudimotOptions&);
    
    OptionParser options; 
    
    static cudimotOptions* gopt;
  };
  
  inline cudimotOptions& cudimotOptions::getInstance(){
    if(gopt == NULL)
      gopt = new cudimotOptions();
    
    return *gopt;
  }
  
  inline cudimotOptions::cudimotOptions():
	verbose(string("-V,--verbose"), false,
		string("switch on diagnostic messages"),
		false, no_argument),
	help(string("-h,--help"), false,
		string("display this message"),
		false, no_argument),
	logdir(string("--ld,--logdir"), string("logdir"),
		string("log directory (default is logdir)"),
		false, requires_argument),
	forcedir(string("--forcedir"),false,
		string("Use the actual directory name given - i.e. don't add + to make a new directory"),
		false,no_argument),
	datafile(string("--data"), string("data"),
		string("data file"),
		true, requires_argument),  
	maskfile(string("--maskfile"), string("nodif_brain_mask"),
		string("mask file"),
		true, requires_argument),
	subjdir(string("--subjdir"), string(""),
		string("Subject directory"),
		true, requires_argument),
	partsdir(string("--partsdir"), string("partsdir"),
		string("Directory where different parts of the data/temporal-results will be stored"),
		true, requires_argument),
	outputdir(string("--outputdir"), string("outputdir"),
		string("Directory where to write the output files"),
		true, requires_argument),
	idPart(string("--idPart"),0,
		string("Number of the part of the dataset to process [0..N-1]"),
		true, requires_argument),
	nParts(string("--nParts"),0,
		string("Total number of parts of the dataset [1..N]"),
		true, requires_argument),
	fudge(string("--fudge"),1,
		string("ARD fudge factor"),
		false,requires_argument),
	njumps(string("--nj,--njumps"),1250,
		string("Num of jumps to be made by MCMC (default is 1250)"),
		false,requires_argument),
	nburn(string("--bi,--burnin"),5000,
		string("Total num of jumps at start of MCMC to be discarded (default is 5000)"),
		false,requires_argument),
	nburn_noard(string("--bn,--burnin_noard"),0,
		string("num of burnin jumps before the ard is imposed (default is 0)"),
		false,requires_argument),
	sampleevery(string("--se,--sampleevery"),25,
		string("Num of jumps for each sample (MCMC) (default is 25)"),
		false,requires_argument),
	updateproposalevery(string("--upe,--updateproposalevery"),40,
		string("Num of jumps for each update to the proposal density std (MCMC) (default is 40)"),
		false,requires_argument),
	seed(string("--seed"),8219,
		string("seed for pseudo random number generator"),
		false,requires_argument),
	runLevMar(string("--runLevMar"),true,
		string("Run Levenberg(-Marquardt if activated) algorithm"),
		false, no_argument),
	useMarquardt(string("--useMarquardt"),true,
		string("Use Marquardt contribution in Levenberg algorithm"),
		false, no_argument),
	iterLevMar(string("iterLevMar"),200,
		string("Maximum number of iterations in Levenberg(-Marquardt) algorithm (default is 200)"),
		false,requires_argument),
   	runMCMC(string("--runMCMC"),false,
		string("Run MCMC algorithm and get distribution of estimates"),
		false,no_argument),
	no_ard(string("--noard"),false,
		string("Turn ARD off on all fibres"),
		false,no_argument),
	all_ard(string("--allard"),false,
		string("Turn ARD on on all fibres"),
		false,no_argument),
	rician(string("--rician"),false,
		string("Use Rician noise modelling in MCMC"),
		false,no_argument),
	CFP(string("--CFP"), string(""),
		string("File with a list of ASCCI files for specifying the common fixed parameters of the model"),
		false, requires_argument),  
	init_params(string("--init_params"), string(""),
		string("File with a list of NIFTI files for the initialization of the model parameters"),
		false, requires_argument),

   options("cudimot","cudimot --help (for list of options)\n")
     {
       try {
	options.add(verbose);
	options.add(help);
	options.add(logdir);
	options.add(forcedir);
	options.add(datafile);
	options.add(maskfile);
	options.add(subjdir);
	options.add(partsdir);
	options.add(outputdir);
	options.add(idPart);
	options.add(nParts);
	options.add(fudge);
	options.add(njumps);
	options.add(nburn);
	options.add(nburn_noard);
	options.add(sampleevery);
	options.add(updateproposalevery);
	options.add(seed);
	options.add(runMCMC); 
	options.add(no_ard);
	options.add(all_ard);
	options.add(runLevMar);
	options.add(useMarquardt);
	options.add(iterLevMar);
	options.add(rician);
	options.add(CFP);
	options.add(init_params);
     }
     catch(X_OptionError& e) {
       options.usage();
       cerr << endl << e.what() << endl;
     } 
     catch(std::exception &e) {
       cerr << e.what() << endl;
     }    
     
   }
}

#endif





