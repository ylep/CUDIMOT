ifndef modelname
$(error Error: modelname has not been set)
endif

include $(FSLCONFDIR)/default.mk

NVCC = ${CUDA}/bin/nvcc
BIN=./bin

MAX_REGISTERS= -maxrregcount 255
TYPE := $(shell cat compileOptions/type.h | grep MyType | cut -f 2 -d ' ')
ifeq ($(TYPE),float)
  MAX_REGISTERS= -maxrregcount 64
else
  MAX_REGISTERS= -maxrregcount 128
endif

MODELPATH= mymodels/$(modelname)

NVCC_FLAGS = -I$(MODELPATH) -O3 -dc $(MAX_REGISTERS)
#-Xptxas -v 
#-G -lineinfo

CUDA_INC=-I${CUDA}/lib -I${CUDA}/lib64

CUDA_INC = -I. -I${FSLDIR}/extras/include/newmat -I${FSLDIR}/include
SM_20 = -gencode arch=compute_20,code=sm_20
SM_21 = -gencode arch=compute_20,code=sm_21
SM_30 = -gencode arch=compute_30,code=sm_30
SM_35 = -gencode arch=compute_35,code=sm_35
SM_37 = -gencode arch=compute_37,code=sm_37
SM_50 = -gencode arch=compute_50,code=sm_50
SM_52 = -gencode arch=compute_52,code=sm_52
SM_60 = -gencode arch=compute_60,code=sm_60
SM_61 = -gencode arch=compute_61,code=sm_61
GPU_CARDS = $(SM_35)

PROJNAME = CUDIMOT

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_CPROB} -I${INC_PROB} -I${INC_BOOST} -I${INC_ZLIB} -I$(MODELPATH)
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

DLIBS = -lwarpfns -lbasisfield -lmeshclass -lnewimage -lutils -lmiscmaths -lnewmat -lnewran -lfslio -lniftiio -lznz -lcprob -lprob -lm -lz

CUDIMOT=bin/${modelname}

CUDIMOT_CUDA_OBJS=bin/modelparameters.o bin/init_gpu.o bin/dMRI_Data.o bin/Model.o bin/Parameters.o bin/GridSearch.o bin/Levenberg_Marquardt.o bin/MCMC.o bin/getPredictedSignal.o

CUDIMOT_OBJS=bin/link_cudimot_gpu.o bin/cudimot.o bin/cudimotoptions.o

SGEBEDPOST = bedpost
SGEBEDPOSTX = bedpostx bedpostx_postproc.sh bedpostx_preproc.sh bedpostx_single_slice.sh bedpostx_datacheck

SCRIPTS = ${modelname}@info
FILES = cart2spherical split_parts_${modelname} ${modelname} merge_parts_${modelname} testFunctions_${modelname}
XFILES=$(addprefix $(BIN)/, $(FILES))

cleanall:
	rm bin/*.o
#find ./ -maxdepth 1 -type f ! -name "*.*" ! -name "Makefile" -delete
cleanbin:
	rm bin/*
debugging:
	rm -f bin/testFunctions_${modelname}
	rm -f bin/Levenberg_Marquardt.o
	rm -f bin/MCMC.o
	make install

all: ${XFILES}

bin/cart2spherical: 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ cart2spherical.cc ${DLIBS} 

bin/cudimotoptions.o:
	${CXX} ${CXXFLAGS} ${LDFLAGS} -c -o $@ cudimotoptions.cc ${DLIBS} 

bin/split_parts_${modelname}: bin/cudimotoptions.o bin/link_cudimot_gpu.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ bin/cudimotoptions.o bin/link_cudimot_gpu.o split_parts.cc ${DLIBS} $(CUDIMOT_CUDA_OBJS) -lcudart -lboost_filesystem -lboost_system -L${CUDA}/lib64 -L${CUDA}/lib

bin/merge_parts_${modelname}: bin/cudimotoptions.o bin/link_cudimot_gpu.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ bin/cudimotoptions.o bin/link_cudimot_gpu.o merge_parts.cc ${DLIBS} $(CUDIMOT_CUDA_OBJS) -lcudart -lboost_filesystem -lboost_system -L${CUDA}/lib64 -L${CUDA}/lib

bin/init_gpu.o: 
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) -o $@ init_gpu.cu $(CUDA_INC)

bin/dMRI_Data.o: 
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) -o $@ dMRI_Data.cu $(CUDA_INC)

bin/Parameters.o: 
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) -o $@ Parameters.cu $(CUDA_INC)

bin/Model.o: 
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) -o $@ Model.cu $(CUDA_INC)

bin/modelparameters.o: 
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) $(MODELPATH)/modelparameters.cc -o bin/modelparameters.o $(CUDA_INC)

bin/GridSearch.o: 	
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) -o $@ GridSearch.cu $(CUDA_INC)

bin/Levenberg_Marquardt.o: 	
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) -o $@ Levenberg_Marquardt.cu $(CUDA_INC)

bin/MCMC.o: 	
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) -o $@ MCMC.cu $(CUDA_INC)

bin/getPredictedSignal.o: 	
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) -o $@ getPredictedSignal.cu $(CUDA_INC)

bin/link_cudimot_gpu.o:	$(CUDIMOT_CUDA_OBJS)
		$(NVCC) $(GPU_CARDS) -dlink $(CUDIMOT_CUDA_OBJS) -o $@ -L${CUDA}/lib64 -L${CUDA}/lib

bin/cudimot.o:
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) -o $@ cudimot.cc $(CUDA_INC)

${CUDIMOT}:	${CUDIMOT_OBJS}
		${CXX} ${CXXFLAGS} ${LDFLAGS} -o bin/${modelname} ${CUDIMOT_OBJS} $(CUDIMOT_CUDA_OBJS) ${DLIBS} -lcudart -L${CUDA}/lib64 -L${CUDA}/lib
		./saveModel_info.sh

bin/testFunctions_${modelname}: 
	$(NVCC) $(GPU_CARDS) -I$(MODELPATH) -O3 $(MAX_REGISTERS) $(MODELPATH)/modelparameters.cc testFunctions.cu -o bin/testFunctions_${modelname} $(CUDA_INC)
