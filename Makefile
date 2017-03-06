ifndef modelname
$(error modelname is not set)
endif

include $(FSLCONFDIR)/default.mk

NVCC = ${CUDA}/bin/nvcc
NVCC_FLAGS = -O3 -dc -Xptxas -v
CUDA_INC=-I${CUDA}/lib -I${CUDA}/lib64
#-maxrregcount=64
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

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_CPROB} -I${INC_PROB} -I${INC_BOOST} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

DLIBS = -lwarpfns -lbasisfield -lmeshclass -lnewimage -lutils -lmiscmaths -lnewmat -lnewran -lfslio -lniftiio -lznz -lcprob -lprob -lm -lz

SPLIT_PARTS=split_parts
MERGE_PARTS=merge_parts
CUDIMOT=${modelname}

SPLIT_PARTS_OBJS=split_parts.o cudimotoptions.o link_cudimot_gpu.o
MERGE_PARTS_OBJS=merge_parts.o cudimotoptions.o link_cudimot_gpu.o 
CUDIMOT_CUDA_OBJS=init_gpu.o dMRI_Data.o Model.o Parameters.o Levenberg_Marquardt.o MCMC.o
CUDIMOT_OBJS=link_cudimot_gpu.o cudimot.o cudimotoptions.o

SGEBEDPOST = bedpost
SGEBEDPOSTX = bedpostx bedpostx_postproc.sh bedpostx_preproc.sh bedpostx_single_slice.sh bedpostx_datacheck

XFILES = split_parts ${modelname} merge_parts

myclean:
	rm -f *.o
	rm -f *.so

all: ${XFILES}

${SPLIT_PARTS}: ${SPLIT_PARTS_OBJS}
		${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SPLIT_PARTS_OBJS} ${DLIBS} $(CUDIMOT_CUDA_OBJS) -lcudart -L${CUDA}/lib64 -L${CUDA}/lib

${MERGE_PARTS}: ${MERGE_PARTS_OBJS}
		${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MERGE_PARTS_OBJS} ${DLIBS} $(CUDIMOT_CUDA_OBJS) -lcudart -L${CUDA}/lib64 -L${CUDA}/lib

init_gpu.o: 
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) init_gpu.cu $(CUDA_INC)

dMRI_Data.o: 
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) dMRI_Data.cu $(CUDA_INC)

Parameters.o: 
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) Parameters.cu $(CUDA_INC)

Model.o: 
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) Model.cu $(CUDA_INC)

Levenberg_Marquardt.o: 	
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) Levenberg_Marquardt.cu $(CUDA_INC)

MCMC.o: 	
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) MCMC.cu $(CUDA_INC)

link_cudimot_gpu.o:	$(CUDIMOT_CUDA_OBJS)
		$(NVCC) $(GPU_CARDS) -dlink $(CUDIMOT_CUDA_OBJS) -o link_cudimot_gpu.o -L${CUDA}/lib64 -L${CUDA}/lib

cudimot.o:
		$(NVCC) $(GPU_CARDS) $(NVCC_FLAGS) cudimot.cc $(CUDA_INC)

${CUDIMOT}:	${CUDIMOT_OBJS}
		${CXX} ${CXXFLAGS} ${LDFLAGS} -o ${modelname} ${CUDIMOT_OBJS} $(CUDIMOT_CUDA_OBJS) ${DLIBS} -lcudart -L${CUDA}/lib64 -L${CUDA}/lib

