CPP=g++
CC=gcc
CFLAGS=-I. -I../../Common -fopenmp -O2 -fPIC
MODEL_OBJECTS=ElaOrthModel.o DFT.o Elastic_Modeling_Job.o Elastic_Shot.o Elastic_Resampled_Trace.o Elastic_Trace_Header.o Elastic_Trace.o Elastic_SEGY_File.o Elastic_SEGY_File_Receiver_Range.o ASCII2EBCDIC.o
PROPAGATOR_OBJECTS=Elastic_Buffer.o Elastic_Pipeline.o Elastic_Propagator.o
CUDA_FILES=Simple_Copy.cu Compute_ABC.cu Earth_Model.cu Propagate_Particle_Velocities.cu Propagate_Stresses_Orthorhombic.cu Elastic_Kernels.cu

all: gpuViscoElasticOrthorhomicFD
	cp ./gpuViscoElasticOrthorhomicFD ../../bin/fd3dElaGpu1shot.exe

# -DGPU-DEBUG
libelaorthogpu.so: $(CUDA_FILES) $(PROPAGATOR_OBJECTS) $(MODEL_OBJECTS) 
	nvcc $(PROPAGATOR_OBJECTS) $(MODEL_OBJECTS) -shared -use_fast_math -gencode arch=compute_30,code=sm_30 -O2 -Xptxas -dlcm=ca -Xcompiler=-fopenmp,-fPIC,-lcrypt -Xptxas -v Elastic_Kernels.cu -o $@
	cp libelaorthogpu.so ../../lib

gpuViscoElasticOrthorhomicFD: libelaorthogpu.so gpuViscoElasticOrthorhomicFD.cpp
	$(CPP) $(CFLAGS) -lelaorthogpu -lvoxet -lcrypt -L../../lib -I$(CUDAPATH)/include gpuViscoElasticOrthorhomicFD.cpp -o $@

testElaOrthModel:: libelaorthogpu.so testElaOrthModel.cpp
	$(CPP) -O2 -fopenmp -lelaorthogpu -lvoxet -lcrypt -L../../lib -I$(CUDAPATH)/include testElaOrthModel.cpp -o $@

%.o: %.c
	$(CC) -c -I$(CUDAPATH)/include $(CFLAGS) $*.c

%.o: %.cpp
	$(CPP) -c -I$(CUDAPATH)/include $(CFLAGS) $*.cpp

