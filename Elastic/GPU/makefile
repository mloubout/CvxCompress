CC=gcc
CFLAGS=-fopenmp -O2
MODEL_OBJECTS=DFT.o Elastic_Modeling_Job.o Elastic_Shot.o Elastic_SEGY_File.o Elastic_SEGY_File_Receiver_Range.o Global_Coordinate_System.o Voxet.o Voxet_Property.o 
PROPAGATOR_OBJECTS=Elastic_Buffer.o Elastic_Pipeline.o Elastic_Propagator.o
CUDA_FILES=Simple_Copy.cu Compute_ABC.cu Earth_Model.cu Propagate_Particle_Velocities.cu Propagate_Stresses_Orthorhombic.cu Elastic_Kernels.cu

# -DGPU-DEBUG
all: $(CUDA_FILES) $(PROPAGATOR_OBJECTS) $(MODEL_OBJECTS) gpuViscoElasticOrthorhomicFD.cpp
	nvcc $(PROPAGATOR_OBJECTS) $(MODEL_OBJECTS) gpuViscoElasticOrthorhomicFD.cpp -use_fast_math -arch=sm_30 -O -Xptxas -dlcm=ca -Xcompiler=-fopenmp -Xptxas -v Elastic_Kernels.cu -o gpuViscoElasticOrthorhomicFD

memcpy: memcpy.cpp
	nvcc memcpy.cpp -arch=sm_20 -O3 -Xptxas -dlcm=ca -Xcompiler=-fopenmp -Xptxas -v -o memcpy

Elastic_Modeling_Job.o: Elastic_Modeling_Job.cpp
	$(CC) -c -I$(CUDAPATH)/include -fopenmp -O Elastic_Modeling_Job.cpp

%.o: %.c
	$(CC) -c -I$(CUDAPATH)/include $(CFLAGS) $*.c

%.o: %.cpp
	$(CC) -c -I$(CUDAPATH)/include $(CFLAGS) $*.cpp

