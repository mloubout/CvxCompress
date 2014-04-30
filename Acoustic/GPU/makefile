CC=gcc
CFLAGS=-fopenmp -O2
OBJECTS=Parmfile_Reader.o Volume_Iterator.o Acoustic_Earth_Model.o Attribute_Reader.o Attribute_Transformer.o Index_Calculator.o UVW_XYZ_Calculator.o

gpuAcouTTI3D: gpuAcouTTI3D.cu $(OBJECTS)
	nvcc $(OBJECTS) -arch=sm_20 -O3 -Xptxas -dlcm=ca -Xcompiler=-fopenmp -Xptxas -v gpuAcouTTI3D.cu -o gpuAcouTTI3D

%.o: %.c
	$(CC) -c -I$(CUDAPATH)/include $(CFLAGS) $*.c

%.o: %.cpp
	$(CC) -c -I$(CUDAPATH)/include $(CFLAGS) $*.cpp

