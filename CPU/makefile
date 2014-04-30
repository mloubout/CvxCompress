#--------------------------------------------------------------
# Define the macro OBJECTS
CC=icc
CFLAGS=-openmp -O2 -inline-forceinline -restrict -DSTANDALONE
OBJECTS=Acoustic_Kernels_Main.o Acoustic_Kernels_SSE.o Acoustic_Kernels_AVX.o

all: tti3mod

tti3mod: $(OBJECTS)
	$(CC) $(OBJECTS) -openmp -lrt -lstdc++ -static -o $@ 

Acoustic_Kernels_AVX.o: Acoustic_Kernels_AVX.cpp
	$(CC) -c -mavx $(CFLAGS) $*.cpp

%.o: %.c
	$(CC) -c $(CUDAPATH)/include $(CFLAGS) $*.c

%.o: %.cpp
	$(CC) -c -I$(CUDAPATH)/include $(CFLAGS) $*.cpp

clean:
	rm -f *.o  

