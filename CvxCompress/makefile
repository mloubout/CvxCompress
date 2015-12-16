CC=icc
F77=ifort
CFLAGS=-openmp -O3 -fPIC -mavx -inline-forceinline -qopt-report=5
FFLAGS=-O3 -fPIC -mavx
OBJECTS=CvxCompress.o Wavelet_Transform_Slow.o Wavelet_Transform_Fast.o Run_Length_Encode_Slow.o Block_Copy.o Read_Raw_Volume.o

CvxCompress_Test: Ds79_Base.cpp CvxCompress_Test.cpp $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -lrt -lstdc++ CvxCompress_Test.cpp -o CvxCompress_Test

CvxCompress_GenCode: CvxCompress_GenCode.cpp CvxCompress.hxx Wavelet_Transform_Slow.cpp
	$(CC) -O2 Wavelet_Transform_Slow.cpp -lrt -lstdc++ CvxCompress_GenCode.cpp -o CvxCompress_GenCode

Test_Compression: Ds79_Base.cpp Test_Compression.cpp $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -lrt -lstdc++ Test_Compression.cpp -o Test_Compression

Quick_Test: Run_Length_Encode_Fast.cpp Run_Length_Encode_Slow.cpp Quick_Test.cpp
	$(CC) $(CFLAGS) -DDEBUG_ENCODE Run_Length_Encode_Fast.cpp Run_Length_Encode_Slow.cpp Quick_Test.cpp -lrt -lstdc++ -o Quick_Test

Ds79_Base.cpp: CvxCompress.hxx CvxCompress_GenCode
	./CvxCompress_GenCode

%.o: %.f
	$(F77) -c $(FFLAGS) $*.f

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

%.o: %.cpp
	$(CC) -c $(CFLAGS) $*.cpp
