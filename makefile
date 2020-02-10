CXX=icpc
CC=icc
CFLAGS=-qopenmp -O3 -fPIC -mavx -inline-forceinline -I.
LDFLAGS=-qopenmp -L. -lcvxcompress

# I think this disabled for now
PAPI_PATH := /devl/geophys/util/papi/papi-5.7.0
PAPI_INCLUDES := -DPAPI -I$(PAPI_PATH)/include
PAPI_LDFLAGS := -L$(PAPI_PATH)/lib -lpapi

OBJECTS=CvxCompress.o Wavelet_Transform_Slow.o Wavelet_Transform_Fast.o Run_Length_Encode_Slow.o Block_Copy.o Read_Raw_Volume.o

all: CvxCompress_Test Test_Compression Compress_Sams_Stuff Compress_Guojians_Stuff

libcvxcompress.so : Ds79_Base.cpp $(OBJECTS)
	$(CXX) -shared -qopenmp -o libcvxcompress.so $(OBJECTS)

CvxCompress_Test: libcvxcompress.so CvxCompress_Test.o
	$(CXX) $(LDFLAGS) $^  -o CvxCompress_Test

CvxCompress_GenCode: CvxCompress_GenCode.o CvxCompress.hxx Wavelet_Transform_Slow.o
	$(CXX) -O2 -qopenmp CvxCompress_GenCode.o Wavelet_Transform_Slow.o -o CvxCompress_GenCode

Test_Compression: libcvxcompress.so Test_Compression.o
	$(CXX) $(LDFLAGS) $^ $(PAPI_INCLUDES) -o Test_Compression

# requires libsegy which is not included in this repo, so don't include in "all" target for now. 
Compress_SEGY: Compress_SEGY.o libcvxcompress.so 
	$(CXX) $(LDFLAGS) $< -lsegy -o Compress_SEGY

Compress_Sams_Stuff: libcvxcompress.so Compress_Sams_Stuff.o
	$(CXX) $(LDFLAGS) $^  -o Compress_Sams_Stuff

Compress_Guojians_Stuff: Compress_Guojians_Stuff.o libcvxcompress.so 
	$(CXX) $(LDFLAGS) $<  -o Compress_Guojians_Stuff

Ds79_Base.cpp: CvxCompress.hxx CvxCompress_GenCode
	./CvxCompress_GenCode

# CvxCompress.o: CvxCompress.cpp
# 	$(CXX) -c -DPAPI $(CFLAGS) $(PAPI_INCLUDES) $*.cpp

Test_With_Generated_Input: Test_With_Generated_Input.o libcvxcompress.so 
	$(CXX) $(LDFLAGS) $<  -L. -lcvxcompress  -o $@

CvxCompress.o: CvxCompress.cpp
	$(CXX) -c $(CFLAGS) $(PAPI_INCLUDES) $< -o $@

%.o: %.f
	$(F77) -c $(FFLAGS) $*.f

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $*.cpp

clean:
	rm -f *.o
	rm -f libcvxcompress.so CvxCompress_Test CvxCompress_GenCode Test_Compression Compress_Sams_Stuff
	rm -f Ds79_Base.cpp Us79_Base.cpp
