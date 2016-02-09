CC=icc
CFLAGS=-openmp -O3 -fPIC -mavx -inline-forceinline
LDFLAGS=-shared -openmp
PAPI=-I/cpfs/lfs01/ESDRD/tjhc/vend/papi-5.4.3/src -I/cpfs/lfs01/ESDRD/tjhc/vend/papi-5.4.3/src/testlib -L/cpfs/lfs01/ESDRD/tjhc/vend/papi-5.4.3/src /cpfs/lfs01/ESDRD/tjhc/vend/papi-5.4.3/src/testlib/libtestlib.a /cpfs/lfs01/ESDRD/tjhc/vend/papi-5.4.3/src/libpapi.a

OBJECTS=CvxCompress.o Wavelet_Transform_Slow.o Wavelet_Transform_Fast.o Run_Length_Encode_Slow.o Block_Copy.o Read_Raw_Volume.o

all: CvxCompress_Test Test_Compression Compress_Sams_Stuff

libcvxcompress.so : Ds79_Base.cpp $(OBJECTS)
	$(CC) $(LDFLAGS) -o libcvxcompress.so $(OBJECTS)
	cp libcvxcompress.so ../lib

CvxCompress_Test: libcvxcompress.so CvxCompress_Test.cpp
	$(CC) $(CFLAGS) $^ -lrt -lstdc++ -o CvxCompress_Test

CvxCompress_GenCode: CvxCompress_GenCode.cpp CvxCompress.hxx Wavelet_Transform_Slow.cpp
	$(CC) -O2 CvxCompress_GenCode.cpp Wavelet_Transform_Slow.cpp -lrt -lstdc++ -o CvxCompress_GenCode

Test_Compression: libcvxcompress.so Test_Compression.cpp
	$(CC) $(CFLAGS) $^ $(PAPI) -lrt -lstdc++ -o Test_Compression

Compress_Sams_Stuff: libcvxcompress.so Compress_Sams_Stuff.cpp
	$(CC) $(CFLAGS) $^ $(PAPI) -lpthread -lrt -lstdc++ -o Compress_Sams_Stuff

Ds79_Base.cpp: CvxCompress.hxx CvxCompress_GenCode
	./CvxCompress_GenCode

%.o: %.f
	$(F77) -c $(FFLAGS) $*.f

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

%.o: %.cpp
	$(CC) -c $(CFLAGS) $*.cpp

clean:
	rm -f *.o
	rm -f libcvxcompress.so CvxCompress_Test CvxCompress_GenCode Test_Compression Compress_Sams_Stuff
	rm -f Ds79_Base.cpp Us79_Base.cpp
