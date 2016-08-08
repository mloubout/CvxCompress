CC=icc
F77=ifort
CFLAGS=-openmp -O2 -fPIC -mavx -inline-forceinline -I.
FFLAGS=-O2 -fPIC -mavx
LDFLAGS=-shared -openmp

all: libvoxet.so libsegy.so
	cp ./libvoxet.so ../lib
	cp ./libsegy.so ../lib

libvoxet.so: Global_Coordinate_System.cpp Voxet.cpp Voxet_Property.cpp swapbytes.c Voxet_Memory_Mapper.cpp
	$(CC) $(CFLAGS) $^ -shared -o $@
	cp libvoxet.so ../lib

libsegy.so: SEGY_File.cpp swapbytes.c
	$(CC) $(CFLAGS) $^ -shared -o $@
	cp libsegy.so ../lib

%.o: %.f
	$(F77) -c $(FFLAGS) $*.f

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

%.o: %.cpp
	$(CC) -c $(CFLAGS) $*.cpp

clean:
	rm -f *.o
	rm -f libvoxet.so
	rm -f libsegy.so

