CC=icc
F77=ifort
CFLAGS=-openmp -O3 -fPIC -mavx -inline-forceinline -I.
FFLAGS=-O3 -fPIC -mavx
LDFLAGS=-shared -openmp

all: libvoxet.so libsegy.so

libvoxet.so: Global_Coordinate_System.cpp Voxet.cpp Voxet_Property.cpp swapbytes.c
	$(CC) $(CFLAGS) $^ -shared -o libvoxet.so
	cp libvoxet.so ../lib

libsegy.so: SEGY_File.cpp swapbytes.c
	$(CC) $(CFLAGS) $^ -shared -o libsegy.so
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

