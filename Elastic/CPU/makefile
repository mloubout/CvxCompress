CC=icc
CFLAGS=-openmp -O3 -mavx -inline-forceinline
OBJECTS=Vanilla_Visco_Elastic_CPU_Propagator.o

# -DGPU-DEBUG
all: $(OBJECTS) main.cpp
	$(CC) $(CFLAGS) $(OBJECTS) main.cpp -lrt -o Vanilla

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c

%.o: %.cpp
	$(CC) $(CFLAGS) -c $*.cpp

