all:
	/vend/intel/impi/4.0.1.007/intel64/bin/mpif90 -fc=ifort -vec-report3 -O3 -inline all \
	-xW -I /vend/intel/impi/4.0.1.007/intel64/include \
	-L /vend/intel/impi/4.0.1.007/intel64/lib -openmp -c \
	bessi0.f90 \
	computeStrainRateElaOnKite.f90 \
	computeStrainRateElaOffKite.f90 \
	fit.f90 \
	generateABCFunction.f90 \
	generateSincWeights.f90 \
	interpolateTrace.f90 \
	prepWavelet.f90 \
	readModelProperties.f90 \
	rotateCIJLocalToGlobal.f90 \
        coordGlobalToLocal.f90 \
        coordLocalToGlobal.f90 \
	smoothSeafloor.f90 \
	updateParticleVelocity.f90 \
        updateStressElaIso.f90 \
        updateStressElaVTI.f90 \
        updateStressElaVOr.f90 \
        updateStressElaTri.f90 \
	writeGocadSnapShot.f90 \
	FDAniQ.f90

	/vend/intel/composerxe-2011.2.137/bin/intel64/icc -c \
	swapBytes.c \
	getSign.c \
	writeSEGY.c 

	/vend/intel/impi/4.0.1.007/intel64/bin/mpif90 -fc=ifort -O3 -inline all \
        -xW -I /vend/intel/impi/4.0.1.007/intel64/include \
	-L /vend/intel/impi/4.0.1.007/intel64/lib -openmp \
	bessi0.o \
	computeStrainRateElaOnKite.o \
	computeStrainRateElaOffKite.o \
	fit.o \
	generateABCFunction.o \
	generateSincWeights.o \
	interpolateTrace.o \
	prepWavelet.o \
	readModelProperties.o \
	rotateCIJLocalToGlobal.o \
        coordGlobalToLocal.o \
        coordLocalToGlobal.o \
	smoothSeafloor.o \
	updateParticleVelocity.o \
        updateStressElaIso.o \
        updateStressElaVTI.o \
        updateStressElaVOr.o \
        updateStressElaTri.o \
	writeGocadSnapShot.o \
	swapBytes.o \
	getSign.o \
	writeSEGY.o \
	FDAniQ.o \
	-o xFDA

clean:
	rm -f xFDA
	rm -f *.o
