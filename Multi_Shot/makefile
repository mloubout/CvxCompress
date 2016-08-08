all: common singleshot multishot utilities

common:
	mkdir -p ../lib
	cd ../Common; make -B -f ./makefile

singleshot:
	mkdir -p ../lib
	cd ../Elastic/GPU; make -B -f ./makefile

multishot: 
	mkdir -p ../bin
	cd ../Multi_Shot/ArrNDElasticForwardModel; make -B -f ./makefile

utilities:
	mkdir -p ../bin
	cd Utils; make -B -f ./makefile

clean:
	rm -f ../lib/*
	rm -f ../bin/*

