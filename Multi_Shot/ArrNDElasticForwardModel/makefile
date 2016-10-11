all: fd3dElaGpu.exe
	cp ./fd3dElaGpu.exe ../../bin

fd3dElaGpu.exe: ArrNDElasticForwardModel.cpp
	icc -O2 -openmp -I.. -I/devl/geophys/atm/3.3/include/ -I../../Common -I../../Elastic/GPU/ -I$(CUDAPATH)/include/ -L/devl/geophys/atm/3.3/lib/intel -L../../lib/ -L$(CUDAPATH)/lib64 $^ -lvoxet -lelaorthogpu -lnuma -lappTaskMgr -lcudart_static -o $@
