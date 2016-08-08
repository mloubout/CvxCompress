/*
 * ArrNDReduceMpi.h
 *
 *  Created on: Aug 26, 2014
 *      Author: uwna
 */

#ifndef ARRNDREDUCEMPI_H_
#define ARRNDREDUCEMPI_H_

#include <mpi.h>
#include "ArrND.h"
#include "ArrNDReduce.h"

template <class T2> class ArrNDReduceMpi:public ArrNDReduce{

	ArrND<T2> &_outloc;
	ArrND<T2> &_out;
	ArrND<T2> *_buf;
	int _myrank;
	int _mode;

public:

	enum mode{copy,sum};

	ArrNDReduceMpi(ArrND<T2> &outloc, ArrND<T2> &out);

	~ArrNDReduceMpi();

	void reduceLocalMaster();
	void reduceLocalSlave(int tag);
	void reduceGlobalMaster(int myrank);

	void setMode(int mode){_mode=mode;}

};

template <class T2> ArrNDReduceMpi<T2>::ArrNDReduceMpi(ArrND<T2> &outloc, ArrND<T2> &out):_outloc(outloc),_out(out),_mode(0){
	MPI_Comm_rank(MPI_COMM_WORLD, &_myrank);
	if(_myrank==0){
		_buf=new ArrND<T2>(out.size());
	}
}
template <class T2> ArrNDReduceMpi<T2>::~ArrNDReduceMpi(){
	if(_myrank==0) delete _buf;
}
template <class T2> void ArrNDReduceMpi<T2>::reduceLocalSlave(int tag){

	// send back results. Need to be careful to send back the proper mpi tag
	cout<<"slave: sending in slave reduce "<<endl;
	MPI_Send(_outloc.datptr(),
			_outloc.len()*sizeof(T2),
			MPI_CHAR,
			0,
			tag,
			MPI_COMM_WORLD);
	cout<<"slave: done sending in slave reduce,  tag="<<tag<<endl;

}

template <class T2> void ArrNDReduceMpi<T2>::reduceLocalMaster(){

	MPI_Status status;

	cout<<"master: receiving in reduce"<<endl;
	MPI_Recv(_outloc.datptr(),
			_outloc.len()*sizeof(T2),
			MPI_CHAR,
			MPI_ANY_SOURCE,
			MPI_ANY_TAG,
			MPI_COMM_WORLD,
			&status);
	cout<<"master: received in reduce, tag="<<status.MPI_TAG<<endl;

	// write output
	cout<<"master: writing to output at "<<status.MPI_TAG<<" "<<*_outloc.datptr()<<endl;
	if(_mode==copy){
		_out[status.MPI_TAG]<<_outloc;// receive data from master
	}
	if(_mode==sum){
		(*_buf)+=_outloc;
	}

}
template <class T2> void ArrNDReduceMpi<T2>::reduceGlobalMaster(int myrank){

	if(_myrank==0){
		if(_mode==sum) _out<<*_buf;
	}

}
#endif /* ArrNDReduceMpi_H_ */
