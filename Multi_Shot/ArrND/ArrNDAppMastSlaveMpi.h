/*
 * ArrNDAppMastSlaveMpi.h
 *
 *  Created on: Apr 18, 2014
 *      Author: uwna
 */

#ifndef ARRNDAPPMASTSLAVEMPI_H_
#define ARRNDAPPMASTSLAVEMPI_H_

#include <mpi.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stack>

#include "ArrND.h"
#include "ArrNDReduce.h"
#include "ArrNDTask.h"

#define DIETAG 10000

template<class T1, class T2> class ArrNDAppMastSlaveMpi {

protected:

	ArrND<T1> &_in;
	ArrNDReduce *_reduce;
	ArrNDTask<T1> &_task;


	// Simple task manager
	std::stack<int> _atms;

	// data lock
	int _proclock;


	void master();
	void slave();

public:
	ArrNDAppMastSlaveMpi(ArrND<T1> &in, ArrNDTask<T1> &task):_in(in),_task(task),_reduce(0){}
	~ArrNDAppMastSlaveMpi();
	void setReduce(ArrNDReduce *reduce){_reduce=reduce;}

	int run();

};


template<class T1, class T2> ArrNDAppMastSlaveMpi<T1,T2>::~ArrNDAppMastSlaveMpi() {

	//delete _atm;

}
template<class T1, class T2> void ArrNDAppMastSlaveMpi<T1, T2>::master() {

	MPI_Status status;

	// set up stack
	int ntasks = _in.size()[0];
	for(int itask=0;itask<ntasks;itask++) _atms.push(itask);
	cout<<"Set up task list ntasks="<<ntasks<<endl;

	// seed the nodes
	int nrank;
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);
	cout<<"master: mpi nrank = "<<nrank<<endl;
	for(int irank=1;irank<nrank;irank++){

			cout<<"master: Doing irank="<<irank<<endl;

			// get the task
			int itask=_atms.top();
			_atms.pop();

			// read data in loc mem on master
			_task._in<<_in[itask];

			// send data to loc mem on node
			cout<<"master: sending in master, seed,  tag="<<itask<<endl;
			MPI_Send(_task._in.datptr(),
					 _task._in.len()*sizeof(T1),
					 MPI_CHAR,
					 irank,
					 itask,
					 MPI_COMM_WORLD);
			cout<<"master: done sending in master, seed,  tag="<<itask<<endl;
	}
	while(!_atms.empty()){

		// get task
		int itask=_atms.top();
		cout<<"master: popped itask="<<itask<<endl;
		_atms.pop();

		if(_reduce) _reduce->reduceLocalMaster();

		// recieve data lock
		cout<<"master: receiving proclock"<<endl;
		MPI_Recv(&_proclock,
				sizeof(int),
				MPI_CHAR,
				MPI_ANY_SOURCE,
				MPI_ANY_TAG,
				MPI_COMM_WORLD,
				&status);
		cout<<"master: received proclock, tag = "<<status.MPI_TAG<<endl;

		/*
		MPI_Recv(_task._out.datptr(),
				_task._out.len()*sizeof(T2),
				MPI_CHAR,
				MPI_ANY_SOURCE,
				MPI_ANY_TAG,
				MPI_COMM_WORLD,
				&status);

		// write output
		_out[status.MPI_TAG]<<_task._out;
		*/

		// read data in loc mem on master
		_task._in<<_in[itask];

		// send data to loc mem on node
		cout<<"master: sending in atms loop, tag="<<itask<<endl;
		MPI_Send(_task._in.datptr(),
				_task._in.len()*sizeof(T1),
				MPI_CHAR,
				status.MPI_SOURCE,
				itask,
				MPI_COMM_WORLD);
	}
	for(int irank=1;irank<nrank;irank++){
		cout<<"master: cleanup, doing receive"<<endl;
		if(_reduce) _reduce->reduceLocalMaster();
		/*
		MPI_Recv(_task._out.datptr(),
				_task._out.len()*sizeof(T2),
				MPI_CHAR,
				MPI_ANY_SOURCE,
				MPI_ANY_TAG,
				MPI_COMM_WORLD,
				&status);

		// write output
		_out[status.MPI_TAG]<<_task._out;
		*/
	}

	for(int irank=1;irank<nrank;irank++){
		cout<<"master: sending die tag"<<endl;
		MPI_Send(0,0,MPI_CHAR,irank,DIETAG,MPI_COMM_WORLD);
	}


}
template<class T1, class T2> void ArrNDAppMastSlaveMpi<T1, T2>::slave() {

	MPI_Status status;



	while(1){

		// receive data from master
		cout<<"slave: receiving"<<endl;
		MPI_Recv(_task._in.datptr(),
				_task._in.len()*sizeof(T1),
				MPI_CHAR,
				0,
				MPI_ANY_TAG,
				MPI_COMM_WORLD,
				&status);
		cout<<"slave: received tag="<<status.MPI_TAG<<endl;

		if(status.MPI_TAG==DIETAG) return;

		// do the task
		cout<<"slave: doing task, tag="<<status.MPI_TAG<<endl;
		_task.doit(status.MPI_TAG);

		if(_reduce) _reduce->reduceLocalSlave(status.MPI_TAG);

		// send back lock. Need to be careful to send back the proper mpi tag
		cout<<"slave:sending back proclock"<<endl;
		MPI_Send(&_proclock,
				sizeof(int),
				MPI_CHAR,
				0,
				status.MPI_TAG,
				MPI_COMM_WORLD);

		cout<<"slave: done sending back proclock, tag ="<<status.MPI_TAG<<endl;

		/*
		// send back results. Need to be careful to send back the proper mpi tag
		MPI_Send(_task._out.datptr(),
				_task._out.len()*sizeof(T2),
				MPI_CHAR,
				0,
				status.MPI_TAG,
				MPI_COMM_WORLD);
		*/
	}


}
template<class T1, class T2> int ArrNDAppMastSlaveMpi<T1, T2>::run() {

		int myrank;

		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		if(myrank==0){
			master();
		}else{
			slave();
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if(_reduce) _reduce->reduceGlobalMaster(myrank);

		return 1;
}

#endif /* ARRNDAPPMASTSlaveMpi_H_ */
