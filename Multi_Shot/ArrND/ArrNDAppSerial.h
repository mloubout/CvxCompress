/*
 * ArrNDAppSerial.h
 *
 *  Created on: Apr 18, 2014
 *      Author: uwna
 */

#ifndef ARRNDAPPSERIAL_H_
#define ARRNDAPPSERIAL_H_

#include <string>
#include <sstream>

#include "ArrND.h"
#include "ArrNDTask.h"


template<class T1, class T2> class ArrNDAppSerial {


protected:

	ArrND<T1> &_in;
	ArrND<T2> &_out;
	ArrNDTask<T1,T2> &_task;



public:
	ArrNDAppSerial(ArrND<T1> &in, ArrND<T2> &out, ArrNDTask<T1,T2> &task):_in(in),_out(out),_task(task){}
	~ArrNDAppSerial();

	void init(ArrND<T1> &in, ArrND<T2> &out, ArrNDTask<T1,T2> &task);
	void run();

};

template<class T1, class T2> void ArrNDAppSerial<T1,T2>::init(ArrND<T1>& in, ArrND<T2>& out, ArrNDTask<T1,T2>& task) {
	_in=in;
	_out=out;
	_task=task;

	// optional - gonna need factories for this
	//initParallel();
}
template<class T1, class T2> ArrNDAppSerial<T1,T2>::~ArrNDAppSerial() {

}


template<class T1, class T2> void ArrNDAppSerial<T1,T2>::run() {



	int ntasks = _in.size()[0];
	cout<<"ntasks="<<ntasks<<endl;
	for (int i = 0; i < ntasks; i++){

			cout<<"Doing task="<<i<<endl;

			//read
			cout<<"doing read"<<endl;
			_task._in=_in[i];
			cout<<"done read"<<endl;

			stringstream ss;
			 ss<<i;
			 string taskid=ss.str();
			_task.doit(taskid);

			//write
			cout<<"doing write"<<endl;
			_out[i]=_task._out;
			cout<<"done write"<<endl;

	}

}

#endif /* ARRNDAPPMASTSLAVE_H_ */
