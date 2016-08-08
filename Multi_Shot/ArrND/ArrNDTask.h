/*
 * ArrNDTask.h
 *
 *  Created on: Apr 19, 2014
 *      Author: uwna
 */

#ifndef ARRNDTASK_H_
#define ARRNDTASK_H_

#include <string>
#include "ArrND.h"

template <class T1, class T2> class ArrNDTask {
public:

	ArrND<T1> &_in;
	ArrND<T2> &_in_aux;
	ArrNDTask(ArrND<T1> &in, ArrND<T2> &in_aux):_in(in),_in_aux(in_aux){}
	~ArrNDTask(){};

	virtual void doit(int itask)=0;
};
#endif /* ARRNDTASK_H_ */
