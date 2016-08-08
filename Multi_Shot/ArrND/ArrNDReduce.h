/*
 * ArrNDReduce.h
 *
 *  Created on: Aug 26, 2014
 *      Author: uwna
 */

#ifndef ARRNDREDUCE_H_
#define ARRNDREDUCE_H_

class ArrNDReduce{

public:


	virtual void reduceLocalMaster()=0;
	virtual void reduceLocalSlave(int tag)=0;
	virtual void reduceGlobalMaster(int myrank)=0;

};

#endif /* ArrNDReduceMpi_H_ */
