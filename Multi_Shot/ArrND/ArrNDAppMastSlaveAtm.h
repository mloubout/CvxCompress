/*
 * ArrNDAppMastSlave.h
 *
 *  Created on: Apr 18, 2014
 *      Author: uwna
 */

#ifndef ARRNDAPPMASTSLAVEATM_H_
#define ARRNDAPPMASTSLAVEATM_H_

#include <string>
#include <fstream>
#include <sstream>
#include "ArrND.h"
#include "ArrNDTask.h"
#include "appTaskMgr.h"

template<class T1, class T2> class ArrNDAppMastSlaveAtm {


protected:

	ArrND<T1> &_in;
	ArrND<T2> &_in_aux;
	ArrNDTask<T1,T2> &_task;

	std::string _atmfile;
	int _taskTimeoutSecs;
	int _sleepTime;

	void doAtm(AppTaskMgr &atm);

public:
	ArrNDAppMastSlaveAtm(ArrND<T1> &in, ArrND<T2> &in_aux, ArrNDTask<T1,T2> &task, std::string atmfile);
	~ArrNDAppMastSlaveAtm();

	int getTaskTimeOutSecs() {return _taskTimeoutSecs;}
	int getSleepTime() {return _sleepTime;}

	void setTaskTimeOutSecs(int t) {_taskTimeoutSecs = t;}
	void setSleepTime(int t) {_sleepTime = t;}

	void run();

};

template<class T1, class T2> ArrNDAppMastSlaveAtm<T1,T2>::ArrNDAppMastSlaveAtm(ArrND<T1> &in, ArrND<T2> &in_aux, ArrNDTask<T1,T2> &task, std::string atmfile):
		_in(in),_in_aux(in_aux),_task(task),_atmfile(atmfile){

	_taskTimeoutSecs=300;
	_sleepTime=30;
}

template<class T1, class T2> void ArrNDAppMastSlaveAtm<T1,T2>::doAtm(AppTaskMgr &atm){

	try {
		for ( AppTaskMgr::taskInfo tInfo = atm.taskOut(_taskTimeoutSecs); !tInfo.complete; // Exit loop when complete
				tInfo = atm.taskOut(_taskTimeoutSecs)) {
			string taskId = tInfo.taskId;

			if (taskId == "initialize") {

				// set up task list
				int ntasks = _in.size()[0];
				std::string tmp = _atmfile + "_e"+".tasks";
				std::ofstream outs(tmp.c_str());
				for (int itask = 0; itask < ntasks; itask++) {
					outs << itask<<endl;
				}
				outs.close();

				// check in task
				atm.taskIn(tInfo.taskId);
				cout<<"Time out="<<_taskTimeoutSecs<<endl;
				cout<<"Sleep time="<<_sleepTime<<endl;

			}
			else if(taskId=="finalize"){

			}
			else if (taskId != "") {

				// do task
				int itask = atoi(tInfo.taskId.taskData.c_str());
				cout<<"***doing itask="<<itask<<endl;
				_task._in<<_in[itask];
				_task._in_aux<<_in_aux[itask];
				_task.doit(itask);

				// check in task
				atm.taskIn(tInfo.taskId);

			}
			else {
				cout<<"taskId=="<<taskId<<endl;
				sleep (_sleepTime);
			}
		}

	} catch (AppTaskMgr::ATMException &e) // Catches both ATMIOException and ATMInvalidTaskException
	{
		cerr << "Caught ATMException in ArrNDAppMastSlaveAtm::doAtm: " << e.msg << endl;
		exit(-1);
	}

}

template<class T1, class T2> ArrNDAppMastSlaveAtm<T1,T2>::~ArrNDAppMastSlaveAtm() {

}

template<class T1, class T2>  void ArrNDAppMastSlaveAtm<T1,T2>::run() {

	{
		std::string tmp = _atmfile + "_i";
		std::string tmp1 = _atmfile + "_i"+".tasks";

		std::ofstream outs(tmp1.c_str());
		outs << "initialize";
		outs.close();

		try {
			AppTaskMgr atm(tmp);
			doAtm(atm);

		} catch (AppTaskMgr::ATMException &e) {
			// Catches both ATMIOException and ATMInvalidTaskException
			cerr << " Caught ATMException: " << e.msg << endl;
			exit(-1);
		} catch (...) {
			cerr << " Unexpected error while creating initializing AppTaskMgr."
					<< endl;
			exit(-1);
		}

	}
	{
		std::string tmp = _atmfile + "_e";

		try {
			AppTaskMgr atm(tmp);
			doAtm(atm);
		} catch (AppTaskMgr::ATMException &e) {
			// Catches both ATMIOException and ATMInvalidTaskException
			cerr << " Caught ATMException: " << e.msg << endl;
			exit(-1);
		} catch (...) {
			cerr << " Unexpected error while creating initializing AppTaskMgr."
					<< endl;
			exit(-1);
		}

	}
}

#endif /* ARRNDAPPMASTSLAVEATM_H_ */
