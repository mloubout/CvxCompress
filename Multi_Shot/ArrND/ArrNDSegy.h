/*
 * SegyND.cpp


 * TODO: Need to add geometry info in output.
 *  Created on: May 29, 2014
 *      Author: atzv
 */
#ifndef ARRNDSEGY_H_
#define ARRNDSEGY_H_

#include "ArrND.h"
#include "GeomTrace.h"
#include "/home/atzv/workspace/SegyWriterCPP/src/SegYFileHeader.h"
#include "/home/atzv/workspace/SegyWriterCPP/src/SegYTraceHeader.h"
#include "/home/atzv/workspace/SegyWriterCPP/src/SegYWriterUtil.h"


class ArrNDSegy: public ArrND<float> {

public:
	SegYFileHeader fhdr;
	ArrNDSegy(vector<long> size, string filename);
	ArrNDSegy(vector<long> size, string filename, ArrND<GeomTrace>);
	~ArrNDSegy();
	/**
	 * Override to write to Seg-Y file.
	 */
	void transfer(ArrND<float> &arr);
	ArrNDSegy &operator<<(ArrND<float> &arr);
};

/**
 *
 * File does not have the .segy extension
 */
ArrNDSegy::ArrNDSegy(std::vector<long> size, std::string filename):ArrND<float>(size, filename){
}

ArrNDSegy::~ArrNDSegy(){
	if(_dat!=0) delete [] _dat;
	_dats.close();

}

void ArrNDSegy::transfer(ArrND<float> &arr) {
	if(_tfrsze!=arr._tfrsze){
		cerr<<"Error in mem transfer: sizes do not agree"<<endl;
		cerr<<"tfrsize="<<_tfrsze<<endl;
		cerr<<"arr.tfrsize="<<arr._tfrsze<<endl;
	}

	float* src;

	if (arr._media == 0) {
		src = arr.datptr();
	}
	if (arr._media == 1) {
		arr.read(src, _tfrsze);
	}
	_dats.seekp(ios::beg);

	//Getting number of traces
	short nTraces =_tfrsze/_size[_tfrdim];
	short nSamp = arr._size[_tfrdim];

	fhdr.Samples_Per_Trace = nSamp;
	fhdr.writeSegYFileHeader(_dats);
	std::cout<<_dats.tellp()<<std::endl;

	for(int i=0;i<nTraces;i++) {
		float* data;
		std::copy(src+(i*nSamp),src+((i+1)*nSamp),data);
		SegYTraceHeader trcHdr;
		trcHdr.cdp = i;
		trcHdr.writeTraceHeader(_dats);
		SegYWriterUtil::swap4bytes((int *)data,nSamp);
		_dats.write((char *)data,nSamp*sizeof(float));
		std::cout<<_dats.tellp()<<std::endl;
	}

return;
}

 ArrNDSegy &ArrNDSegy::operator<<(ArrND<float> &arr){

	if(_tfrsze!=arr._tfrsze){
		cerr<<"Error in mem transfer: sizes do not agree"<<endl;
	}

	transfer(arr);
	 std::cout<<"Got here"<<std::endl;

	_tfrsze=_maxsze;
	_tfrloc=0;
	_tfrdim=0;
	arr._tfrsze=arr._maxsze;
	arr._tfrloc=0;
	arr._tfrdim=0;

	return *this;
}
#endif
