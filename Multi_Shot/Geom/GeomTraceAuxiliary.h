/*
* GeomTrace.h
 *
 *  Created on: Aug 02, 2016
 *      Author: tjhc
 *
 */

#ifndef GEOMTRACEAUXILIARY_H_
#define GEOMTRACEAUXILIARY_H_

class GeomTraceAuxiliary {
private:
	int _sffid;	// source ffid
	int _seqno;	// source sequence number
	int _sil;	// source inline
	int _sxl;	// source xline
	int _rffid;	// receiver ffid
	int _ril;	// receiver inline
	int _rxl;	// receiver xline

public:
	GeomTraceAuxiliary();
	GeomTraceAuxiliary(GeomTraceAuxiliary& gt);
	~GeomTraceAuxiliary();

	// copy operator
	GeomTraceAuxiliary &operator=(const GeomTraceAuxiliary &gt);

	//Access functions

	int getSourceFFID() {return _sffid;}
	void setSourceFFID(int ffid) {_sffid=ffid;}

	int getSeqNo() {return _seqno;}
	void setSeqNo(int seqno) {_seqno=seqno;}

	int getSourceInline() {return _sil;}
	void setSourceInline(int il) {_sil=il;}

	int getSourceXline() {return _sxl;}
	void setSourceXline(int xl) {_sxl=xl;}

	int getReceiverFFID() {return _rffid;}
	void setReceiverFFID(int ffid) {_rffid=ffid;}
	
	int getReceiverInline() {return _ril;}
	void setReceiverInline(int il) {_ril=il;}

	int getReceiverXline() {return _rxl;}
	void setReceiverXline(int xl) {_rxl=xl;}
};

GeomTraceAuxiliary::GeomTraceAuxiliary(){

	_sffid = 0;
	_seqno = 0;
	_sil = 0;
	_sxl = 0;
	_rffid = 0;
	_ril = 0;
	_rxl = 0;
}
GeomTraceAuxiliary::GeomTraceAuxiliary(GeomTraceAuxiliary& gt)
{
	_sffid = gt._sffid;
	_seqno = gt._seqno;
	_sil = gt._sil;
	_sxl = gt._sxl;
	_rffid = gt._rffid;
	_ril = gt._ril;
	_rxl = gt._rxl;
}
GeomTraceAuxiliary::~GeomTraceAuxiliary(){}
GeomTraceAuxiliary& GeomTraceAuxiliary::operator=(const GeomTraceAuxiliary &gt){

	_sffid = gt._sffid;
	_seqno = gt._seqno;
	_sil = gt._sil;
	_sxl = gt._sxl;
	_rffid = gt._rffid;
	_ril = gt._ril;
	_rxl = gt._rxl;
  
  return *this;

}

#endif /* GEOMTRACEAUXILIARY_H_ */
