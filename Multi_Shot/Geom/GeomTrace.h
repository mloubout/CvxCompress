/*
 * GeomTrace.h
 *
 *  Created on: Apr 17, 2014
 *      Author: uwna
 *
 */

#ifndef GEOMTRACE_H_
#define GEOMTRACE_H_

#include <time.h>

class GeomTrace {
private:
	double _sx;
	double _sy;
	double _sz;
	double _rx;
	double _ry;
	double _rz;
	unsigned long _shotTime;
	long   _sortindex;

public:
	GeomTrace();
	GeomTrace(double sx, double sy, double sz, double rx, double ry, double rz, bool live, int sortindex);
	GeomTrace(double sx, double sy, double rx, double ry, bool live, int sortindex);
	~GeomTrace();

	// copy operator
	GeomTrace &operator=(const GeomTrace &gt);

	//Access functions

	void Encode_Bits(unsigned long& acc, int val, int first_bit, int nbits)
	{
		unsigned long mask = ((1 << nbits)-1);
		unsigned long lval = val;
		lval = (lval & mask) << first_bit;
		acc = acc & ~(mask << first_bit);
		acc = acc | lval;
	}

	int Decode_Bits(unsigned long acc, int first_bit, int nbits) const 
	{
		unsigned long mask = ((1 << nbits)-1);
		return (acc >> first_bit) & mask;
	}

	bool isLive() const {
		return Decode_Bits(_shotTime,0,1) != 0 ? true : false;
	}

	void setLive(bool live) {
		if (live)
		{
			Encode_Bits(_shotTime,1,0,1);
		}
		else
		{
			Encode_Bits(_shotTime,0,0,1);
		}
	}

	time_t getShotTime() const {
		struct tm tms;
		tms.tm_isdst = 0;
		tms.tm_sec = Decode_Bits(_shotTime,1,6);
		tms.tm_min = Decode_Bits(_shotTime,7,6);
		tms.tm_hour = Decode_Bits(_shotTime,13,5);
		tms.tm_mday = Decode_Bits(_shotTime,18,5);
		tms.tm_mon = Decode_Bits(_shotTime,23,4);
		tms.tm_year = Decode_Bits(_shotTime,27,16);
		time_t timestamp = mktime(&tms);
		//printf("timestamp = %s\n",asctime(gmtime(&timestamp)));
		return timestamp;
	}

	void setShotTime(const time_t shotTime) {
		time_t UTCshotTime = shotTime;
		struct tm* ptm = localtime(&UTCshotTime);
		Encode_Bits(_shotTime,ptm->tm_sec,1,6);
		Encode_Bits(_shotTime,ptm->tm_min,7,6);
		Encode_Bits(_shotTime,ptm->tm_hour,13,5);
		Encode_Bits(_shotTime,ptm->tm_mday,18,5);
		Encode_Bits(_shotTime,ptm->tm_mon,23,4);
		Encode_Bits(_shotTime,ptm->tm_year,27,16);
	}

	int getShotTimeUSec() const {
		return Decode_Bits(_shotTime,43,20);
	}

	void setShotTimeUSec(const int usec) {
		Encode_Bits(_shotTime,usec,43,20);
	}

	double getRx() const {
		return _rx;
	}

	void setRx(double rx) {
		_rx = rx;
	}

	double getRy() const {
		return _ry;
	}

	void setRy(double ry) {
		_ry = ry;
	}

	double getRz() const {
		return _rz;
	}

	void setRz(double rz) {
		_rz = rz;
	}

	long getSortindex() const {
		return _sortindex;
	}

	void setSortindex(long sortindex) {
		_sortindex = sortindex;
	}

	double getSx() const {
		return _sx;
	}

	void setSx(double sx) {
		_sx = sx;
	}

	double getSy() const {
		return _sy;
	}

	void setSy(double sy) {
		_sy = sy;
	}

	double getSz() const {
		return _sz;
	}

	void setSz(double sz) {
		_sz = sz;
	}
};

GeomTrace::GeomTrace(){

  _sx=0.l;
  _sy=0.l;
  _sz=0.l;
  _rx=0.l;
  _ry=0.l;
  _rz=0.l;
  _shotTime=0;
  _sortindex=0;

}
GeomTrace::GeomTrace(double sx, double sy, double sz, double rx, double ry, double rz, bool live, int sortindex){

  _sx=sx;
  _sy=sy;
  _sz=sz;
  _rx=rx;
  _ry=ry;
  _rz=rz;
  _shotTime=0;
  setLive(live);
  _sortindex=sortindex;

}
GeomTrace::GeomTrace(double sx, double sy, double rx, double ry, bool live, int sortindex){

  _sx=sx;
  _sy=sy;
  _sz=0.l;
  _rx=rx;
  _ry=ry;
  _rz=0.l;
  _shotTime=0;
  setLive(live);
  _sortindex=sortindex;

}
GeomTrace::~GeomTrace(){}
GeomTrace &GeomTrace::operator=(const GeomTrace &gt){

  _sx=gt._sx;
  _sy=gt._sy;
  _sz=gt._sz;
  _rx=gt._rx;
  _ry=gt._ry;
  _rz=gt._rz;
  _shotTime=gt._shotTime;
  _sortindex=gt._sortindex;
  
  return *this;

}

#endif /* GEOMTRACE_H_ */
