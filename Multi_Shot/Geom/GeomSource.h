/*
 * GeomSource.h
 *
 *  Created on: Apr 17, 2014
 *      Author: uwna
 *
 */

#ifndef GEOMSOURCE_H_
#define GEOMSOURCE_H_

class GeomSource {
private:
	double _sx;
	double _sy;
	double _sz;
	bool _live;
	long   _sortindex;

public:
	GeomSource();
	GeomSource(double sx, double sy, double sz, bool live, int sortindex);
	~GeomSource();

	double sx(){return _sx;}
	double sy(){return _sy;}
	double sz(){return _sz;}
	bool live(){return _live;}
	long sortindex(){return _sortindex;}
	void setSx(double sx){_sx= sx;}
	void setSy(double sy){_sy= sy;}
	void setSz(double sz){_sz= sz;}
	void setLive(bool live){_live= live;}
	void setSortindex(long sortindex){_sortindex= sortindex;}


	// copy operator
	GeomSource &operator=(const GeomSource &gt);
};

GeomSource::GeomSource(){

  _sx=0.l;
  _sy=0.l;
  _sz=0.l;
  _live=true;
  _sortindex=0;

}
GeomSource::GeomSource(double sx, double sy, double sz, bool live, int sortindex){

  _sx=sx;
  _sy=sy;
  _sz=sz;
  _live=live;
  _sortindex=sortindex;

}

GeomSource::~GeomSource(){}
GeomSource &GeomSource::operator=(const GeomSource &gt){

  _sx=gt._sx;
  _sy=gt._sy;
  _sz=gt._sz;
  _live=gt._live;
  _sortindex=gt._sortindex;
  
  return *this;

}

#endif /* GeomSource_H_ */
