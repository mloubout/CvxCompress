/*
 * ArrND.h
 *
 *  Created on: Apr 18, 2014
 *      Author: uwna
 */

#ifndef ARRND_H_
#define ARRND_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex.h>
#include <list>
#include <map>
#include <sys/types.h>
#include <stdexcept>
#include <cstdlib>
#include <cstring>


using namespace std;

template<class T> class ArrND {
friend class ArrNDSegy;
protected:

	// media type
	bool _media;

	// data
	T *_dat;

	// vars
	int _dim;
	std::vector<long> _size;

	// transfer params
	long _tfrdim;
	long _tfrloc;
	long _tfrsze;
	long _maxsze;

	// file data
	std::string	_filename;
	std::fstream _dats;

	void openDat();
	void transfer(ArrND &arr);
	void transferVarSize(ArrND &arr);
	void init();

	//Header file stuff
	std::vector<double> _delta;
	std::vector<double> _origin;
	std::vector<std::string> _units;
	std::vector<std::string> _axis;
	std::string _format;
	std::fstream _heds;
	void openHed();
	void closeHed();
	void getRegExp(regmatch_t* m, regex_t& reg, std::string& line, std::string& name, std::string& value);

	// Option to output variable ArrND sizes across the first dimension - option varsize
	bool _mode_varsize;
	long _varsize_index;
	size_t _varsize_offset_global;
	size_t *_varsize_offset;
	size_t *_varsize_offset_ptr;



public:

	// general functions
	ArrND(std::vector<long> size);
	ArrND(std::vector<long> size, T* dat);
	ArrND(std::vector<long> size, std::string filename);
	ArrND(std::string filename);
	~ArrND();

	T* datptr();
	T& dat();
	std::vector<long>  &size(){return _size;}
	long len(){return _maxsze;}

	ArrND &operator[](long i);
	ArrND &operator<<(ArrND &arr);
	ArrND &operator+=(ArrND &arr);

	bool operator==(ArrND &arr){
		if(this==arr) return 1;
		return 0;
	}


	void read(T *dat, long n);
	void write(T *dat, long n);

	void readHed(std::vector<long>& size, std::vector<double>& delta, std::vector<double>& origin, std::vector<std::string>& axis, std::vector<std::string>& units,
			std::string& format, int& dim);
	void writeHed(std::vector<double> delta, std::vector<double> origin, std::vector<std::string> axis, std::vector<std::string> units,
			std::string format);

	// functions related to varsize
	void setModeVarSize();
	void initVarSize();
	void cleanVarSize();
	vector<long> sizeVarSize(int i);

	// mem allocation switch
	bool _memswitch;

};
template<class T> ArrND<T>::ArrND(std::vector<long> size){
	_memswitch=1;
	_size=size;
	init();
}
template<class T> ArrND<T>::ArrND(std::vector<long> size, T* dat){
	_memswitch=0;
	_dat=dat;
	_size=size;
	init();
}/**
 * This constructor creates a file object from an existing file
 * Transfer variables are initialized, data stream is opened.
 * Header file is not read, header variables are not initialized!
 */
template<class T> ArrND<T>::ArrND(std::vector<long> size, std::string filename){
	_dat=0;
	_size=size;
	_dim=size.size();
	_tfrsze=1;
	_tfrloc=0;
	for(int i=0;i<_dim;i++) _tfrsze*=_size[i];
	_maxsze=_tfrsze;
	_tfrdim=0;

	// open file
	_filename=filename;
	openDat();
	_media=1;

	// var size
	initVarSize();

}
template<class T> void ArrND<T>::init(){
	_dim=_size.size();
	_tfrsze=1;
	_tfrloc=0;
	for(int i=0;i<_dim;i++) _tfrsze*=_size[i];
	_maxsze=_tfrsze;
	_tfrdim=0;

	// allocate mem
	if(_memswitch) {
		_dat= new T[_maxsze];

		for (int i=0;i<_maxsze;i++) {
			_dat[i] = T();
		}
	}
	_media=0;

	// var size
	initVarSize();
}
template<class T> ArrND<T>::~ArrND(){
	if(_dat!=0&&_memswitch) delete [] _dat;
	_dats.close();

	// var size
	cleanVarSize();
}

template<class T> ArrND<T>::ArrND(std::string filename){

	_filename=filename;

	readHed(_size, _delta, _origin, _axis, _units, _format, _dim);

	_dat=0;
	_tfrsze=1;
	_tfrloc=0;
	for(int i=0;i<_dim;i++) _tfrsze*=_size[i];
	_maxsze=_tfrsze;
	_tfrdim=0;

	// open file
	_filename=filename;
	openDat();
	_media=1;

	// var size
	initVarSize();

}

template<class T> void ArrND<T>::initVarSize(){
	_mode_varsize=0;
	_varsize_offset=0L;
}
template<class T> void ArrND<T>::cleanVarSize(){
	if(_varsize_offset!=0L) delete _varsize_offset;
}
template<class T> void ArrND<T>::setModeVarSize(){
	_mode_varsize=1;
	_varsize_offset_global=0;
	_varsize_offset = new size_t[_size[0]];
	_varsize_offset_ptr= _varsize_offset;
}
template<class T> void ArrND<T>::write(T *dat, long n){
	_dats.seekp(_tfrloc*sizeof(T), ios_base::beg);
	_dats.write((char *) dat,n*sizeof(T));
}
template<class T> void ArrND<T>::read(T *dat, long n){
	_dats.seekg(_tfrloc*sizeof(T), ios_base::beg);
	_dats.read((char *) dat,n*sizeof(T));
}
template<class T> void ArrND<T>::transfer(ArrND<T> &arr){

	//cout<<"Doing transfer, tfrloc="<<_tfrloc<<" "<<arr._tfrloc<<endl;

	if(_tfrsze!=arr._tfrsze){
		cerr<<"Error in mem transfer: sizes do not agree"<<endl;
	}

	if (this->_media == 0 && arr._media == 0) {
		T* dest = datptr();
		T* src = arr.datptr();
		memcpy(dest, src, _tfrsze * sizeof(T));
		return;
	}
	if (_media == 1 && arr._media == 0) {
		write(arr.datptr(),arr._maxsze);
		return;
	}
	if (_media == 0 && arr._media == 1) {
		arr.read(datptr(), _tfrsze);
		return;
	}

}
template<class T> vector<long> ArrND<T>::sizeVarSize(int i){

	int dim;

	_dats.seekg(_varsize_offset[i], ios_base::beg);
	_dats.read((char *)&dim, sizeof(int));

	vector<long> size(dim);
	for(int i=0;i<dim;i++){
		long tmp;
		_dats.read((char *)&tmp, sizeof(long));
		size[i]=tmp;
	}
	return size;

}
template<class T> void ArrND<T>::transferVarSize(ArrND<T> &arr){

	//cout<<"Doing transfer, tfrloc="<<_tfrloc<<" "<<arr._tfrloc<<endl;

	if (this->_media == 0 && arr._media == 0) {
		if(_tfrsze!=arr._tfrsze){
				cerr<<"Error in mem transfer: sizes do not agree"<<endl;
			}
		T* dest = datptr();
		T* src = arr.datptr();
		memcpy(dest, src, _tfrsze * sizeof(T));
		return;
	}
	if (_media == 1 && arr._media == 0) {
		*_varsize_offset_ptr=_varsize_offset_global;
		_varsize_offset_ptr++;
		int dim = arr._dim;
		_dats.seekp(_varsize_offset_global, ios_base::beg);
		_dats.write((char *)&dim, sizeof(int));
		_varsize_offset_global+=sizeof(int);
		for (int i = 0; i < dim; i++) {
			_dats.write((char *)arr.size()[i], sizeof(long));
			_varsize_offset_global+=sizeof(long);
		}
		_dats.write((char *)arr.datptr(),arr._maxsze*sizeof(T));
		_varsize_offset_global+=arr._maxsze*sizeof(T);
		return;
	}
	if (_media == 0 && arr._media == 1) {
		size_t offset = arr._varsize_offset[arr._varsize_index];
		offset+=sizeof(int);
		for(int i=0;i<_dim;i++) offset+=sizeof(long);
		_dats.seekg(offset, ios_base::beg);
		arr.read(datptr(), _maxsze);
		return;
	}

}
template<class T> T* ArrND<T>::datptr(){
	T* tmp=_dat+_tfrloc;
	_tfrsze=_maxsze;
	_tfrdim=0;
	_tfrloc=0;
	return (tmp);
}
template<class T> T& ArrND<T>::dat(){
	T* tmp=_dat+_tfrloc;
	_tfrsze=_maxsze;
	_tfrdim=0;
	_tfrloc=0;
	return (*tmp);
}

/**
 * Create and write to a header file (_size and _filname are from the object, other information is provided by user)
 */
template<class T> void ArrND<T>::writeHed(std::vector<double> delta, std::vector<double> origin,
		std::vector<std::string> axis, std::vector<std::string> units, std::string format){

	// open header file
	openHed();
	cout<<"opened _heds"<<_heds.good()<<endl;
	_heds<<"axis= ";
	for(int i=0;i<ArrND<T>::_dim;i++)_heds<<" "<<axis[i];
	_heds<<endl;

	_heds<<"size=";
	for(int i=0;i<ArrND<T>::_dim;i++)_heds<<" "<<_size[i];
	_heds<<endl;

	_heds<<"delta=";
	for(int i=0;i<ArrND<T>::_dim;i++)_heds<<" "<<delta[i];
	_heds<<endl;

	_heds<<"origin=";
	for(int i=0;i<ArrND<T>::_dim;i++)_heds<<" "<<origin[i];
	_heds<<endl;

	_heds<<"units= ";
	for(int i=0;i<ArrND<T>::_dim;i++)_heds<<" "<<units[i];
	_heds<<endl;
	_heds<<"format= "<<format<<endl;
	_heds<<"data= "<<_filename<<endl;

	closeHed();
}


template<class T> ArrND<T> &ArrND<T>::operator+=(ArrND<T> &arr){
	if(_tfrsze!=arr._tfrsze){
		cerr<<"Error in mem transfer: sizes do not agree"<<endl;
	}

	if (this->_media == 0 && arr._media == 0) {
		T* dest = datptr();
		T* src = arr.datptr();
		for(long i=0;i<_tfrsze;i++)dest[i]+=src[i];
	}

}

template<class T> ArrND<T> &ArrND<T>::operator[](long i){

	_tfrsze/=_size[_tfrdim];
	_tfrloc+=_tfrsze*i;
	_tfrdim++;

	// varsize variables
	_varsize_index=i;

	return *this;
}

template<class T> ArrND<T> &ArrND<T>::operator<<(ArrND<T> &arr){

	if(_tfrsze!=arr._tfrsze){
		cerr<<"Error in mem transfer: sizes do not agree"<<endl;
		cerr<<"source size = "<<_tfrsze<<endl;
		cerr<<"dest size = "<<arr._tfrsze<<endl;
	}

	if(_mode_varsize){
		transferVarSize(arr);
	}
	else{
		transfer(arr);
	}

	_tfrsze=_maxsze;
	_tfrloc=0;
	_tfrdim=0;
	arr._tfrsze=arr._maxsze;
	arr._tfrdim=0;
	arr._tfrloc=0;

	return *this;
}
template<class T> void ArrND<T>::openDat(){

	std::string tmp=_filename;
	ifstream tmps(tmp.c_str());
	if(tmps) {
		try {
			_dats.open(tmp.c_str(), ios::in|ios::out);
		}
		catch(fstream::failure e) {
			cerr<<"Exception in input data file open, file="<<tmp<<endl;
		}
	}
	else {
		try {
			_dats.open(tmp.c_str(), ios::out);
		}
		catch(fstream::failure e) {
			cerr<<"Exception in input data file open, file="<<tmp<<endl;
		}
	}
}

//Same as openDat(), except that it opens the header file
template<class T> void ArrND<T>::openHed(){

	std::string tmp=_filename+".arr";
	ifstream tmps(tmp.c_str());
	if(tmps) {
		try {
			_heds.open(tmp.c_str(), ios::in|ios::out);
		}
		catch(fstream::failure e) {
			cerr<<"Exception in input header file open, file="<<tmp<<endl;
		}
	}
	else {
		try {
			_heds.open(tmp.c_str(), ios::out);
		}
		catch(fstream::failure e) {
			cerr<<"Exception in input header file open, file="<<tmp<<endl;
		}
	}
}

template<class T> void ArrND<T>::closeHed(){
_heds.close();
}

template<class T> void ArrND<T>::readHed(std::vector<long>& size, std::vector<double>& delta, std::vector<double>& origin,
		std::vector<std::string>& axis, std::vector<std::string>& units, std::string& format, int& dim){

    openHed();

	std::string line;
	std::string name;
	std::string value;
	std::stringstream iss;

	regmatch_t m[3];      // Identifies matched substrings
	regex_t    reg;       // Compiled regular expression

	// generate reg expression
	regcomp(&reg, "^\\s*([^=]+)\\s*=\\s*(.*)\\s*$" , REG_EXTENDED);

	list<string> llist;

	while(std::getline(_heds,line)){
		llist.push_back(line);
	}

	// Get number of dims
	for(list<string>::iterator iter=llist.begin();iter!=llist.end();iter++){
		line=(*iter);
		getRegExp(m,reg,line,name,value);
		//cout<<"name="<<name<<endl;
		//cout<<"value="<<value<<endl;
		if(name=="axis"){
			iss.clear();
			iss<<value;
			dim=0;
			string tmp;
			while(iss>>tmp){
				if(name=="axis"){
					dim++;
				}
			}
		}
	}

	//Sanity tests
	if (_dim!=dim){
		cout<<"WARNING:Dimensions do not match"<<endl;
		cout<<"Dimension in object="<<_dim<<endl;
		cout<<"Dimension in file="<<dim<<endl;
	}

	// allocate
	//cout<<"about to realloc"<<endl;
	size.resize(dim);
	delta.resize(dim);
	origin.resize(dim);
	axis.resize(dim);
	units.resize(dim);
	//cout<<"done realloc"<<endl;

	// Fetch params
	//cout<<"fetching params"<<endl;
	for(list<string>::iterator iter=llist.begin();iter!=llist.end();iter++){
		line=(*iter);
		getRegExp(m,reg,line,name,value);
		//cout<<"name="<<name<<endl;
		//cout<<"value="<<value<<endl;
		iss.clear();
		iss<<value;
		int cnt=0;
		string tmp;
		while(iss>>tmp){
			if(name=="axis")axis[cnt]=tmp;
			if(name=="size")size[cnt]=atoi(tmp.c_str());
			if(name=="delta")delta[cnt]=atof(tmp.c_str());
			if(name=="origin")origin[cnt]=atof(tmp.c_str());
			if(name=="units")  units[cnt]=tmp;
			if(name=="format") format=tmp;
			cnt++;
		}
	}

	//Sanity check
	for (int i=0;i<dim;i++) {
		if(size[i]!=_size[i]) {
			cout<<"WARNING:Sizes of dimension "<<i<<" do not match"<<endl;
			cout<<"Size in object="<<_size[i]<<endl;
			cout<<"Size in file="<<size[i]<<endl;
		}
	}

//	cout<< "closing header"<<endl;
	closeHed();
}

template<class T> void ArrND<T>::getRegExp(regmatch_t* m, regex_t& reg, string& line, string& name, string& value){

	if (regexec(&reg, line.c_str(), 3, m, 0) == 0){
		name  = line.substr(m[1].rm_so, m[1].rm_eo - m[1].rm_so);
		value = line.substr(m[2].rm_so, m[2].rm_eo - m[2].rm_so);
	}
}


#endif /* ARRND_H_ */
