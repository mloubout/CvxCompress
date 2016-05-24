#ifndef SEGY_FILE_H
#define SEGY_FILE_H

#include <stdio.h>
#include <string.h>

#include <SEGY_Header_Structs.h>

class SEGY_Reel_Id_Header_1
{
public:
	SEGY_Reel_Id_Header_1(FILE* fp, bool ascii);
	SEGY_Reel_Id_Header_1();
	~SEGY_Reel_Id_Header_1();

	bool Is_Valid() {return _Is_Valid;}
	void Copy(SEGY_Reel_Id_Header_1* reel_id_hdr_1);

	void Write(FILE* fp);

private:
	bool _Is_Valid;
	char* _hdr;
};

class SEGY_Reel_Id_Header_2
{
public:
	SEGY_Reel_Id_Header_2(FILE* fp, bool ascii);
	SEGY_Reel_Id_Header_2();
	~SEGY_Reel_Id_Header_2();

	bool Is_Valid() {return _Is_Valid;}
	void Copy(SEGY_Reel_Id_Header_2* reel_id_hdr_2);
	
	void Write(FILE* fp);

private:
	bool _Is_Valid;
	char* _hdr;
};

class SEGY_Trace_Header
{
public:
	SEGY_Trace_Header(FILE* fp, bool ascii);
	SEGY_Trace_Header(int nsamp);
	~SEGY_Trace_Header();

	bool Is_Valid() {return _Is_Valid;}
	void Copy(SEGY_Trace_Header* hdr);

	void Print();

	int Get_NSAMP();
	double Get_Sample_Rate();
	
	int Get_FFID();

	double Get_Src_X();
	double Get_SrcD_X();
	double Get_Src_Y();
	double Get_SrcD_Y();
	double Get_Src_Z();

	double Get_Rec_X();
	double Get_RecD_X();
	double Get_Rec_Y(); 
	double Get_RecD_Y();
	double Get_Rec_Z();

	double Get_Src_H2O();
	double Get_Rec_H2O();

	double Get_Src_Water_Vp();
	double Get_Rec_Water_Vp();

	int Get_ILine();
	void Set_ILine(int IL);

	int Get_XLine();
	void Set_XLine(int XL);

        void Write(FILE* fp);

private:
	bool _Is_Valid;
	char* _hdr;
};

class SEGY_Trace
{
public:
	SEGY_Trace(FILE* fp, bool ascii);
	SEGY_Trace(int nsamp);
	~SEGY_Trace();

	bool Is_Valid() {return _Is_Valid;}	
	void Copy(SEGY_Trace* trace);

	void Print();
	
	class SEGY_Trace_Header* Get_Trace_Header() {return _hdr;}
	float* Get_Samples();

	void Write(FILE* fp);

private:
	bool _Is_Valid;
	FILE* _fp;
	
	class SEGY_Trace_Header* _hdr;

	size_t _sample_offset;
	float* _samples;
};

class SEGY_File
{
public:
	SEGY_File(const char* path, bool ascii);
	SEGY_File(int num_traces, int nsamp);
	~SEGY_File();
	
	void Print();
	bool Is_Valid() {return _Is_Valid;}	

	SEGY_Reel_Id_Header_1* Get_Reel_Id_Header_1() {return _reel_id_hdr1;}
	SEGY_Reel_Id_Header_2* Get_Reel_Id_Header_2() {return _reel_id_hdr2;}

	int Get_Number_Of_Traces() {return _num_traces;}

	SEGY_Trace* Get_Trace(int idx)
	{
		if (idx >= 0 && idx < _num_traces)
		{
			return _traces[idx];
		}
		else
		{
			return 0L;
		}
	}

	void Write(const char* path);

private:
	bool _Is_Valid;
	FILE* _fp;
	
	class SEGY_Reel_Id_Header_1* _reel_id_hdr1;
	class SEGY_Reel_Id_Header_2* _reel_id_hdr2;

	class SEGY_Trace** _traces;
	int _max_num_traces;
	int _num_traces;
};

#endif

