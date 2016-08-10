#include <math.h>
#include <stdio.h>
#include <string.h>

#include <SEGY_Header_Structs.h>
#include <SEGY_File.h>
#include <swapbytes.h>

SEGY_Reel_Id_Header_1::SEGY_Reel_Id_Header_1(FILE* fp, bool ascii)
{
	_hdr = new char[3200];
	if (ascii)
	{
		for (int i = 0;  i < 3200;  ++i) _hdr[i] = (char)0;
		_Is_Valid = true;
	}
	else
	{
		int nread = fread((void*)_hdr, 1, 3200, fp);
		_Is_Valid = nread == 3200 ? true : false;
	}
}

SEGY_Reel_Id_Header_1::SEGY_Reel_Id_Header_1()
{
	_hdr = new char[3200];
	for (int i = 0;  i < 3200;  ++i) _hdr[i] = 0;
	_Is_Valid = true;
}

SEGY_Reel_Id_Header_1::~SEGY_Reel_Id_Header_1()
{
	delete [] _hdr;
}

void SEGY_Reel_Id_Header_1::Write(FILE* fp)
{
	fwrite((void*)_hdr, 1, 3200, fp);
}

void SEGY_Reel_Id_Header_1::Copy(SEGY_Reel_Id_Header_1* reel_id_hdr_1)
{
	for (int i = 0;  i < 3200;  ++i) _hdr[i] = reel_id_hdr_1->_hdr[i];
}

SEGY_Reel_Id_Header_2::SEGY_Reel_Id_Header_2(FILE* fp, bool ascii)
{
	_Is_Valid = false;
	_hdr = new char[400];
	if (ascii)
	{
		for (int i = 0;  i < 400;  ++i) _hdr[i] = (char)0;
		_Is_Valid = true;
	}
	else
	{
		int nread = fread((void*)_hdr, 1, 400, fp);
		_Is_Valid = nread == 400 ? true : false;
	}
}

SEGY_Reel_Id_Header_2::SEGY_Reel_Id_Header_2()
{
	_hdr = new char[400];
	for (int i = 0;  i < 400;  ++i) _hdr[i] = 0;
	_Is_Valid = true;
}

SEGY_Reel_Id_Header_2::~SEGY_Reel_Id_Header_2()
{
	delete [] _hdr;
}

void SEGY_Reel_Id_Header_2::Write(FILE* fp)
{
	fwrite((void*)_hdr, 1, 400, fp);
}

void SEGY_Reel_Id_Header_2::Copy(SEGY_Reel_Id_Header_2* reel_id_hdr_2)
{
	for (int i = 0;  i < 400;  ++i) _hdr[i] = reel_id_hdr_2->_hdr[i];
}

SEGY_Trace_Header::SEGY_Trace_Header(FILE* fp, bool ascii)
{
	_Is_Valid = false;
	_hdr = new char[240];
	if (ascii)
	{
		for (int i = 0;  i < 240;  ++i) _hdr[i] = (char)0;
		struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;
		short insamp = 0;
		swap2bytes(&insamp,1);
		hdr->nsamp = insamp;
		double srcx, srcy, srcz, recx, recy, recz;
		if (fscanf(fp, "%lf %lf %lf %lf %lf %lf", &srcx, &srcy, &srcz, &recx, &recy, &recz) == 6)
		{
			_Is_Valid = true;
			int isrcx = (int)round(srcx*100.0);
			int isrcy = (int)round(srcy*100.0);
			int isrcz = (int)round(srcz*100.0);
			int irecx = (int)round(recx*100.0);
			int irecy = (int)round(recy*100.0);
			int irecz = -(int)round(recz*100.0);
			swap4bytes(&isrcx,1);
			swap4bytes(&isrcy,1);
			swap4bytes(&isrcz,1);
			swap4bytes(&irecx,1);
			swap4bytes(&irecy,1);
			swap4bytes(&irecz,1);
			hdr->srcx = isrcx;
			hdr->srcy = isrcy;
			hdr->srcdepth = isrcz;
			hdr->recx = irecx;
			hdr->recy = irecy;
			hdr->recelev = irecz;
		}
	}
	else
	{
		int nread = fread((void*)_hdr, 1, 240, fp);
		_Is_Valid = nread == 240 ? true : false;
	}
}

SEGY_Trace_Header::SEGY_Trace_Header(int nsamp)
{
	_Is_Valid = true;
	_hdr = new char[240];
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;
	short insamp = nsamp;
	swap2bytes(&insamp,1);
        hdr->nsamp = insamp;
}

void SEGY_Trace_Header::Copy(SEGY_Trace_Header* hdr)
{
	for (int i = 0;  i < 240;  ++i) _hdr[i] = hdr->_hdr[i];
}

SEGY_Trace_Header::~SEGY_Trace_Header()
{
	delete [] _hdr;
}

void SEGY_Trace_Header::Print()
{
	printf("Src=[%lf,%lf,%lf], Rec=[%lf,%lf,%lf]\n",Get_Src_X(),Get_Src_Y(),Get_Src_Z(),Get_Rec_X(),Get_Rec_Y(),Get_Rec_Z());
}

short SEGY_Trace_Header::_Get_Short(int bytepos)
{
	short rval = 0;
	char* dst = (char*)&rval;
	char* src = _hdr + bytepos;
	dst[0] = src[1];
	dst[1] = src[0];
	return rval;
}

int SEGY_Trace_Header::_Get_Int(int bytepos)
{
	int rval = 0;
	char* dst = (char*)&rval;
	char* src = _hdr + bytepos;
	dst[0] = src[3];
	dst[1] = src[2];
	dst[2] = src[1];
	dst[3] = src[0];
	return rval;
}

float SEGY_Trace_Header::_Get_Float(int bytepos)
{
	float rval = 0.0f;
	char* dst = (char*)&rval;
	char* src = _hdr + bytepos;
	dst[0] = src[3];
	dst[1] = src[2];
	dst[2] = src[1];
	dst[3] = src[0];
	return rval;
}

double SEGY_Trace_Header::_Get_Double(int bytepos)
{
	double rval = 0.0;
	char* dst = (char*)&rval;
	char* src = _hdr + bytepos;
	dst[7] = src[0];
	dst[6] = src[1];
	dst[5] = src[2];
	dst[4] = src[3];
	dst[3] = src[4];
	dst[2] = src[5];
	dst[1] = src[6];
	dst[0] = src[7];
	return rval;
}

int SEGY_Trace_Header::Get_NSAMP()
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;
	short insamp = hdr->nsamp;
	swap2bytes(&insamp,1);  
	return (int)insamp;
}

int SEGY_Trace_Header::Get_FFID()
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;
	int isrc = hdr->isrc;
	swap4bytes(&isrc,1);
	return isrc;
}

int SEGY_Trace_Header::Get_ILine()
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;
	int iline_no = hdr->iline_no;
	swap4bytes(&iline_no,1);
	return iline_no;
}

void SEGY_Trace_Header::Set_ILine(int IL)
{
	int iline_no = IL;
	swap4bytes(&iline_no,1);
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;
	hdr->iline_no = iline_no;
}

int SEGY_Trace_Header::Get_XLine()
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;
	int xline_no = hdr->xline_no;
	swap4bytes(&xline_no,1);
	return xline_no;
}

void SEGY_Trace_Header::Set_XLine(int XL)
{
	int xline_no = XL;
	swap4bytes(&xline_no,1);
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;
	hdr->xline_no = xline_no;
}

double SEGY_Trace_Header::Get_Sample_Rate()
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;
	short dtmicro = hdr->dtmicro;
	swap2bytes(&dtmicro,1);
	return (double)dtmicro * 1e-6;
}

double SEGY_Trace_Header::Get_Src_X() 
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;  
	int isrcx = hdr->srcx;  
	swap4bytes(&isrcx,1);  
	return (double)isrcx/100.0;
}

double SEGY_Trace_Header::Get_SrcD_X()
{
	struct custom_trc_id_hdr_t* hdr = (struct custom_trc_id_hdr_t*)_hdr;
        double sou_xd = hdr->sou_xd;
        swap8bytes((long*)&sou_xd,1);
	return sou_xd;
}

double SEGY_Trace_Header::Get_Src_Y() 
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;  
	int isrcy = hdr->srcy;  
	swap4bytes(&isrcy,1);  
	return (double)isrcy/100.0;
}

double SEGY_Trace_Header::Get_SrcD_Y()
{
	struct custom_trc_id_hdr_t* hdr = (struct custom_trc_id_hdr_t*)_hdr;
        double sou_yd = hdr->sou_yd;
        swap8bytes((long*)&sou_yd,1);
	return sou_yd;
}

double SEGY_Trace_Header::Get_Src_Z() 
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;  
	int isrcz = hdr->srcdepth;  
	swap4bytes(&isrcz,1);  
	return (double)isrcz/100.0;
}

double SEGY_Trace_Header::Get_Rec_X() 
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;  
	int irecx = hdr->recx;  
	swap4bytes(&irecx,1);  
	return (double)irecx/100.0;
}

double SEGY_Trace_Header::Get_RecD_X()
{
	struct custom_trc_id_hdr_t* hdr = (struct custom_trc_id_hdr_t*)_hdr;
        double rec_xd = hdr->rec_xd;
        swap8bytes((long*)&rec_xd,1);
	return rec_xd;
}

double SEGY_Trace_Header::Get_Rec_Y() 
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;  
	int irecy = hdr->recy;  
	swap4bytes(&irecy,1);  
	return (double)irecy/100.0;
}

double SEGY_Trace_Header::Get_RecD_Y()
{
	struct custom_trc_id_hdr_t* hdr = (struct custom_trc_id_hdr_t*)_hdr;
        double rec_yd = hdr->rec_yd;
        swap8bytes((long*)&rec_yd,1);
	return rec_yd;
}

double SEGY_Trace_Header::Get_Rec_Z() 
{
	struct trc_id_hdr_t* hdr = (struct trc_id_hdr_t*)_hdr;  
	int irecz = hdr->recelev;  
	swap4bytes(&irecz,1);  
	return -(double)irecz/100.0;
}

double SEGY_Trace_Header::Get_Src_H2O()
{
	struct custom_trc_id_hdr_t* hdr = (struct custom_trc_id_hdr_t*)_hdr;
	float d = hdr->src_water_model_depth;
	swap4bytes((int*)&d,1);
	return (double)d;
}

double SEGY_Trace_Header::Get_Rec_H2O()
{
	struct custom_trc_id_hdr_t* hdr = (struct custom_trc_id_hdr_t*)_hdr;
	float d = hdr->rec_water_model_depth;
	swap4bytes((int*)&d,1);
	return (double)d;
}

double SEGY_Trace_Header::Get_Src_Water_Vp()
{
	struct custom_trc_id_hdr_t* hdr = (struct custom_trc_id_hdr_t*)_hdr;
        float d = hdr->src_water_Vp;
	swap4bytes((int*)&d,1);
	return (double)d;
}

double SEGY_Trace_Header::Get_Rec_Water_Vp()
{
	struct custom_trc_id_hdr_t* hdr = (struct custom_trc_id_hdr_t*)_hdr;
	float d = hdr->rec_water_Vp;
	swap4bytes((int*)&d,1);
	return (double)d;
}

int SEGY_Trace_Header::Get_Custom1_SEQ_NO()
{
	return _Get_Int(0);
}

short SEGY_Trace_Header::Get_Custom1_GUN_SEQ()
{
	return _Get_Short(4);
}

short SEGY_Trace_Header::Get_Custom1_COMPON()
{
	return _Get_Short(6);
}

int SEGY_Trace_Header::Get_Custom1_FFID()
{
	return _Get_Int(8);
}

int SEGY_Trace_Header::Get_Custom1_OFFSET()
{
	return _Get_Int(36);
}

int SEGY_Trace_Header::Get_Custom1_RCV_ELEV()
{
	return _Get_Int(40);
}

int SEGY_Trace_Header::Get_Custom1_DEPTH()
{
	return _Get_Int(48);
}

int SEGY_Trace_Header::Get_Custom1_SOU_H2OD()
{
	return _Get_Int(60);
}

int SEGY_Trace_Header::Get_Custom1_REC_H2OD()
{
	return _Get_Int(64);
}

double SEGY_Trace_Header::Get_Custom1_AOFFSET()
{
	return _Get_Double(90);
}

short SEGY_Trace_Header::Get_Custom1_FLAG_VWXYZT()
{
	return _Get_Short(118);
}

double SEGY_Trace_Header::Get_Custom1_SOU_XD()
{
	return _Get_Double(120);
}

double SEGY_Trace_Header::Get_Custom1_SOU_YD()
{
	return _Get_Double(128);
}

double SEGY_Trace_Header::Get_Custom1_REC_XD()
{
	return _Get_Double(136);
}

double SEGY_Trace_Header::Get_Custom1_REC_YD()
{
	return _Get_Double(144);
}

int SEGY_Trace_Header::Get_Custom1_RCV_STAT()
{
	return _Get_Int(152);
}

short SEGY_Trace_Header::Get_Custom1_YEAR()
{
	return _Get_Short(156);
}

short SEGY_Trace_Header::Get_Custom1_DAY_OF_YEAR()
{
	return _Get_Short(158);
}

short SEGY_Trace_Header::Get_Custom1_HOUR_OF_DAY()
{
	return _Get_Short(160);
}

short SEGY_Trace_Header::Get_Custom1_MINUTE_OF_HOUR()
{
	return _Get_Short(162);
}

short SEGY_Trace_Header::Get_Custom1_SECOND_OF_MINUTE()
{
	return _Get_Short(164);
}

int SEGY_Trace_Header::Get_Custom1_USEC_OF_SECOND()
{
	return _Get_Int(168);
}

int SEGY_Trace_Header::Get_Custom1_SOU_LINE()
{
	return _Get_Int(172);
}

int SEGY_Trace_Header::Get_Custom1_SHOT_POINT()
{
	return _Get_Int(176);
}

int SEGY_Trace_Header::Get_Custom1_RCV_LINE()
{
	return _Get_Int(188);
}

int SEGY_Trace_Header::Get_Custom1_RCV_POINT()
{
	return _Get_Int(192);
}

void SEGY_Trace_Header::Write(FILE* fp)
{
	fwrite((void*)_hdr, 1, 240, fp);
}

SEGY_Trace::SEGY_Trace(FILE* fp, bool ascii)
{
	_Is_Valid = false;
	_fp = 0L;
	_samples = 0L;
	_sample_offset = 0L;
	_hdr = new SEGY_Trace_Header(fp,ascii);
	if (_hdr->Is_Valid())
	{
		_fp = fp;
		_sample_offset = ftell(fp);
		if (_hdr->Get_NSAMP() > 0) fseek(fp, _hdr->Get_NSAMP() * sizeof(float), SEEK_CUR);
		_Is_Valid = true;
	}
}

SEGY_Trace::SEGY_Trace(int nsamp)
{
	_Is_Valid = true;
	_fp = 0L;
	_samples = new float[nsamp];
	for (int i = 0;  i < nsamp;  ++i) _samples[i] = 0.0f;
	_hdr = new SEGY_Trace_Header(nsamp);
}

SEGY_Trace::~SEGY_Trace()
{
	delete _hdr;
	delete [] _samples;
}

void SEGY_Trace::Copy(SEGY_Trace* trace)
{
	_hdr->Copy(trace->Get_Trace_Header());
	if (_samples == 0L) _samples = new float[_hdr->Get_NSAMP()];
	for (int i = 0;  i < _hdr->Get_NSAMP();  ++i) _samples[i] = trace->Get_Samples()[i];
}

void SEGY_Trace::Print()
{
	_hdr->Print();
}

float* SEGY_Trace::Get_Samples()
{
	if (_samples == 0L)
	{
		_samples = new float[_hdr->Get_NSAMP()];
		memset((void*)_samples, 0, _hdr->Get_NSAMP()*sizeof(float));
//#pragma omp critical
//{
		fseek(_fp, _sample_offset, SEEK_SET);
		int nread = fread((void*)_samples, sizeof(float), _hdr->Get_NSAMP(), _fp);
		swap4bytes((int*)_samples, nread);
		if (nread != _hdr->Get_NSAMP()) printf("SEGY_Trace::Get_Samples - Error! Read %d samples, expected %d\n",nread,_hdr->Get_NSAMP());
//}
	}
	return _samples;
}

SEGY_File::SEGY_File(const char* path, bool ascii)
{
	_Is_Valid = false;
	_fp = fopen(path, "rb");
	_traces = 0L;
	_max_num_traces = 0;
	_num_traces = 0;
	if (_fp != 0L)
	{
		_Is_Valid = true;
		_reel_id_hdr1 = new class SEGY_Reel_Id_Header_1(_fp,ascii);
		_reel_id_hdr2 = new class SEGY_Reel_Id_Header_2(_fp,ascii);
		bool done = false;
		do
		{
			class SEGY_Trace* trace = new SEGY_Trace(_fp,ascii);
			if (feof(_fp) || !trace->Is_Valid())
			{
				done = true;
			}
			else
			{
				if (_num_traces+1 > _max_num_traces)
				{
					if (_max_num_traces == 0)
					{
						_traces = new SEGY_Trace*[1024];
						_max_num_traces = 1024;
						for (int i = 0;  i < _max_num_traces;  ++i) _traces[i] = 0L;
					}
					else
					{
						_max_num_traces = _max_num_traces * 2;
						class SEGY_Trace** old_traces = _traces;
						_traces = new SEGY_Trace*[_max_num_traces];
						for (int i = 0;  i < _num_traces;  ++i) _traces[i] = old_traces[i];
						for (int i = _num_traces;  i < _max_num_traces;  ++i) _traces[i] = 0L;
						delete [] old_traces;
					}
				}
				_traces[_num_traces] = trace;
				++_num_traces;
			}
		} while (!done);
	}
}

SEGY_File::~SEGY_File()
{
	if (_fp != 0L) fclose(_fp);
	if (_num_traces > 0)
	{
		for (int i = 0;  i < _num_traces;  ++i) delete _traces[i];
		delete [] _traces;
	}
	delete _reel_id_hdr1;
	delete _reel_id_hdr2;
}

void SEGY_File::Print()
{
	for (int i = 0;  i < _num_traces;  ++i) _traces[i]->Print();
}
	
SEGY_File::SEGY_File(int num_traces, int nsamp)
{
	_reel_id_hdr1 = new SEGY_Reel_Id_Header_1();
	_reel_id_hdr2 = new SEGY_Reel_Id_Header_2();
	_num_traces = num_traces;
	_traces = new SEGY_Trace*[_num_traces];
	for (int i = 0;  i < _num_traces;  ++i) _traces[i] = new SEGY_Trace(nsamp);
}

void SEGY_Trace::Write(FILE* fp)
{
	_hdr->Write(fp);
	int nsamp = _hdr->Get_NSAMP();
	float* src = Get_Samples();
	float* dst = new float[nsamp];
	for (int i = 0;  i < nsamp;  ++i) dst[i] = src[i];
	swap4bytes((int*)dst,nsamp);
	fwrite(dst,sizeof(float),nsamp,fp);
	delete [] dst;
}

void SEGY_File::Write(const char* path)
{
	FILE* fp = fopen(path, "wb");
	if (fp != 0L)
	{
		_reel_id_hdr1->Write(fp);
		_reel_id_hdr2->Write(fp);
		for (int i = 0;  i < _num_traces;  ++i) _traces[i]->Write(fp);
		fclose(fp);
	}
	else
	{
		printf("Error! Unable to open %s for writing.\n",path);
	}
}

