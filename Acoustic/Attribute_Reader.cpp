#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Attribute_Reader.hxx"

void Attribute_Reader::_open(const char* filename, float const_val, char*& dupfilename, FILE*& fp, float& val_Constant)
{
	printf("Attribute_Reader::_open(filename=%s, const_val=%f)\n",filename!=0L?filename:"nil",const_val);
	val_Constant = const_val;
	dupfilename = _strdup(filename);
	if (dupfilename != 0L)
	{
		fp = fopen(dupfilename, "r");
		if (fp != 0L)
			printf("Success : opened %s for reading\n",dupfilename);
		else
			printf("FAILED to open %s for reading\n",dupfilename);
	}
	else
	{
		fp = 0L;
	}
}

char* Attribute_Reader::_strdup(const char* str)
{
	if (str != 0L)
	{
		int len = strlen(str);
		if (len > 0)
		{
			char* dupstr = new char[len+1];
			memcpy((void*)dupstr, (void*)str, len);
			dupstr[len] = (char)0;
			return dupstr;
		}
		else
		{
			return 0L;
		}
	}
	else
	{
		return 0L;
	}
}

void Attribute_Reader::_Read(FILE* fp, long offset, int n, float* inpbuf, float const_val)
{
	if (fp != 0L)
	{
		//printf("Attempting to read %d floats from offset %ld\n",n,offset);
		fseek(fp, offset, SEEK_SET);
		fread(inpbuf, sizeof(float), n, fp);
		if(_swap_Flag) swap4bytes((int*)inpbuf, n);
	}
	else
	{
		for (int i = 0;  i < n;  ++i) inpbuf[i] = const_val;
	}
	for (int i = 0;  i < n;  ++i) inpbuf[i] = Transform(inpbuf[i]);
}

/***** swap4bytes:  dcba --> abcd *****/
void Attribute_Reader::swap4bytes(int *i4, int n)
{
	const int bmask = 16711680, cmask = 65280, dmask = 255;
	for(int k=0;  k<n;  ++k)
	{ 
		int i = i4[k];
		int a =  i << 24;          
		int b = (i << 8)  & bmask;
		int c = (i >> 8) & cmask;  
		int d = (i >> 24) & dmask;
		i4[k] = a | b | c | d ;
	}
}

Attribute_Loader::Attribute_Loader(
		const char* arg_name,
		int swap_Flag,
		Index_Calculator* icalc,
		UVW_XYZ_Calculator* ucalc
		)
		: Attribute_Reader(arg_name,swap_Flag)
{
	_icalc = icalc;
	_ucalc = ucalc;
}

Attribute_Loader::~Attribute_Loader()
{
}

//
// compress value into bits.
//
int 
Attribute_Loader::Compress_Attribute(float val, float min, float range, int mask, int shift)
{
	int val_idx = (int)(((val - min) * range * (float)mask) + 0.5f);
        if (val_idx < 0) val_idx = 0;
        if (val_idx > mask) val_idx = mask;
        if (shift > 0) val_idx = val_idx << shift;
	//printf("val_idx=%d, val=%e, min=%e, range=%e, mask=%d, shift=%d\n",val_idx,val,min,range,mask,shift);
        return val_idx;
}

//
// decompress bits into value.
float 
Attribute_Loader::Decompress_Attribute(int attr, float min, float binsize, int mask, int shift)
{
	int ival = shift > 0 ? ((attr >> shift) & mask) : (attr & mask);
	float val = (float)ival * binsize + min;
	return val;
}

//
// u,v,w coordinates in input file maps to this index in compressed earth model.
//
unsigned long Attribute_Loader::_compressedIDX(
		int u,
		int v,
		int w
		)
{
	int iX,iY,iZ;
	_ucalc->Compute_XYZ_From_UVW(u,v,w,iX,iY,iZ);
	return _icalc->Calculate_Index(iX,iY,iZ);
}

Single_Attribute_Loader::Single_Attribute_Loader(
		int* attr,
		Index_Calculator* icalc,
		UVW_XYZ_Calculator* ucalc,
		const char* name,
		const char* filename,
		float const_val,
		float min,
		float range,
		int swap_Flag,
		int prop_mask,
		int prop_shift
		)
		: Attribute_Loader(name,swap_Flag,icalc,ucalc)
{
	_attr = attr;		// no duplication.
	_open(filename, const_val, _filename, _fp, _val_Constant);
	_inpbuf = 0L;
        _inpbuflen = 0;
	_min = min;
	_range = range;
	_prop_mask = prop_mask;
	_prop_shift = prop_shift;
}

Single_Attribute_Loader::~Single_Attribute_Loader()
{
	if (_filename != 0L) delete [] _filename;
	if (_inpbuf != 0L) delete [] _inpbuf;
	if (_fp != 0L) fclose(_fp);
	// name handled by super class
	// attr, icalc, ucalc are managed by calling class and should not be deallocated.
}

void Single_Attribute_Loader::Read(long offset, int ubeg, int uend, int v, int w)
{
	int actual_dimu = uend - ubeg + 1;
	if (actual_dimu > _inpbuflen)
	{
		if (_inpbuf != 0L) delete [] _inpbuf;
		_inpbuf = new float[actual_dimu];
		_inpbuflen = actual_dimu;
	}

	_Read(_fp, offset, actual_dimu, _inpbuf, _val_Constant);
	for (int i = 0;  i < actual_dimu;  ++i)
	{
		int attr_val = Compress_Attribute(_inpbuf[i],_min,_range,_prop_mask,_prop_shift);
		unsigned long idx = _compressedIDX(ubeg+i,v,w);
		_attr[idx] |= attr_val;
		//if (attr_val != 0) printf("attr[%ld] |= %d\n",idx,attr_val);
	}
}

Vp_Loader::Vp_Loader(
		int* attr,
		Index_Calculator* icalc,
		UVW_XYZ_Calculator* ucalc,
		const char* name,
		const char* filename,
		float const_val,
		float min,
		float range,
		int swap_Flag,
		int prop_mask,
		int prop_shift,
		float dt
		)
		: Single_Attribute_Loader(attr,icalc,ucalc,name,filename,const_val,min,range,swap_Flag,prop_mask,prop_shift)
{
	_dt = dt;
}

float Vp_Loader::Transform(float val)
{
	float rval = val * val * _dt * _dt;
	//printf("Vp_Loader::Transform - val = %e, rval = %e\n",val,rval);
	return rval;
}

Q_Loader::Q_Loader(
		int* attr,
		Index_Calculator* icalc,
		UVW_XYZ_Calculator* ucalc,
		const char* name,
		const char* filename,
		float const_val,
		float min,
		float range,
		int swap_Flag,
		int prop_mask,
		int prop_shift,
		float fq,
		float dt
		)
		: Single_Attribute_Loader(attr,icalc,ucalc,name,filename,const_val,min,range,swap_Flag,prop_mask,prop_shift)
{
	_trans = new Transform_Q(fq,dt);
}

Q_Loader::~Q_Loader()
{
	if (_trans != 0L) delete _trans;
}

float Q_Loader::Transform(float val)
{
	float inv_Q;
	_trans->Transform(val,inv_Q);
	return inv_Q;
}

Eps_Dta_Loader::Eps_Dta_Loader(
		int* attr,
		Index_Calculator* icalc,
		UVW_XYZ_Calculator* ucalc,
		const char* name,
		const char* eps_or_eta_File,
		float eps_or_eta_Constant,
		int eta_Flag,
		float eps_min,
		float eps_range,
		int eps_mask,
		int eps_shift,
		const char* dta_File,
		float dta_Constant,
		float dta_min,
		float dta_range,
		int dta_mask,
		int dta_shift,
		int swap_Flag
		)
		: Attribute_Loader(name, swap_Flag, icalc, ucalc)
{
	_attr = attr;
	_open(dta_File, dta_Constant, _dta_File, _dta_fp, _dta_Constant);
	_open(eps_or_eta_File, eps_or_eta_Constant, _eps_or_eta_File, _eps_or_eta_fp, _eps_or_eta_Constant);
	_trans = new Transform_Eps(eta_Flag);
	_inpbuf_dta = 0L;
	_inpbuf_eps_or_eta = 0L;
	_inpbuflen = 0;
	_eps_min = eps_min;
	_eps_range = eps_range;
	_eps_mask = eps_mask;
	_eps_shift = eps_shift;
	_dta_min = dta_min;
	_dta_range = dta_range;
	_dta_mask = dta_mask;
	_dta_shift = dta_shift;
}

Eps_Dta_Loader::~Eps_Dta_Loader()
{
	if (_dta_File != 0L) delete [] _dta_File;
	if (_eps_or_eta_File != 0L) delete [] _eps_or_eta_File;
	if (_dta_fp != 0L) fclose(_dta_fp);
	if (_eps_or_eta_fp != 0L) fclose(_eps_or_eta_fp);
	if (_inpbuf_dta != 0L) delete [] _inpbuf_dta;
	if (_inpbuf_eps_or_eta != 0L) delete [] _inpbuf_eps_or_eta;
	if (_trans != 0L) delete _trans;
}

void Eps_Dta_Loader::Read(long offset, int ubeg, int uend, int v, int w)
{
	int actual_dimu = uend - ubeg + 1;
	if (actual_dimu > 0)
	{
		if (actual_dimu > _inpbuflen)
		{
			if (_inpbuf_dta != 0L) delete [] _inpbuf_dta;
			_inpbuf_dta = new float[actual_dimu];
			if (_inpbuf_eps_or_eta != 0L) delete [] _inpbuf_eps_or_eta;
			_inpbuf_eps_or_eta = new float[actual_dimu];
			_inpbuflen = actual_dimu;
		}

		_Read(_dta_fp, offset, actual_dimu, _inpbuf_dta, _dta_Constant);
		_Read(_eps_or_eta_fp, offset, actual_dimu, _inpbuf_eps_or_eta, _eps_or_eta_Constant);
		for (int i = 0;  i < actual_dimu;  ++i)
		{
			float dta = _inpbuf_dta[i];

			float eps;
			_trans->Transform(_inpbuf_eps_or_eta[i],dta,eps);

			int dta_idx = Compress_Attribute(dta,_dta_min,_dta_range,_dta_mask,_dta_shift);
			int eps_idx = Compress_Attribute(eps,_eps_min,_eps_range,_eps_mask,_eps_shift);

			_attr[_compressedIDX(ubeg+i,v,w)] |= (dta_idx | eps_idx);
		}
	}
}

Dip_Azm_Loader::Dip_Azm_Loader(
		int* attr,
		Index_Calculator* icalc,
		UVW_XYZ_Calculator* ucalc,
		const char* name,
		const char* dip_or_dx_File,
		float dip_or_dx_Constant,
		float dip_min,
		float dip_range,
		int dip_mask,
		int dip_shift,
		const char* azm_or_dy_File,
		float azm_or_dy_Constant,
		float azm_min,
		float azm_range,
		int azm_mask,
		int azm_shift,
		int swap_Flag,
		int degrees_Flag,
		int dipxdipy_Flag
		)
: Attribute_Loader(name, swap_Flag, icalc, ucalc)
{
	_attr = attr;
	_open(dip_or_dx_File, dip_or_dx_Constant, _dip_or_dx_File, _dip_or_dx_fp, _dip_or_dx_Constant);
	_open(azm_or_dy_File, azm_or_dy_Constant, _azm_or_dy_File, _azm_or_dy_fp, _azm_or_dy_Constant);
	_inpbuf_dip_or_dx = 0L;
	_inpbuf_azm_or_dy = 0L;
	_inpbuflen = 0;
	_trans = new Transform_Dip_Azm(degrees_Flag,dipxdipy_Flag);
	_dip_min = dip_min;
	_dip_range = dip_range;
	_dip_mask = dip_mask;
	_dip_shift = dip_shift;
	_azm_min = azm_min;
	_azm_range = azm_range;
	_azm_mask = azm_mask;
	_azm_shift = azm_shift;
}

Dip_Azm_Loader::~Dip_Azm_Loader()
{
	if (_dip_or_dx_File != 0L) delete [] _dip_or_dx_File;
	if (_azm_or_dy_File != 0L) delete [] _azm_or_dy_File;
	if (_dip_or_dx_fp != 0L) fclose(_dip_or_dx_fp);
	if (_azm_or_dy_fp != 0L) fclose(_azm_or_dy_fp);
	if (_inpbuf_dip_or_dx != 0L) delete [] _inpbuf_dip_or_dx;
	if (_inpbuf_azm_or_dy != 0L) delete [] _inpbuf_azm_or_dy;
	if (_trans != 0L) delete _trans;
}

void Dip_Azm_Loader::Read(long offset, int ubeg, int uend, int v, int w)
{
	int actual_dimu = uend - ubeg + 1;
	if (actual_dimu > 0)
	{
		if (actual_dimu > _inpbuflen)
		{
			if (_inpbuf_dip_or_dx != 0L) delete [] _inpbuf_dip_or_dx;
			_inpbuf_dip_or_dx = new float[actual_dimu];
			if (_inpbuf_azm_or_dy != 0L) delete [] _inpbuf_azm_or_dy;
			_inpbuf_azm_or_dy = new float[actual_dimu];
			_inpbuflen = actual_dimu;
		}

		_Read(_dip_or_dx_fp, offset, actual_dimu, _inpbuf_dip_or_dx, _dip_or_dx_Constant);
		_Read(_azm_or_dy_fp, offset, actual_dimu, _inpbuf_azm_or_dy, _azm_or_dy_Constant);
		for (int i = 0;  i < actual_dimu;  ++i)
		{
			float dip, azm;
			_trans->Transform(_inpbuf_dip_or_dx[i],_inpbuf_azm_or_dy[i],dip,azm);
			int dip_idx = Compress_Attribute(dip,_dip_min,_dip_range,_dip_mask,_dip_shift);
			int azm_idx = Compress_Attribute(azm,_azm_min,_azm_range,_azm_mask,_azm_shift);
			_attr[_compressedIDX(ubeg+i,v,w)] |= (dip_idx | azm_idx);
		}
	}
}

Single_Attribute_Scanner::Single_Attribute_Scanner(
		const char* name,
		const char* filename,
		float const_val,
		int swap_Flag
		)
		: Attribute_Reader(name,swap_Flag)
{
	_open(filename, const_val, _filename, _fp, _val_Constant);
	_inpbuf = 0L;
	_inpbuflen = 0;
	_min = 1e37f;
	_max = -1e37f;
}

Single_Attribute_Scanner::~Single_Attribute_Scanner()
{
	if (_filename != 0L) delete [] _filename;
	if (_fp != 0L) fclose(_fp);
	if (_inpbuf != 0L) delete [] _inpbuf;
}

float Single_Attribute_Scanner::Get_Min()
{
	return _min;
}

float Single_Attribute_Scanner::Get_Max()
{
	return _max;
}

void Single_Attribute_Scanner::Read(long offset, int ubeg, int uend, int v, int w)
{
	int actual_dimu = uend - ubeg + 1;
	if (actual_dimu > _inpbuflen)
	{
		if (_inpbuf != 0L) delete [] _inpbuf;
		_inpbuf = new float[actual_dimu];
		_inpbuflen = actual_dimu;
	}

	_Read(_fp, offset, actual_dimu, _inpbuf, _val_Constant);
	for (int i = 0;  i < actual_dimu;  ++i)
	{
		float val = _inpbuf[i];
		if (val > _max) _max = val;
		if (val < _min) _min = val;
	}
}

Q_Scanner::Q_Scanner(
		const char* name,
		const char* filename,
		float const_val,
		int swap_Flag,
		float fq,
		float dt
		)
: Single_Attribute_Scanner(name,filename,const_val,swap_Flag)
{
	_trans = new Transform_Q(fq,dt);
}

Q_Scanner::~Q_Scanner()
{
	if (_trans) delete _trans;
}

float Q_Scanner::Transform(float Q)
{	
	float inv_Q;
	_trans->Transform(Q, inv_Q);
	return inv_Q;
}

Dip_Azm_Scanner::Dip_Azm_Scanner(
		const char* name,
		const char* dip_or_dx_File,
		const char* azm_or_dy_File,
		float dip_or_dx_Constant,
		float azm_or_dy_Constant,
		int swap_Flag,
		int degrees_Flag,
		int dipxdipy_Flag
		)
: Attribute_Reader(name,swap_Flag)
{
	_open(dip_or_dx_File, dip_or_dx_Constant, _dip_or_dx_File, _dip_or_dx_fp, _dip_or_dx_Constant);
	_open(azm_or_dy_File, azm_or_dy_Constant, _azm_or_dy_File, _azm_or_dy_fp, _azm_or_dy_Constant);
	_inpbuf_dip_or_dx = 0L;
	_inpbuf_azm_or_dy = 0L;
	_inpbuflen = 0;
	_dip_Min = 1e37f;
	_dip_Max = -1e37f;
	_azm_Min = 1e37f;
	_azm_Max = -1e37f;
	_trans = new Transform_Dip_Azm(degrees_Flag,dipxdipy_Flag);
}

Dip_Azm_Scanner::~Dip_Azm_Scanner()
{
	if (_dip_or_dx_File != 0L) delete [] _dip_or_dx_File;
	if (_azm_or_dy_File != 0L) delete [] _azm_or_dy_File;
	if (_dip_or_dx_fp != 0L) fclose(_dip_or_dx_fp);
	if (_azm_or_dy_fp != 0L) fclose(_azm_or_dy_fp);
	if (_inpbuf_dip_or_dx != 0L) delete [] _inpbuf_dip_or_dx;
	if (_inpbuf_azm_or_dy != 0L) delete [] _inpbuf_azm_or_dy;
	if (_trans != 0L) delete _trans;
}

float Dip_Azm_Scanner::Get_Dip_Min()
{
	return _dip_Min;
}

float Dip_Azm_Scanner::Get_Dip_Max()
{
	return _dip_Max;
}

float Dip_Azm_Scanner::Get_Azm_Min()
{
	return _azm_Min;
}

float Dip_Azm_Scanner::Get_Azm_Max()
{
	return _azm_Max;
}

void Dip_Azm_Scanner::Read(long offset, int ubeg, int uend, int v, int w)
{
	int actual_dimu = uend - ubeg + 1;
	if (actual_dimu > _inpbuflen)
	{
		if (_inpbuf_dip_or_dx != 0L) delete [] _inpbuf_dip_or_dx;
		_inpbuf_dip_or_dx = new float[actual_dimu];
		if (_inpbuf_azm_or_dy != 0L) delete [] _inpbuf_azm_or_dy;
		_inpbuf_azm_or_dy = new float[actual_dimu];
		_inpbuflen = actual_dimu;
	}

	_Read(_dip_or_dx_fp, offset, actual_dimu, _inpbuf_dip_or_dx, _dip_or_dx_Constant);
	_Read(_azm_or_dy_fp, offset, actual_dimu, _inpbuf_azm_or_dy, _azm_or_dy_Constant);
	for (int i = 0;  i < actual_dimu;  ++i)
	{
		float dip, azm;
		_trans->Transform(_inpbuf_dip_or_dx[i],_inpbuf_azm_or_dy[i],dip,azm);
		if (dip > _dip_Max) _dip_Max = dip;
		if (dip < _dip_Min) _dip_Min = dip;
		if (azm > _azm_Max) _azm_Max = azm;
		if (azm < _azm_Min) _azm_Min = azm;
	}
}

Vp_Eps_Dta_Scanner::Vp_Eps_Dta_Scanner(
		const char* name,
		const char* vp_File,
		const char* dta_File,
		const char* eps_or_eta_File,
		float vp_Constant,
		float dta_Constant,
		float eps_or_eta_Constant,
		int eta_Flag,
		int swap_Flag
		)
: Attribute_Reader(name,swap_Flag)
{
	_open(vp_File, vp_Constant, _vp_File, _vp_fp, _vp_Constant);
	_open(dta_File, dta_Constant, _dta_File, _dta_fp, _dta_Constant);
	_open(eps_or_eta_File, eps_or_eta_Constant, _eps_or_eta_File, _eps_or_eta_fp, _eps_or_eta_Constant);
	_trans = new Transform_Eps(eta_Flag);
	_inpbuf_vp = 0L;
	_inpbuf_dta = 0L;
	_inpbuf_eps_or_eta = 0L;
	_inpbuflen = 0;
	_vp_X_min = 1e37f;
	_vp_X_max = -1e37f;
	_vp_Z_min = 1e37f;
	_vp_Z_max = -1e37f;
	_dta_min = 1e37f;
	_dta_max = -1e37f;
	_eps_min = 1e37f;
	_eps_max = -1e37f;
}

Vp_Eps_Dta_Scanner::~Vp_Eps_Dta_Scanner()
{
	if (_vp_File != 0L) delete [] _vp_File;
	if (_dta_File != 0L) delete [] _dta_File;
	if (_eps_or_eta_File != 0L) delete [] _eps_or_eta_File;
	if (_vp_fp != 0L) fclose(_vp_fp);
	if (_dta_fp != 0L) fclose(_dta_fp);
	if (_eps_or_eta_fp != 0L) fclose(_eps_or_eta_fp);
	if (_inpbuf_vp != 0L) delete [] _inpbuf_vp;
	if (_inpbuf_dta != 0L) delete [] _inpbuf_dta;
	if (_inpbuf_eps_or_eta != 0L) delete [] _inpbuf_eps_or_eta;
	if (_trans != 0L) delete _trans;
}

float Vp_Eps_Dta_Scanner::Get_vp_X_Min()
{
	return _vp_X_min;
}

float Vp_Eps_Dta_Scanner::Get_vp_X_Max()
{
	return _vp_X_max;
}

float Vp_Eps_Dta_Scanner::Get_vp_Z_Min()
{
	return _vp_Z_min;
}

float Vp_Eps_Dta_Scanner::Get_vp_Z_Max()
{
	return _vp_Z_max;
}

float Vp_Eps_Dta_Scanner::Get_Eps_Min()
{
	return _eps_min;
}

float Vp_Eps_Dta_Scanner::Get_Eps_Max()
{
	return _eps_max;
}

float Vp_Eps_Dta_Scanner::Get_Dta_Min()
{
	return _dta_min;
}

float Vp_Eps_Dta_Scanner::Get_Dta_Max()
{
	return _dta_max;
}

void Vp_Eps_Dta_Scanner::Read(long offset, int ubeg, int uend, int v, int w)
{
	int actual_dimu = uend - ubeg + 1;
	if (actual_dimu > _inpbuflen)
	{
		if (_inpbuf_vp != 0L) delete [] _inpbuf_vp;
		_inpbuf_vp = new float[actual_dimu];
		if (_inpbuf_dta != 0L) delete [] _inpbuf_dta;
		_inpbuf_dta = new float[actual_dimu];
		if (_inpbuf_eps_or_eta != 0L) delete [] _inpbuf_eps_or_eta;
		_inpbuf_eps_or_eta = new float[actual_dimu];
		_inpbuflen = actual_dimu;
	}

	//printf("Vp_Eps_Dta_Scanner::Read(offset=%ld, ubeg=%d, uend=%d, v=%d, w=%d)\n",offset,ubeg,uend,v,w);

	_Read(_vp_fp, offset, actual_dimu, _inpbuf_vp, _vp_Constant);
	_Read(_dta_fp, offset, actual_dimu, _inpbuf_dta, _dta_Constant);
	_Read(_eps_or_eta_fp, offset, actual_dimu, _inpbuf_eps_or_eta, _eps_or_eta_Constant);
	for (int i = 0;  i < actual_dimu;  ++i)
	{
		float vp_Z = _inpbuf_vp[i];
		float dta = _inpbuf_dta[i];

		float eps;
		_trans->Transform(_inpbuf_eps_or_eta[i],dta,eps);

		//printf("vp_Z = %f, dta = %f, eps = %f\n",vp_Z,dta,eps);

		float vp_X = sqrtf(1.0f + 2.0f * eps) * vp_Z;

		if (vp_Z > _vp_Z_max) _vp_Z_max = vp_Z;
		if (vp_Z < _vp_Z_min) _vp_Z_min = vp_Z;
		if (vp_X > _vp_X_max) _vp_X_max = vp_X;
		if (vp_X < _vp_X_min) _vp_X_min = vp_X;
		if (eps > _eps_max) _eps_max = eps;
		if (eps < _eps_min) _eps_min = eps;
		if (dta > _dta_max) _dta_max = dta;
		if (dta < _dta_min) _dta_min = dta;
	}
}
