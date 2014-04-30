#ifndef ATTRIBUTE_READER_HXX
#define ATTRIBUTE_READER_HXX

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Index_Calculator.hxx"
#include "UVW_XYZ_Calculator.hxx"
#include "Attribute_Transformer.hxx"

class Attribute_Reader
{
	public:
		Attribute_Reader(const char* arg_name, int swap_Flag) {_swap_Flag = swap_Flag; name = _strdup(arg_name);}
		virtual ~Attribute_Reader() {if (name != 0L) delete [] name;}

		virtual float Transform(float val) {return val;}				// transform value, useful for Q, dip, azm etc.
		virtual void Read(long offset, int ubeg, int uend, int v, int w) = 0;		// read one stripe [ubeg,uend] at v,w

		char* name;

	protected:
		int _swap_Flag;

		void _open(const char* filename, float const_val, char*& dupfilename, FILE*& fp, float& val_Constant);
		char* _strdup(const char* str);
		void _Read(FILE* fp, long offset, int n, float* inpbuf, float const_val);
		void swap4bytes(int *i4, int n);
};

class Attribute_Loader : public Attribute_Reader
{
	public:
		Attribute_Loader(
				const char* name,
				int swap_Flag,
				Index_Calculator* icalc,
                                UVW_XYZ_Calculator* ucalc
				);
		virtual ~Attribute_Loader();

		static int Compress_Attribute(float val, float min, float range, int mask, int shift);
		static float Decompress_Attribute(int attr, float min, float binsize, int mask, int shift);

	protected:
		Index_Calculator* _icalc;
		UVW_XYZ_Calculator* _ucalc;

		unsigned long _compressedIDX(int u, int v, int w);
};

class Single_Attribute_Loader : public Attribute_Loader
{
	public:
		Single_Attribute_Loader(
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
				);

		virtual ~Single_Attribute_Loader();

		void Read(long offset, int ubeg, int uend, int v, int w);

	protected:
		int* _attr;
		Index_Calculator* _icalc;
		UVW_XYZ_Calculator* _ucalc;

		char* _filename;
		FILE* _fp;
		float _val_Constant;

		float _min;
		float _range;

		int _prop_mask;
		int _prop_shift;

                int _inpbuflen;
                float* _inpbuf;
};

class Vp_Loader : public Single_Attribute_Loader
{
	public:
		Vp_Loader(
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
				);

		float Transform(float val);

	private:
		float _dt;
};

class Q_Loader : public Single_Attribute_Loader
{
	public:
		Q_Loader(
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
				);
		virtual ~Q_Loader();

		float Transform(float val);

	private:
		Transform_Q* _trans;
};

class Eps_Dta_Loader : public Attribute_Loader
{
	public:
		Eps_Dta_Loader(
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
				);
		virtual ~Eps_Dta_Loader();

		void Read(long offset, int ubeg, int uend, int v, int w);

	private:
		int* _attr;
		Index_Calculator* icalc;
		UVW_XYZ_Calculator* ucalc;
		Transform_Eps* _trans;

		float _eps_min;
		float _eps_range;
		int _eps_mask;
		int _eps_shift;

		float _dta_min;
		float _dta_range;
		int _dta_mask;
		int _dta_shift;

		char* _dta_File;
		char* _eps_or_eta_File;

		float _dta_Constant;
		float _eps_or_eta_Constant;

		FILE* _dta_fp;
		FILE* _eps_or_eta_fp;

		int _inpbuflen;
		float* _inpbuf_dta;
		float* _inpbuf_eps_or_eta;
};

class Dip_Azm_Loader : public Attribute_Loader
{
	public:
		Dip_Azm_Loader(
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
				);
		virtual ~Dip_Azm_Loader();

		void Read(long offset, int ubeg, int uend, int v, int w);

	private:
		int* _attr;
		Index_Calculator* icalc;
		UVW_XYZ_Calculator* ucalc;
		Transform_Dip_Azm* _trans;

		float _dip_min;
		float _dip_range;
		int _dip_mask;
		int _dip_shift;

		float _azm_min;
		float _azm_range;
		int _azm_mask;
		int _azm_shift;

		char* _dip_or_dx_File;
		char* _azm_or_dy_File;

		FILE* _dip_or_dx_fp;
		FILE* _azm_or_dy_fp;

		float _dip_or_dx_Constant;
		float _azm_or_dy_Constant;

		int _inpbuflen;
		float* _inpbuf_dip_or_dx;
		float* _inpbuf_azm_or_dy;
};

class Single_Attribute_Scanner : public Attribute_Reader
{
	public:
		Single_Attribute_Scanner(
				const char* name,
				const char* filename,
				float const_val,
				int swap_Flag
				);

		~Single_Attribute_Scanner();

		float Get_Min();
		float Get_Max();
		void Read(long offset, int ubeg, int uend, int v, int w);

	protected:
		char* _filename;
		FILE* _fp;
		float _val_Constant;

		int _inpbuflen;
		float* _inpbuf;

		float _min;
		float _max;		
};

class Q_Scanner : public Single_Attribute_Scanner
{
	public:
		Q_Scanner(
				const char* name,
				const char* filename,
				float const_val,
				int swap_Flag,
				float fq,
				float dt
				);

		~Q_Scanner();

		float Transform(float Q);

	private:
		Transform_Q* _trans;
};

class Dip_Azm_Scanner : public Attribute_Reader
{
	public:
		Dip_Azm_Scanner(
				const char* name,
				const char* dip_or_dx_File,
				const char* azm_or_dy_File,
				float dip_or_dx_Constant,
				float azm_or_dy_Constant,
				int swap_Flag,
				int degrees_Flag,
				int dipxdipy_Flag
				);

		~Dip_Azm_Scanner();

		float Get_Dip_Min();
		float Get_Dip_Max();
		float Get_Azm_Min();
		float Get_Azm_Max();

		void Read(long offset, int ubeg, int uend, int v, int w);

	private:
		Transform_Dip_Azm* _trans;

		char* _dip_or_dx_File;
		char* _azm_or_dy_File;

		FILE* _dip_or_dx_fp;
		FILE* _azm_or_dy_fp;

		float _dip_or_dx_Constant;
		float _azm_or_dy_Constant;

		int _inpbuflen;
		float* _inpbuf_dip_or_dx;
		float* _inpbuf_azm_or_dy;

		float _dip_Min;
		float _dip_Max;
		float _azm_Min;
		float _azm_Max;
};

class Vp_Eps_Dta_Scanner : public Attribute_Reader
{
	public:
		Vp_Eps_Dta_Scanner(
				const char* name,
				const char* vp_File,
				const char* dta_File,
				const char* eps_or_eta_File,
				float vp_Constant,
				float dta_Constant,
				float eps_or_eta_Constant,
				int eta_Flag,
				int swap_Flag
				);

		~Vp_Eps_Dta_Scanner();

		float Get_vp_X_Min();
		float Get_vp_X_Max();
		float Get_vp_Z_Min();
		float Get_vp_Z_Max();
		float Get_Eps_Min();
		float Get_Eps_Max();
		float Get_Dta_Min();
		float Get_Dta_Max();

		void Read(long offset, int ubeg, int uend, int v, int w);

	private:
		Transform_Eps* _trans;

		char* _vp_File;
		char* _dta_File;
		char* _eps_or_eta_File;

		float _vp_Constant;
		float _dta_Constant;
		float _eps_or_eta_Constant;

		FILE* _vp_fp;
		FILE* _dta_fp;
		FILE* _eps_or_eta_fp;

		int _inpbuflen;
		float* _inpbuf_vp;
		float* _inpbuf_dta;
		float* _inpbuf_eps_or_eta;

		float _vp_Z_min;
		float _vp_Z_max;
		float _vp_X_min;
		float _vp_X_max;
		float _eps_min;
		float _eps_max;
		float _dta_min;
		float _dta_max;
};
#endif

