#ifndef ACOUSTIC_EARTH_MODEL_HXX
#define ACOUSTIC_EARTH_MODEL_HXX

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "Index_Calculator.hxx"
#include "Attribute_Reader.hxx"
#include "Volume_Iterator.hxx"

#define VELMASK 32767
#define EPSMASK 255
#define DELMASK 255
#define C44C33MASK 1
#define QMASK 255
#define DENMASK 255
#define DIPMASK 255
#define AZMMASK 255

#define SHIFTEps 15
#define SHIFTDel 23
#define SHIFTC44C33 31

#define SHIFTQ 24
#define SHIFTDen 16
#define SHIFTAzm  8

class Acoustic_Earth_Model
{
public:
	Acoustic_Earth_Model(
		int Log_Level,
		int Kernel_Type,                        // 0->ISO, 1->VTI, 2->TTI
		int Target_Platform,                    // 0->CPU, 1->GPU
		const char* vp_File,
		const char* den_File,
		const char* eps_or_eta_File,
		const char* dta_File,
		const char* dip_or_dx_File,
		const char* azm_or_dy_File,
		const char* Q_File,
		int eta_Flag,
		int dipxdipy_Flag,
		int degrees_Flag,
		int swapin_Flag,
		int fast_Axis,
		int med_Axis,
		int slow_Axis,
		int dimx,
		int dimy,
		int dimz,
		int x0,
		int x1,
		int y0,
		int y1,
		int z0,
		int z1,
		int sky_z,                              // number of ZZzz'z in sky layer
		float sky_scale_factor,
		float dh,
		float dz
		);

	virtual ~Acoustic_Earth_Model();

	void Set_Log_Level(int log_level);

	//
	// Returns 1 if earth model and sub volume parameters are valid, 0 otherwise.
	// 
	int Is_Valid();

	//
	// Create compressed earth model.
	// Argument fq provides center frequency for Q calculations.
	// If fq == 0, center frequency is calculated as fmax/3.
	//
	void Create_Compressed_Earth_Model(
		int OTflag,
		float gamfac,
		float fq
		);

	//
	// Scan earth model to determine some values required to compute run time parameters.
	// This method should only be called for dryruns.
	// Calling Create_Compressed_Earth_Model will achieve the same thing.
	//
	void Scan_Earth_Model(
		int OTflag,
		float gamfac
		);

	int Get_Compressed_Earth_Model(
		int*& PadVelAnis,
		int*& PadDenAng,
		int***& VelAnis,
		int***& DenAng
		);

	void Compute_FMAX(int OTflag);
	void Compute_DT(
		int OTflag,
		float gamfac
		);

	int FMAX_Is_Valid();
	float Get_FMAX();

	int DT_Is_Valid();
	float Get_DT();

	int Get_NX();
	int Get_Actual_NX();
	int Get_NY();
	int Get_NZ();

	float Get_Dip_Scaler() {return _dip_binsize;}
	float Get_Dip_Min() {return _dip_Min;}
	float Get_Azm_Scaler() {return _azm_binsize;}
	float Get_Azm_Min() {return _azm_Min;}
	float Get_Den_Scaler() {return _den_binsize;}
	float Get_Den_Min() {return _den_Min;}
	float Get_Q_Scaler() {return _Q_binsize;}
	float Get_Q_Min() {return _Q_Min;}
	float Get_Vp_Scaler() {return _vp_binsize;}
	float Get_Vp_Min() {return _vp_min;}
	float Get_Dta_Scaler() {return _dta_binsize;}
	float Get_Dta_Min() {return _dta_Min;}
	float Get_Eps_Scaler() {return _eps_binsize;}
	float Get_Eps_Min() {return _eps_Min;}

	//
	// Write out X-Z cross-section of earth model.
	// Filename arguments can be nil if you don't want to output all attributes.
	//
	void Write_XZ_Slice_Gnuplot(
		int iY,
		const char* Vp_filename,
		const char* Eps_filename,
		const char* Dta_filename,
		const char* C44C33_filename,
		const char* Den_filename,
		const char* Dip_filename,
		const char* Azm_filename,
		const char* Q_filename
		);

	void Write_XY_Slice_Gnuplot(
		int iZ,
		const char* Vp_filename,
		const char* Eps_filename,
		const char* Dta_filename,
		const char* C44C33_filename,
		const char* Den_filename,
		const char* Dip_filename,
		const char* Azm_filename,
		const char* Q_filename
		);

private:
	int _Kernel_Type;
	int _Target_Platform;

	int _Valid_PropParm_Ranges;
	float _vp_X_Min;
	float _vp_X_Max;
	float _vp_Z_Min;
	float _vp_Z_Max;
	float _den_Min;
	float _den_Max;
	float _eps_Min;
	float _eps_Max;
	float _dta_Min;
	float _dta_Max;
	float _dip_Min;
	float _dip_Max;
	float _azm_Min;
	float _azm_Max;
	float _Q_Min;
	float _Q_Max;

	float _vp_min;
	float _vp_max;
	float _vp_range;
	float _vp_binsize;

	float _eps_range;
	float _eps_binsize;

	float _dta_range;
	float _dta_binsize;

	float _den_range;
	float _den_binsize;

	float _dip_range;
	float _dip_binsize;

	float _azm_range;
	float _azm_binsize;
	
	float _Q_range;
	float _Q_binsize;

	float _dh;
	float _dz;

	// propagation parameters
	int _Valid_FMAX;
	float _fmax;

	int _Valid_DT;
	float _dt;
	float _fq;

	int _Valid_Compressed_Earth_Model;
	int* _PadVelAnis;
	int* _PadDenAng;
	int*** _VelAnis;
	int*** _DenAng;

	int _Valid_Earth_Model_Parameters;
	int _log_level;

	char* _vp_File;
	char* _den_File;
	char* _eps_or_eta_File;
	char* _dta_File;
	char* _dip_or_dx_File;
	char* _azm_or_dy_File;
	char* _Q_File;

	int _vp_Is_Constant;
	int _den_Is_Constant;
	int _eps_or_eta_Is_Constant;
	int _dta_Is_Constant;
	int _dip_or_dx_Is_Constant;
	int _azm_or_dy_Is_Constant;
	int _Q_Is_Constant;

	float _vp_Constant;
	float _den_Constant;
	float _eps_or_eta_Constant;
	float _dta_Constant;
	float _dip_or_dx_Constant;
	float _azm_or_dy_Constant;
	float _Q_Constant;

	int _eta_Flag;
	int _dipxdipy_Flag;
	int _degrees_Flag;
	int _swap_Flag;

	Index_Calculator* _icalc;	
	UVW_XYZ_Calculator* _ucalc;
	Volume_Iterator* _volume_iterator;

	int _dimx;
	int _dimy;
	int _dimz;

	int _Valid_Sub_Volume;
	int _x0;
	int _x1;
	int _y0;
	int _y1;
	int _z0;
	int _z1;

	int _actual_nx;
	int _nx;
	int _ny;
	int _nz;
	int _xh;
	int _yh;
	int _zh;

	unsigned long _dim;

	int _sky_z;
	float _sky_scale_factor;

	int _GPU_HostAlloc;

	void omp_memset(void* p, unsigned long wflen);
	void Check_File_Arg(const char* lbl, const char* inp_File, char*& File, float& Const_Val, int& Is_Const_Val);
	void Check_Earth_Model_Parameters();
	void _Add_Q_To_Earth_Model(float fq);
	void _Create_Compressed_Earth_Model_No_Q(int OTflag, float gamfac);
	void Find_Min_Max_Q(
		const char* name,
		const char* filename,
		float const_val,
		float fq,
		float dt,
		float& min,
		float& max
		);
	void Find_Min_Max_DIP_AZM(
		const char* name,
		const char* dip_or_dx_File,
		const char* azm_or_dy_File,
		float dip_or_dx_Constant,
		float azm_or_dy_Constant,
		int degrees_Flag,
		int dipxdipy_Flag,
		float& dip_Min,
		float& dip_Max,
		float& azm_Min,
		float& azm_Max
		);
	void Find_Min_Max_Single_Attribute(
		const char* name,
		const char* filename,
		float const_val,
		float& min,
		float& max
		);
	void Find_Min_Max_VP_DTA_EPS(
		const char* name,
		const char* vp_File,
		const char* dta_File,
		const char* eps_or_eta_File,
		float vp_Constant,
		float dta_Constant,
		float eps_or_eta_Constant,
		int eta_Flag,
		float& vpminX,
		float& vpmaxX,
		float& vpminZ,
		float& vpmaxZ,
		float& Epsmin,
		float& Epsmax,
		float& Delmin,
		float& Delmax
		);

	void Add_Vp(
		int* PadVelAnis,
		const char* name,
		const char* vp_File,
		float vp_Constant,
		float vp_min,
		float vp_range,
		int vel_mask,
		int vel_shift
		);
	void Add_Eps_Dta(
		int* PadVelAnis,
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
		int dta_shift
		);
	void Add_Property(
		int* PadDenAng,
		const char* name,
		const char* filename,
		float const_val,
		float min,
		float range,
		int prop_mask,
		int prop_shift
		);
	void Add_Dip_Azm(
		int* PadDenAng,
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
		int degrees_Flag,
		int dipxdipy_Flag
		);
	void Add_Q(
		int* PadVelAnis,
		const char* name,
		const char* Q_File,
		float Q_Constant,
		float Q_min,
		float Q_range,
		int Q_mask,
		int Q_shift,
		float fq,
		float dt
		);

	void Destroy_Compressed_Earth_Model();
};
#endif

