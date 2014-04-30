#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <cuda_runtime.h>

#include "Index_Calculator.hxx"
#include "UVW_XYZ_Calculator.hxx"
#include "Attribute_Reader.hxx"
#include "Acoustic_Earth_Model.hxx"
#include "Attribute_Transformer.hxx"

Acoustic_Earth_Model::Acoustic_Earth_Model(
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
		float VsoVp0,
		float newfmax,
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
		int sky_z,				// number of ZZzz'z in sky layer
		float sky_scale_factor,
		int absorb_z,
		float dh,
		float dz
		)
{
	_log_level = Log_Level;

	if (_log_level >= 4)
	{
		printf("Kernel_Type = %d\n",Kernel_Type);
		printf("Target_Platform = %d\n",Target_Platform);
		printf("vp_File = %s\n",vp_File!=0L?vp_File:"nil");
		printf("den_File = %s\n",den_File!=0L?den_File:"nil");
		printf("eps_or_eta_File = %s\n",eps_or_eta_File!=0L?eps_or_eta_File:"nil");
		printf("dta_File = %s\n",dta_File!=0L?dta_File:"nil");
		printf("dip_or_dx_File = %s\n",dip_or_dx_File!=0L?dip_or_dx_File:"nil");
		printf("azm_or_dy_File = %s\n",azm_or_dy_File!=0L?azm_or_dy_File:"nil");
		printf("Q_File = %s\n",Q_File!=0L?Q_File:"nil");
		printf("eta_Flag = %d\n",eta_Flag);
		printf("dipxdipy_Flag = %d\n",dipxdipy_Flag);
		printf("degrees_Flag = %d\n",degrees_Flag);
		printf("swapin_Flag = %d\n",swapin_Flag);
		printf("VsoVp0 = %f\n",VsoVp0);
		printf("newfmax = %f\n",newfmax);
		printf("fast_Axis = %d\n",fast_Axis);
		printf("med_Axis = %d\n",med_Axis);
		printf("slow_Axis = %d\n",slow_Axis);
		printf("dimx = %d\n",dimx);
		printf("dimy = %d\n",dimy);
		printf("dimz = %d\n",dimz);
		printf("x0 = %d\n",x0);
		printf("x1 = %d\n",x1);
		printf("y0 = %d\n",y0);
		printf("y1 = %d\n",y1);
		printf("z0 = %d\n",z0);
		printf("z1 = %d\n",z1);
		printf("sky_z = %d\n",sky_z);
		printf("sky_scale_factor = %e\n",sky_scale_factor);
		printf("absorb_z = %d\n",absorb_z);
		printf("dh = %f\n",dh);
		printf("dz = %f\n",dz);
	}

	_Valid_Compressed_Earth_Model = 0;
	_Valid_Earth_Model_Parameters = 1;
	_newfmax = newfmax;
	_computed_fmax = 0.0f;
	_fmax = 0.0f;
	_Valid_FMAX = 0;
	_Valid_DT = 0;

	_Kernel_Type = Kernel_Type;	
	_Target_Platform = Target_Platform;

	_vp_File = 0L;
	_den_File = 0L;
	_eps_or_eta_File = 0L;
	_dta_File = 0L;
	_dip_or_dx_File = 0L;
	_azm_or_dy_File = 0L;
	_Q_File = 0L;

	Check_File_Arg("vp", vp_File, _vp_File, _vp_Constant, _vp_Is_Constant);
	Check_File_Arg("den", den_File, _den_File, _den_Constant, _den_Is_Constant);
	Check_File_Arg("eps_or_eta", eps_or_eta_File, _eps_or_eta_File, _eps_or_eta_Constant, _eps_or_eta_Is_Constant);
	Check_File_Arg("dta", dta_File, _dta_File, _dta_Constant, _dta_Is_Constant);
	Check_File_Arg("dip_or_dx", dip_or_dx_File, _dip_or_dx_File, _dip_or_dx_Constant, _dip_or_dx_Is_Constant);
	Check_File_Arg("azm_or_dy", azm_or_dy_File, _azm_or_dy_File, _azm_or_dy_Constant, _azm_or_dy_Is_Constant);
	Check_File_Arg("Q", Q_File, _Q_File, _Q_Constant, _Q_Is_Constant);

	_eta_Flag = eta_Flag;
	_dipxdipy_Flag = dipxdipy_Flag;
	_degrees_Flag = degrees_Flag;
	_swap_Flag = swapin_Flag;

	_VsoVp0 = VsoVp0;

	_sky_z = sky_z;
	_sky_scale_factor = sky_scale_factor;
	_absorb_z = absorb_z;

	_dimx = dimx;
	_dimy = dimy;
	_dimz = dimz;
	if ((x0 <= x1 && x0 >= 0 && x1 < _dimx) && (y0 <= y1 && y0 >= 0 && y1 < _dimy) && (z0 <= z1 && z0 >= 0 && z1 < _dimz)) 
		_Valid_Sub_Volume = 1;
	else
		_Valid_Sub_Volume = 0;
	_x0 = x0;
	_x1 = x1;
	_y0 = y0;
	_y1 = y1;
	_z0 = z0;
	_z1 = z1;

	_actual_nx = _x1 - _x0 + 1;
	_actual_ny = _y1 - _y0 + 1;
	_actual_nz = _z1 - _z0 + 1;
	if (_Target_Platform == 0) // CPU
	{
		_nx = ((_actual_nx + 7) >> 3) << 3;
		_ny = _actual_ny;
		_nz = _actual_nz;
		_xh = 16;
		_yh = 9;
		_zh = 0;
	}
	else if (_Target_Platform == 1) // GPU
	{
		_nx = ((_actual_nx + 15) >> 4) << 4;
		_ny = ((_actual_ny + 3) >> 2) << 2;
		_nz = (((_actual_nz + _sky_z + _absorb_z + 3) >> 2) << 2) - _sky_z - _absorb_z;
		_xh = 0;
		_yh = 0;
		_zh = 0;
	}
	_dim = (unsigned long)(_nx+2*_xh) * (unsigned long)(_ny+2*_yh) * (unsigned long)(_nz+_sky_z+_absorb_z+2*_zh) * (unsigned long)4;

	_volume_iterator = new Volume_Iterator(_log_level,fast_Axis,med_Axis,slow_Axis,dimx,dimy,dimz,x0,x1,y0,y1,z0,z1);

	_dh = dh;
	_dz = dz;

	_icalc = new Index_Calculator(_Target_Platform,_Kernel_Type==0?1:2,_nx,_ny,_nz,_sky_z,_absorb_z,_xh,_yh,_zh);
	_ucalc = new UVW_XYZ_Calculator(fast_Axis,med_Axis,slow_Axis);

	_GPU_HostAlloc = 0;
	_PadVelAnis = 0L;
	_PadDenAng = 0L;

	Check_Earth_Model_Parameters();

	_Valid_PropParm_Ranges = 0;
	_Valid_Compressed_Earth_Model = 0;
	_VelAnis = 0L;
	_DenAng = 0L;
	_PadVelAnis = 0L;
	_PadDenAng = 0L;
}

Acoustic_Earth_Model::~Acoustic_Earth_Model()
{
	// NB smart-ass!
	// I know delete operator handles nil operand correctly,
	// but it disturbs me to see any memory de-allocation without check for nil.
	if (_vp_File != 0L) delete [] _vp_File;
	if (_den_File != 0L) delete [] _den_File;
	if (_eps_or_eta_File != 0L) delete [] _eps_or_eta_File;
	if (_dta_File != 0L) delete [] _dta_File;
	if (_dip_or_dx_File != 0L) delete [] _dip_or_dx_File;
	if (_azm_or_dy_File != 0L) delete [] _azm_or_dy_File;
	if (_Q_File != 0L) delete [] _Q_File;

	if (_ucalc != 0L) delete _ucalc;
	if (_icalc != 0L) delete _icalc;
	if (_volume_iterator != 0L) delete _volume_iterator;

	Destroy_Compressed_Earth_Model();
}

void Acoustic_Earth_Model::Set_Log_Level(int log_level)
{
	_log_level = log_level;
}

//
// Returns 1 if earth model and sub volume parameters are valid, 0 otherwise.// 
int Acoustic_Earth_Model::Is_Valid()
{
	return _Valid_Earth_Model_Parameters && _Valid_Sub_Volume ? 1 : 0;
}

//
// Create compressed earth model.
// Argument fq provides center frequency for Q calculations.
// If fq == 0, center frequency is calculated as fmax/3.
//
void Acoustic_Earth_Model::Create_Compressed_Earth_Model(
		int OTflag,
		float gamfac,
		float fq,
		float dtout
		)
{
	_Create_Compressed_Earth_Model_No_Q(OTflag,gamfac,dtout);
	_Add_Q_To_Earth_Model(fq!=0.0f?fq:_fmax/3.0f);

	// extend earth model to accomodate sky layer and/or absorbing layer
	int z_free_surface = 1;
	int z_top_of_sky_layer = -_sky_z;
	int z_top_of_absorb_layer = z_top_of_sky_layer - _absorb_z;
	if (z_top_of_absorb_layer < z_free_surface)
	{
		// To-Do: Multiply density by sky_factor.
		// Current code simply replicates top layer.
		for (int iz = z_top_of_absorb_layer;  iz < z_free_surface;  ++iz)
		{
			for (int iy = -_yh;  iy < _ny + _yh;  ++iy)
			{
				for (int ix = -_xh;  ix < _nx + _xh;  ++ix)
				{
					unsigned long dst_idx = _icalc->Calculate_Index(ix,iy,iz);
					unsigned long src_idx = _icalc->Calculate_Index(ix,iy,z_free_surface);
					if (_Kernel_Type > 0) _PadDenAng[dst_idx] = _PadDenAng[src_idx];
					_PadVelAnis[dst_idx] = _PadVelAnis[src_idx];
				}
			}
		}
	}
}

//
// Scan earth model to determine some values required to compute run time parameters.
// This method should only be called for dryruns.
// Calling Create_Compressed_Earth_Model will achieve the same thing.
//
void Acoustic_Earth_Model::Scan_Earth_Model(
	int OTflag,
	float gamfac,
	float dtout
	)
{
	Find_Min_Max_VP_DTA_EPS("Scanning VP, DTA & EPS",_vp_File,_dta_File,_eps_or_eta_File,_vp_Constant,_dta_Constant,_eps_or_eta_Constant,_eta_Flag,_vp_X_Min,_vp_X_Max,_vp_Z_Min,_vp_Z_Max,_eps_Min,_eps_Max,_dta_Min,_dta_Max);
	Compute_FMAX(OTflag);
	Compute_DT(OTflag,gamfac,dtout);
}

int Acoustic_Earth_Model::Get_Compressed_Earth_Model(
		int*& PadVelAnis,
		int*& PadDenAng,
		int***& VelAnis,
		int***& DenAng
		)
{
	if (_Valid_Compressed_Earth_Model)
	{
		PadVelAnis = _PadVelAnis;
		PadDenAng = _PadDenAng;
		VelAnis = _VelAnis;
		DenAng = _DenAng;
		return -1;
	}
	else
	{
		PadVelAnis = 0L;
		PadDenAng = 0L;
		VelAnis = 0L;
		DenAng = 0L;
		return 0;
	}
}

void Acoustic_Earth_Model::Compute_FMAX(int OTflag)
{
	float khmax;
	if(OTflag==4)
	{
		khmax = 2.41f;
	}
	else /* OT2 */
	{
		khmax = 1.76f;
	}
	float minfac = (_vp_X_Min / _dh < _vp_Z_Min / _dz) ? _vp_X_Min / _dh : _vp_Z_Min / _dz;
	_computed_fmax = khmax * minfac / (2.0f * M_PI);
	if (_newfmax < 0.0f)
	{
		_fmax = _computed_fmax * (-_newfmax);
	}
	else if (_newfmax == 0.0f)
	{
		_fmax = _computed_fmax;
	}
	else
	{
		_fmax = _newfmax;
	}
	_Valid_FMAX = 1;
	printf("_fmax = %f, _computed_fmax = %f, _vp_X_Min = %f, _vp_X_Max = %f, _vp_Z_Min = %f, _vp_Z_Max = %f\n",_fmax,_computed_fmax,_vp_X_Min,_vp_X_Max,_vp_Z_Min,_vp_Z_Max);
}

void Acoustic_Earth_Model::Compute_DT(
		int OTflag,
		float gamfac,
		float dtout
		)
{
	float gam;
	if(OTflag==4)  
	{ 
		gam = 0.6601f*gamfac; 
	}
	else /* OT2 */ 
	{
		gam = 0.3801f*gamfac;
	}
	float vpmax = _vp_X_Max > _vp_Z_Max ? _vp_X_Max : _vp_Z_Max;
	_dt = gam * ( (_dz / vpmax < _dh / vpmax) ? _dz / vpmax : _dh / vpmax );
	if (_dt > dtout) _dt = dtout;
	_Valid_DT = 1;
	printf("_dt = %f, _vp_X_Min = %f, _vp_X_Max = %f, _vp_Z_Min = %f, _vp_Z_Max = %f\n",_dt,_vp_X_Min,_vp_X_Max,_vp_Z_Min,_vp_Z_Max);
}

int Acoustic_Earth_Model::FMAX_Is_Valid()
{
	return _Valid_FMAX;
}

float Acoustic_Earth_Model::Get_FMAX()
{
	return _fmax;
}

int Acoustic_Earth_Model::DT_Is_Valid()
{
	return _Valid_DT;
}

float Acoustic_Earth_Model::Get_DT()
{
	return _dt;
}

int Acoustic_Earth_Model::Get_Actual_NX()
{
	return _actual_nx;
}

int Acoustic_Earth_Model::Get_Actual_NY()
{
	return _actual_ny;
}

int Acoustic_Earth_Model::Get_Actual_NZ()
{
	return _actual_nz;
}

int Acoustic_Earth_Model::Get_NX()
{
	return _nx;
}

int Acoustic_Earth_Model::Get_NY()
{
	return _ny;
}

int Acoustic_Earth_Model::Get_NZ()
{
	return _nz;
}

void Acoustic_Earth_Model::Write_XZ_Slice_Gnuplot(
		int iY,
		const char* Vp_filename,
		const char* Eps_filename,
		const char* Dta_filename,
		const char* C44C33_filename,
		const char* Den_filename,
		const char* Dip_filename,
		const char* Azm_filename,
		const char* Q_filename
		)
{
	printf("Acoustic_Earth_Model::Write_XZ_Slice_Gnuplot()\n");

	FILE* Vp_fp = 0L;
	FILE* Eps_fp = 0L;
	FILE* Dta_fp = 0L;
	FILE* C44C33_fp = 0L;
	FILE* Den_fp = 0L;
	FILE* Dip_fp = 0L;
	FILE* Azm_fp = 0L;
	FILE* Q_fp = 0L;
	if (_Kernel_Type > 1)
	{
		if (Dip_filename != 0L) Dip_fp = fopen(Dip_filename, "w");
		if (Azm_filename != 0L) Azm_fp = fopen(Azm_filename, "w");
	}
	if (_Kernel_Type > 0)
	{
		if (Eps_filename != 0L) Eps_fp = fopen(Eps_filename, "w");
		if (Dta_filename != 0L) Dta_fp = fopen(Dta_filename, "w");
		if (C44C33_filename != 0L) C44C33_fp = fopen(C44C33_filename, "w");
	}
	if (Vp_filename != 0L) Vp_fp = fopen(Vp_filename, "w");
	if (Den_filename != 0L) Den_fp = fopen(Den_filename, "w");
	if (Q_filename != 0L) Q_fp = fopen(Q_filename, "w");

	for (int iZ = -(_sky_z+_absorb_z);  iZ < _nz;  ++iZ)
	{
		for (int iX = 0;  iX < _nx;  ++iX)
		{
			unsigned long idx = _icalc->Calculate_Index(iX,iY,iZ);
			if (_Kernel_Type == 0)
			{
				int VelAnis = _PadVelAnis[idx];
				float Vp = Attribute_Loader::Decompress_Attribute(VelAnis,_vp_min,_vp_binsize,VELMASK,0);
				float Den = Attribute_Loader::Decompress_Attribute(VelAnis,_den_Min,_den_binsize,DENMASK,SHIFTDen);
				float inv_Q = Attribute_Loader::Decompress_Attribute(VelAnis,_Q_Min,_Q_binsize,QMASK,SHIFTQ);
				if (Vp_fp != 0L) fprintf(Vp_fp, "%d %d %e\n",iX,iZ,Vp);
				if (Den_fp != 0L) fprintf(Den_fp, "%d %d %e\n",iX,iZ,Den);
				if (Q_fp != 0L) fprintf(Q_fp, "%d %d %e\n",iX,iZ,inv_Q);
			}
			else if (_Kernel_Type == 1)
			{
				int VelAnis = _PadVelAnis[idx];
				float Vp = Attribute_Loader::Decompress_Attribute(VelAnis,_vp_min,_vp_binsize,VELMASK,0);
				float Eps = Attribute_Loader::Decompress_Attribute(VelAnis,_eps_Min,_eps_binsize,EPSMASK,SHIFTEps);
				float Dta = Attribute_Loader::Decompress_Attribute(VelAnis,_dta_Min,_dta_binsize,DELMASK,SHIFTDel);
				float C44C33 = Attribute_Loader::Decompress_Attribute(VelAnis,_C44C33_Min,_C44C33_binsize,C44C33MASK,SHIFTC44C33);
				if (Vp_fp != 0L) fprintf(Vp_fp, "%d %d %e\n",iX,iZ,Vp);
				if (Eps_fp != 0L) fprintf(Eps_fp, "%d %d %e\n",iX,iZ,Eps);
				if (Dta_fp != 0L) fprintf(Dta_fp, "%d %d %e\n",iX,iZ,Dta);
				if (C44C33_fp != 0L) fprintf(C44C33_fp, "%d %d %e\n",iX,iZ,C44C33);

				int DenAng = _PadDenAng[idx];
				float Den = Attribute_Loader::Decompress_Attribute(DenAng,_den_Min,_den_binsize,DENMASK,SHIFTDen);
				float inv_Q = Attribute_Loader::Decompress_Attribute(DenAng,_Q_Min,_Q_binsize,QMASK,SHIFTQ);
				if (Den_fp != 0L) fprintf(Den_fp, "%d %d %e\n",iX,iZ,Den);
                                if (Q_fp != 0L) fprintf(Q_fp, "%d %d %e\n",iX,iZ,inv_Q);
			}
			else if (_Kernel_Type == 2)
			{
				int VelAnis = _PadVelAnis[idx];
				float Vp = Attribute_Loader::Decompress_Attribute(VelAnis,_vp_min,_vp_binsize,VELMASK,0);
				float Eps = Attribute_Loader::Decompress_Attribute(VelAnis,_eps_Min,_eps_binsize,EPSMASK,SHIFTEps);
				float Dta = Attribute_Loader::Decompress_Attribute(VelAnis,_dta_Min,_dta_binsize,DELMASK,SHIFTDel);
				float C44C33 = Attribute_Loader::Decompress_Attribute(VelAnis,_C44C33_Min,_C44C33_binsize,C44C33MASK,SHIFTC44C33);
				if (Vp_fp != 0L) fprintf(Vp_fp, "%d %d %e\n",iX,iZ,Vp);
				if (Eps_fp != 0L) fprintf(Eps_fp, "%d %d %e\n",iX,iZ,Eps);
				if (Dta_fp != 0L) fprintf(Dta_fp, "%d %d %e\n",iX,iZ,Dta);
				if (C44C33_fp != 0L) fprintf(C44C33_fp, "%d %d %e\n",iX,iZ,C44C33);

				int DenAng = _PadDenAng[idx];
				float Den = Attribute_Loader::Decompress_Attribute(DenAng,_den_Min,_den_binsize,DENMASK,SHIFTDen);
				float Dip = Attribute_Loader::Decompress_Attribute(DenAng,_dip_Min,_dip_binsize,DIPMASK,0);
				float Azm = Attribute_Loader::Decompress_Attribute(DenAng,_azm_Min,_azm_binsize,AZMMASK,SHIFTAzm);
				float inv_Q = Attribute_Loader::Decompress_Attribute(DenAng,_Q_Min,_Q_binsize,QMASK,SHIFTQ);
				if (Den_fp != 0L) fprintf(Den_fp, "%d %d %e\n",iX,iZ,Den);
				if (Dip_fp != 0L) fprintf(Dip_fp, "%d %d %e\n",iX,iZ,Dip);
				if (Azm_fp != 0L) fprintf(Azm_fp, "%d %d %e\n",iX,iZ,Azm);
                                if (Q_fp != 0L) fprintf(Q_fp, "%d %d %e\n",iX,iZ,inv_Q);
			}
		}
		if (Vp_fp != 0L) fprintf(Vp_fp, "\n");
		if (Den_fp != 0L) fprintf(Den_fp, "\n");
		if (Eps_fp != 0L) fprintf(Eps_fp, "\n");
		if (Dta_fp != 0L) fprintf(Dta_fp, "\n");
		if (C44C33_fp != 0L) fprintf(C44C33_fp, "\n");
		if (Dip_fp != 0L) fprintf(Dip_fp, "\n");
		if (Azm_fp != 0L) fprintf(Azm_fp, "\n");
		if (Q_fp != 0L) fprintf(Q_fp, "\n");
	}

	if (Dip_fp != 0L) fclose(Dip_fp);
	if (Azm_fp != 0L) fclose(Azm_fp);
	if (Eps_fp != 0L) fclose(Eps_fp);
	if (Dta_fp != 0L) fclose(Dta_fp);
	if (C44C33_fp != 0L) fclose(C44C33_fp);
	if (Vp_fp != 0L) fclose(Vp_fp);
	if (Den_fp != 0L) fclose(Den_fp);
	if (Q_fp != 0L) fclose(Q_fp);
}

void Acoustic_Earth_Model::Write_XY_Slice_Gnuplot(
		int iZ,
		const char* Vp_filename,
		const char* Eps_filename,
		const char* Dta_filename,
		const char* C44C33_filename,
		const char* Den_filename,
		const char* Dip_filename,
		const char* Azm_filename,
		const char* Q_filename
		)
{
	FILE* Vp_fp = 0L;
	FILE* Eps_fp = 0L;
	FILE* Dta_fp = 0L;
	FILE* C44C33_fp = 0L;
	FILE* Den_fp = 0L;
	FILE* Dip_fp = 0L;
	FILE* Azm_fp = 0L;
	FILE* Q_fp = 0L;
	if (_Kernel_Type > 1)
	{
		if (Dip_filename != 0L) Dip_fp = fopen(Dip_filename, "w");
		if (Azm_filename != 0L) Azm_fp = fopen(Azm_filename, "w");
	}
	else if (_Kernel_Type > 0)
	{
		if (Eps_filename != 0L) Eps_fp = fopen(Eps_filename, "w");
		if (Dta_filename != 0L) Dta_fp = fopen(Dta_filename, "w");
		if (C44C33_filename != 0L) C44C33_fp = fopen(C44C33_filename, "w");
	}
	if (Vp_filename != 0L) Vp_fp = fopen(Vp_filename, "w");
	if (Den_filename != 0L) Den_fp = fopen(Den_filename, "w");
	if (Q_filename != 0L) Q_fp = fopen(Q_filename, "w");

	for (int iY = 0;  iY < _ny;  ++iY)
	{
		for (int iX = 0;  iX < _nx;  ++iX)
		{
			unsigned long idx = _icalc->Calculate_Index(iX,iY,iZ);
			if (_Kernel_Type == 0)
			{
				int VelAnis = _PadVelAnis[idx];
				float Vp = Attribute_Loader::Decompress_Attribute(VelAnis,_vp_min,_vp_binsize,VELMASK,0);
				float Den = Attribute_Loader::Decompress_Attribute(VelAnis,_den_Min,_den_binsize,DENMASK,SHIFTDen);
				float inv_Q = Attribute_Loader::Decompress_Attribute(VelAnis,_Q_Min,_Q_binsize,QMASK,SHIFTQ);
				if (Vp_fp != 0L) fprintf(Vp_fp, "%d %d %e\n",iX,iZ,Vp);
				if (Den_fp != 0L) fprintf(Den_fp, "%d %d %e\n",iX,iZ,Den);
				if (Q_fp != 0L) fprintf(Q_fp, "%d %d %e\n",iX,iZ,inv_Q);
			}
			else if (_Kernel_Type == 1)
			{
				int VelAnis = _PadVelAnis[idx];
				float Vp = Attribute_Loader::Decompress_Attribute(VelAnis,_vp_min,_vp_binsize,VELMASK,0);
				float Eps = Attribute_Loader::Decompress_Attribute(VelAnis,_eps_Min,_eps_binsize,EPSMASK,SHIFTEps);
				float Dta = Attribute_Loader::Decompress_Attribute(VelAnis,_dta_Min,_dta_binsize,DELMASK,SHIFTDel);
				float C44C33 = Attribute_Loader::Decompress_Attribute(VelAnis,_C44C33_Min,_C44C33_binsize,C44C33MASK,SHIFTC44C33);
				if (Vp_fp != 0L) fprintf(Vp_fp, "%d %d %e\n",iX,iZ,Vp);
				if (Eps_fp != 0L) fprintf(Eps_fp, "%d %d %e\n",iX,iZ,Eps);
				if (Dta_fp != 0L) fprintf(Dta_fp, "%d %d %e\n",iX,iZ,Dta);
				if (C44C33_fp != 0L) fprintf(C44C33_fp, "%d %d %e\n",iX,iZ,C44C33);

				int DenAng = _PadDenAng[idx];
				float Den = Attribute_Loader::Decompress_Attribute(DenAng,_den_Min,_den_binsize,DENMASK,SHIFTDen);
				float inv_Q = Attribute_Loader::Decompress_Attribute(DenAng,_Q_Min,_Q_binsize,QMASK,SHIFTQ);
				if (Den_fp != 0L) fprintf(Den_fp, "%d %d %e\n",iX,iZ,Den);
                                if (Q_fp != 0L) fprintf(Q_fp, "%d %d %e\n",iX,iZ,inv_Q);
			}
			else if (_Kernel_Type == 2)
			{
				int VelAnis = _PadVelAnis[idx];
				float Vp = Attribute_Loader::Decompress_Attribute(VelAnis,_vp_min,_vp_binsize,VELMASK,0);
				float Eps = Attribute_Loader::Decompress_Attribute(VelAnis,_eps_Min,_eps_binsize,EPSMASK,SHIFTEps);
				float Dta = Attribute_Loader::Decompress_Attribute(VelAnis,_dta_Min,_dta_binsize,DELMASK,SHIFTDel);
				float C44C33 = Attribute_Loader::Decompress_Attribute(VelAnis,_C44C33_Min,_C44C33_binsize,C44C33MASK,SHIFTC44C33);
				if (Vp_fp != 0L) fprintf(Vp_fp, "%d %d %e\n",iX,iZ,Vp);
				if (Eps_fp != 0L) fprintf(Eps_fp, "%d %d %e\n",iX,iZ,Eps);
				if (Dta_fp != 0L) fprintf(Dta_fp, "%d %d %e\n",iX,iZ,Dta);
				if (C44C33_fp != 0L) fprintf(C44C33_fp, "%d %d %e\n",iX,iZ,C44C33);

				int DenAng = _PadDenAng[idx];
				float Den = Attribute_Loader::Decompress_Attribute(DenAng,_den_Min,_den_binsize,DENMASK,SHIFTDen);
				float Dip = Attribute_Loader::Decompress_Attribute(DenAng,_dip_Min,_dip_binsize,DIPMASK,0);
				float Azm = Attribute_Loader::Decompress_Attribute(DenAng,_azm_Min,_azm_binsize,AZMMASK,SHIFTAzm);
				float inv_Q = Attribute_Loader::Decompress_Attribute(DenAng,_Q_Min,_Q_binsize,QMASK,SHIFTQ);
				if (Den_fp != 0L) fprintf(Den_fp, "%d %d %e\n",iX,iZ,Den);
				if (Dip_fp != 0L) fprintf(Dip_fp, "%d %d %e\n",iX,iZ,Dip);
				if (Azm_fp != 0L) fprintf(Azm_fp, "%d %d %e\n",iX,iZ,Azm);
                                if (Q_fp != 0L) fprintf(Q_fp, "%d %d %e\n",iX,iZ,inv_Q);
			}
		}
		if (Vp_fp != 0L) fprintf(Vp_fp, "\n");
		if (Den_fp != 0L) fprintf(Den_fp, "\n");
		if (Eps_fp != 0L) fprintf(Eps_fp, "\n");
		if (Dta_fp != 0L) fprintf(Dta_fp, "\n");
		if (C44C33_fp != 0L) fprintf(C44C33_fp, "\n");
		if (Dip_fp != 0L) fprintf(Dip_fp, "\n");
		if (Azm_fp != 0L) fprintf(Azm_fp, "\n");
		if (Q_fp != 0L) fprintf(Q_fp, "\n");
	}

	if (Dip_fp != 0L) fclose(Dip_fp);
	if (Azm_fp != 0L) fclose(Azm_fp);
	if (Eps_fp != 0L) fclose(Eps_fp);
	if (Dta_fp != 0L) fclose(Dta_fp);
	if (C44C33_fp != 0L) fclose(C44C33_fp);
	if (Vp_fp != 0L) fclose(Vp_fp);
	if (Den_fp != 0L) fclose(Den_fp);
	if (Q_fp != 0L) fclose(Q_fp);
}

void Acoustic_Earth_Model::omp_memset(void* p, unsigned long wflen)
{
	unsigned long wflen_32 = wflen >> 5;
#pragma omp parallel for
	for (int i = 0;  i < 32;  ++i)
	{
		unsigned long pos = wflen_32 * (unsigned long)i;
		memset((void*)(((char*)p)+pos), 0, wflen_32);
	}
}

void Acoustic_Earth_Model::Check_File_Arg(const char* lbl, const char* inp_File, char*& File, float& Const_Val, int& Is_Const_Val)
{
	if (File != 0L) delete [] File;
	File = 0L;

	char num_str[4096];
	int nmatched = sscanf(inp_File, "const:%s", num_str);
	if (nmatched == 1)
	{
		Const_Val = atof(num_str);
		Is_Const_Val = 1;
		printf("%s : using constant value of %f\n",lbl,Const_Val);
	}
	else
	{
		char scratch_filename[4096];
		sprintf(scratch_filename, "/scratch/aTTI_Cached/%s", inp_File);
		FILE* fp1 = fopen(scratch_filename, "r");
		if (fp1 != 0L)
		{
			fclose(fp1);
			Is_Const_Val = 0;
			int len = strlen(scratch_filename);
			File = new char[len+1];
			memcpy((void*)File, (void*)scratch_filename, len);
			File[len] = 0;
			printf("%s : Using cached file %s\n",lbl,scratch_filename);
		}
		else
		{
			FILE* fp2 = fopen(inp_File, "r");
			if (fp2 != 0L)
			{
				fclose(fp2);
				Is_Const_Val = 0;
				int len = strlen(inp_File);
				File = new char[len+1];
				memcpy((void*)File, (void*)inp_File, len);
				File[len] = 0;
				printf("%s : %s exists\n",lbl,inp_File);
			}
			else
			{
				printf("%s : BAD EARTH MODEL FILENAME (%s)!\n",lbl,inp_File);
				_Valid_Earth_Model_Parameters = 0;
			}
		}
	}
}

void Acoustic_Earth_Model::Check_Earth_Model_Parameters()
{
	if (
			(_vp_Is_Constant         && _vp_Constant         <= 0) ||
			(_den_Is_Constant        && _den_Constant        <= 0) ||
			(_eps_or_eta_Is_Constant && _eps_or_eta_Constant <  0) ||
			(_dta_Is_Constant        && _dta_Constant        <  0) ||
			(_dip_or_dx_Is_Constant  && _dip_or_dx_Constant  <  0) ||
			(_azm_or_dy_Is_Constant  && _azm_or_dy_Constant  <  0) ||
			(_Q_Is_Constant          && _Q_Constant          <= 0) ||
			_dimx < 1 ||
			_dimy < 1 ||
			_dimz < 1
	   )
	{
		printf("bookie: _vp_Constant=%f, _den_Constant=%f, _eps_or_eta_Constant=%f, _dta_Constant=%f, _dip_or_dx_Constant=%f, _azm_or_dy_Constant=%f, _Q_Constant=%f\n",_vp_Constant,_den_Constant,_eps_or_eta_Constant,_dta_Constant,_dip_or_dx_Constant,_azm_or_dy_Constant,_Q_Constant);
		_Valid_Earth_Model_Parameters = 0;
	}
	else
	{
		if (!_ucalc->Is_Valid()) 
		{
			printf("cookie\n");
			_Valid_Earth_Model_Parameters = 0;
		}
	}
}

void Acoustic_Earth_Model::_Add_Q_To_Earth_Model(float fq)
{
	_fq = fq;

	Find_Min_Max_Q("Scanning Q",_Q_File,_Q_Constant,_fq,_dt,_Q_Min,_Q_Max);

	_Q_range = (_Q_Min == _Q_Max) ? 0.0f : 1.0f / (_Q_Max - _Q_Min);
	_Q_binsize = (_Q_Min == _Q_Max) ? 0.0f : (_Q_Max - _Q_Min) / QMASK;

	if (_Kernel_Type == 0) // ISO
	{
		Add_Q(_PadVelAnis,"Loading Q",_Q_File,_Q_Constant,_Q_Min,_Q_range,QMASK,SHIFTQ,_fq,_dt);
	}
	else // VTI and TTI
	{
		Add_Q(_PadDenAng,"Loading Q",_Q_File,_Q_Constant,_Q_Min,_Q_range,QMASK,SHIFTQ,_fq,_dt);
	}
}

void Acoustic_Earth_Model::_Create_Compressed_Earth_Model_No_Q(
		int OTflag,
		float gamfac,
		float dtout
		)
{
	Destroy_Compressed_Earth_Model();

	// Scan vp, eps, dta properties.
	Scan_Earth_Model(OTflag,gamfac,dtout);

	Find_Min_Max_Single_Attribute("Scanning DEN",_den_File,_den_Constant,_den_Min,_den_Max);
	Find_Min_Max_DIP_AZM("Scanning DIP & AZM",_dip_or_dx_File,_azm_or_dy_File,_dip_or_dx_Constant,_azm_or_dy_Constant,_degrees_Flag,_dipxdipy_Flag,_dip_Min,_dip_Max,_azm_Min,_azm_Max);

	_vp_min = _vp_Z_Min;
	_vp_max = _vp_Z_Max;
	_vp_range = (_vp_min == _vp_max) ? 0.0f : 1.0f / (_vp_max - _vp_min);
	_vp_binsize = (_vp_min == _vp_max) ? 0.0f : (_vp_max - _vp_min) / VELMASK;

	_eps_range = (_eps_Min == _eps_Max) ? 0.0f : 1.0f / (_eps_Max - _eps_Min);
	_eps_binsize = (_eps_Min == _eps_Max) ? 0.0f : (_eps_Max - _eps_Min) / EPSMASK;

	_dta_range = (_dta_Min == _dta_Max) ? 0.0f : 1.0f / (_dta_Max - _dta_Min);
	_dta_binsize = (_dta_Min == _dta_Max) ? 0.0f : (_dta_Max - _dta_Min) / DELMASK;

	_C44C33_Min = 0.0f;
	_C44C33_Max = _VsoVp0 * _VsoVp0;
	_C44C33_range = 1.0f / (_C44C33_Max - _C44C33_Min);
	_C44C33_binsize = _C44C33_Max - _C44C33_Min;

	_den_range = (_den_Min == _den_Max) ? 0.0f : 1.0f / (_den_Max - _den_Min);
	_den_binsize = (_den_Min == _den_Max) ? 0.0f : (_den_Max - _den_Min) / DENMASK;

	_dip_range = (_dip_Min == _dip_Max) ? 0.0f : 1.0f / (_dip_Max - _dip_Min);
	_dip_binsize = (_dip_Min == _dip_Max) ? 0.0f : (_dip_Max - _dip_Min) / DIPMASK;

	_azm_range = (_azm_Min == _azm_Max) ? 0.0f : 1.0f / (_azm_Max - _azm_Min);
	_azm_binsize = (_azm_Min == _azm_Max) ? 0.0f : (_azm_Max - _azm_Min) / AZMMASK;
	
	if (_log_level >= 3)
	{
		printf("Vp=[%.0f,%.0f]\n",_vp_Z_Min,_vp_Z_Max);
		printf("Eps=[%.4f,%.4f]\n",_eps_Min,_eps_Max);
		printf("Dta=[%.4f,%.4f]\n",_dta_Min,_dta_Max);
		printf("Den=[%.4f,%.4f]\n",_den_Min,_den_Max);
		printf("Dip=[%.4f,%.4f]\n",_dip_Min,_dip_Max);
		printf("Azm=[%.4f,%.4f]\n",_azm_Min,_azm_Max);
	}

	// Allocate memory for compressed earth model
	if (_Target_Platform == 0) // CPU
	{
		_GPU_HostAlloc = 0;
		if (_Kernel_Type > 0)
		{
			int err1 = posix_memalign((void**)&_PadVelAnis, 64, _dim);
			int err2 = posix_memalign((void**)&_PadDenAng, 64, _dim);
			if (err1 != 0 || err2 != 0) return;
			_VelAnis = new int**[_nz];
			_DenAng = new int**[_nz];
			for (int iZ = 0;  iZ < _nz;  ++iZ)
			{
				_VelAnis[iZ] = new int*[_ny];
				_DenAng[iZ] = new int*[_ny];
				for (int iY = 0;  iY < _ny;  ++iY)
				{
					unsigned long idx = _icalc->Calculate_Index(0,iY,iZ);
					_VelAnis[iZ][iY] = _PadVelAnis + idx;
					_DenAng[iZ][iY] = _PadDenAng + idx;
				}
			}
		}
		else
		{
			int err1 = posix_memalign((void**)&_PadVelAnis, 64, _dim);
			if (err1 != 0) return;
			_VelAnis = new int**[_nz];
			for (int iZ = 0;  iZ < _nz;  ++iZ)
			{
				_VelAnis[iZ] = new int*[_ny];
				for (int iY = 0;  iY < _ny;  ++iY)
				{
					unsigned long idx = _icalc->Calculate_Index(0,iY,iZ);
					_VelAnis[iZ][iY] = _PadVelAnis + idx;
				}
			}
		}
		if (_PadVelAnis != 0L) omp_memset(_PadVelAnis, _dim);
		if (_PadDenAng != 0L) omp_memset(_PadDenAng, _dim);
	}
	else if (_Target_Platform == 1) // GPU
	{
		_GPU_HostAlloc = 1;	
		_VelAnis = 0L;
		_DenAng = 0L;
		if (_Kernel_Type > 0)
		{
			cudaError_t err1 = cudaHostAlloc((void**)&_PadDenAng, _dim * (long)2, cudaHostAllocDefault);
			_PadVelAnis = _PadDenAng + _icalc->Get_Field_Offset();
			if (err1 != cudaSuccess) return;
			omp_memset(_PadDenAng, _dim * (long)2);
		}
		else
		{
			cudaError_t err1 = cudaHostAlloc((void**)&_PadVelAnis, _dim, cudaHostAllocDefault);
			if (err1 != cudaSuccess) return;
			omp_memset(_PadVelAnis, _dim);
		}
	}

	if (_Kernel_Type == 0) // ISO
	{
		Add_Vp(_PadVelAnis,"Loading VP",_vp_File,_vp_Constant,_vp_min,_vp_range,VELMASK,0);
		Add_Property(_PadVelAnis,"Loading DEN",_den_File,_den_Constant,_den_Min,_den_range,DENMASK,SHIFTDen);
	}
	else // VTI and TTI
	{
		// VelAnis
		Add_Vp(_PadVelAnis,"Loading VP",_vp_File,_vp_Constant,_vp_min,_vp_range,VELMASK,0);
		Add_Eps_Dta(_PadVelAnis,"Loading EPS & DTA",_eps_or_eta_File,_eps_or_eta_Constant,_eta_Flag,_eps_Min,_eps_range,EPSMASK,SHIFTEps,_dta_File,_dta_Constant,_dta_Min,_dta_range,DELMASK,SHIFTDel);
		for (int iZ = 0;  iZ < _nz;  ++iZ)
		{
			for (int iY = 0;  iY < _ny;  ++iY)
			{
				for (int iX = 0;  iX < _nx;  ++iX)
				{
					unsigned long idx = _icalc->Calculate_Index(iX,iY,iZ);
					int mask = (EPSMASK << SHIFTEps) | (DELMASK << SHIFTDel);
					if ((_PadVelAnis[idx] & mask) != 0)
					{
						_PadVelAnis[idx] |= (1 << SHIFTC44C33);
					}
				}
			}
		}

		// DenAng
		Add_Property(_PadDenAng,"Loading DEN",_den_File,_den_Constant,_den_Min,_den_range,DENMASK,SHIFTDen);
		if (_Kernel_Type > 1)
		{
			Add_Dip_Azm(_PadDenAng,"Loading DIP & AZM",_dip_or_dx_File,_dip_or_dx_Constant,_dip_Min,_dip_range,DIPMASK,0,_azm_or_dy_File,_azm_or_dy_Constant,_azm_Min,_azm_range,AZMMASK,SHIFTAzm,_degrees_Flag,_dipxdipy_Flag);
		}
	}

	_Valid_Compressed_Earth_Model = 1;
}

void Acoustic_Earth_Model::Find_Min_Max_Q(
		const char* name,
		const char* filename,
		float const_val,
		float fq,
		float dt,
		float& min,
		float& max
		)
{
	Q_Scanner* reader = new Q_Scanner(name,filename,const_val,_swap_Flag,fq,dt);
	_volume_iterator->Iterate(reader);
	min = reader->Get_Min();
	max = reader->Get_Max();
	delete reader;
}

void Acoustic_Earth_Model::Find_Min_Max_DIP_AZM(
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
)
{
	if (_Kernel_Type == 2)
	{
		Dip_Azm_Scanner* reader = new Dip_Azm_Scanner(name,dip_or_dx_File,azm_or_dy_File,dip_or_dx_Constant,azm_or_dy_Constant,_swap_Flag,degrees_Flag,dipxdipy_Flag);
		_volume_iterator->Iterate(reader);
		dip_Min = reader->Get_Dip_Min();
		dip_Max = reader->Get_Dip_Max();
		azm_Min = reader->Get_Azm_Min();
		azm_Max = reader->Get_Azm_Max();
		delete reader;
	}
	else
	{
		dip_Min = 0.0f;
		dip_Max = 0.0f;
		azm_Min = 0.0f;
		azm_Max = 0.0f;
	}
}

void Acoustic_Earth_Model::Find_Min_Max_Single_Attribute(
		const char* name,
		const char* filename,
		float const_val,
		float& min,
		float& max
		)
{
	Single_Attribute_Scanner* reader = new Single_Attribute_Scanner(name,filename,const_val,_swap_Flag);
	_volume_iterator->Iterate(reader);
	min = reader->Get_Min();
	max = reader->Get_Max();
	delete reader;
}

void Acoustic_Earth_Model::Find_Min_Max_VP_DTA_EPS(
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
		)
{
	Vp_Eps_Dta_Scanner* vp_min_max = 0L;
	if (_Kernel_Type == 0)
	{
		vp_min_max = new Vp_Eps_Dta_Scanner(name,vp_File,0L,0L,vp_Constant,0.0f,0.0f,0,_swap_Flag);
	}
	else
	{
		vp_min_max = new Vp_Eps_Dta_Scanner(name,vp_File,dta_File,eps_or_eta_File,vp_Constant,dta_Constant,eps_or_eta_Constant,eta_Flag,_swap_Flag);
	}
	_volume_iterator->Iterate(vp_min_max);
	vpminX = vp_min_max->Get_vp_X_Min();
	vpmaxX = vp_min_max->Get_vp_X_Max();
	vpminZ = vp_min_max->Get_vp_Z_Min();
	vpmaxZ = vp_min_max->Get_vp_Z_Max();
	Epsmin = vp_min_max->Get_Eps_Min();
	Epsmax = vp_min_max->Get_Eps_Max();
	Delmin = vp_min_max->Get_Dta_Min();
	Delmax = vp_min_max->Get_Dta_Max();
	delete vp_min_max;
	_Valid_PropParm_Ranges = 1;
}

void Acoustic_Earth_Model::Add_Vp(
		int* PadVelAnis,	
		const char* name,
		const char* vp_File,
		float vp_Constant,
		float vp_min,
		float vp_range,
		int vel_mask,
		int vel_shift
		)
{
	Single_Attribute_Loader* reader = new Single_Attribute_Loader(PadVelAnis,_icalc,_ucalc,name,vp_File,vp_Constant,vp_min,vp_range,_swap_Flag,vel_mask,vel_shift);
	_volume_iterator->Iterate(reader);
	delete reader;
}

void Acoustic_Earth_Model::Add_Eps_Dta(
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
		)
{
	Eps_Dta_Loader* reader = new Eps_Dta_Loader(PadVelAnis,_icalc,_ucalc,name,eps_or_eta_File,eps_or_eta_Constant,eta_Flag,eps_min,eps_range,eps_mask,eps_shift,dta_File,dta_Constant,dta_min,dta_range,dta_mask,dta_shift,_swap_Flag);
	_volume_iterator->Iterate(reader);
	delete reader;
}

void Acoustic_Earth_Model::Add_Property(
		int* PadDenAng,
		const char* name,
		const char* filename,
		float const_val,
		float min,
		float range,
		int prop_mask,
		int prop_shift
		)
{
	Single_Attribute_Loader* reader = new Single_Attribute_Loader(PadDenAng,_icalc,_ucalc,name,filename,const_val,min,range,_swap_Flag,prop_mask,prop_shift);
	_volume_iterator->Iterate(reader);
	delete reader;
}

void Acoustic_Earth_Model::Add_Dip_Azm(
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
		)
{
	Dip_Azm_Loader* reader = new Dip_Azm_Loader(PadDenAng,_icalc,_ucalc,name,dip_or_dx_File,dip_or_dx_Constant,dip_min,dip_range,dip_mask,dip_shift,azm_or_dy_File,azm_or_dy_Constant,azm_min,azm_range,azm_mask,azm_shift,_swap_Flag,degrees_Flag,dipxdipy_Flag);
	_volume_iterator->Iterate(reader);
        delete reader;
}

void Acoustic_Earth_Model::Add_Q(
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
		)
{
	Q_Loader* reader = new Q_Loader(PadVelAnis,_icalc,_ucalc,name,Q_File,Q_Constant,Q_min,Q_range,_swap_Flag,Q_mask,Q_shift,fq,dt);
	_volume_iterator->Iterate(reader);
        delete reader;
}

void Acoustic_Earth_Model::Destroy_Compressed_Earth_Model()
{
	if (_GPU_HostAlloc)
	{
		if (_PadDenAng != 0L)
		{
			cudaFreeHost(_PadDenAng);
			_PadDenAng = 0L;
		}
	}
	else
	{
		if (_PadVelAnis != 0L)
		{
			free(_PadVelAnis);
			_PadVelAnis = 0L;
		}
		if (_PadDenAng != 0L)
		{
			free(_PadDenAng);
			_PadDenAng = 0L;
		}
		if (_VelAnis != 0L)
		{
			for (int iZ = 0;  iZ < _nz;  ++iZ)
			{
				if (_VelAnis[iZ] != 0l)
				{
					delete [] _VelAnis[iZ];
				}
			}
			delete _VelAnis;
			_VelAnis = 0L;
		}
		if (_DenAng != 0L)
		{
			for (int iZ = 0;  iZ < _nz;  ++iZ)
			{
				if (_DenAng[iZ] != 0l)
				{
					delete [] _DenAng[iZ];
				}
			}
			delete _DenAng;
			_DenAng = 0L;
		}
	}

	_Valid_PropParm_Ranges = 0;
	_Valid_Earth_Model_Parameters = 1;
	_Valid_Compressed_Earth_Model = 0;

	_vp_X_Min = 1e37f;
	_vp_X_Max = -1e37f;
	_vp_Z_Min = 1e37f;
	_vp_Z_Max = -1e37f;
	_den_Min = 1e37f;
	_den_Max = -1e37f;
	_eps_Min = 1e37f;
	_eps_Max = -1e37f;
	_dta_Min = 1e37f;
	_dta_Max = -1e37f;
	_dip_Min = 1e37f;
	_dip_Max = -1e37f;
	_azm_Min = 1e37f;
	_azm_Max = -1e37f;
	_Q_Min = 1e37f;
	_Q_Max = -1e37f;
}

float Acoustic_Earth_Model::Get_Density(
		int xsrc,
		int ysrc,
		int zsrc
		)
{
	unsigned long idx = _icalc->Calculate_Index(xsrc,ysrc,zsrc);
        if (_Kernel_Type == 0)
        {
                int iVelAnis = _PadVelAnis[idx];
                return (float)((iVelAnis >> SHIFTDen) & DENMASK) * _den_binsize + _den_Min;
        }
        else
        {
                int iDenAng = _PadDenAng[idx];
                return (float)((iDenAng >> SHIFTDen) & DENMASK) * _den_binsize + _den_Min;
        }
}

float Acoustic_Earth_Model::Get_Vp(
		int xsrc,
		int ysrc,
		int zsrc
	    )
{
	unsigned long idx = _icalc->Calculate_Index(xsrc,ysrc,zsrc);
        int iVelAnis = _PadVelAnis[idx];
        return (float)(iVelAnis & VELMASK) * _vp_binsize + _vp_Z_Min;
}

float Acoustic_Earth_Model::Get_Bulk_Modulus(
		int xsrc,
		int ysrc,
		int zsrc
		)
{
	float Vp = Get_Vp(xsrc,ysrc,zsrc);
        float Dn = Get_Density(xsrc,ysrc,zsrc);
        return Dn * Vp * Vp;
}

