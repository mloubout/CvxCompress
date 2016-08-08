#ifndef ELASTIC_ORTHORHOMBIC_VARIABLE_WATER_VELOCITY_HXX
#define ELASTIC_ORTHORHOMBIC_VARIABLE_WATER_VELOCITY_HXX

#include <Global_Coordinate_System.hxx>
#include <Voxet.hxx>
#include <Voxet_Property.hxx>
#include <Voxet_Memory_Mapper.hxx>
#include <swapbytes.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <map>
#include <omp.h>

class Variable_Water_Velocity
{
public:
	Variable_Water_Velocity()
	{
		_has_been_initialized = false;
		_xoff = _yoff = _zoff = 0;
		_nx = _ny = _nz = 0;
		_stride_x = _stride_y = _stride_z = 0;
		_wbIdx = 0L;
		_wbIdxMax = 0;
		_VwIdx = -1;
		_ratio = 1.0f;
		_Vwxyzt_voxet = 0L;
		_numVwInterp = 0;
		_VwInterp_start_time = 0;
		_VwInterp = 0;
		_numPIES = 0;
		_PIES_time = 0L;
		_PIES_perturbation = 0L;
	}

	virtual ~Variable_Water_Velocity()
	{
		delete _Vwxyzt_voxet;
		delete [] _wbIdx;
		delete [] _VwInterp_start_time;
		delete [] _VwInterp;
		delete [] _PIES_time;
		delete [] _PIES_perturbation;
	}

	//!
	//! Create voxet object from the .vo file in the given path
	//! Returns integer indicating status.
	//! 0 -> success
	//! 1 -> voxet already exists
	//! 2 -> voxet has no global coordinate system
	int Create_Voxet(const char* path, const char* parmfile_path, int line_num, int log_level)
	{
		if (_Vwxyzt_voxet != 0L)
		{
			printf("%s (line %d) : Error - Vwxyzt_Voxet cannot appear more than once in file.\n",parmfile_path,line_num);
			return 1;
		}
		else
		{
			_Vwxyzt_voxet = new Voxet(log_level,(char*)path);
			if (_Vwxyzt_voxet->Get_Global_Coordinate_System() == 0L)
			{
				printf("%s (line %d) : Error - Voxet contains no global coordinate system information.\n",parmfile_path,line_num);
				return 2;
			}
			else
			{
				if (log_level >= 3) printf("Vwxyzt :: water volumes will be read from %s.\n",path);
				return 0;
			}
		}
	}

	//!
	//! Set UVW transpose.
	//! Returns integer indicating status.
	//! 0 -> success
	//! 1 -> missing voxet
	//! 2 -> Set_Transpose failed
	int Set_Transpose(const char* transpose, const char* parmfile_path, int line_num, int log_level)
	{
		if (_Vwxyzt_voxet == 0L)
		{
			printf("%s (line %d): Error - Vwxyzt_Transpose UVW cannot appear before Vwxyzt_Voxet.\n",parmfile_path,line_num);
			return 1;
		}
		else
		{
			if (_Vwxyzt_voxet->Get_Global_Coordinate_System()->Set_Transpose(transpose))
			{
				if (log_level > 3) printf("Vwxyzt :: Transpose set to uvw -> %s\n",transpose);
				return 0;
			}
			else
			{
				printf("%s (line %d) : Error - Vwxyzt :: Set transpose to uvw -> %s failed.\n",parmfile_path,line_num,transpose);
				return 1;
			}
		}
	}

	//!
	//! Add water volume.
	//! Returns integer indicating status.
	//! 0 -> success
	//! 1 -> missing voxet
	//! 2 -> unknown attribute
	int Add_Water_Volume(const char* moniker, time_t start_time, const char* parmfile_path, int line_num, int log_level)
	{
		if (_Vwxyzt_voxet == 0)
		{
			printf("%s (line %d): Error - PROPERTY cannot appear before Vwxyzt_Voxet.\n",parmfile_path,line_num);
			return 1;
		}
		else
		{
			Voxet_Property* prop = _Vwxyzt_voxet->Get_Property_By_Moniker(moniker);
			if (prop == 0L)
			{
				printf("%s (line %d): Error - Voxet file does not have property %s.\n",parmfile_path,line_num,moniker);
				printf("The voxet has the following properties:\n");
				for (int i = 0;  i < _Vwxyzt_voxet->Get_Number_Of_Properties();  ++i)
				{
					printf("  %s\n",_Vwxyzt_voxet->Get_Property_By_Index(i)->Get_Moniker());
				}
				return 2;
			}
			else
			{
				_Vwxyzt_volumes[start_time] = moniker;
				if (log_level >= 3)
				{
					printf("Vwxyzt interpolation volume will be read from %s, starting at time %s",prop->Get_Full_Path(),asctime(localtime(&start_time)));
				}
				return 0;
			}
		}
	}

	//!
	//! Add PIES file path.
	//! Returns integer status indicator:
	//! 0->success
	//! 1->empty filename
	//! 2->unable to open file for reading
	int Add_PIES_File(const char* filepath, const char* parmfile_path, int line_num, int log_level)
	{
		if (strlen(filepath) > 0)
		{
			std::fstream fs;
			fs.open(filepath,std::fstream::in);
			if (fs.good())
			{
				fs.close();
				_PIES_File = (std::string)filepath;
				if (log_level >= 3) printf("PIES data will be read from %s.\n",_PIES_File.c_str());
				return 0;
			}
			else
			{
				printf("%s (line %d): Unable to open PIES file %d for reading.\n",parmfile_path,line_num,_PIES_File.c_str());
				return 2;
			}
		}
		else
		{
			printf("%s (line %d): PIES file path cannot be an empty string.\n",parmfile_path,line_num);
			return 1;
		}
	}
	
	bool Has_Been_Initialized()
	{
		return _has_been_initialized;
	}

	bool Ready_For_Initialization()
	{
		return (_Vwxyzt_voxet != 0L) ? true : false;
	}

	void Initialize(
		Voxet_Memory_Mapper* mapper,		// we will share memory mapper object with the earth model
		Global_Coordinate_System* gcs,		// earth model gcs
		float* Vs				// Vs field memory mapped. NB! Need byte swap
		)
	{
		assert(mapper != 0L);
		assert(gcs != 0L);
		assert(Vs != 0L);
		assert(_Vwxyzt_volumes.size() > 0);
		Global_Coordinate_System* Vw_gcs = _Vwxyzt_voxet->Get_Global_Coordinate_System();
		_nx = Vw_gcs->Get_NX();
		_ny = Vw_gcs->Get_NY();
		_nz = Vw_gcs->Get_NZ();
		_numVwInterp = _Vwxyzt_volumes.size();
		_VwInterp_start_time = new time_t[_numVwInterp];
		_VwInterp = new float*[_numVwInterp];
		_VwInterp_min = 1.0f;
		_VwInterp_max = 0.0f;
		int iVw = 0;
		int num_threads;
#pragma omp parallel
		{
			num_threads = omp_get_num_threads();
		}
		for (std::map<time_t,std::string>::iterator it = _Vwxyzt_volumes.begin();  it != _Vwxyzt_volumes.end();  ++it)
		{
			_VwInterp_start_time[iVw] = it->first;
			Voxet_Property* prop = _Vwxyzt_voxet->Get_Property_By_Moniker(it->second.c_str());
			float* p = mapper->Get_Memory_Mapped_File(prop->Get_Full_Path());
			_VwInterp[iVw] = p;
			long nn = (long)_nx * (long)_ny * (long)_nz;
#pragma omp parallel for
			for (long iThr = 0;  iThr < num_threads;  ++iThr)
			{
				long i0 = iThr * nn / (long)num_threads;
				long i1 = (iThr+1) * nn / (long)num_threads;
				float min = p[0];
				float max = min;
				for (long i = i0;  i < i1;  ++i)
				{
					float v = p[i];
					swap4bytes((int*)&v,1);
					if (v < min) min = v;
					if (v > max) max = v;
				}
#pragma omp critical
				{
					if (_VwInterp_min > _VwInterp_max)
					{
						_VwInterp_min = min;
						_VwInterp_max = max;
					}
					else
					{
						if (min < _VwInterp_min) _VwInterp_min = min;
						if (max > _VwInterp_max) _VwInterp_max = max;
					}
				}
			}
			++iVw;
		}
		// info needed to traverse em fields
		int em_nx = gcs->Get_NX();
		int em_ny = gcs->Get_NY();
		int em_nz = gcs->Get_NZ();
		int em_stride_x, em_stride_y, em_stride_z;
		gcs->Convert_Local_Index_To_Transposed_Index(1,gcs->Get_NU(),gcs->Get_NU()*gcs->Get_NV(),em_stride_x,em_stride_y,em_stride_z);
		// info needed to traverse water volumes and merge those with em fields
		double g0,g1,g2;
		Vw_gcs->Convert_Local_To_Global(0,0,0,g0,g1,g2);
		gcs->Convert_Global_To_Transposed_Index(g0,g1,g2,_xoff,_yoff,_zoff);
		Vw_gcs->Convert_Local_Index_To_Transposed_Index(1,Vw_gcs->Get_NU(),Vw_gcs->Get_NU()*Vw_gcs->Get_NV(),_stride_x,_stride_y,_stride_z);
		_wbIdxMax = 0;
		_wbIdx = new int[_nx*_ny];
#pragma omp parallel for
		for (int y = 0;  y < _ny;  ++y)
		{
			for (int x = 0;  x < _nx;  ++x)
			{
				int em_x = x + _xoff;
				int em_y = y + _yoff;
				if (em_x < 0 || em_x >= em_nx || em_y < 0 || em_y >= em_ny)
				{
					_wbIdx[y*_nx+x] = 0;
				}
				else
				{
					bool done = false;
					for (int em_z = 0;  em_z < em_nz && !done;  ++em_z)
					{
						float val = Vs[(long)em_z*(long)em_stride_z+(long)em_y*(long)em_stride_y+(long)em_x*(long)em_stride_x];
						swap4bytes((int*)&val,1);
						if (val > 0.0f)
						{
							done = true;
							int idx = y*gcs->Get_NX() + x;
							_wbIdx[y*_nx+x] = em_z;
							if (em_z > _wbIdxMax) _wbIdxMax = em_z;
						}
					}
				}
			}
		}
		printf("_wbIdxMax = %d\n",_wbIdxMax);
		// extrapolate water depths in the same manner the earth model is extrapolated
		// x first
#pragma omp parallel for
		for (int y = 0;  y < _ny;  ++y)
		{
			int xl = 0;
			for (;  xl < _nx && _wbIdx[y*_nx+xl] == 0;  ++xl);
			if (xl > 0 && xl < _nx)
			{
				int wb = _wbIdx[y*_nx+xl];
				for (int x = 0;  x < xl;  ++x) _wbIdx[y*_nx+x] = wb;
			}
			int xr = _nx-1;
			for (;  xr >= 0 && _wbIdx[y*_nx+xr] == 0;  --xr);
			if (xr < (_nx-1) && xr >= 0)
			{
				int wb = _wbIdx[y*_nx+xr];
				for (int x = xr+1;  x < _nx;  ++x) _wbIdx[y*_nx+x] = wb;
			}
		}
		// then y
#pragma omp parallel for
		for (int x = 0;  x < _nx;  ++x)
		{
			int yl = 0;
			for (;  yl < _ny && _wbIdx[yl*_nx+x] == 0;  ++yl);
			if (yl > 0 && yl < _ny)
			{
				int wb = _wbIdx[yl*_nx+x];
				for (int y = 0;  y < yl;  ++y) _wbIdx[y*_nx+x] = wb;
			}
			int yr = _ny-1;
			for (;  yr >= 0 && _wbIdx[yr*_nx+x] == 0;  --yr);
			if (yr < (_ny-1) && yr >= 0)
			{
				int wb = _wbIdx[yr*_nx+x];
				for (int y = yr+1;  y < _ny;  ++y) _wbIdx[y*_nx+x] = wb;
			}
		}
		// verify that no cell have wbIdx == 0
		for (int y = 0;  y < _ny;  ++y)
		{
			for (int x = 0;  x < _nx;  ++x)
			{
				assert(_wbIdx[y*_nx+x] > 0);
			}
		}
		// read in PIES data
		std::map<time_t,float> pies_data;
		std::fstream fs;
		fs.open(_PIES_File.c_str(),std::fstream::in);
		assert(fs.good());
		float min_PIES = 1.0f;
		float max_PIES = 0.0f;
		while (!fs.eof())
		{
			char str[4096];
			fs.getline(str,4096);
			if (strlen(str) > 0 && str[0] != 'H')
			{
				char date_str[128];
				char time_str[128];
				double fractional_days;
				double PIES;
				if (sscanf(str,"%s %s %lf %lf",date_str,time_str,&fractional_days,&PIES) == 4)
				{
					char date_time_str[256];
					sprintf(date_time_str,"%s %s",date_str,time_str);
					struct tm tms;
					strptime(date_time_str,"%d-%b-%y %H:%M:%S",&tms);
					time_t rec_time = mktime(&tms);
					char* s = asctime(localtime(&rec_time));
					s[strlen(s)-1] = 0; // remove trailing newline character
					//printf("%s %f %.6f\n",s,fractional_days,PIES);
					pies_data[rec_time] = PIES;
					if (min_PIES > max_PIES)
					{
						min_PIES = max_PIES = PIES;
					}
					else
					{
						if (PIES < min_PIES) min_PIES = PIES;
						if (PIES > max_PIES) max_PIES = PIES;
					}
				}
			}
		}
		fs.close();
		assert(pies_data.size() > 0);
		_numPIES = pies_data.size();
		_PIES_time = new time_t[_numPIES];
		_PIES_perturbation = new float[_numPIES];
		int pie_cnt = 0;
		for (std::map<time_t,float>::iterator it = pies_data.begin();  it != pies_data.end();  ++it)
		{
			_PIES_time[pie_cnt] = it->first;
			_PIES_perturbation[pie_cnt] = it->second;
			++pie_cnt;
		}
		// include PIES data in min,max
		_VwInterp_min += min_PIES;
		_VwInterp_max += max_PIES;
		_has_been_initialized = true;
	}

	inline int Get_Maximum_Water_Depth() {_wbIdxMax;}

	void Set_Shot_Time(time_t shot_time)
	{
		if (_numVwInterp == 1)
		{
			_VwIdx = -1;
			_ratio = 1.0f;
		}
		else
		{
			if (shot_time < _VwInterp_start_time[0])
			{
				_VwIdx = 0;
				_ratio = 0.0f;
			}
			else if (shot_time >= _VwInterp_start_time[_numVwInterp-1])
			{
				_VwIdx = _numVwInterp-2;
				_ratio = 1.0f;
			}
			else
			{
				for (int i = 0;  i < _numVwInterp-1;  ++i)
				{
					if (shot_time >= _VwInterp_start_time[i] && shot_time < _VwInterp_start_time[i+1])
					{
						_VwIdx = i;
						_ratio = difftime(shot_time,_VwInterp_start_time[i]) / difftime(_VwInterp_start_time[i+1],_VwInterp_start_time[i]);
					}
				}
			}
		}
		printf("_VwIdx=%d, _ratio=%e\n",_VwIdx,_ratio);
		if (_numPIES <= 0)
		{
			_VwCurPIES = 0.0f;
		}
		else if (shot_time <= _PIES_time[0])
		{
			_VwCurPIES = _PIES_perturbation[0];
		}
		else if (shot_time >= _PIES_time[_numPIES-1])
		{
			_VwCurPIES = _PIES_perturbation[_numPIES-1];
		}
		else
		{
			for (int i = 0;  i < (_numPIES-1);  ++i)
			{
				if (shot_time >= _PIES_time[i] && shot_time < _PIES_time[i+1])
				{
					double PIES_ratio = difftime(shot_time,_PIES_time[i]) / difftime(_PIES_time[i+1],_PIES_time[i]);
					assert(PIES_ratio >= 0.0 && PIES_ratio <= 1.0);
					_VwCurPIES = PIES_ratio * _PIES_perturbation[i+1] + (1.0 - PIES_ratio) * _PIES_perturbation[i];
				}
			}
		}
	}

	inline float Compute_Velocity_Increment(
			int em_x,
			int em_y,
			int em_z
			)
	{
		int x = em_x - _xoff;
		int y = em_y - _yoff;
		int z = em_z - _zoff;
		int inside_range = (x >= 0 && x < _nx && y >= 0 && y < _ny) ? true : false;
		if (!inside_range || z >= _wbIdx[y*_nx+x])
		{
			return 0.0f;
		}
		else
		{
			long idx = (long)z * (long)_stride_z + (long)y * (long)_stride_y + (long)x * (long)_stride_x;
			if (_VwIdx < 0)
			{
				float v = _VwInterp[0][idx];
				swap4bytes((int*)&v,1);
				return v + _VwCurPIES;
			}
			else
			{
				float v0 = _VwInterp[_VwIdx][idx];
				float v1 = _VwInterp[_VwIdx+1][idx];
				swap4bytes((int*)&v0,1);
				swap4bytes((int*)&v1,1);
				return (v1 * _ratio + v0 * (1.0 - _ratio)) + _VwCurPIES;
			}
		}
	}

	int Get_Water_Bottom(int em_x, int em_y)
	{
		int x = em_x - _xoff;
		int y = em_y - _yoff;
		assert(x >= 0 && x < _nx && y >= 0 && y < _ny);
		return _wbIdx[y*_nx+x];
	}

	float Get_Min() {return _VwInterp_min;}
	float Get_Max() {return _VwInterp_max;}

private:
	bool _has_been_initialized;
	Voxet* _Vwxyzt_voxet;
	std::map<time_t,std::string> _Vwxyzt_volumes;
	int _xoff,_yoff,_zoff;   		// origin of water volumes falls at this local index in the earth model
	int _nx, _ny, _nz;			// water volumes dimensions
	int _stride_x,_stride_y,_stride_z;	// strides for memory mapped water volumes

	int* _wbIdx;
	int _wbIdxMax;
	int _VwIdx;
	float _ratio;
	float _VwCurPIES;

	int _numVwInterp;
	time_t* _VwInterp_start_time;
	float** _VwInterp;
	float _VwInterp_min,_VwInterp_max;

	std::string _PIES_File;
	int _numPIES;
	time_t* _PIES_time;
	float* _PIES_perturbation;
};
#endif
