#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <swapbytes.h>
#include <Voxet_Property.hxx>
#include <Voxet_Memory_Mapper.hxx>

Voxet_Property::Voxet_Property(int logLevel, const char* moniker, int id)
{
	_log_level = logLevel;
	_moniker = _strdup(moniker);
	_id = id;
	_path = 0L;
	_fullpath = 0L;
	_has_min_max = false;
	_min = 0.0;
	_max = 0.0;
	_min_max_from_scan = false;
}

Voxet_Property::~Voxet_Property()
{
	if (_moniker != 0L) delete [] _moniker;
	if (_path != 0L) delete [] _path;
	if (_fullpath != 0L) delete [] _fullpath;
}

int Voxet_Property::Get_ID()
{
	return _id;
}

const char* Voxet_Property::Get_Moniker()
{
	return _moniker;
}

void Voxet_Property::Set_MinMax(double min, double max)
{
	_has_min_max = true;
	_min = min;
	_max = max;
}

double Voxet_Property::Get_Min()
{
	return _min;
}

double Voxet_Property::Get_Max()
{
	return _max;
}

void Voxet_Property::Set_Path(const char* path, const char* fullpath)
{
	_path = _strdup(path);
	_fullpath = _strdup(fullpath);
}

const char* Voxet_Property::Get_Path()
{
	return _path;
}

const char* Voxet_Property::Get_Full_Path()
{
	return _fullpath;
}

void Voxet_Property::Get_MinMax_From_File(Voxet_Memory_Mapper* mapper)
{
	_min_max_from_scan = true;
	if (mapper != 0L)
	{
		float* memfile = mapper->Get_Memory_Mapped_File((std::string)_fullpath);
		if (_log_level >= 3)
		{
			printf("Scanning %s for [MIN,MAX]...",_fullpath);
			fflush(stdout);
		}
		int num_threads;
#pragma omp parallel
		{
			num_threads = omp_get_num_threads();
		}
		size_t nn = mapper->Get_Length((std::string)_fullpath) / sizeof(float);
		float min = 1e38f;
		float max = -1e38f;
#pragma omp parallel for schedule(static,1)
		for (size_t thr = 0;  thr < num_threads;  ++thr)
		{
			size_t beg = (thr * nn) / num_threads;
			size_t end = ((thr+1)* nn) / num_threads;
			float loc_min = 1e38f;
			float loc_max = -1e38f;
			for (size_t idx = beg;  idx < end;  ++idx)
			{
				float v = memfile[idx];
				swap4bytes((int*)&v,1);
				if (v < loc_min) loc_min = v;
				if (v > loc_max) loc_max = v;
			}
#pragma omp critical
			{
				if (loc_min < min) min = loc_min;
				if (loc_max > max) max = loc_max;
			}
		}
		Set_MinMax(min,max);
		printf(" MIN=%e, MAX=%e\n",min,max);
	}
	else
	{
		FILE* fp = fopen(_fullpath, "rb");
		if (fp != 0L)
		{
			if (_log_level >= 3)
			{
				printf("Scanning %s for [MIN,MAX]...",_fullpath);
				fflush(stdout);
			}
			float min = 1e38f;
			float max = -1e38f;
			float f[1024];
			for (
					size_t nread = fread(f, sizeof(float), 1024, fp);
					nread > 0;
					nread = fread(f, sizeof(float), 1024, fp))
			{
				swap4bytes((int*)f, nread);
				for (size_t i = 0;  i < nread;  ++i)
				{
					float v = f[i];
					if (v < min) min = v;
					if (v > max) max = v;
				}
			}
			fclose(fp);
			Set_MinMax(min,max);
			printf(" MIN=%e, MAX=%e\n",min,max);
		}
	}
}

char* Voxet_Property::_strdup(const char* src)
{
	if (src != 0L)
	{
		int n = strlen(src);
		if (n > 0)
		{
			char* buf = new char[n+1];
			for (int i = 0;  i < n;  ++i) buf[i] = src[i];
			buf[n] = 0;
			return buf;
		}
	}
	return 0L;
}

void Voxet_Property::Dump()
{
	printf("Voxet_Property instance %p\n",this);
	printf("_id=%d, _moniker=%s, _path=%s, _fullpath=%s\n",_id,_moniker,_path!=0L?_path:"<nil>",_fullpath!=0L?_fullpath:"<nil>");
}

