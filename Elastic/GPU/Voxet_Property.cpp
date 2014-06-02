#include <stdio.h>
#include <string.h>
#include "Voxet_Property.hxx"

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

void Voxet_Property::_swap_endian(float* v)
{
        int* iv = (int*)v;
        *iv =
                (((*iv) >> 24) & 255) |
                (((*iv) >>  8) & 65280) |
                (((*iv) & 65280) <<  8) |
                (((*iv) & 255) << 24);
}

void Voxet_Property::Get_MinMax_From_File()
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
                        for (int i = 0;  i < nread;  ++i)
                        {
                                _swap_endian(f+i);
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

