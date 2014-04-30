#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include "Voxet_Property.hxx"
#include "Voxet.hxx"
#include "Global_Coordinate_System.hxx"

Voxet::Voxet(int log_level, char* path)
{
	_log_level = log_level;
	_path = strdup(path);
	_num_properties = 0;
	_properties = new Voxet_Property*[MAX_PROPERTIES];
	for (int i = 0;  i < MAX_PROPERTIES;  ++i) _properties[i] = 0L;
	_Read(_path);
}

Voxet::~Voxet()
{
	if (_path != 0L) delete [] _path;
	if (_global_coordinate_system != 0L) delete _global_coordinate_system;
	for (int i = 0;  i < MAX_PROPERTIES;  ++i)
		if (_properties[i] != 0L) delete _properties[i];
	delete [] _properties;
}

void Voxet::_Read(const char* voxet_path)
{
	_global_coordinate_system = new Global_Coordinate_System(voxet_path);

	// scan for properties
	_num_properties = 0;
	FILE* fp = fopen(voxet_path, "r");
	if (fp != 0L)
	{
		char str[4096];
		for (char* s = fgets(str,4096,fp);  s != 0L;  s = fgets(str,4096,fp))
		{
			int id;
			char moniker[4096];
			if (sscanf(s, "PROPERTY %d %s", &id, moniker) == 2)
			{
				if (_log_level > 3) printf("PROPERTY %d %s\n",id,moniker);
				_properties[_num_properties] = new Voxet_Property(moniker,id);
				++_num_properties;
			}
			char prop_path[4096];
			if (sscanf(s, "PROP_FILE %d %s", &id, prop_path) == 2)
			{
				char scratch[4096];
				char fullpath[4096];
				sprintf(scratch, "%s", voxet_path);  // dirname() overwrites its input string
				sprintf(fullpath, "%s/%s", dirname(scratch),prop_path);
				if (_log_level > 3) printf("PROP_FILE %d %s\n",id,fullpath);
				Voxet_Property* prop = Get_Property(id);
				if (prop != 0L)
				{
					prop->Set_Path(prop_path,fullpath);
				}
			}
		}
		fclose(fp);
	}
}

Global_Coordinate_System* Voxet::Get_Global_Coordinate_System()
{
	return _global_coordinate_system;
}

int Voxet::Get_Number_Of_Properties()
{
	return _num_properties;
}

Voxet_Property* Voxet::Get_Property(int id)
{
	for (int i = 0;  i < _num_properties;  ++i)
	{
		if (_properties[i]->Get_ID() == id)
		{
			return _properties[i];
		}
	}
	return 0L;
}

Voxet_Property* Voxet::Get_Property(const char* moniker)
{
	for (int i = 0;  i < _num_properties;  ++i)
        {
                if (strcmp(_properties[i]->Get_Moniker(), moniker) == 0)
                {
                        return _properties[i];
                }
        }
        return 0L;
}

char* Voxet::_strdup(const char* src)
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

void Voxet::Dump()
{
	printf("Voxet instance %p\n",this);
	if (_global_coordinate_system != 0L) _global_coordinate_system->Dump();
	for (int i = 0;  i < _num_properties;  ++i)
	{
		_properties[i]->Dump();
	}
}

