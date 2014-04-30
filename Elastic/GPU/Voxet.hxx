#ifndef CVX_ESDRD_MI_TMJ_VOXET_HXX
#define CVX_ESDRD_MI_TMJ_VOXET_HXX

class Voxet_Property;
class Global_Coordinate_System;

class Voxet
{
public:
	Voxet(int log_level, char* path);
	~Voxet();

	void Dump();  // dump debug info to stdout

	Global_Coordinate_System* Get_Global_Coordinate_System();

	int Get_Number_Of_Properties();
	Voxet_Property* Get_Property(int index);
	Voxet_Property* Get_Property(const char* moniker);

private:
	static const int MAX_PROPERTIES = 2048;
	int _log_level;

	char* _path;
	Global_Coordinate_System* _global_coordinate_system;

	char* _strdup(const char* str);
	void _Read(const char* path);

	int _num_properties;
	Voxet_Property** _properties;
};

#endif

