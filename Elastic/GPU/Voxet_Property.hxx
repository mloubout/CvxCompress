#ifndef CVX_ESDRD_MI_TMJ_VOXET_PROPERTY_HXX
#define CVX_ESDRD_MI_TMJ_VOXET_PROPERTY_HXX

class Voxet_Property
{
public:
	Voxet_Property(const char* moniker, int id);
	~Voxet_Property();

	void Dump();  // dump debug info to stdout

	int Get_ID();
	const char* Get_Moniker();

	void Set_MinMax(double min, double max);
	double Get_Min();
	double Get_Max();

	void Set_Path(const char* path, const char* fullpath);

	// get path as stated in .vo file
	const char* Get_Path();

	// get full path, i.e. dirname of .vo file + path
	const char* Get_Full_Path();

private:
	int _id;
	char* _moniker;
	char* _path;
	char* _fullpath;

	double _min;
	double _max;

	char* _strdup(const char* src);
};

#endif

