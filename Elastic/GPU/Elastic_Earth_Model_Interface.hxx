#include "../../Common/Global_Coordinate_System.hxx"

//
// Base class that defines earth model interface for elastic propagator.
// Elastic code uses this interface to extract and compress a sub volume from the earth model.
//

enum Elastic_Attribute_Index_t
{
	Attr_Idx_Vp = 0,
	Attr_Idx_Vs = 1,
	Attr_Idx_Density = 2,
	Attr_Idx_Q = 3,
	Attr_Idx_Dip = 4,
	Attr_Idx_Azimuth = 5,
	Attr_Idx_Rake = 6,
	Attr_Idx_Delta1 = 7,
	Attr_Idx_Delta2 = 8,
	Attr_Idx_Delta3 = 9,
	Attr_Idx_Epsilon1 = 10,
	Attr_Idx_Epsilon2 = 11,
	Attr_Idx_Gamma1 = 12,
	Attr_Idx_Gamma2 = 13
};

//
// This enum tells the propagator what type of earth model you have supplied.
// This will affect what the propagator does of course.
// 
//
enum Elastic_Earth_Model_Type_t
{
	Elastic_ISO_Model,
	Elastic_VTI_Model,
	Elastic_VOR_Model,
	Elastic_TTI_Model,
	Elastic_TOR_Model
};

class Elastic_Earth_Model_Interface
{
public:
	virtual Global_Coordinate_System* Get_Global_Coordinate_System() = 0;

	virtual Elastic_Earth_Model_Type_t Get_Earth_Model_Type() = 0;

	// Get dimensions of earth model.
	virtual int Get_NX() = 0;
	virtual int Get_NY() = 0;
	virtual int Get_NZ() = 0;

	// Get cell sizes. Can be either in ft or m, but must match whatever was used for the velocity properties.
	virtual float Get_DX() = 0;
	virtual float Get_DY() = 0;
	virtual float Get_DZ() = 0;

	// Indicates that reading consecutive values along a given axis is a fast operation.
	// If the earth model is stored in a file, reading along one axis will be faster than the others.
	// In this case, the object should return true only for the fast axis.
	// If the earth model comes from an in-memory array, reading along any axis is fast,
	// so all methods should return true.
	virtual bool Reading_Along_X_Is_Fast() = 0;
	virtual bool Reading_Along_Y_Is_Fast() = 0;
	virtual bool Reading_Along_Z_Is_Fast() = 0;

	virtual bool Get_Min_Max(Elastic_Attribute_Index_t attr_idx, int x0, int x1, int y0, int y1, int z0, int z1, float& min, float& max) = 0;
	virtual bool Read_Sub_Volume(Elastic_Attribute_Index_t attr_idx, int x0, int x1, int y0, int y1, int z0, int z1, float* vol) = 0;
};

