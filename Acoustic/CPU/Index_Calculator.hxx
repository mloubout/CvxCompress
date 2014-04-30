#ifndef INDEX_CALCULATOR_HXX
#define INDEX_CALCULATOR_HXX

class Index_Calculator
{
	public:
		Index_Calculator(int Target_Platform, int num_fields, int nx, int ny, int nz, int sky_z, int xh, int yh, int zh);
		virtual ~Index_Calculator();

		unsigned long Get_Field_Offset();

		unsigned long Calculate_Index(int iX, int iY, int iZ);
	
	private:
		int _Target_Platform;
		unsigned long _num_fields;
		unsigned long _nx;
		unsigned long _ny;
		unsigned long _nz;
		unsigned long _sky_z;
		unsigned long _xh;
		unsigned long _yh;
		unsigned long _zh;
};

#endif

