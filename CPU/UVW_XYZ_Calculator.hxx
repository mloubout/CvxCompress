#ifndef UVW_XYZ_CALCULATOR_HXX
#define UVW_XYZ_CALCULATOR_HXX

class UVW_XYZ_Calculator
{
public:
	UVW_XYZ_Calculator(int fast_Axis, int med_Axis, int slow_Axis);
	virtual ~UVW_XYZ_Calculator();

	int Is_Valid();

	void Compute_XYZ_From_UVW(int u, int v, int w, int& x, int& y, int& z);
	void Compute_UVW_From_XYZ(int x, int y, int z, int& u, int& v, int& w);

private:
	int _fast_Axis;
	int _med_Axis;
	int _slow_Axis;	
};

#endif

