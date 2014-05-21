#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_SEGY_FILE_RECEIVER_RANGE_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_SEGY_FILE_RECEIVER_RANGE_HXX

class Elastic_SEGY_File_Receiver_Range
{
public:
	Elastic_SEGY_File_Receiver_Range(int range_idx);
	~Elastic_SEGY_File_Receiver_Range();

	bool Is_Valid() {return _has_x && _has_y && _has_z;}

	int Get_Range_Idx() {return _range_idx;}

	void Add_X(double start, double end, double interval);
	void Add_Y(double start, double end, double interval);
	void Add_Z(double start, double end, double interval);

	int Compute_Receiver_Locations(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z,
		int*& iline,
		int*& xline,
		int*& trcens  // trace number within ensemble
		);

private:
	int _range_idx;
	bool _has_x, _has_y, _has_z;
	double _x_start, _x_end, _x_interval;
	double _y_start, _y_end, _y_interval;
	double _z_start, _z_end, _z_interval;
};

#endif

