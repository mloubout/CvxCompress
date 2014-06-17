#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_SEGY_FILE_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_SEGY_FILE_HXX

#include "Elastic_Interpolation.hxx"

class Elastic_SEGY_File_Receiver_Range;
class Elastic_Propagator;

class Elastic_SEGY_File
{
public:
	Elastic_SEGY_File(
		int fileidx,
		const char* base_filename,
		double sample_rate,
		double tshift,
		double reclen,
		bool do_P,
		bool do_Vx,
		bool do_Vy,
		bool do_Vz
		);
	~Elastic_SEGY_File();

	bool Is_Valid() {return _Is_Valid;}

	int Get_File_Index() {return _fileidx;}

	Elastic_Interpolation_t Get_Interpolation_Method() {return _interpolation_method;}
	void Set_Interpolation_Method(Elastic_Interpolation_t interpolation_method) {_interpolation_method = interpolation_method;}
	
	double Get_Sample_Rate() {return _sample_rate;}
	double Get_Timeshift() {return _tshift;}
	double Get_Record_Length() {return _reclen;}

        int Get_Selection_Flags();

	void Add_Receiver_Range_X(
		int range_idx,
		double start,
		double end,
		double interval
		);
	void Add_Receiver_Range_Y(
		int range_idx,
		double start,
		double end,
		double interval
		);
	void Add_Receiver_Range_Z(
		int range_idx,
		double start,
		double end,
		double interval
		);

	int Compute_Receiver_Locations(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z
		);
	int Compute_Receiver_Locations(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z,
		int*& iline,
		int*& xline,
		int*& trcens
		);
	
	void Write_SEGY_File(
		float** traces,
		double srcx,
		double srcy,
		double srcz,
		double* recx,
		double* recy,
		double* recz,
		int* iline,
		int* xline,
		int* trcens,
		int num_traces,
		int nsamp,
		int flag
		);

private:
	bool _Is_Valid;
	int _fileidx;
	char* _base_filename;
	double _sample_rate;
	double _tshift;
	double _reclen;
	bool _do_P;
	bool _do_Vx;
	bool _do_Vy;
	bool _do_Vz;
	Elastic_Interpolation_t _interpolation_method;

	void swap2bytes(short *i2, int n);
	void swap4bytes(int *i4, int n);

	Elastic_SEGY_File_Receiver_Range** _rcv_ranges;
	int _num_rcv_ranges;

	Elastic_SEGY_File_Receiver_Range* _Get_Receiver_Range(int range_idx);
};

#endif

