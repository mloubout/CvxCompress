#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_SHOT_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_SHOT_HXX

#include "Elastic_Interpolation.hxx"

class Elastic_Modeling_Job;
class Elastic_SEGY_File;
class Elastic_Propagator;
class Elastic_Pipeline;

class Elastic_Shot
{
public:
	Elastic_Shot(int log_level, Elastic_Modeling_Job* job, int souidx, double x, double y, double z);
	~Elastic_Shot();

	double Get_Propagation_Time();

	int Get_Source_Index() {return _souidx;}
	void Set_Source_Index(int souidx){_souidx =souidx;}

	int Get_Ordertime() {return _ordertime;}

	Elastic_Interpolation_t Get_Source_Interpolation_Method() {return _souintrp;}
	void Set_Source_Interpolation_Method(Elastic_Interpolation_t intrp) {_souintrp = intrp;}

	double Get_Source_X() {return _x;}
	double Get_Source_Y() {return _y;}
	double Get_Source_Z() {return _z;}

	void Set_Source_X(double x){_x =x;}
	void Set_Source_Y(double y){_y =y;}
	void Set_Source_Z(double z){_z =z;}

	void Set_Source(double x, double y, double z){_x=x;_y=y;_z=z;}

	double Get_Propagation_Source_X();
	double Get_Propagation_Source_Y();
	double Get_Propagation_Source_Z();

	bool Inject_Source(int discrete_timestep) {return discrete_timestep >= 0 && discrete_timestep < _tsrc ? true : false;}
	double Get_Source_Wavelet_Sample(int discrete_timestep, bool Is_Particle_Velocity)
	{
		if (discrete_timestep >= 0 && discrete_timestep < _tsrc)
		{
			if (Is_Particle_Velocity && _soutype == Source_Type_Force)
			{
				return _stf[discrete_timestep];
			}
			else
			{
				return _stf_int[discrete_timestep];
			}
		}
		else
		{
			return 0.0;
		}
	}

	static const int Source_Type_Unknown = 0;
	static const int Source_Type_Force = 1;
	static const int Source_Type_Velocity = 2;
	static const int Source_Type_Pressure = 3;

	void Set_Source_Type(int soutype) {_soutype = soutype;}
	int Get_Source_Type() {return _soutype;}
	const char* Get_Source_Type_String();
	bool Get_Source_Is_Force() {return _soutype == Source_Type_Force ? true : false;}
	bool Get_Source_Is_Velocity() {return _soutype == Source_Type_Velocity ? true : false;}
	bool Get_Source_Is_Pressure() {return _soutype == Source_Type_Pressure ? true : false;}

	void Set_Amplitudes(double ampl1, double ampl2, double ampl3) {_ampl1 = ampl1; _ampl2 = ampl2; _ampl3 = ampl3;}
	double Get_Amplitude1() {return _ampl1;}
	double Get_Amplitude2() {return _ampl2;}
	double Get_Amplitude3() {return _ampl3;}

	bool Use_Builtin_Source_Wavelet(const char* wavetype, double max_freq, const char* parmfile_path, int line_num);
	bool Read_Source_Wavelet_From_File(const char* wavelet_path, double max_freq, int filter_order);

	void Prepare_Source_Wavelet(double dt, bool debug_output_source_wavelet);

	double Find_Deepest_Receiver();


	void Add_SEGY_File(Elastic_SEGY_File* segy_file);

	void Add_Receiver_Array(int nrec, 
		double* rec_x,
		double* rec_y,
		double* rec_z,
		int* iline,
		int* xline,
		int* trcens);

	Elastic_SEGY_File* Get_SEGY_File(int segy_file_idx);

	int Get_Number_Of_SEGY_Files() {return _num_segy_files;}
	Elastic_SEGY_File* Get_SEGY_File_by_Index(int idx)
	{
		if (idx >= 0 && idx < _num_segy_files)
			return _segy_files[idx];
		else
			return 0L;
	}

	void Allocate_Pinned_Host_Memory(Elastic_Propagator* prop);
	void Free_Pinned_Host_Memory(Elastic_Propagator* prop);

	void Create_Trace_Resample_Buffers(Elastic_Propagator* prop);
	void Free_Trace_Resample_Buffers();

	void Calculate_RX_Locations_And_Results_Size(Elastic_Propagator* prop, int device_id, int& rxloc_size, int& rxres_size, int& num_full_compute);
	void Start_Extract_Receiver_Values_From_Device(
			Elastic_Propagator* prop, 
			Elastic_Pipeline* pipe,
			int device_id, 
			int* block_offsets, 
			int* timesteps, 
			int* num_rx,
			int* flags,
			int num_blocks, 
			float** d_rxloc_block, 
			float* d_rxres,
			float* h_rxres
			);
	void Extract_Receiver_Values_From_Device(
			Elastic_Propagator* prop, 
			Elastic_Pipeline* pipe,
			int device_id, 
			int* block_offsets, 
			int* timesteps, 
			int* num_rx,
			int* flags,
			int num_blocks, 
			float** d_rxloc_block, 
			float* d_rxres,
			float* h_rxres
			);
	void DEMUX_Receiver_Values(
			Elastic_Pipeline* pipe,
			int* block_offsets, 
			int* timesteps, 
			int* num_rx,
			int* flags,
			int num_blocks, 
			float* h_rxres
			);
	void Resample_Receiver_Traces();

	void Write_SEGY_Files();

private:
	int _log_level;
	Elastic_Modeling_Job* _job;

	int _souidx;
	double _x;
	double _y;
	double _z;
	int _soutype;
	Elastic_Interpolation_t _souintrp;
	double _ampl1, _ampl2, _ampl3;
	int _ordertime;

	int _wavetype;
	double _max_freq;
	char* _wavelet_path;
	int _filter_order;

	// from Joe Stefani
	int _tsrc;
	double* _stf;
	double* _stf_int;
	void _src(double dt, double fmax, int type, char* stfname, int* tsrc, double* stf);
	// Ricker wavelet
	void _generate_ricker_wavelet(double dt, double fmax, int* tsrc, double* stf);
	
	Elastic_SEGY_File** _segy_files;
	int _num_segy_files;

        bool _Range_Intersects(int x0, int x1, int x_lo, int x_hi);
        bool _Receiver_Intersects(Elastic_Interpolation_t interpolation_method, int x0, int x1, int y0, int y1, float recx, float recy);
	void _Create_Receiver_Transfer_Buffers(Elastic_Propagator* prop);

        int _totSteps;
        int _nBlks;
        int _num_pipes;
        float** _h_rcv_loc;
        float*** _h_rcv_binned;
	float** _h_pinned_rcv_loc;
        int*** _h_rcv_trcidx;
        int*** _h_rcv_trcflag;
        int* _h_rcv_loc_size_f;

	int _num_traces;
	double* _h_trace_rcv_x;
	double* _h_trace_rcv_y;
	double* _h_trace_rcv_z;
	int* _h_trace_iline;
	int* _h_trace_xline;
	int* _h_trace_trcens;
	double* _h_trace_tshift;
	double* _h_trace_sample_rate;
	float** _h_trace_in;
	float** _h_trace_out;
	int* _h_trace_flag;
	int* _h_trace_idx_in;
	int* _h_trace_idx_in_nn;
	int* _h_trace_idx_out;
	int* _h_trace_iFile;  // iFile no for this trace
	unsigned int* _h_trace_touched;
	//bool* _h_trace_touched;

	// indexed on _num_segy_files
	int* _h_trace_nsamp_in;
	int* _h_trace_nsamp_out;
	int** _h_trace_out_idxM;
	float*** _h_trace_out_sinc_coeffs;

	unsigned long _Comp_RxLoc_Length(float* rxloc);
	unsigned long _Comp_RxRes_Length(float* rxloc, int flags);
};

#endif

