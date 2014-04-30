#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_SHOT_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_SHOT_HXX

class Elastic_Modeling_Job;
class Elastic_SEGY_File;
class Elastic_Propagator;

class Elastic_Shot
{
public:
	Elastic_Shot(int log_level, Elastic_Modeling_Job* job, int souidx, double x, double y, double z);
	~Elastic_Shot();

	bool Set_Propagation_Time(double PropTime, const char* parmfile_path, int line_num);
	double Get_Propagation_Time() {return _propagation_time;}

	int Get_Source_Index() {return _souidx;}

	int Get_Ordertime() {return _ordertime;}

	double Get_Source_X() {return _x;}
	double Get_Source_Y() {return _y;}
	double Get_Source_Z() {return _z;}

	double Get_Propagation_Source_X();
	double Get_Propagation_Source_Y();
	double Get_Propagation_Source_Z();

	bool Inject_Source(int discrete_timestep) {return discrete_timestep >= 0 && discrete_timestep < _tsrc ? true : false;}
	double Get_Source_Wavelet_Sample(int discrete_timestep) {return discrete_timestep >= 0 && discrete_timestep < _tsrc ? _stf[discrete_timestep] : 0.0f;}

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

	void Prepare_Source_Wavelet(double dt);

	void Add_SEGY_File(Elastic_SEGY_File* segy_file);
	Elastic_SEGY_File* Get_SEGY_File(int segy_file_idx);
	void Create_Receiver_Transfer_Buffers(Elastic_Propagator* prop);

	void Shift_Receiver_Transfer_Buffers();

private:
	int _log_level;
	Elastic_Modeling_Job* _job;

	int _souidx;
	double _x;
	double _y;
	double _z;
	int _soutype;
	double _ampl1, _ampl2, _ampl3;
	int _ordertime;

	int _wavetype;
	double _max_freq;

	// from Joe Stefani
	int _tsrc;
	double* _stf;
	void _src(double dt, double fmax, int type, char* stfname, int* tsrc, double* stf);

	double _propagation_time;

	Elastic_SEGY_File** _segy_files;
	int _num_segy_files;
};

#endif

