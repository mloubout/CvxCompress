#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_MODELING_JOB_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_MODELING_JOB_HXX

class Voxet;
class Voxet_Property;
class Elastic_Propagator;
class Elastic_Shot;

class Elastic_Modeling_Job
{
public:
	Elastic_Modeling_Job(int log_level, const char* parmfile_path);
	~Elastic_Modeling_Job();

	bool Is_Valid() {return _Is_Valid;}
	int Get_Log_Level() {return _log_level;}

	//
	// Get dimension and origin of propagation volume.
	// Origin can be negative if any of the axes are extended.
	//

	int Get_Propagation_NX();
	int Get_Propagation_NY();
	int Get_Propagation_NZ();

	int Get_Propagation_X0();
	int Get_Propagation_Y0();
	int Get_Propagation_Z0();

	double Get_DX();
	double Get_DY();
	double Get_DZ();

	bool Freesurface_Enabled();	
	bool Source_Ghost_Enabled();
	bool Receiver_Ghost_Enabled();

	bool Use_Isotropic_Sphere_During_Source_Injection() {return _use_isotropic_sphere_during_source_injection;}

	// z index in propagation volume corresponding to sea surface level.
	// this will be zero if freesurface is enabled and a positive value otherwise.

	float Get_Vpvert_Avg_Top();
	float Get_Vpvert_Avg_Bot();

	int Get_NABC_SDX();
	int Get_NABC_SDY();
	int Get_NABC_TOP();
	int Get_NABC_BOT();
	bool Get_NABC_SDX_Extend();
	bool Get_NABC_SDY_Extend();
	bool Get_NABC_TOP_Extend();
	bool Get_NABC_BOT_Extend();

	float Get_Courant_Factor() {return _Courant_Factor;}

	void Add_Shot(Elastic_Shot* shot);
	Elastic_Shot* Get_Shot(int souidx);
	Elastic_Shot* Get_Shot_By_Index(int idx);
	int Get_Number_Of_Shots() {return _num_shots;}

	float Get_FQ() {return _fq;}

	int Get_Number_Of_Earth_Model_Attributes();
	int Get_Earth_Model_Attribute_Index(const char* moniker);
	const char* Get_Earth_Model_Attribute_Moniker(int attr_idx);

	float Get_Earth_Model_Attribute(int attr_idx, int ix, int iy, int iz, bool& error);
	void Set_Earth_Model_Attribute(int attr_idx, int ix, int iy, int iz, float new_value, bool& error);

	float Get_Earth_Model_Attribute_Min(int attr_idx, bool& error);
	float Get_Earth_Model_Attribute_Max(int attr_idx, bool& error);
	float Get_Earth_Model_Attribute_Range(int attr_idx, bool& error);

	float Get_Earth_Model_Attribute(int attr_idx, int ix, int iy, int iz);
	float Get_Earth_Model_Attribute_Min(int attr_idx);
	float Get_Earth_Model_Attribute_Max(int attr_idx);
	float Get_Earth_Model_Attribute_Range(int attr_idx);

	float Get_IsoOrEarth_Model_Attribute_Min(int attr_idx, bool isosphere);
	float Get_IsoOrEarth_Model_Attribute_Range(int attr_idx, bool isosphere);

	void Write_Earth_Model_Attribute_XZ_Slice(const char* path, int attr_idx, int iy);
	void Write_Earth_Model_XZ_Slice(const char* path, int iy);

	void Write_Earth_Model_Attribute_XY_Slice(const char* path, int attr_idx, int iz);
	void Write_Earth_Model_XY_Slice(const char* path, int iz);

	// as name implies, a hack to test what happens if we manipulate Qp to mute Scholte waves
	// traveling along the sea floor.
	void HACK_Mute_Sea_Floor();

	//
	// wf_type : 0->Vx, 1->Vy, 2->Vz, 3->P
	//
	void Write_XZ_Slice(const char* path, int wf_type, int iy);
	void Write_XY_Slice(const char* path, int wf_type, int iz);

	Voxet* Get_Voxet() {return _voxet;}

	int Get_Number_Of_GPU_Pipes();
	void Set_Number_Of_GPU_Pipes(int num_pipes);

	int Get_Steps_Per_GPU();
	void Set_Steps_Per_GPU(int num_timesteps);

	const int* Get_GPU_Devices();
	void Set_GPU_Devices(const int* device_ids, int num_devices);
	int Get_Number_Of_GPU_Devices();

	static const int Attr_Idx_Vp = 0;
	static const int Attr_Idx_Vs = 1;
	static const int Attr_Idx_Density = 2;
	static const int Attr_Idx_Q = 3;
	static const int Attr_Idx_Dip = 4;
	static const int Attr_Idx_Azimuth = 5;
	static const int Attr_Idx_Rake = 6;
	static const int Attr_Idx_Delta1 = 7;
	static const int Attr_Idx_Delta2 = 8;
	static const int Attr_Idx_Delta3 = 9;
	static const int Attr_Idx_Epsilon1 = 10;
	static const int Attr_Idx_Epsilon2 = 11;
	static const int Attr_Idx_Gamma1 = 12;
	static const int Attr_Idx_Gamma2 = 13;

private:
	friend class Elastic_Propagator;
	Elastic_Propagator* _propagator;

	bool _Is_Valid;
	int _log_level;

	int _prop_nx;
	int _prop_ny;
	int _prop_nz;
	
	int _prop_x0;
	int _prop_y0;
	int _prop_z0;

	Voxet* _voxet;

	int _num_em_props;
	Voxet_Property** _props;
	float* _const_vals;

	char** _pck_moniker;
	int* _pck_mask;
	int* _pck_shft;
	int* _pck_widx;
	float* _pck_min;
	float* _pck_max;
	float* _pck_range;
	float* _pck_iso;

	bool _use_isotropic_sphere_during_source_injection;

	void _Pack_Earth_Model_Attribute(unsigned int& word, int attr_idx, float val);
	float _Unpack_Earth_Model_Attribute(unsigned int word, int attr_idx);

	double _fq;

	float _vpvert_avgtop;
	float _vpvert_avgbot;
	
	int _sub_ix0;
	int _sub_ix1;
	int _sub_iy0;
	int _sub_iy1;
	int _sub_iz0;
	int _sub_iz1;

	int _nabc_sdx;
	int _nabc_sdy;
	int _nabc_top;
	int _nabc_bot;

	bool _nabc_sdx_extend;
	bool _nabc_sdy_extend;
	bool _nabc_top_extend;
	bool _nabc_bot_extend;

	bool _freesurface_enabled;
	bool _source_ghost_enabled;
	bool _receiver_ghost_enabled;
	int _sea_surface_z;

	Elastic_Shot** _shots;
	int _num_shots;

	int* _GPU_Devices;
	int _num_GPU_Devices;
	int _GPU_Pipes;
	int _Steps_Per_GPU;

	float _Courant_Factor;

	char* _tolower(char* str);
	void _swap_endian(float* v);
	bool _Check_Property(
			const char* property_name,
			Voxet_Property* prop,
			double const_val,
			size_t expected_file_size
			);
	bool _Calculate_Sub_Volume(
			const char* name,
			const char* parmfile_path,
			int line_num,
			int dim,
			double cell_size,
			double sub_min,
			double sub_max,
			char* sub_unit,
			int& ilu0,
			int& ilu1
			);
	bool _Calculate_ABC_Sponge(
			const char* name,
			const char* parmfile_path,
			int line_num,
			double abc_size,
			char* abc_unit,
			char* abc_flag,
			int dim,
			double cell_size,
			int& nabc_size,
			bool& nabc_flag
			);

	void _Read_Earth_Model(Elastic_Propagator* propagator);
};

#endif

