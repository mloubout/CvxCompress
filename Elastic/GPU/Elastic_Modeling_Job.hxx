#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_MODELING_JOB_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_MODELING_JOB_HXX

#include <cstdio>
#include <cassert>
#include <istream>
#include <map>
#include <list>
#include <string>

class Voxet;
class Voxet_Property;
class Voxet_Memory_Mapper;
class Elastic_Propagator;
class Elastic_Shot;
class Variable_Water_Velocity;

class Elastic_Modeling_Job
{
public:
	Elastic_Modeling_Job(int log_level, const char* parmfile_path, Voxet_Memory_Mapper* mapper, Variable_Water_Velocity* Vwxyzt);
	Elastic_Modeling_Job(int log_level, const char* parmfile_path, Voxet_Memory_Mapper* mapper, Variable_Water_Velocity* Vwxyzt, std::istream& fs);
	~Elastic_Modeling_Job();

	static void Print_Version_Information()
	{
		printf("\nCVX 3D Elastic Orthorhombic Modeling - v2.22 - 08/17/16\n\n");
	}

	bool Is_Valid() {return _Is_Valid;}
	int Get_Log_Level() {return _log_level;}

	int Get_Spatial_Order() {return _spatial_order;}

	//!
	//! Return true if this model includes a variable (in time and space) water column.
	//!
	bool Is_Vwxyzt();

	//!
	//! Return true if parmfile requested debug output of Vp.
	//!
	bool Vp_QC_Output_Enabled() {return _Vp_QC_Output;}
	
	//!
	//! Return true if parmfile requested debug output of a 1D Vp profile at source location.
	//!
	bool Vp_QC_1D_Profile_Enabled() {return _Vp_QC_1D_Profile;}

	const char* Get_EBCDIC_Header_Filename() {return _ebcdic_header_filename;}

	bool Anchor_Is_Source() {return _sub_origin == 0 ? true : false;}
	bool Anchor_Is_Model_Origin() {return _sub_origin == 1 ? true : false;}
	bool Anchor_Is_Source_And_Receivers() {return _sub_origin == 2 ? true : false;}
	bool Anchor_Is_Receivers() {return _sub_origin == 3 ? true : false;}
	void Compute_Subvolume();

	bool Web_Allowed() {return _web_allowed;}

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

	bool Lower_Q_Seafloor_Enabled() {return _lower_Q_seafloor_enabled;}
	bool Scholte_Only() {return _scholte_only;}
	bool Extend_Model_If_Necessary() {return _extend_model_if_necessary;}

	float Get_Courant_Factor() {return _Courant_Factor;}

	void Add_Shot(Elastic_Shot* shot);
	Elastic_Shot* Get_Shot(int souidx);
	Elastic_Shot* Get_Shot_By_Index(int idx);
	int Get_Number_Of_Shots() {return _num_shots;}

	float Get_FQ() {return _fq;}

	int Get_Number_Of_Earth_Model_Attributes();
	int Get_Earth_Model_Attribute_Index(const char* moniker);
	const char* Get_Earth_Model_Attribute_Moniker(int attr_idx);

	//
	// If silent is false and the coordinates are outside of the modeling domain, the program will be terminated.
	// If silent is true and the coordinates are outisde of the modeling domain, an error code is returned instead.
	//
	float Get_Earth_Model_Attribute(int attr_idx, int ix, int iy, int iz, bool silent, bool& error);
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

	bool Compute_Model_Water_Depth_And_Avg_Vp(int ix, int iy, float& model_water_depth, float& model_water_Vp);

	void Write_Earth_Model_Attribute_XZ_Slice(const char* path, int attr_idx, int iy);
	void Write_Earth_Model_XZ_Slice(const char* path, int iy);

	void Write_Earth_Model_Attribute_YZ_Slice(const char* path, int attr_idx, int ix);
	void Write_Earth_Model_YZ_Slice(const char* path, int ix);

	void Write_Earth_Model_Attribute_XY_Slice(const char* path, int attr_idx, int iz);
	void Write_Earth_Model_XY_Slice(const char* path, int iz);

	// a hack to test what happens if we manipulate Qp to mute Scholte waves
	// traveling along the sea floor.
	void Lower_Q_Seafloor(bool scholte_only);

	//
	// wf_type : 0->Vx, 1->Vy, 2->Vz, 3->P
	//
	void Write_YZ_Slice(const char* path, int wf_type, int ix);
	void Write_XZ_Slice(const char* path, int wf_type, int iy);
	void Write_XY_Slice(const char* path, int wf_type, int iz);

	//!
	//! Returns TRUE if voxet memory mapper was enabled in parmfile.
	//!
	bool Memory_Mapper_Enabled() {return _mapper_enabled;}
	//!
	//! Provide an optional memory mapper object for the voxet.
	//! If provided, all files in the voxet will be read through mmap()'ed files inside this object.
	//! This object is managed outside of this class so that it can be shared among multiple instances.
	//!
	void Set_Memory_Mapper(Voxet_Memory_Mapper* mapper) {_mapper=mapper;}
	Voxet_Memory_Mapper* Get_Memory_Mapper() {return _mapper;}

	Voxet* Get_Voxet() {return _voxet;}
	int Get_Number_of_Voxet_Properties() {return _num_em_props;}
	Voxet_Property* Get_Voxet_Property(int prop_idx)
	{
		if (prop_idx >= 0 && prop_idx < _num_em_props)
		{
			return _props[prop_idx];
		}
		return 0L;
	}

	int Get_Number_Of_GPU_Pipes();
	void Set_Number_Of_GPU_Pipes(int num_pipes);

	int Get_Steps_Per_GPU();
	void Set_Steps_Per_GPU(int num_timesteps);

	const int* Get_GPU_Devices();
	void Set_GPU_Devices(const int* device_ids, int num_devices);
	int Get_Number_Of_GPU_Devices();

	int Get_Number_Of_Parallel_Shots() {return _Num_Parallel_Shots;}

	//!
	//! Write selected fields from propagation earth model to a voxet.
	//! The propagation earth model is the earth model the propagator sees, after compression, variable water column etc.
	//!
	void Write_Propagation_Earth_Model_To_Voxet(const char* base_filename, std::list<int> fields);

	//!
	//! Write 1D vertical profile at source location.
	//!
	void Write_Leis_Debug_Trace(const char* base_filename, long ffid, double sx, double sy, std::list<int> fields);

	const int Attr_Idx_Vp;
	const int Attr_Idx_Vs;
	const int Attr_Idx_Density;
	const int Attr_Idx_Q;
	const int Attr_Idx_Dip;
	const int Attr_Idx_Azimuth;
	const int Attr_Idx_Rake;
	const int Attr_Idx_Delta1;
	const int Attr_Idx_Delta2;
	const int Attr_Idx_Delta3;
	const int Attr_Idx_Epsilon1;
	const int Attr_Idx_Epsilon2;
	const int Attr_Idx_Gamma1;
	const int Attr_Idx_Gamma2;

	std::string Get_Attribute_String(const int Attr_Idx)
	{
		if (Attr_Idx == Attr_Idx_Vp)
		{
			return "Vp";
		}
		else if (Attr_Idx == Attr_Idx_Vs)
		{
			return "Vs";
		}
		else if (Attr_Idx == Attr_Idx_Density)
		{
			return "Density";
		}
		else if (Attr_Idx == Attr_Idx_Q)
		{
			return "Q";
		}
		else if (Attr_Idx == Attr_Idx_Dip)
		{
			return "Dip";
		}
		else if (Attr_Idx == Attr_Idx_Azimuth)
		{
			return "Azimuth";
		}
		else if (Attr_Idx == Attr_Idx_Rake)
		{
			return "Rake";
		}
		else if (Attr_Idx == Attr_Idx_Delta1)
		{
			return "Delta1";
		}
		else if (Attr_Idx == Attr_Idx_Delta2)
		{
			return "Delta2";
		}
		else if (Attr_Idx == Attr_Idx_Delta3)
		{
			return "Delta3";
		}
		else if (Attr_Idx == Attr_Idx_Epsilon1)
		{
			return "Epsilon1";
		}
		else if (Attr_Idx == Attr_Idx_Epsilon2)
		{
			return "Epsilon2";
		}
		else if (Attr_Idx == Attr_Idx_Gamma1)
		{
			return "Gamma1";
		}
		else if (Attr_Idx == Attr_Idx_Gamma2)
		{
			return "Gamma2";
		}
		else
		{
			assert(false);
		}
	}

private:
	friend class Elastic_Propagator;
	Elastic_Propagator* _propagator;

	void _initialize(int log_level, const char* parmfile_path, std::istream& fs);
	bool _read_parmfile(int log_level, const char* parmfile_path, std::istream& fs);

	bool _Is_Valid;
	int _log_level;

	int _spatial_order;

	char* _ebcdic_header_filename;

	int _prop_nx;
	int _prop_ny;
	int _prop_nz;
	
	int _prop_x0;
	int _prop_y0;
	int _prop_z0;

	bool _mapper_enabled;
	Voxet_Memory_Mapper* _mapper;
	Voxet* _voxet;
	bool _Vp_QC_1D_Profile;
	bool _Vp_QC_Output;

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
	
	// 0->Source, 1->Volume, 2->Source & Receivers, 3->Receivers
	int _sub_origin;

	bool _sub_x_set;
	bool _sub_y_set;
	bool _sub_z_set;

	int _parm_sub_ix0;
	int _parm_sub_ix1;
	int _parm_sub_iy0;
	int _parm_sub_iy1;
	int _parm_sub_iz0;
	int _parm_sub_iz1;

	int _parm_nabc_sdx;
	int _parm_nabc_sdy;
	int _parm_nabc_top;
	int _parm_nabc_bot;

	bool _parm_nabc_sdx_extend;
	bool _parm_nabc_sdy_extend;
	bool _parm_nabc_top_extend;
	bool _parm_nabc_bot_extend;

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
	
	bool _lower_Q_seafloor_enabled;
	bool _extend_model_if_necessary;
	bool _scholte_only;

	Elastic_Shot** _shots;
	int _num_shots;

	int* _GPU_Devices;
	int _num_GPU_Devices;
	int _GPU_Pipes;
	int _Steps_Per_GPU;
	int _Num_Parallel_Shots;

	bool _web_allowed;

	Variable_Water_Velocity* _Vwxyzt_Computer;

	float _Courant_Factor;

	char* _tolower(char* str);
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
